#pragma once

//
// The fluid update is a conservative finite-volume Euler step with Rusanov
// interface fluxes. It is deliberately first order: some numerical diffusion
// is preferable here to the unbounded centered advection used by the original
// course-work experiment. Gravity comes from a mean-subtracted periodic
// Poisson solve, which avoids the old 1/r singularity and O(N^4) force loop.

namespace gc_iter {

int64_t fsize = 100;
double time_step = 0.04;
double adiabat = 1.67;
double thermal_radiation = 0.01;
double grav_const = 0.005;
double initial_internal_energy = 0.05;
constexpr double density_floor = 1e-6;
constexpr double energy_floor = 1e-4;
constexpr double courant_number = 0.35;
constexpr int poisson_iterations = 80;
constexpr double poisson_relaxation = 2. / 3.;

std::atomic<bool> iter_pause{ true };
std::atomic<bool> simulation_running{ false };
greqit grei_base(static_cast<size_t>(fsize), initial_internal_energy);
std::mutex local_lock;
int step_counter = 0;
double simulated_time = 0.;
double last_time_step = 0.;

struct conserved {
	double density;
	double x_momentum;
	double y_momentum;
	double total_energy;
};

struct solver_workspace {
	std::vector<conserved> current;
	std::vector<conserved> next;
	std::vector<conserved> x_flux;
	std::vector<conserved> y_flux;
	std::vector<double> potential;
	std::vector<double> potential_buffer;
	std::vector<double> poisson_source;
	std::vector<double> x_acceleration;
	std::vector<double> y_acceleration;
	std::vector<double> internal_energy;

	explicit solver_workspace(size_t side = 0) {
		resize(side);
	}

	void resize(size_t side) {
		const size_t cells = side * side;
		current.resize(cells);
		next.resize(cells);
		x_flux.resize(cells);
		y_flux.resize(cells);
		potential.resize(cells, 0.);
		potential_buffer.resize(cells, 0.);
		poisson_source.resize(cells, 0.);
		x_acceleration.resize(cells, 0.);
		y_acceleration.resize(cells, 0.);
		internal_energy.resize(cells, initial_internal_energy);
	}
};

struct field_diagnostics {
	bool finite = true;
	double mass = 0.;
	double min_density = (std::numeric_limits<double>::max)();
	double max_density = 0.;
	double max_speed = 0.;
	double second_moment = 0.;
};

struct stability_report {
	bool success = false;
	int completed_steps = 0;
	double elapsed_time = 0.;
	field_diagnostics initial;
	field_diagnostics final;
};

inline size_t grid_index(int64_t x, int64_t y, size_t side) {
	const int64_t n = static_cast<int64_t>(side);
	x %= n;
	y %= n;
	if (x < 0)
		x += n;
	if (y < 0)
		y += n;
	return static_cast<size_t>(y) * side + static_cast<size_t>(x);
}

inline conserved add(const conserved& left, const conserved& right) {
	return {
		left.density + right.density,
		left.x_momentum + right.x_momentum,
		left.y_momentum + right.y_momentum,
		left.total_energy + right.total_energy
	};
}

inline conserved subtract(const conserved& left, const conserved& right) {
	return {
		left.density - right.density,
		left.x_momentum - right.x_momentum,
		left.y_momentum - right.y_momentum,
		left.total_energy - right.total_energy
	};
}

inline conserved multiply(double factor, const conserved& value) {
	return {
		factor * value.density,
		factor * value.x_momentum,
		factor * value.y_momentum,
		factor * value.total_energy
	};
}

inline void primitive_values(const conserved& value, double& density, double& x_speed,
	double& y_speed, double& pressure, double& sound_speed, double& internal_energy) {
	density = (std::max)(value.density, density_floor);
	x_speed = value.x_momentum / density;
	y_speed = value.y_momentum / density;
	const double kinetic_energy = 0.5 * density * (x_speed * x_speed + y_speed * y_speed);
	const double internal_energy_density = (std::max)(
		value.total_energy - kinetic_energy,
		density * energy_floor);
	internal_energy = internal_energy_density / density;
	pressure = (adiabat - 1.) * internal_energy_density;
	sound_speed = sqrt((std::max)(adiabat * pressure / density, 0.));
}

inline conserved physical_flux(const conserved& value, bool x_direction, double& signal_speed) {
	double density, x_speed, y_speed, pressure, sound_speed, internal_energy;
	primitive_values(value, density, x_speed, y_speed, pressure, sound_speed, internal_energy);
	const double normal_speed = x_direction ? x_speed : y_speed;
	signal_speed = fabs(normal_speed) + sound_speed;
	if (x_direction) {
		return {
			value.x_momentum,
			value.x_momentum * x_speed + pressure,
			value.y_momentum * x_speed,
			(value.total_energy + pressure) * x_speed
		};
	}
	return {
		value.y_momentum,
		value.x_momentum * y_speed,
		value.y_momentum * y_speed + pressure,
		(value.total_energy + pressure) * y_speed
	};
}

inline conserved rusanov_flux(const conserved& left, const conserved& right, bool x_direction) {
	double left_signal, right_signal;
	const conserved left_flux = physical_flux(left, x_direction, left_signal);
	const conserved right_flux = physical_flux(right, x_direction, right_signal);
	const double signal = (std::max)(left_signal, right_signal);
	return multiply(0.5, subtract(
		add(left_flux, right_flux),
		multiply(signal, subtract(right, left))));
}

inline void prepare_conserved(const greqit& state, solver_workspace& workspace) {
	const size_t side = state.size();
	if (workspace.current.size() != side * side)
		workspace.resize(side);
	for (size_t y = 0; y < side; ++y) {
		for (size_t x = 0; x < side; ++x) {
			const size_t at = y * side + x;
			const double raw_density = state.density.at(x, y);
			const double density = std::isfinite(raw_density) ?
				(std::max)(raw_density, density_floor) : density_floor;
			const double raw_x_speed = state.x_speed.at(x, y);
			const double raw_y_speed = state.y_speed.at(x, y);
			const double x_speed = std::isfinite(raw_x_speed) ? raw_x_speed : 0.;
			const double y_speed = std::isfinite(raw_y_speed) ? raw_y_speed : 0.;
			const double raw_energy = state.energy.at(x, y);
			const double energy = std::isfinite(raw_energy) ?
				(std::max)(raw_energy, energy_floor) : initial_internal_energy;
			workspace.current[at] = {
				density,
				density * x_speed,
				density * y_speed,
				density * energy + 0.5 * density * (x_speed * x_speed + y_speed * y_speed)
			};
		}
	}
}

inline void solve_periodic_gravity(const greqit& state, solver_workspace& workspace) {
	const size_t side = state.size();
	const size_t cells = side * side;
	double mean_density = 0.;
	for (const conserved& value : workspace.current)
		mean_density += value.density;
	mean_density /= static_cast<double>(cells);

	for (size_t at = 0; at < cells; ++at) {
		workspace.poisson_source[at] =
			4. * pi * grav_const * (workspace.current[at].density - mean_density);
	}

	for (int iteration = 0; iteration < poisson_iterations; ++iteration) {
		for (size_t y = 0; y < side; ++y) {
			for (size_t x = 0; x < side; ++x) {
				const size_t at = y * side + x;
				const double jacobi_value = 0.25 * (
					workspace.potential[grid_index(static_cast<int64_t>(x) - 1, y, side)] +
					workspace.potential[grid_index(static_cast<int64_t>(x) + 1, y, side)] +
					workspace.potential[grid_index(x, static_cast<int64_t>(y) - 1, side)] +
					workspace.potential[grid_index(x, static_cast<int64_t>(y) + 1, side)] -
					workspace.poisson_source[at]);
				workspace.potential_buffer[at] =
					(1. - poisson_relaxation) * workspace.potential[at] +
					poisson_relaxation * jacobi_value;
			}
		}
		workspace.potential.swap(workspace.potential_buffer);
	}

	for (size_t y = 0; y < side; ++y) {
		for (size_t x = 0; x < side; ++x) {
			const size_t at = y * side + x;
			workspace.x_acceleration[at] = -0.5 * (
				workspace.potential[grid_index(static_cast<int64_t>(x) + 1, y, side)] -
				workspace.potential[grid_index(static_cast<int64_t>(x) - 1, y, side)]);
			workspace.y_acceleration[at] = -0.5 * (
				workspace.potential[grid_index(x, static_cast<int64_t>(y) + 1, side)] -
				workspace.potential[grid_index(x, static_cast<int64_t>(y) - 1, side)]);
		}
	}
}

inline double stable_time_step(const solver_workspace& workspace) {
	double max_signal_speed = 0.;
	double max_acceleration = 0.;
	for (size_t at = 0; at < workspace.current.size(); ++at) {
		double density, x_speed, y_speed, pressure, sound_speed, internal_energy;
		primitive_values(workspace.current[at], density, x_speed, y_speed,
			pressure, sound_speed, internal_energy);
		max_signal_speed = (std::max)(max_signal_speed,
			(std::max)(fabs(x_speed), fabs(y_speed)) + sound_speed);
		max_acceleration = (std::max)(max_acceleration, sqrt(
			workspace.x_acceleration[at] * workspace.x_acceleration[at] +
			workspace.y_acceleration[at] * workspace.y_acceleration[at]));
	}
	const double fluid_limit = courant_number / (std::max)(max_signal_speed, 1e-12);
	const double gravity_limit = courant_number * sqrt(
		1. / (std::max)(max_acceleration, 1e-12));
	return (std::min)(time_step, (std::min)(fluid_limit, gravity_limit));
}

inline double advance_one_step(greqit& state, solver_workspace& workspace) {
	const size_t side = state.size();
	prepare_conserved(state, workspace);
	solve_periodic_gravity(state, workspace);
	const double dt = stable_time_step(workspace);

	for (size_t y = 0; y < side; ++y) {
		for (size_t x = 0; x < side; ++x) {
			const size_t at = y * side + x;
			workspace.x_flux[at] = rusanov_flux(
				workspace.current[at],
				workspace.current[grid_index(static_cast<int64_t>(x) + 1, y, side)],
				true);
			workspace.y_flux[at] = rusanov_flux(
				workspace.current[at],
				workspace.current[grid_index(x, static_cast<int64_t>(y) + 1, side)],
				false);
		}
	}

	for (size_t y = 0; y < side; ++y) {
		for (size_t x = 0; x < side; ++x) {
			const size_t at = y * side + x;
			const conserved flux_difference = add(
				subtract(workspace.x_flux[at],
					workspace.x_flux[grid_index(static_cast<int64_t>(x) - 1, y, side)]),
				subtract(workspace.y_flux[at],
					workspace.y_flux[grid_index(x, static_cast<int64_t>(y) - 1, side)]));
			conserved updated = subtract(
				workspace.current[at],
				multiply(dt, flux_difference));
			if (!std::isfinite(updated.density) || !std::isfinite(updated.x_momentum) ||
				!std::isfinite(updated.y_momentum) || !std::isfinite(updated.total_energy)) {
				updated = workspace.current[at];
			}

			updated.density = (std::max)(updated.density, density_floor);
			double x_speed = updated.x_momentum / updated.density;
			double y_speed = updated.y_momentum / updated.density;
			const double kinetic_energy = 0.5 * updated.density *
				(x_speed * x_speed + y_speed * y_speed);
			const double internal_energy = (std::max)(
				(updated.total_energy - kinetic_energy) / updated.density,
				energy_floor);

			x_speed += dt * workspace.x_acceleration[at];
			y_speed += dt * workspace.y_acceleration[at];
			updated.x_momentum = updated.density * x_speed;
			updated.y_momentum = updated.density * y_speed;
			updated.total_energy = updated.density * internal_energy +
				0.5 * updated.density * (x_speed * x_speed + y_speed * y_speed);
			workspace.next[at] = updated;
			workspace.internal_energy[at] = internal_energy;
		}
	}

	greqit next_state(side, initial_internal_energy);
	for (size_t y = 0; y < side; ++y) {
		for (size_t x = 0; x < side; ++x) {
			const size_t at = y * side + x;
			const conserved& updated = workspace.next[at];
			const double laplacian_energy =
				workspace.internal_energy[grid_index(static_cast<int64_t>(x) - 1, y, side)] +
				workspace.internal_energy[grid_index(static_cast<int64_t>(x) + 1, y, side)] +
				workspace.internal_energy[grid_index(x, static_cast<int64_t>(y) - 1, side)] +
				workspace.internal_energy[grid_index(x, static_cast<int64_t>(y) + 1, side)] -
				4. * workspace.internal_energy[at];
			const double internal_energy = (std::max)(
				workspace.internal_energy[at] +
				thermal_radiation * dt * laplacian_energy,
				energy_floor);
			next_state.density.at(x, y) = updated.density;
			next_state.x_speed.at(x, y) = updated.x_momentum / updated.density;
			next_state.y_speed.at(x, y) = updated.y_momentum / updated.density;
			next_state.energy.at(x, y) = internal_energy;
			next_state.x_force.at(x, y) = workspace.x_acceleration[at];
			next_state.y_force.at(x, y) = workspace.y_acceleration[at];
		}
	}
	state.swap(next_state);
	return dt;
}

inline void initialize_two_clouds(greqit& state) {
	const double side = static_cast<double>(state.size());
	const double width = side / 24.;
	for (size_t y = 0; y < state.size(); ++y) {
		for (size_t x = 0; x < state.size(); ++x) {
			const double first_x = static_cast<double>(x) - side * 0.35;
			const double first_y = static_cast<double>(y) - side * 0.35;
			const double second_x = static_cast<double>(x) - side * 0.65;
			const double second_y = static_cast<double>(y) - side * 0.65;
			state.density.at(x, y) = density_floor +
				exp(-(first_x * first_x + first_y * first_y) / (width * width)) +
				exp(-(second_x * second_x + second_y * second_y) / (width * width));
			state.energy.at(x, y) = initial_internal_energy;
			state.x_speed.at(x, y) = 0.;
			state.y_speed.at(x, y) = 0.;
			state.x_force.at(x, y) = 0.;
			state.y_force.at(x, y) = 0.;
		}
	}
}

inline field_diagnostics diagnose(const greqit& state) {
	field_diagnostics answer;
	const double center = (static_cast<double>(state.size()) - 1.) * 0.5;
	for (size_t y = 0; y < state.size(); ++y) {
		for (size_t x = 0; x < state.size(); ++x) {
			const double density = state.density.at(x, y);
			const double x_speed = state.x_speed.at(x, y);
			const double y_speed = state.y_speed.at(x, y);
			const double energy = state.energy.at(x, y);
			answer.finite = answer.finite && std::isfinite(density) &&
				std::isfinite(x_speed) && std::isfinite(y_speed) && std::isfinite(energy);
			answer.mass += density;
			answer.min_density = (std::min)(answer.min_density, density);
			answer.max_density = (std::max)(answer.max_density, density);
			answer.max_speed = (std::max)(answer.max_speed,
				sqrt(x_speed * x_speed + y_speed * y_speed));
			const double dx = static_cast<double>(x) - center;
			const double dy = static_cast<double>(y) - center;
			answer.second_moment += density * (dx * dx + dy * dy);
		}
	}
	answer.second_moment /= (std::max)(answer.mass, density_floor);
	return answer;
}

inline stability_report run_stability_check(int requested_steps = 200) {
	stability_report report;
	greqit state(static_cast<size_t>(fsize), initial_internal_energy);
	initialize_two_clouds(state);
	solver_workspace workspace(state.size());
	report.initial = diagnose(state);
	for (int step = 0; step < requested_steps; ++step) {
		report.elapsed_time += advance_one_step(state, workspace);
		report.completed_steps = step + 1;
		if (!diagnose(state).finite)
			break;
	}
	report.final = diagnose(state);
	const double relative_mass_error = fabs(report.final.mass - report.initial.mass) /
		(std::max)(report.initial.mass, density_floor);
	report.success = report.completed_steps == requested_steps && report.final.finite &&
		report.final.min_density >= density_floor && relative_mass_error < 0.01 &&
		report.final.max_speed < 20.;
	return report;
}

inline void create_iter_threads() {
	bool expected = false;
	if (!simulation_running.compare_exchange_strong(expected, true))
		return;
	std::thread([]() {
		solver_workspace workspace(static_cast<size_t>(fsize));
		while (simulation_running.load()) {
			if (iter_pause.load()) {
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				continue;
			}
			greqit local_state(static_cast<size_t>(fsize), initial_internal_energy);
			{
				std::lock_guard<std::mutex> lock(local_lock);
				local_state = grei_base;
			}
			const double dt = advance_one_step(local_state, workspace);
			{
				std::lock_guard<std::mutex> lock(local_lock);
				grei_base.swap(local_state);
				last_time_step = dt;
				simulated_time += dt;
				++step_counter;
			}
			if ((step_counter % 25) == 0) {
				printf("step:%i dt:%.6f time:%.3f\n",
					step_counter, dt, simulated_time);
			}
		}
	}).detach();
}

}
