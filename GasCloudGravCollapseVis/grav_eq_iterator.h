#pragma once
#include <thread>
#include <complex>
#include <functional>
#include "field_vis.h"
#include "weird_hacks.h"
#include "consts.h"
#include "pooled_thread.h"

using namespace std;

enum class d {
	dx, dy, dxy
};

struct greqit {
	dsfield density,
		x_speed,
		y_speed,
		x_force,
		y_force,
		energy;
	greqit(size_t n, double initial_temp = 1.) {
		density = dsfield(n, 0, continue_edge, c_continue_edge);
		x_speed = dsfield(n, 0, continue_edge, c_continue_edge);
		y_speed = dsfield(n, 0, continue_edge, c_continue_edge);
		x_force = dsfield(n, 0, continue_edge, c_continue_edge);
		y_force = dsfield(n, 0, continue_edge, c_continue_edge);
		energy = dsfield(n, initial_temp, continue_edge, c_continue_edge);
		for (auto &y : energy.fd) {
			for (auto &x : y) {
				x = initial_temp;
			}
		}
	}
	greqit(const dsfield &dsf_p, double initial_energy = 1.) {
		density = dsf_p;
		x_speed = dsfield(dsf_p.size(), 0, continue_edge, c_continue_edge);
		y_speed = dsfield(dsf_p.size(), 0, continue_edge, c_continue_edge);
		x_force = dsfield(dsf_p.size(), 0, continue_edge, c_continue_edge);
		y_force = dsfield(dsf_p.size(), 0, continue_edge, c_continue_edge);
		energy = dsfield(dsf_p.size(), initial_energy, continue_edge, c_continue_edge);
		for (auto &y : energy.fd) {
			for (auto &x : y) {
				x = initial_energy;
			}
		}
	}
	greqit(const dsfield &dsf_p, const dsfield &dsf_vx, const dsfield &dsf_vy, double initial_temp = 1.) {
		density = dsf_p;
		x_speed = dsf_vx;
		y_speed = dsf_vy;
		x_force = dsfield(dsf_p.size(), 0, continue_edge, c_continue_edge);
		y_force = dsfield(dsf_p.size(), 0, continue_edge, c_continue_edge);
		energy = dsfield(dsf_p.size(), initial_temp, continue_edge, c_continue_edge);
		for (auto &y : energy.fd) {
			for (auto &x : y) {
				x = initial_temp;
			}
		}
	}
	void swap(greqit &grei) {
		density.swap(grei.density);
		x_speed.swap(grei.x_speed);
		y_speed.swap(grei.y_speed);
		x_force.swap(grei.x_force);
		y_force.swap(grei.y_force);
		energy.swap(grei.energy);
	}
	inline size_t size() {
		return density.size();
	}
};

namespace d_h2 {
	inline double _first_finite_difference(double left1, double right1) {
		return (right1 - left1)*0.5;
	}
	inline double _first_finite_difference_x(const dsfield &dsf, int64_t x, int64_t y) {
		return _first_finite_difference(
			dsf.at(x - 1, y),
			dsf.at(x + 1, y)
		);
	}
	inline double _first_finite_difference_y(const dsfield &dsf, int64_t x, int64_t y) {
		return _first_finite_difference(
			dsf.at(x, y - 1),
			dsf.at(x, y + 1)
		);
	}
	inline double _second_finite_difference(double left1, double center, double right1) {
		return (left1 - 2.*center + right1);
	}
	inline double _second_finite_difference_xx(const dsfield &dsf, int64_t x, int64_t y) {
		return _second_finite_difference(
			dsf.at(x - 1, y),
			dsf.at(x, y),
			dsf.at(x + 1, y)
		);
	}
	inline double _second_finite_difference_yy(const dsfield &dsf, int64_t x, int64_t y) {
		return _second_finite_difference(
			dsf.at(x, y - 1),
			dsf.at(x, y),
			dsf.at(x, y + 1)
		);
	}
	inline double second_difference(const dsfield &dsf, int64_t x, int64_t y, d diff) {
		switch (diff) {
		case d::dx:
			return _second_finite_difference_xx(dsf, x, y);
		case d::dy:
			return _second_finite_difference_yy(dsf, x, y);
		}
		return 0;
	}
	inline double first_difference(const dsfield &dsf, int64_t x, int64_t y, d diff) {
		switch (diff) {
		case d::dx:
			return _first_finite_difference_x(dsf, x, y);
		case d::dy:
			return _first_finite_difference_y(dsf, x, y);
		}
	}
}

namespace d_h4 {
	inline double _first_finite_difference(double left2, double left1, double right1, const double& right2) {
		return (.25*left2 - 2 * left1 + 2 * right1 - .25*right2) / 3.;
	}
	inline double _first_finite_difference_x(const dsfield &dsf, int64_t x, int64_t y) {
		return _first_finite_difference(
			dsf.at(x - 2, y),
			dsf.at(x - 1, y),
			dsf.at(x + 1, y),
			dsf.at(x + 2, y)
		);
	}
	inline double _first_finite_difference_y(const dsfield &dsf, int64_t x, int64_t y) {
		return _first_finite_difference(
			dsf.at(x, y - 2),
			dsf.at(x, y - 1),
			dsf.at(x, y + 1),
			dsf.at(x, y + 2)
		);
	}
	inline double _second_finite_difference(double left2, double left1, double center, double right1, double right2) {
		return (-.25*left2 + 4. * left1 - 7.5*center + 4. * right1 - .25*right2) / 3.;
	}
	inline double _second_finite_difference_xx(const dsfield &dsf, int64_t x, int64_t y) {
		return _second_finite_difference(
			dsf.at(x - 2, y),
			dsf.at(x - 1, y),
			dsf.at(x, y),
			dsf.at(x + 1, y),
			dsf.at(x + 2, y)
		);
	}
	inline double _second_finite_difference_yy(const dsfield &dsf, int64_t x, int64_t y) {
		return _second_finite_difference(
			dsf.at(x, y - 2),
			dsf.at(x, y - 1),
			dsf.at(x, y),
			dsf.at(x, y + 1),
			dsf.at(x, y + 2)
		);
	}
	inline double second_difference(const dsfield &dsf, int64_t x, int64_t y, d diff) {
		switch (diff) {
		case d::dx:
			return _second_finite_difference_xx(dsf, x, y);
		case d::dy:
			return _second_finite_difference_yy(dsf, x, y);
		}
		return 0;
	}
	inline double first_difference(const dsfield &dsf, int64_t x, int64_t y, d diff) {
		if (diff == d::dx)
			return _first_finite_difference_x(dsf, x, y);
		else
			return _first_finite_difference_y(dsf, x, y);
	}
}

namespace d_op {
	inline double divergence_h4(const dsfield &dsf, int64_t x, int64_t y) {
		return d_h4::first_difference(dsf, x, y, d::dx) + d_h4::first_difference(dsf, x, y, d::dy);
	}
	inline double divergence_h2(const dsfield &dsf, int64_t x, int64_t y) {
		return d_h2::first_difference(dsf, x, y, d::dx) + d_h2::first_difference(dsf, x, y, d::dy);
	}
	inline double divergence_h4(const dsfield &x_dsf, const dsfield &y_dsf, int64_t x, int64_t y) {
		return d_h4::first_difference(x_dsf, x, y, d::dx) + d_h4::first_difference(y_dsf, x, y, d::dy);
	}
	inline double divergence_h2(const dsfield &x_dsf, const dsfield &y_dsf, int64_t x, int64_t y) {
		return d_h2::first_difference(x_dsf, x, y, d::dx) + d_h2::first_difference(y_dsf, x, y, d::dy);
	}
	inline void cross_product(double a1, double a2, double a3, double b1, double b2, double b3, double &out1, double &out2, double &out3) {
		out1 = a2 * b3 - a3 * b2;
		out2 = a1 * b3 - a3 * b1;
		out3 = a1 * b2 - a2 * b1;
	}
	inline double DF_h4_curl_2d_operator(const dsfield &x_dsf, const dsfield &y_dsf, int64_t x, int64_t y) {
		return d_h2::first_difference(y_dsf, x, y, d::dy) - d_h4::first_difference(x_dsf, x, y, d::dx);
	}
	inline double DF_h2_curl_2d_operator(const dsfield &x_dsf, const dsfield &y_dsf, int64_t x, int64_t y) {
		return d_h2::first_difference(y_dsf, x, y, d::dy) - d_h2::first_difference(x_dsf, x, y, d::dx);
	}
	inline void DF_h4_lamb_operator_at(const dsfield &x_dsf, const dsfield &y_dsf, int64_t x, int64_t y, double &out_x, double &out_y) {
		double curl_z = DF_h4_curl_2d_operator(x_dsf, y_dsf, x, y);
		cross_product(0, 0, curl_z, x_dsf.at(x, y), y_dsf.at(x, y), 0, out_x, out_y, curl_z);
	}
	inline void DF_h2_lamb_operator_at(const dsfield &x_dsf, const dsfield &y_dsf, int64_t x, int64_t y, double &out_x, double &out_y) {
		double curl_z = DF_h2_curl_2d_operator(x_dsf, y_dsf, x, y);
		cross_product(0, 0, curl_z, x_dsf.at(x, y), y_dsf.at(x, y), 0, out_x, out_y, curl_z);
	}
}

namespace gc_iter {
	typedef struct {
		double denergy, ddensity, dxspeed, dyspeed, xforce, yforce;
	} step_ans;
	step_ans operator*(double mul, step_ans sa) {
		sa.ddensity *= mul;
		sa.denergy *= mul;
		sa.dxspeed *= mul;
		sa.dyspeed *= mul;
		sa.xforce *= mul;
		sa.yforce *= mul;
		return sa;
	}
	step_ans operator+(step_ans fa, step_ans sa) {
		sa.ddensity += fa.ddensity;
		sa.denergy += fa.denergy;
		sa.dxspeed += fa.dxspeed;
		sa.dyspeed += fa.dyspeed;
		sa.xforce += fa.xforce;
		sa.yforce += fa.yforce;
		return sa;
	}
	void apply_step_ans(greqit& gr, const step_ans &sa, int64_t x, int64_t y, bool add_flag) {
		gr.density.at(x, y) = (add_flag * gr.density.at(x, y)) + sa.ddensity;
		gr.energy.at(x, y) = add_flag * gr.energy.at(x, y) + sa.denergy;
		gr.x_speed.at(x, y) = add_flag * gr.x_speed.at(x, y) + sa.dxspeed;
		gr.y_speed.at(x, y) = add_flag * gr.y_speed.at(x, y) + sa.dyspeed;
		gr.x_force.at(x, y) = add_flag * gr.x_force.at(x, y) + sa.xforce;
		gr.y_force.at(x, y) = add_flag * gr.y_force.at(x, y) + sa.yforce;
	}
	step_ans sum_step_ans(const greqit& gr, step_ans sa, int64_t x, int64_t y, bool add_flag) {
		sa.ddensity = add_flag * gr.density.at(x, y) + sa.ddensity;
		sa.denergy = add_flag * gr.energy.at(x, y) + sa.denergy;
		sa.dxspeed = add_flag * gr.x_speed.at(x, y) + sa.dxspeed;
		sa.dyspeed = add_flag * gr.y_speed.at(x, y) + sa.dyspeed;
		sa.xforce = add_flag * gr.x_force.at(x, y) + sa.xforce;
		sa.yforce = add_flag * gr.y_force.at(x, y) + sa.yforce;
		return sa;
	}

	int64_t fsize = 65;
	double time_step = 0.05;
	double adiabat = 1.67;
	double grav_const = 1;
	////inter-thread-ey variables////
	volatile bool iter_pause = true, iter_break = false;
	volatile int threads_count = max(thread::hardware_concurrency() - 1, 1u);
	greqit grei_base(fsize, 2);
	greqit grei_buffer(fsize, 2);
	greqit grei_i_buffer(fsize, 2);
	greqit grei_f0buffer(fsize, 2);
	std::recursive_mutex global_pause_lock;
	std::mutex local_lock;
	vector<pooled_thread*> threads;
	int __step_counter = 0;
	inline double rad_3(int64_t x, int64_t y) {
		if (!x && !y)
			return 0.;
		else
			return 1. / pow(x * x + y * y, 1.5);
	}
	inline double Fx(int64_t at_x, int64_t at_y) {
		double sum = 0;
		for (int64_t x = 0; x < fsize; x++) {
			for (int64_t y = 0; y < fsize; y++) {
				sum += rad_3(x - at_x, y - at_y)* (x - at_x)*grei_base.density.at(x, y);
			}
		}
		return sum;
	}
	inline double Fy(int64_t at_x, int64_t at_y) {
		double sum = 0;
		for (int64_t x = 0; x < fsize; x++) {
			for (int64_t y = 0; y < fsize; y++) {
				sum += rad_3(x - at_x, y - at_y)* (y - at_y)*grei_base.density.at(x, y);
			}
		}
		return sum;
	}
	inline step_ans iter_grei_at(int64_t x, int64_t y, const greqit& gfield) {
		step_ans ans{0};
		constexpr bool is_test = true;
		if (is_test) {
			ans.denergy =
				d_h4::second_difference(gfield.density, x, y, d::dx) + d_h4::second_difference(gfield.density, x, y, d::dy);
			ans.ddensity = gfield.energy.at(x, y);
			ans.dxspeed = 0;
			ans.dyspeed = 0;
			ans.xforce = 0;
			ans.yforce = 0;
		}
		else{
			ans.xforce = 0;// grav_const * Fx(x, y);
			ans.yforce = -0.01;// grav_const* Fy(x, y);

			ans.ddensity = (
				-(
					gfield.density.at(x, y) * (d_op::divergence_h2(gfield.x_speed, gfield.y_speed, x, y)) +
					gfield.x_speed.at(x, y) * d_h2::first_difference(gfield.density, x, y, d::dx) +
					gfield.y_speed.at(x, y) * d_h2::first_difference(gfield.density, x, y, d::dy)
				)
			);
			ans.denergy = (
				-(
					(adiabat - 1) * gfield.energy.at(x, y) * d_op::divergence_h2(gfield.x_speed, gfield.y_speed, x, y) + (
						gfield.x_speed.at(x, y) * d_h2::first_difference(gfield.energy, x, y, d::dx) +
						gfield.y_speed.at(x, y) * d_h2::first_difference(gfield.energy, x, y, d::dy)
					)
				)
			);
			ans.dxspeed = (
				-(
					(adiabat - 1) * (
						d_h2::first_difference(gfield.density, x, y, d::dx) * gfield.energy.at(x, y) / gfield.density.at(x, y) +
						d_h2::first_difference(gfield.energy, x, y, d::dx)
					) +	(
						gfield.x_speed.at(x, y) * d_h2::first_difference(gfield.x_speed, x, y, d::dx) +
						gfield.y_speed.at(x, y) * d_h2::first_difference(gfield.x_speed, x, y, d::dy)
					)
				)
				+
				gfield.x_force.at(x, y)
				);
			ans.dyspeed = (
				-(
					(adiabat - 1) * (
						d_h2::first_difference(gfield.density, x, y, d::dy) * gfield.energy.at(x, y) / gfield.density.at(x, y) +
						d_h2::first_difference(gfield.energy, x, y, d::dy)
					) + (
						gfield.x_speed.at(x, y) * d_h2::first_difference(gfield.y_speed, x, y, d::dx) +
						gfield.y_speed.at(x, y) * d_h2::first_difference(gfield.y_speed, x, y, d::dy)
					)
				)
				+
				gfield.y_force.at(x, y)
			);
		}
		return ans;
	}

	void create_iter_threads() {
		gc_iter::iter_pause = false;
		gc_iter::iter_break = false;

		typedef struct {
			int id;
			int64_t start, end, fsize;
			int *step_counter;
		} thread_info;

		local_lock.lock();
		for (int thid = 0; thid < threads_count; thid++) {
			const int64_t start = (thid)*fsize / (threads_count);
			const int64_t end = (thid + 1) * fsize / (threads_count);

			threads.push_back(new pooled_thread()); // executors
			auto t = threads.back()->__void_ptr_accsess();
			*t = (void*)(new thread_info{thid, start, end, fsize, &__step_counter});
			threads.back()->set_new_function([](void** void_ptr) {
				thread_info** pptr = (thread_info**)void_ptr;
				const int s_cnt = *(*pptr)->step_counter;
				constexpr bool is_dbg_printf = false;
				constexpr step_ans zero_sa{ 0,0,0,0,0,0 };
				global_pause_lock.lock();
				for (int64_t x = (*pptr)->start; x < (*pptr)->end; x++) {
					for (int64_t y = 0; y < (*pptr)->fsize; y++) {
						if (!s_cnt) {//prognosis
							auto dgrei = time_step * iter_grei_at(x, y, grei_base);
							apply_step_ans(grei_f0buffer, dgrei, x, y, false);
							dgrei = sum_step_ans(grei_base, dgrei, x, y, true);
							apply_step_ans(grei_i_buffer, 
									dgrei, 
								x, y, false);

							if(is_dbg_printf) 
								printf("%i: (%i:%i) D:%lf E:%lf X:%lf Y:%lf FX:%lf FY:%lf\n", s_cnt, x, y, dgrei.ddensity, dgrei.denergy, dgrei.dxspeed, dgrei.dyspeed, dgrei.xforce, dgrei.yforce);
						}
						else {//correction
							auto dbuffer = time_step * iter_grei_at(x, y, grei_buffer);
							auto f_0 = sum_step_ans(grei_f0buffer, zero_sa, x, y, true);
							dbuffer = 0.5 * (f_0 + dbuffer);
							apply_step_ans(grei_i_buffer,
								sum_step_ans(grei_base, dbuffer, x, y, true),
								x, y, false);

							if (is_dbg_printf) {
								printf("%i: (%i:%i) D:%lf E:%lf X:%lf Y:%lf FX:%lf FY:%lf\n", s_cnt, x, y, dbuffer.ddensity, dbuffer.denergy, dbuffer.dxspeed, dbuffer.dyspeed, dbuffer.xforce, dbuffer.yforce);
							}
						}
					}
				}
				global_pause_lock.unlock();
			});
			threads.back()->sign_awaiting();
		}
		threads.push_back(new pooled_thread());//observer
		threads.back()->set_new_awaiting_time(10);
		threads.back()->set_new_default_state(pooled_thread::state::waiting);
		threads.back()->set_new_function([](void** ptr) {
			for (auto ptr : threads)
				if (ptr != threads.back() && ptr->get_state() != pooled_thread::state::idle) {
					threads.back()->sign_awaiting();
					return;
				}

			global_pause_lock.lock();

			if (__step_counter == 5) {
				grei_base.swap(grei_i_buffer);
				__step_counter = 0;
				printf("prediction\n");
			}
			else {
				grei_i_buffer.swap(grei_buffer);
				__step_counter++;
				printf("correction\n");
			}

			global_pause_lock.unlock();
			for (auto ptr : threads)
				ptr->sign_awaiting();
		});
		threads.back()->sign_awaiting();
		local_lock.unlock();
	}
}

