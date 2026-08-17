#pragma once

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <mutex>
#include <thread>
#include <vector>

#include "consts.h"
#include "field_vis.h"

struct greqit {
	dsfield density;
	dsfield x_speed;
	dsfield y_speed;
	dsfield x_force;
	dsfield y_force;
	dsfield energy;

	explicit greqit(size_t side, double initial_energy = 1.) :
		density(side, 0., continue_edge, c_continue_edge),
		x_speed(side, 0., continue_edge, c_continue_edge),
		y_speed(side, 0., continue_edge, c_continue_edge),
		x_force(side, 0., continue_edge, c_continue_edge),
		y_force(side, 0., continue_edge, c_continue_edge),
		energy(side, initial_energy, continue_edge, c_continue_edge) {
		for (auto& row : energy.fd) {
			for (double& value : row)
				value = initial_energy;
		}
	}

	explicit greqit(const dsfield& initial_density, double initial_energy = 1.) :
		greqit(initial_density.size(), initial_energy) {
		density = initial_density;
	}

	greqit(const dsfield& initial_density, const dsfield& initial_x_speed,
		const dsfield& initial_y_speed, double initial_energy = 1.) :
		greqit(initial_density.size(), initial_energy) {
		density = initial_density;
		x_speed = initial_x_speed;
		y_speed = initial_y_speed;
	}

	void swap(greqit& other) {
		density.swap(other.density);
		x_speed.swap(other.x_speed);
		y_speed.swap(other.y_speed);
		x_force.swap(other.x_force);
		y_force.swap(other.y_force);
		energy.swap(other.energy);
	}

	size_t size() const {
		return density.size();
	}
};

#include "stable_grid_solver.h"
