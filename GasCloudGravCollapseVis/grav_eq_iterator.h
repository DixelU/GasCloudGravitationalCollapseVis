#pragma once
#include <thread>
#include <complex>
#include "field_vis.h"
#include "weird_hacks.h"
#include "consts.h"

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
	inline double _finite_difference_1st_order(double left1, double right1) {
		return (right1 - left1)*0.5;
	}
	inline double _finite_difference_1st_order_x(dsfield &dsf, int64_t x, int64_t y) {
		return _finite_difference_1st_order(
			dsf.at(x - 1, y),
			dsf.at(x + 1, y)
		);
	}
	inline double _finite_difference_1st_order_y(dsfield &dsf, int64_t x, int64_t y) {
		return _finite_difference_1st_order(
			dsf.at(x, y - 1),
			dsf.at(x, y + 1)
		);
	}
	inline double _finite_difference_2nd_order(double left1, double center, double right1) {
		return (left1 - 2.*center + right1);
	}
	inline double _finite_difference_2nd_order_xx(dsfield &dsf, int64_t x, int64_t y) {
		return _finite_difference_2nd_order(
			dsf.at(x - 1, y),
			dsf.at(x, y),
			dsf.at(x + 1, y)
		);
	}
	inline double _finite_difference_2nd_order_yy(dsfield &dsf, int64_t x, int64_t y) {
		return _finite_difference_2nd_order(
			dsf.at(x, y - 1),
			dsf.at(x, y),
			dsf.at(x, y + 1)
		);
	}
	inline double DF_2_order(dsfield &dsf, int64_t x, int64_t y, d diff) {
		switch (diff) {
		case d::dx:
			return _finite_difference_2nd_order_xx(dsf, x, y);
		case d::dy:
			return _finite_difference_2nd_order_yy(dsf, x, y);
		}
		return 0;
	}
	inline double DF_1_order(dsfield &dsf, int64_t x, int64_t y, d diff) {
		switch (diff) {
		case d::dx:
			return _finite_difference_1st_order_x(dsf, x, y);
		case d::dy:
			return _finite_difference_1st_order_y(dsf, x, y);
		}
	}
}

namespace d_h4 {
	inline double _finite_difference_1st_order(double left2, double left1, double right1, const double& right2) {
		return (.25*left2 - 2 * left1 + 2 * right1 - .25*right2) / 3.;
	}
	inline double _finite_difference_1st_order_x(const dsfield &dsf, int64_t x, int64_t y) {
		return _finite_difference_1st_order(
			dsf.at(x - 2, y),
			dsf.at(x - 1, y),
			dsf.at(x + 1, y),
			dsf.at(x + 2, y)
		);
	}
	inline double _finite_difference_1st_order_y(dsfield &dsf, int64_t x, int64_t y) {
		return _finite_difference_1st_order(
			dsf.at(x, y - 2),
			dsf.at(x, y - 1),
			dsf.at(x, y + 1),
			dsf.at(x, y + 2)
		);
	}
	inline double _finite_difference_2nd_order(double left2, double left1, double center, double right1, double right2) {
		return (-.25*left2 + 4. * left1 - 7.5*center + 4. * right1 - .25*right2) / 3.;
	}
	inline double _finite_difference_2nd_order_xx(dsfield &dsf, int64_t x, int64_t y) {
		return _finite_difference_2nd_order(
			dsf.at(x - 2, y),
			dsf.at(x - 1, y),
			dsf.at(x, y),
			dsf.at(x + 1, y),
			dsf.at(x + 2, y)
		);
	}
	inline double _finite_difference_2nd_order_yy(dsfield &dsf, int64_t x, int64_t y) {
		return _finite_difference_2nd_order(
			dsf.at(x, y - 2),
			dsf.at(x, y - 1),
			dsf.at(x, y),
			dsf.at(x, y + 1),
			dsf.at(x, y + 2)
		);
	}
	inline double DF_2_order(dsfield &dsf, int64_t x, int64_t y, d diff) {
		switch (diff) {
		case d::dx:
			return _finite_difference_2nd_order_xx(dsf, x, y);
		case d::dy:
			return _finite_difference_2nd_order_yy(dsf, x, y);
		}
		return 0;
	}
	inline double DF_1_order(dsfield &dsf, int64_t x, int64_t y, d diff) {
		if (diff == d::dx)
			return _finite_difference_1st_order_x(dsf, x, y);
		else
			return _finite_difference_1st_order_y(dsf, x, y);
	}
}

namespace d_op {
	inline double divergence_h4(dsfield &dsf, int64_t x, int64_t y) {
		return d_h4::DF_1_order(dsf, x, y, d::dx) + d_h4::DF_1_order(dsf, x, y, d::dy);
	}
	inline double divergence_h2(dsfield &dsf, int64_t x, int64_t y) {
		return d_h2::DF_1_order(dsf, x, y, d::dx) + d_h2::DF_1_order(dsf, x, y, d::dy);
	}
	inline double divergence_h4(dsfield &x_dsf, dsfield &y_dsf, int64_t x, int64_t y) {
		return d_h4::DF_1_order(x_dsf, x, y, d::dx) + d_h4::DF_1_order(y_dsf, x, y, d::dy);
	}
	inline double divergence_h2(dsfield &x_dsf, dsfield &y_dsf, int64_t x, int64_t y) {
		return d_h2::DF_1_order(x_dsf, x, y, d::dx) + d_h2::DF_1_order(y_dsf, x, y, d::dy);
	}
	inline void cross_product(double a1, double a2, double a3, double b1, double b2, double b3, double &out1, double &out2, double &out3) {
		out1 = a2 * b3 - a3 * b2;
		out2 = a1 * b3 - a3 * b1;
		out3 = a1 * b2 - a2 * b1;
	}
	inline double DF_h4_curl_2d_operator(dsfield &x_dsf, dsfield &y_dsf, int64_t x, int64_t y) {
		return d_h2::DF_1_order(y_dsf, x, y, d::dy) - d_h4::DF_1_order(x_dsf, x, y, d::dx);
	}
	inline double DF_h2_curl_2d_operator(dsfield &x_dsf, dsfield &y_dsf, int64_t x, int64_t y) {
		return d_h2::DF_1_order(y_dsf, x, y, d::dy) - d_h2::DF_1_order(x_dsf, x, y, d::dx);
	}
	inline void DF_h4_lamb_operator_at(dsfield &x_dsf, dsfield &y_dsf, int64_t x, int64_t y, double &out_x, double &out_y) {
		double curl_z = DF_h4_curl_2d_operator(x_dsf, y_dsf, x, y);
		cross_product(0, 0, curl_z, x_dsf.at(x, y), y_dsf.at(x, y), 0, out_x, out_y, curl_z);
	}
	inline void DF_h2_lamb_operator_at(dsfield &x_dsf, dsfield &y_dsf, int64_t x, int64_t y, double &out_x, double &out_y) {
		double curl_z = DF_h2_curl_2d_operator(x_dsf, y_dsf, x, y);
		cross_product(0, 0, curl_z, x_dsf.at(x, y), y_dsf.at(x, y), 0, out_x, out_y, curl_z);
	}
}

namespace gc_iter {
	int64_t fsize = 32;
	double iter_speed = 0.01;
	double adiabat = 1.67;
	double grav_const = 1;
	////inter-thread-ey variables////
	volatile bool iter_pause = true, iter_break = false;
	volatile int threads_count = max(thread::hardware_concurrency() - 1, 3u);
	greqit grei_base(fsize, 2);
	greqit grei_buffer(fsize, 2);
	volatile bool* finised_flags = nullptr;
	std::recursive_mutex* rec_mutexes = nullptr;
	std::recursive_mutex global_pause_lock;
	std::mutex local_lock;
	double rad_3(int64_t x, int64_t y) {
		if (!x && !y)
			return 0.;
		else
			return 1. / pow(x * x + y * y, 1.5);
	}
	double Fx(int64_t at_x, int64_t at_y) {
		double sum = 0;
		for (int64_t x = 0; x < fsize; x++) {
			for (int64_t y = 0; y < fsize; y++) {
				sum += rad_3(x - at_x, y - at_y)* (x - at_x)*grei_base.density.at(x, y);
			}
		}
		return sum;
	}
	double Fy(int64_t at_x, int64_t at_y) {
		double sum = 0;
		for (int64_t x = 0; x < fsize; x++) {
			for (int64_t y = 0; y < fsize; y++) {
				sum += rad_3(x - at_x, y - at_y)* (y - at_y)*grei_base.density.at(x, y);
			}
		}
		return sum;
	}

	void iter_grei_at(int64_t x, int64_t y) {
		double t_T = 0;

		grei_buffer.x_force.at(x, y) = grav_const * Fx(x, y) / grei_base.density.at(x, y);
		grei_buffer.y_force.at(x, y) = grav_const * Fy(x, y) / grei_base.density.at(x, y);

		grei_buffer.density.at(x, y) = grei_base.density.at(x, y) + iter_speed * (
			grei_base.density.at(x, y) * (d_op::divergence_h2(grei_base.x_speed, grei_base.y_speed, x, y)) +
			grei_base.x_speed.at(x, y) * d_h2::DF_1_order(grei_base.density, x, y, d::dx) +
			grei_base.y_speed.at(x, y) * d_h2::DF_1_order(grei_base.density, x, y, d::dy)
		);
		grei_buffer.energy.at(x, y) = grei_base.energy.at(x, y) + iter_speed * (
			(adiabat - 1) * grei_base.energy.at(x, y) * d_op::divergence_h2(grei_base.x_speed, grei_base.y_speed, x, y) + (
				grei_base.x_speed.at(x, y) * d_h2::DF_1_order(grei_base.energy, x, y, d::dx) +
				grei_base.y_speed.at(x, y) * d_h2::DF_1_order(grei_base.energy, x, y, d::dy) 
			)
		);
 		grei_buffer.x_speed.at(x, y) = grei_base.x_speed.at(x, y) + iter_speed * (
			(adiabat - 1) * (
				d_h2::DF_1_order(grei_base.density, x, y, d::dx) * grei_base.energy.at(x, y) +
				d_h2::DF_1_order(grei_base.energy, x, y, d::dx) * grei_base.density.at(x, y)
				) / grei_base.density.at(x, y) +
				(
					grei_buffer.x_speed.at(x, y) * d_h2::DF_1_order(grei_buffer.x_speed, x, y, d::dx) +
					grei_buffer.y_speed.at(x, y) * d_h2::DF_1_order(grei_buffer.x_speed, x, y, d::dy)
				) 
			+
			grei_base.x_force.at(x, y)
			);
		grei_buffer.y_speed.at(x, y) = grei_base.y_speed.at(x, y) + iter_speed * (
			(adiabat - 1) * (
				d_h2::DF_1_order(grei_base.density, x, y, d::dy) * grei_base.energy.at(x, y) +
				d_h2::DF_1_order(grei_base.energy, x, y, d::dy) * grei_base.density.at(x, y)
				) / grei_base.density.at(x, y) +
				(
					grei_buffer.x_speed.at(x, y) * d_h2::DF_1_order(grei_buffer.y_speed, x, y, d::dx) +
					grei_buffer.y_speed.at(x, y) * d_h2::DF_1_order(grei_buffer.y_speed, x, y, d::dy)
				)
			+
			grei_base.y_force.at(x, y)
		);
	}

	void create_iter_threads() {
		gc_iter::iter_pause = false;
		gc_iter::iter_break = false;
		if (rec_mutexes)
			delete[] rec_mutexes;
		if (finised_flags)
			delete[] finised_flags;
		finised_flags = new bool[threads_count];
		rec_mutexes = new std::recursive_mutex[threads_count];
		for (int thi = 0; thi < threads_count; thi++)
			finised_flags[thi] = false;
		
		local_lock.lock();
		for (int thid = 0; thid < threads_count; thid++) {
			std::thread th([&](const int thi) {
				finised_flags[thi] = false;
				const int64_t start = (thi) * fsize / (threads_count);
				const int64_t end = (thi + 1) * fsize / (threads_count);
				while (!iter_break) {
					local_lock.lock();
					local_lock.unlock();
					finised_flags[thi] = false;
					rec_mutexes[thi].lock();
					for (int64_t x = start; x < end; x++) {
						for (int64_t y = 0; y < fsize; y++) {
							iter_grei_at(x, y);
						}
					}
					finised_flags[thi] = true;
					rec_mutexes[thi].unlock();
					//printf("unlocked\n");
					global_pause_lock.lock();
					global_pause_lock.unlock();
					//printf("released\n");
				}
				}, thid);
			th.detach();
		}
		std::thread th_checker([&]() {
			while (!iter_break) {
				local_lock.lock();
				local_lock.unlock();
				global_pause_lock.lock();
				bool flag = true;
				while (flag) {
					Sleep(5);
					flag = false;
					for (int thi = 0; thi < threads_count; thi++) {
						if (!finised_flags[thi]) {
							flag = true;
							break;
						}
					}
				}
				//printf("G: unlocked\n");
				grei_base.swap(grei_buffer);
				global_pause_lock.unlock();
				//printf("G: released\n");
				local_lock.lock();
				Sleep(13);
				local_lock.unlock();
			}
			});
		th_checker.detach();

		local_lock.unlock();
	}
	void create_single_thread() {
		thread th([&]() {
			while (true) {
				global_pause_lock.lock();
				for (int64_t x = 0; x < fsize; x++) {
					for (int64_t y = 0; y < fsize; y++) {
						iter_grei_at(x, y);
					}
				}
				grei_base.swap(grei_buffer);
				global_pause_lock.unlock();
			}
			});
		th.detach();
	}

}

