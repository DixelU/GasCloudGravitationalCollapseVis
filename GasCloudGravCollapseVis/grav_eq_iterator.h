#pragma once
#include <thread>
#include "field_vis.h"
#include "weird_hacks.h"
#include "fftw3.h"

#pragma comment (lib, "libfftw3-3.lib")

using namespace std;


enum class d{
	dx, dy, dxy
};

struct greqit {
	dsfield density,
		x_speed,
		y_speed,
		x_force,
		y_force,
		temperature;
	greqit(size_t n, double initial_temp = 1.) {
		density = dsfield(n, 0, false);
		x_speed = dsfield(n, 0, false);
		y_speed = dsfield(n, 0, false);
		x_force = dsfield(n, 0);
		y_force = dsfield(n, 0);
		temperature = dsfield(n, initial_temp);
		for (auto &y : temperature.fd) {
			for (auto &x : y) {
				x = initial_temp;
			}
		}
	}
	greqit(const dsfield &dsf_p, double initial_temp = 1.) {
		density = dsf_p;
		x_speed = dsfield(dsf_p.size(), 0, false);
		y_speed = dsfield(dsf_p.size(), 0, false);
		x_force = dsfield(dsf_p.size(), 0);
		y_force = dsfield(dsf_p.size(), 0);
		temperature = dsfield(dsf_p.size(), initial_temp);
		for (auto &y : temperature.fd) {
			for (auto &x : y) {
				x = initial_temp;
			}
		}
	}
	greqit(const dsfield &dsf_p, const dsfield &dsf_vx, const dsfield &dsf_vy, double initial_temp = 1.) {
		density = dsf_p;
		x_speed = dsf_vx;
		y_speed = dsf_vy;
		x_force = dsfield(dsf_p.size(), 0);
		y_force = dsfield(dsf_p.size(), 0);
		temperature = dsfield(dsf_p.size(), initial_temp);
		for (auto &y : temperature.fd) {
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
		temperature.swap(grei.temperature);
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
	inline double _finite_difference_1st_order_x(dsfield &dsf, int64_t x, int64_t y) {
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
	inline double inv_rad(int64_t x, int64_t y) {
		if (!x && !y)
			return 2.;
		else
			return Q_rsqrt((double)(x*x+y*y));
	}
}

namespace gc_iter {
	int64_t fsize = 100;
	double iter_speed=0.01;
	double pressure_transform_coef = 1.37;
	double gas_viscosity = 0.01;
	double grav_const = 0.1;
	////inter-thread-ey variables////
	volatile bool iter_pause = true, iter_break = false;
	volatile int threads_count = max(thread::hardware_concurrency() - 1, 3u);
	greqit grei_base(fsize);
	greqit grei_buffer(fsize);
	dsfield fourier_coef, fourier_coef_buffer;
	volatile int* f_flags = nullptr;//if the thread completed cycle, it sets f_flags[thid] to one and waits until it will be set back to zero
	

	void iter_grei_at(int64_t x, int64_t y) {
		grei_buffer.density.at(x, y) = grei_base.density.at(x, y) + iter_speed * (
			grei_base.density.at(x, y)*(d_op::divergence_h2(grei_base.x_speed, grei_base.y_speed, x, y)) +
			grei_base.x_speed.at(x, y)*d_h2::DF_1_order(grei_base.density, x, y, d::dx) +
			grei_base.y_speed.at(x, y)*d_h2::DF_1_order(grei_base.density, x, y, d::dy)
			////ok
		);

		grei_buffer.temperature.at(x, y) = grei_base.temperature.at(x, y) + iter_speed * (
			d_op::divergence_h2(grei_buffer.x_speed, grei_buffer.y_speed, x, y)/grei_base.density.at(x,y)
		);
		grei_buffer.x_speed.at(x, y) = grei_base.x_speed.at(x, y) + iter_speed * (
			(
				grei_buffer.x_speed.at(x, y)*d_h2::DF_1_order(grei_buffer.x_speed, x, y, d::dx) +
				grei_buffer.y_speed.at(x, y)*d_h2::DF_1_order(grei_buffer.x_speed, x, y, d::dy)
			) +
			pressure_transform_coef * (
				d_h2::DF_1_order(grei_base.density, x, y, d::dx) * grei_base.temperature.at(x, y) +
				d_h2::DF_1_order(grei_base.temperature, x, y, d::dx) * grei_base.density.at(x, y)
			) +
			grei_base.x_force.at(x, y)
		);
		grei_buffer.y_speed.at(x, y) = grei_base.y_speed.at(x, y) + iter_speed * (
			(
				grei_buffer.x_speed.at(x, y)*d_h2::DF_1_order(grei_buffer.y_speed, x, y, d::dx) +
				grei_buffer.y_speed.at(x, y)*d_h2::DF_1_order(grei_buffer.y_speed, x, y, d::dy)
			) +
			pressure_transform_coef * (	
				d_h2::DF_1_order(grei_base.density, x, y, d::dy) * grei_base.temperature.at(x, y) +
				d_h2::DF_1_order(grei_base.temperature, x, y, d::dy) * grei_base.density.at(x, y)
			) + 
			grei_base.y_force.at(x, y)
		);
		////recalculate integral of force here...
		//grei_buffer.x_force.at(x, y) = 1;
		//grei_buffer.y_force.at(x, y) = -grav_const * (y - fsize / 2.);
		//grei_buffer.x_force.at(x, y) = -grav_const * (x - fsize / 2.);
	}
	void do_fast_fourier_x(int64_t x) {

	}
	void do_fast_fourier_y(int64_t y) {

	}
	void do_fourier_by_thread_id(int64_t thid) {

	}


	void create_iter_threads() {
		f_flags = new int[threads_count];
		for (int thid = 0; thid < threads_count; thid++) {
			f_flags[thid] = 1;
			thread th([&](const int thi) {
				const int64_t start = (thid - (thid >> 2) - 1)* grei_base.size() / (threads_count - (threads_count >> 2) - 1);
				const int64_t end = (thid - (thid >> 2)) * grei_base.size() / (threads_count - (threads_count >> 2) - 1);
				const int64_t fourier_id = ((thi & 3)) ? 0 : ((thi >> 2) + 1);
				while (!iter_break) {
					if (fourier_id) {//do a fourier transform for a part of pressure func

					}
					else {
						for (int64_t x = start; x < end; x++) {
							for (int64_t y = 0; y < grei_base.size(); y++) {
								
								iter_grei_at(x, y);
							}
						}
					}
					//Sleep(1000);
					f_flags[thi] = 0;
					while (!f_flags[thi] || iter_pause) {
						Sleep(5);
					}
				}
			}, thid);
			th.detach();
		}
		thread th_checker([&]() {
			bool flag;
			while (!iter_break) {
				flag = true;
				Sleep(10);
				for (int thi = 0; thi < threads_count; thi++) {
					if (f_flags[thi]) {
						flag = false;
						break;
					}
				}
				if (flag) {
					fourier_coef.swap(fourier_coef_buffer);
					grei_base.swap(grei_buffer);
					for (int thi = 0; thi < threads_count; thi++) {
						f_flags[thi] = 1;
					}
				}
			}
		});
		th_checker.detach();
	}


}

