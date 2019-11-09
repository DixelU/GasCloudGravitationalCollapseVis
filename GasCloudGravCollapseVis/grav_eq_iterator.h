#pragma once
#include "field_vis.h"
#include "weird_hacks.h"

enum class d{
	dx, dy, dxy
};

struct greqit {
	dsfield density,x_speed,y_speed,grav_potential,entropy;
	greqit(size_t n) {
		density = dsfield(n);
		x_speed = dsfield(n);
		y_speed = dsfield(n);
		grav_potential = dsfield(n);
		entropy = dsfield(n);
	}
	greqit(const dsfield &dsf_p) {
		density = dsf_p;
		x_speed = dsfield(dsf_p.size());
		y_speed = dsfield(dsf_p.size());
		grav_potential = dsfield(dsf_p.size());
		entropy = dsfield(dsf_p.size());
	}
	greqit(const dsfield &dsf_p, const dsfield &dsf_vx, const dsfield &dsf_vy) {
		density = dsf_p;
		x_speed = dsf_vx;
		y_speed = dsf_vy;
		grav_potential = dsfield(dsf_p.size());
		entropy = dsfield(dsf_p.size());
	}
	void init() {

	}
	void swap(greqit &grei) {
		density.swap(grei.density);
		x_speed.swap(grei.x_speed);
		y_speed.swap(grei.y_speed);
		grav_potential.swap(grei.grav_potential);
		entropy.swap(grei.entropy);
	}
	inline size_t size() {
		return density.size();
	}
};

namespace d_h2 {
	inline double _finite_difference_1st_order(const double& left1, const double& right1) {
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
	inline double _finite_difference_2nd_order(const double& left1, const double& center, const double& right1) {
		return (left1 - 2.*center + right1);
	}
	inline double _finite_mixed_difference_2nd_order(
		const double& x_p45deg/*x+1 y+1*/,
		const double& x_p135deg/*x+1 y-1*/, 
		const double& x_m135deg/*x-1 y-1*/,
		const double& x_m45deg/*x-1 y+1*/) {
		return ((x_p45deg - x_p135deg) - (x_m45deg - x_m135deg))*0.25;
	}
	inline double _finite_difference_2nd_order_xx(dsfield &dsf, int64_t x, int64_t y) {
		return _finite_difference_2nd_order(
			dsf.at(x - 1, y),
			dsf.at(x, y),
			dsf.at(x + 1, y)
		);
	}
	inline double _finite_difference_2nd_order_xy(dsfield &dsf, int64_t x, int64_t y) {
		return _finite_mixed_difference_2nd_order(
			dsf.at(x + 1, y + 1),
			dsf.at(x + 1, y - 1),
			dsf.at(x - 1, y - 1),
			dsf.at(x - 1, y + 1)
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
		case d::dxy:
			return _finite_difference_2nd_order_xy(dsf, x, y);
		}
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
	inline double _finite_difference_1st_order(const double& left2, const double& left1, const double& right1, const double& right2) {
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
	inline double _finite_difference_2nd_order(const double& left2, const double& left1, const double& center, const double& right1, const double& right2) {
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
	}
	inline double DF_1_order(dsfield &dsf, int64_t x, int64_t y, d diff) {
		if (diff == d::dx)
			return _finite_difference_1st_order_x(dsf, x, y);
		else
			return _finite_difference_1st_order_y(dsf, x, y);
	}
}

namespace gc_iter {
	double iter_speed=0.05;
	double pressure_transform_coef = 5;
	double grav_const = 0.001;
	inline void _1stline_h4(greqit &base, greqit &buffer, int64_t x, int64_t y) {
		buffer.density.at(x, y) =
			base.density.at(x, y) - iter_speed * (
				d_h4::_finite_difference_1st_order(
					base.density.at(x - 2, y)*(base.x_speed.at(x - 2, y) + base.y_speed.at(x - 2, y)),
					base.density.at(x - 1, y)*(base.x_speed.at(x - 1, y) + base.y_speed.at(x - 1, y)),
					base.density.at(x + 1, y)*(base.x_speed.at(x + 1, y) + base.y_speed.at(x + 1, y)),
					base.density.at(x + 2, y)*(base.x_speed.at(x + 2, y) + base.y_speed.at(x + 2, y))
				)
				+
				d_h4::_finite_difference_1st_order(
					base.density.at(x, y - 2)*(base.x_speed.at(x, y - 2) + base.y_speed.at(x, y - 2)),
					base.density.at(x, y - 1)*(base.x_speed.at(x, y - 1) + base.y_speed.at(x, y - 1)),
					base.density.at(x, y + 1)*(base.x_speed.at(x, y + 1) + base.y_speed.at(x, y + 1)),
					base.density.at(x, y + 2)*(base.x_speed.at(x, y + 2) + base.y_speed.at(x, y + 2))
				)
			);
	}
	inline void _1stline_h2(greqit &base, greqit &buffer, int64_t x, int64_t y) {
		buffer.density.at(x, y) =
			base.density.at(x, y) - iter_speed * (
				d_h2::_finite_difference_1st_order(
					base.density.at(x - 1, y)*(base.x_speed.at(x - 1, y) + base.y_speed.at(x - 1, y)),
					base.density.at(x + 1, y)*(base.x_speed.at(x + 1, y) + base.y_speed.at(x + 1, y))
				)
				+
				d_h2::_finite_difference_1st_order(
					base.density.at(x, y - 1)*(base.x_speed.at(x, y - 1) + base.y_speed.at(x, y - 1)),
					base.density.at(x, y + 1)*(base.x_speed.at(x, y + 1) + base.y_speed.at(x, y + 1))
				)
			);
	}
	inline void _2ndline_h4(greqit &base, greqit &buffer, int64_t x, int64_t y) {
		///  
		///  P = ro*T; T = ro*c;
		///
		buffer.x_speed.at(x, y) =
			base.x_speed.at(x, y) - iter_speed * (
				base.x_speed.at(x, y)*(d_h4::DF_1_order(base.x_speed, x, y, d::dx) + d_h4::DF_1_order(base.x_speed, x, y, d::dy)) +
				2.*pressure_transform_coef* d_h4::DF_1_order(base.x_speed, x, y, d::dx) + // 1/ro * gradP -> 1/ro * grad(ro^2) -> 1/ro * 2ro*grad(ro) -> 2*grad(ro)
				d_h4::DF_1_order(base.grav_potential, x, y, d::dx)
			);

		buffer.y_speed.at(x, y) =
			base.y_speed.at(x, y) - iter_speed * (
				base.x_speed.at(x, y)*( d_h4::DF_1_order(base.x_speed, x, y, d::dy)) +
				2.*pressure_transform_coef* d_h4::DF_1_order(base.y_speed, x, y, d::dy) + 
				d_h4::DF_1_order(base.grav_potential, x, y, d::dy)
			);
	}
	inline void _2ndline_h2(greqit &base, greqit &buffer, int64_t x, int64_t y) {
		buffer.x_speed.at(x, y) =
			base.x_speed.at(x, y) - iter_speed * (
				base.x_speed.at(x, y)*(d_h2::DF_1_order(base.x_speed, x, y, d::dx)) +
				2.*pressure_transform_coef* d_h2::DF_1_order(base.x_speed, x, y, d::dx) + 
				d_h2::DF_1_order(base.grav_potential, x, y, d::dx)
			);

		buffer.y_speed.at(x, y) =
			base.y_speed.at(x, y) - iter_speed * (
				base.x_speed.at(x, y)*(d_h2::DF_1_order(base.x_speed, x, y, d::dy)) +
				2.*pressure_transform_coef* d_h2::DF_1_order(base.y_speed, x, y, d::dy) + 
				d_h2::DF_1_order(base.grav_potential, x, y, d::dy)
			);
	}
	inline void _3rdline(greqit &base, greqit &buffer, int64_t x, int64_t y) {
		double sum = 0.;
		for (int64_t i = 0; i < base.density.size(); i++) {
			for (int64_t j = 0; j < base.density.size(); j++) {
				if (!((x - i) | (y - j)))
					continue;
				sum += base.density.at(i, j)/sqrt((float)((x - i)*(x - i) + (y - j)*(y - j)));
			}
		}
		buffer.grav_potential.at(x, y) = sum * grav_const;
	}
	inline void _4thline_h4(greqit &base, greqit &buffer, int64_t x, int64_t y) {

	}
	inline void _4thline_h2(greqit &base, greqit &buffer, int64_t x, int64_t y) {

	}
	inline void process_single_point_h2(greqit &base, greqit &buffer, int64_t x, int64_t y) {
		_1stline_h2(base, buffer, x, y);
		_2ndline_h2(base, buffer, x, y);
		_3rdline(base, buffer, x, y);
		_4thline_h2(base, buffer, x, y);
	}
	inline void process_single_point_h4(greqit &base, greqit &buffer, int64_t x, int64_t y) {
		_1stline_h4(base, buffer, x, y);
		_2ndline_h4(base, buffer, x, y);
		_3rdline(base, buffer, x, y);
		_4thline_h4(base, buffer, x, y);
	}
}