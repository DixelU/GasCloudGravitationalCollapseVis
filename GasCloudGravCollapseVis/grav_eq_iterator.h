#pragma once
#include "field_vis.h"

enum class d{
	dx, dy
};

struct greqit {
	dsfield p,vx,vy,fi,s;
	greqit(size_t n) {
		p = dsfield(n);
		vx = dsfield(n);
		vy = dsfield(n);
		fi = dsfield(n);
		s = dsfield(n);
	}
	greqit(const dsfield &dsf_p) {
		p = dsf_p;
		vx = dsfield(dsf_p.size());
		vy = dsfield(dsf_p.size());
		fi = dsfield(dsf_p.size());
		s = dsfield(dsf_p.size());
	}
	greqit(const dsfield &dsf_p, const dsfield &dsf_vx, const dsfield &dsf_vy) {
		p = dsf_p;
		vx = dsf_vx;
		vy = dsf_vy;
		fi = dsfield(dsf_p.size());
		s = dsfield(dsf_p.size());
	}
};

namespace d_h2 {
	inline double _finite_difference_1st_order_x(dsfield &dsf, size_t x, size_t y) {
		return (
			dsf.at(x + 1, y)
			- dsf.at(x - 1, y)
		)*0.5;
	}
	inline double _finite_difference_1st_order_y(dsfield &dsf, size_t x, size_t y) {
		return (
			dsf.at(x + 1, y)
			- dsf.at(x - 1, y)
		)*0.5;
	}
	inline double _finite_difference_2nd_order_xx(dsfield &dsf, size_t x, size_t y) {
		return (
			dsf.at(x - 1, y)
			- 2.*dsf.at(x, y)
			+ dsf.at(x + 1, y)
		);
	}
	inline double _finite_difference_2nd_order_yy(dsfield &dsf, size_t x, size_t y) {
		return (
			dsf.at(x, y - 1)
			- 2.*dsf.at(x, y)
			+ dsf.at(x, y + 1)
		);
	}
	inline double DF_2_order(dsfield &dsf, size_t x, size_t y, d diff) {
		switch (diff) {
		case d::dx:
			return _finite_difference_2nd_order_xx(dsf, x, y);
		case d::dy:
			return _finite_difference_2nd_order_yy(dsf, x, y);
		}
	}
	inline double DF_1_order(dsfield &dsf, size_t x, size_t y, d diff) {
		if (diff == d::dx)
			return _finite_difference_1st_order_x(dsf, x, y);
		else
			return _finite_difference_1st_order_y(dsf, x, y);
	}
}

namespace d_h4 {
	inline double _finite_difference_1st_order_x(dsfield &dsf, size_t x, size_t y) {
		return (
			.25 * dsf.at(x - 2, y)
			- 2 * dsf.at(x - 1, y)
			+ 2 * dsf.at(x + 1, y)
			-.25* dsf.at(x + 2, y)
		) / 3.;
	}
	inline double _finite_difference_1st_order_y(dsfield &dsf, size_t x, size_t y) {
		return (
			.25 * dsf.at(x, y- 2)
			- 2 * dsf.at(x, y - 1)
			+ 2 * dsf.at(x, y + 1)
			-.25* dsf.at(x, y + 2)
		) / 3.;
	}
	inline double _finite_difference_2nd_order_xx(dsfield &dsf, size_t x, size_t y) {
		return (
			-.25 * dsf.at(x - 2, y)
			+ 4. * dsf.at(x - 1, y)
			- 7.5 * dsf.at(x, y)
			+ 4. * dsf.at(x + 1, y)
			- .25 * dsf.at(x + 2, y)
			) / 3.;
	}
	inline double _finite_difference_2nd_order_yy(dsfield &dsf, size_t x, size_t y) {
		return (
			- .25 * dsf.at(x, y - 2)
			+ 4. * dsf.at(x, y - 1)
			- 7.5 * dsf.at(x, y)
			+ 4. * dsf.at(x, y + 1)
			- .25 * dsf.at(x, y + 2)
		) / 3.;
	}
	inline double DF_2_order(dsfield &dsf, size_t x, size_t y, d diff) {
		switch (diff) {
		case d::dx:
			return _finite_difference_2nd_order_xx(dsf, x, y);
		case d::dy:
			return _finite_difference_2nd_order_yy(dsf, x, y);
		}
	}
	inline double DF_1_order(dsfield &dsf, size_t x, size_t y, d diff) {
		if (diff == d::dx)
			return _finite_difference_1st_order_x(dsf, x, y);
		else
			return _finite_difference_1st_order_y(dsf, x, y);
	}
}