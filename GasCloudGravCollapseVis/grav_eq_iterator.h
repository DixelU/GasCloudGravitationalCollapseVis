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
		size_t dsf_edge = dsf.size() - 1;
		if (x > 0 && x < dsf_edge)
			return (
				dsf[y][x + 1] 
				- dsf[y][x - 1]
			)*0.5;
		else {
			int32_t pdx = (!x ? 1 : -1);
			return pdx*(
				-0.5 * dsf[y][x + 2*pdx] 
				+ 2.*dsf[y][x + pdx] 
				- 1.5*dsf[y][x]
			);
		}
	}
	inline double _finite_difference_1st_order_y(dsfield &dsf, size_t x, size_t y) {
		size_t dsf_edge = dsf.size() - 1;
		if (y > 0 && y < dsf_edge)
			return (
				dsf[y + 1][x] 
				- dsf[y - 1][x]
			)*0.5;
		else {
			int32_t pdy = (!y ? 1 : -1);
			return pdy*(
				-0.5*dsf[y + 2 * pdy][x] 
				+ dsf[y + pdy][x] * 2. 
				- 1.5*dsf[y][x]
			);
		}
	}
	inline double _finite_difference_2nd_order_xx(dsfield &dsf, size_t x, size_t y) {
		size_t dsf_edge = dsf.size() - 1;
		if (x > 0 && x < dsf_edge)
			return (
				dsf[y][x - 1] 
				- 2.*dsf[y][x] 
				+ dsf[y][x + 1]
			);
		else {
			int32_t pdx = (!x ? 1 : -1);
			return (
				2 * dsf[y][x] 
				- 5 * dsf[y][x + pdx] 
				+ 4 * dsf[y][x + pdx * 2] 
				- dsf[y][x + pdx * 3]
			);
		}
	}
	inline double _finite_difference_2nd_order_yy(dsfield &dsf, size_t x, size_t y) {
		size_t dsf_edge = dsf.size() - 1;
		if (y > 0 && y < dsf_edge)
			return (
				dsf[y - 1][x] 
				- 2.*dsf[y][x] 
				+ dsf[y + 1][x]
			);
		else {
			int32_t pdy = (!y ? 1 : -1);
			return (
				2 * dsf[y][x] 
				- 5 * dsf[y + pdy][x] 
				+ 4 * dsf[y + pdy * 2][x] 
				- dsf[y + pdy * 3][x]
			);
		}
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
		size_t dsf_edge = dsf.size() - 2;
		if (x > 1 && x < dsf_edge)
			return (
				.25 * dsf[y][x - 2] 
				- 2 * dsf[y][x - 1] 
				+ 2 * dsf[y][x + 1] 
				-.25* dsf[y][x + 2]
				) / 3.;
		else {
			return d_h2::_finite_difference_1st_order_x(dsf, x, y);
		}
	}
	inline double _finite_difference_1st_order_y(dsfield &dsf, size_t x, size_t y) {
		size_t dsf_edge = dsf.size() - 2;
		if (y > 1 && y < dsf_edge)
			return (
				.25 * dsf[y - 2][x]  
				- 2 * dsf[y - 1][x] 
				+ 2 * dsf[y + 1][x] 
				-.25* dsf[y + 2][x]
			) / 3.;
		else {
			return d_h2::_finite_difference_1st_order_y(dsf, x, y);
		}
	}
	inline double _finite_difference_2nd_order_xx(dsfield &dsf, size_t x, size_t y) {
		size_t dsf_edge = dsf.size() - 2;
		if (x > 1 && x < dsf_edge)
			return (
				- .25 * dsf[y][x - 2]  
				+ 16. * dsf[y][x - 1] 
				- 30. * dsf[y][x    ] 
				+ 16. * dsf[y][x + 1] 
				- .25 * dsf[y][x + 2] 
				) / 3.;
		else {
			return d_h2::_finite_difference_2nd_order_xx(dsf, x, y);
		}
	}
	inline double _finite_difference_2nd_order_yy(dsfield &dsf, size_t x, size_t y) {
		size_t dsf_edge = dsf.size() - 2;
		if (y > 1 && y < dsf_edge)
			return (
				- .25 * dsf[y - 2][x]
				+ 16. * dsf[y - 1][x]
				- 30. * dsf[y    ][x]
				+ 16. * dsf[y + 1][x]
				- .25 * dsf[y + 2][x]
				) / 3.;
		else {
			return d_h2::_finite_difference_2nd_order_yy(dsf, x, y);
		}
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