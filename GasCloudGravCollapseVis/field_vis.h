#pragma once
#include <vector>
#include <utility>
#include <random>
#include "multidimentional_point.h"

namespace fv_utils {
	std::default_random_engine gen;
	std::uniform_real_distribution<double> distr(0., 1.);
	inline double rdrand() {
		return distr(gen);
	}
}

typedef std::vector<double> line;
typedef std::vector<Point<2>> pline;

typedef std::vector<line> field;
typedef std::vector<pline> pfield;

using fv_utils::rdrand;

typedef struct drawable_square_field {
	field fd;
	double cell_size;
	drawable_square_field(size_t n=100) {
		line ld(n,0);
		fd.assign(n, ld);
	}
	line& operator[](size_t N) {
		return fd[N];
	}
	void swap(drawable_square_field& dsf) {
		fd.swap(dsf.fd);
	}
	size_t size() const {
		return fd.size();
	}
} dsfield;

void randomise_dsfield(dsfield& dsf, int64_t rastr_rad, double offset = 0.5, double mul = 1., int64_t smoothing_factor=1) {
	dsfield tdsf=dsf;
	for (int64_t y = 0; y < dsf.size(); y++) {
		for (int64_t x = 0; x < dsf.size(); x++) {
			dsf[y][x] = (rdrand() - offset);
		}
	}
	while (smoothing_factor > 0) {
		for (int64_t y = 0; y < dsf.size(); y++) {
			for (int64_t x = 0; x < dsf.size(); x++) {
				double avg = 0;
				int64_t counter = 0;
				for (int64_t ox = -rastr_rad; ox <= rastr_rad; ox++) {
					for (int64_t oy = -rastr_rad; oy <= rastr_rad; oy++) {
						if (y + oy < 0 || y + oy >= dsf.size() || x + ox < 0 || x + ox >= dsf.size())
							continue;
						avg += dsf[y + oy][x + ox];
						counter++;
					}
				}
				tdsf[y][x] = mul * avg / counter;
			}
		}
		dsf.swap(tdsf);
		smoothing_factor--;
	}
}


inline void draw_dsfield(const dsfield& dsf, float center_xpos, float center_ypos, float range, float pixel_size, const float decrement = 1.) {
	float ym = range + center_ypos, xm = range + center_xpos, inverse, cell_size = 2 * range / (dsf.size());
	glPointSize(cell_size/pixel_size);
	glBegin(GL_POINTS);
	float y = -range + center_ypos; 
	float x = 0;
	for (auto &&it_y : dsf.fd) {
		x = -range + center_xpos;
		for (auto &&it_x : it_y) {
			inverse = -it_x;
			glColor3f(it_x*decrement,-it_x*inverse*decrement,inverse*decrement);
			glVertex2f(x, y);
			x += cell_size;
		}
		y += cell_size;
	}
	glEnd();
}