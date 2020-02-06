#pragma once
#include <vector>
#include <utility>
#include <random>
#include <complex>
#include "multidimentional_point.h"

namespace fv_utils {
	std::default_random_engine gen;
	std::mt19937 mtrand(gen);
	std::uniform_real_distribution<double> distr(0., 1.);
	inline double rdrand() {
		return distr(mtrand);
	}
}

using line = std::vector<double>;
using pline = std::vector<Point<2>>;
using cline = std::vector<std::complex<double>>;

using field = std::vector<line>;
using pfield = std::vector<pline>;
using cfield = std::vector<cline>;

using fv_utils::rdrand;

auto constant_edge = [](field &f, int64_t x, int64_t y, double &val) -> double& {
	if (x >= 0 && x < f.size() && y >= 0 && y < f.size())
		return f[y][x];
	else 
		return val;
};
auto c_constant_edge = [](const field &f, int64_t x, int64_t y, const double &val) -> const double& {
	if (x >= 0 && x < f.size() && y >= 0 && y < f.size())
		return f[y][x];
	else 
		return val;
};

auto continue_edge = [](field &f, int64_t x, int64_t y, double &val) -> double& {
	int64_t size = f.size();
	return f[(y + size)% size][(x + size) % size];
};
auto c_continue_edge = [](const field &f, int64_t x, int64_t y, const double &val) -> const double& {
	int64_t size = f.size();
	return f[(y + size) % size][(x + size) % size];
};

typedef struct complex_square_field {
	cfield fd;
	std::complex<double> garbage_var, outside_val;
	bool ivdev;
	int64_t field_size; 
	complex_square_field(size_t n = 100, std::complex<double> outside_val = 0., bool is_value_defined_edge_value = true) : field_size(n), outside_val(outside_val) {
		cline cld(n, 0);
		fd.assign(n, cld);
		ivdev = is_value_defined_edge_value;
	}
	const std::complex<double>& at(int64_t x, int64_t y) const {
		if (x >= 0 && x < field_size && y >= 0 && y < field_size)
			return fd[y][x];
		else {
			if (ivdev)
				return outside_val;
			else {
				if (x < 0)x = 0;
				else x = field_size - 1;
				if (y < 0)y = 0;
				else y = field_size - 1;
				return fd[y][x];
			}
		}
	}
	std::complex<double>& at(int64_t x, int64_t y) {
		if (x >= 0 && x < field_size && y >= 0 && y < field_size)
			return fd[y][x];
		else {
			if (ivdev) {
				garbage_var = outside_val;
				return garbage_var;
			}
			else {
				if (x < 0)x = 0;
				else x = field_size - 1;
				if (y < 0)y = 0;
				else y = field_size - 1;
				garbage_var = fd[y][x];
				return garbage_var;
			}
		}
	}
	cline& operator[](int64_t N) {
		return fd[N];
	}
	const cline& operator[](int64_t N) const {
		return fd[N];
	}
	void swap(complex_square_field& dsf) {
		fd.swap(dsf.fd);
	}
	size_t size() const {
		return fd.size();
	}
} csfield;

typedef struct drawable_square_field {
	field fd;
	double garbage_var, outside_val;
	double&(*edge_handler)(field&, int64_t, int64_t, double&);
	const double&(*const_edge_handler)(const field&, int64_t, int64_t, const double&);
	int64_t field_size;
	double cell_size;
	drawable_square_field(size_t n = 100, double outside_val = 0.,
		double&(*edge_handler)(field&, int64_t, int64_t, double&) = continue_edge,
		const double&(*const_edge_handler)(const field&, int64_t, int64_t, const double&) = c_continue_edge
	) : field_size(n), outside_val(outside_val), edge_handler(edge_handler), const_edge_handler(const_edge_handler){
		line ld(n, 0);
		fd.assign(n, ld);
		garbage_var = 0;
	}
	const double& at(int64_t x, int64_t y) const {
		return const_edge_handler(fd, x, y, outside_val);
	}
	double& at(int64_t x, int64_t y) {
		garbage_var = outside_val;
		return edge_handler(fd, x, y, garbage_var);
	}
	line& operator[](int64_t N) {
		return fd[N];
	}
	const line& operator[](int64_t N) const {
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

inline float extended_edge(float V) {
	constexpr float r = 4.5, q = -1.3;
	return 1 - (r + q * V) / (V * V + r);
}
inline float extended_center(float V) {
	constexpr float r = 9.8, q = 0;
	return 1 - (r + q * V) / (V * V + r);
}

inline void draw_dsfield(const dsfield& dsf, float center_xpos, float center_ypos, float range, float pixel_size, float decrement = 1.) {
	float ym = range + center_ypos, xm = range + center_xpos, inverse, orig, cell_size = 2 * range / (dsf.size());
	glPointSize(cell_size / pixel_size);
	glBegin(GL_POINTS);
	float y = -range + center_ypos + cell_size / 2.;
	float x = 0;
	for (auto&& it_y : dsf.fd) {
		x = -range + center_xpos + cell_size / 2.;
		for (auto&& it_x : it_y) {
			float local_val = (it_x>0?1:-1)*(abs(it_x)) * decrement;
			orig = pow(extended_edge(local_val), 2);
			inverse = pow(extended_edge(-local_val), 2);
			glColor3f(orig*1.2, extended_center(local_val), inverse*1.2);
			glVertex2f(x, y);
			x += cell_size;
		}
		y += cell_size;
	}
	glEnd();
}