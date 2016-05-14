#pragma once
#include "types.h"
namespace dgn {
	namespace utilities {
		quadrant quadrant_of_angle(num_t mu, num_t xi, num_t eta);
		bool is_boundary(interface_direction idir, size_t np, size_t ix, size_t iy, size_t iz);
		bool is_boundary(interface_direction idir, size_t np, size_t ix, size_t iy);
		bool is_boundary(interface_direction idir, size_t np, size_t ix);
	}
}
