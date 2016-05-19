#pragma once
#include "types.h"
namespace dgn {
	class utilities {
	public:
		static quadrant quadrant_of_angle(num_t mu, num_t xi, num_t eta);
		static bool quadrant_x_positive_flag(quadrant quad);
		static bool quadrant_y_positive_flag(quadrant quad);
		static bool quadrant_z_positive_flag(quadrant quad);		
		static bool is_boundary(interface_direction idir, size_t np, size_t ix, size_t iy, size_t iz);
		static bool is_boundary(interface_direction idir, size_t np, size_t ix, size_t iy);
		static bool is_boundary(interface_direction idir, size_t np, size_t ix);
	};
}
