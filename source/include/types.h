#pragma once
#include <array>
#include <xlib.h>

namespace dgn {
#ifdef USE_DOUBLE_PRECISION
	typedef double num_t;
#else
	typedef float num_t;
#endif
	typedef num_t* num_p;
	typedef const num_t* const_num_p;
	typedef size_t* size_p;

	typedef xlib::mkl_ext::matrix<num_t> matrix_t;
	typedef xlib::mkl_ext::vector<num_t> vector_t;
	
	struct vector3 {
		num_t x, y, z;
	};

	//quadrants and list for it
	enum class quadrant {
		quadrant1 = 0,
		quadrant2 = 1,
		quadrant3 = 2,
		quadrant4 = 3,
		quadrant5 = 4,
		quadrant6 = 5,
		quadrant7 = 6,
		quadrant8 = 7,
	};
	const std::array<quadrant, 8> quadrant_list = {
		quadrant::quadrant1,
		quadrant::quadrant2,
		quadrant::quadrant3,
		quadrant::quadrant4,
		quadrant::quadrant5,
		quadrant::quadrant6,
		quadrant::quadrant7,
		quadrant::quadrant8,
	};

	//interface directions and list for it
	enum class interface_direction {
		B = 0,
		F = 1,
		L = 2,
		R = 3,
		D = 4,
		U = 5
	};

	const std::array<interface_direction, 6> interface_direction_list = {
		interface_direction::B,
		interface_direction::F,
		interface_direction::L,
		interface_direction::R,
		interface_direction::D,
		interface_direction::U
	};
}