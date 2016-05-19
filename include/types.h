#pragma once
#include <array>
#include <xlib.h>

namespace dgn {
	enum class dgn_statues {
		dgn_success = 0,
		dgn_failure = 1
	};

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

	class vector3 {
	public:
		vector3() {}
		vector3(num_t x, num_t y, num_t z) 
		: x_(x), y_(y), z_(z)
		{}
		num_t x() const { return x_; }
		num_t y() const { return y_; }
		num_t z() const { return z_; }
	private:
		num_t x_, y_, z_;
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
	const std::array<std::string, 8> quadrant_name_list = {
		"quadrant1",
		"quadrant2",
		"quadrant3",
		"quadrant4",
		"quadrant5",
		"quadrant6",
		"quadrant7",
		"quadrant8"
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

	const std::array<std::string, 6> interface_direction_name_list = {
		"Back",
		"Front",
		"Left",
		"Right",
		"Down",
		"Up"
	};


	static const std::vector<size_t> INTERFACE_TOTAL = { 1, 2, 4, 6 };
	static const std::vector<size_t> QUADRANT_TOTAL = { 1, 2, 4, 8 };

	//normal directions 
	enum class surface_direction {
		X = 0,
		Y = 1,
		Z = 2
	};

	enum class element_direction {
		pre = 0,
		inc = 1
	};
	
}