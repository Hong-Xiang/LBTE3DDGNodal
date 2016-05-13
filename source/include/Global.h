#pragma once
#include "Types.hpp"
#include <vector>
#include <array>

#define DEBUG
#define DOUBLEFLOAT

namespace DG3DNodal {

#ifdef DOUBLEFLOAT
	typedef double num_t;
//#define mkl_omatadd mkl_domatadd

#else
	typedef float num_t;
#endif

	typedef num_t* num_p;
	typedef const num_t* const_num_p;
	typedef size_t* size_p;

	typedef MKLEXT::Matrix<num_t> matrix_t;
	typedef MKLEXT::Vector<num_t> vector_t;

	inline num_t* calloc(size_t num) {
		return static_cast<num_t*>(mkl_calloc(num, sizeof(num_t), 16));
	}

	enum Quadrant {
		Quadrant1 = 0,
		Quadrant2 = 1,
		Quadrant3 = 2,
		Quadrant4 = 3,
		Quadrant5 = 4,
		Quadrant6 = 5,
		Quadrant7 = 6,
		Quadrant8 = 7
	};

	enum InterfaceDirection {
		InterfaceDirectionBack = 0,
		InterfaceDirectionFront = 1,
		InterfaceDirectionLeft = 2,
		InterfaceDirectionRight = 3,
		InterfaceDirectionUnder = 4,
		InterfaceDirectionTop = 5
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