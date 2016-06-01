#pragma once
#include <xlib.h>
#include "types.h"
#include "utilities.h"
#include <exception>

using xlib::operator<<;

namespace dgn {
	class unknown_test_id : public std::exception {};

	class analytical_solution {
	public:
		static num_t solution(num_t x, num_t y, num_t z, num_t mu, num_t xi, num_t eta, num_t sigma);



		static num_t source(num_t x, num_t y, num_t z, num_t mu, num_t xi, num_t eta, num_t sigma);


		static void set_test_id(size_t id) {
			test_id_ = id;
		}

		
	private:
		static size_t test_id_;
	};
}