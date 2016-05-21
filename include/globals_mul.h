#pragma once
#include "types_mul.h"
namespace linear_boltzmann_transport_equation_solver {

	class UnknownTestId : public std::exception {
	public:
		UnknownTestId(size_t id);
		virtual const char* what() const noexcept;
	private:
		size_t id_;		
	};

	class Global {
	public:
		static size_t verbose() { return verbose_; }
		static void verbose_set(size_t veb) { verbose_ = veb; }



	private:
		static size_t verbose_;		
	};

	class AnalyticalSolution {
	public:
		static num_t solution(num_t x, num_t y, num_t z, num_t mu, num_t xi, num_t eta, num_t sigma);
		static num_t source(num_t x, num_t y, num_t z, num_t mu, num_t xi, num_t eta, num_t sigma);

		static size_t test_id();
		static void test_id_set(size_t id);
	private:
		static size_t test_id_;
	};
}