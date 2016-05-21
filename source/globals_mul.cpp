#include "globals_mul.h"

namespace linear_boltzmann_transport_equation_solver {
	size_t Global::verbose_ = 0;
	size_t AnalyticalSolution::test_id_ = 0;

	num_t AnalyticalSolution::solution(num_t x, num_t y, num_t z, num_t mu, num_t xi, num_t eta, num_t sigma)
	{
		switch (test_id_)
		{
		case 0:
			return 1.0;
			break;
		default:
			throw(UnknownTestId(test_id_));
			break;
		}
	}

	num_t AnalyticalSolution::source(num_t x, num_t y, num_t z, num_t mu, num_t xi, num_t eta, num_t sigma)
	{
		switch (test_id_)
		{
		case 0:
			return sigma;
			break;
		default:
			throw(UnknownTestId(test_id_));
			break;
		}
	}

	size_t AnalyticalSolution::test_id()
	{
		return test_id_;
	}

	void AnalyticalSolution::test_id_set(size_t id)
	{
		test_id_ = id;
	}

	UnknownTestId::UnknownTestId(size_t id)
		:	id_(id)
	{		
	}

	const char* UnknownTestId::what() const noexcept
	{
		std::string msg = "Unknown test id: " + xlib::i2s(id_) + ".";
		return msg.c_str();
	}

}