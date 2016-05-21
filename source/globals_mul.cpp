#include "globals_mul.h"

namespace linear_boltzmann_transport_equation_solver {
	size_t Global::verbose_ = 0;
	size_t AnalyticalSolution::test_id_ = 0;

	num_t AnalyticalSolution::solution(num_t x, num_t y, num_t z, num_t mu, num_t xi, num_t eta, num_t sigma)
	{
		switch (test_id_)
		{
			num_t r2, r, ans;
		case 0:
			return 1.0;
			break;
		case 1:
			
			r2 = x*x + y*y + z*z;
			r = std::pow(r2, 0.5);
			ans = std::exp(-sigma*r);
			return ans;
		case 2:
			return std::exp(-sigma*(x*x + y*y + z*z));
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
			num_t r2, r, ans;
		case 0:
			return sigma-0.5;
			break;
		case 1:
			
			r2 = x*x + y*y + z*z;
			r = std::pow(r2, 0.5);
			if (r == 0)
				ans = sigma*std::exp(-sigma*r);
			else
				ans = mu*sigma*x / r*std::exp(-sigma*r) + xi*sigma*y / r*std::exp(-sigma*r) + eta*sigma*z / r*std::exp(-sigma*r) + sigma*std::exp(-sigma*r);
			return ans;
			break;
		case 2:
			r2 = x*x + y*y + z*z;
			ans = mu*(-2)*sigma*x*std::exp(-sigma*r2) + xi*(-2)*sigma*y*std::exp(-sigma*r2) + eta*(-2)*sigma*z*std::exp(-sigma*r2) + (sigma-0.5)*std::exp(-sigma*r2);
			return ans;
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