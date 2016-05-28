#include "global.h"
#include <fstream>

namespace discontinues_galerkin_nodal_solver {
  num_t AnalyticalSolution::solution(num_t x, num_t y, num_t z, num_t mu, num_t xi, num_t eta, num_t sigma, size_t id)
  {
    num_t ans = 0;
    num_t r2 = 0.0;
    num_t r = 0.0;
    switch (id)
    {
    case 0:
      ans = 1.0;
      break;
    case 1:
      ans = x + y + z;
      break;
    case 2:
      ans = x * x + y * y + z * z;
      break;
    case 3:
      ans = (x + 1)*(x - 1)*(y + 1)*(y - 1)*(z + 1)*(z - 1);
      break;
    case 4:
      ans = std::exp(-sigma*(x*x + y*y + z*z));
      break;
    case 5:
      r2 = x*x + y*y + z*z;
      r = std::pow(r2, 0.5);
      ans = std::exp(-sigma*r);
      break;
    default:
      BOOST_THROW_EXCEPTION(ExceptionUnknownTestId());
      break;
    }
    return ans;
  }

  num_t AnalyticalSolution::source(num_t x, num_t y, num_t z, num_t mu, num_t xi, num_t eta, num_t sigma, size_t id)
  {
    num_t ans = 0.0;
    num_t r2 = 0.0;
    num_t r = 0.0;
    switch (id)
    {
    case 0:
      ans = sigma;
      break;
    case 1:
      ans = mu + xi + eta + sigma*(x + y + z);
      break;
    case 2:
      ans = 2 * mu*x + 2 * xi*y + 2 * eta*z + sigma*(x * x + y * y + z * z);
      break;
    case 3:
      ans = 2 * mu*x*(y + 1)*(y - 1)*(z + 1)*(z - 1)
        + 2 * xi*(x + 1)*(x - 1)*y*(z + 1)*(z - 1)
        + 2 * eta*(x + 1)*(x - 1)*(y + 1)*(y - 1)*z
        + sigma*(x + 1)*(x - 1)*(y + 1)*(y - 1)*(z + 1)*(z - 1);
      break;
    case 4:
      r2 = x*x + y*y + z*z;
      ans = mu*(-2)*sigma*x*std::exp(-sigma*r2) + xi*(-2)*sigma*y*std::exp(-sigma*r2) + eta*(-2)*sigma*z*std::exp(-sigma*r2) + sigma*std::exp(-sigma*r2);
      break;
    case 5:
      r2 = x*x + y*y + z*z;
      r = std::pow(r2, 0.5);
      if (r == 0)
        ans = sigma*std::exp(-sigma*r);
      else
        ans = mu*sigma*x / r*std::exp(-sigma*r) + xi*sigma*y / r*std::exp(-sigma*r) + eta*sigma*z / r*std::exp(-sigma*r) + sigma*std::exp(-sigma*r);
      break;
    default:
      BOOST_THROW_EXCEPTION(ExceptionUnknownTestId());
      break;
    }
    return ans;
  }


}