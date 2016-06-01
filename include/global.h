#pragma once
#define  MEATIME
#include "xlib.h"
#include "types.h"
#include "utilities.h"
#include "exception_types.h"
#include "problem_definition.h"

namespace discontinues_galerkin_nodal_solver {
  using xlib::kInfSize;
  using xlib::operator<<;
  using xlib::mkl::operator<<;  
  using xlib::kInfSize;
  using xlib::mkl::Transport;

  class AnalyticalSolution {
  public:
    static num_t solution(num_t x, num_t y, num_t z, 
        num_t mu, num_t xi, num_t eta, num_t sigma, size_t id);
    static num_t source(num_t x, num_t y, num_t z,
      num_t mu, num_t xi, num_t eta, num_t sigma, size_t id);
  };


}


