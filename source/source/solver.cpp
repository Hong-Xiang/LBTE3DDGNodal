#include "solver.h"

dgn::source_generator::source_generator(const system_matrix_angle& s, const mesh& m)
{
	size_t szb = s.basis_element();
	size_t n = m.memory_element_total(szb);
	xp_ = xlib::mkl_ext::xcalloc<num_t>(n);
	yp_ = xlib::mkl_ext::xcalloc<num_t>(n);
	zp_ = xlib::mkl_ext::xcalloc<num_t>(n);

}

dgn::source_generator::~source_generator()
{

}
