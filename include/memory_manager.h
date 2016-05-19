#pragma once

namespace linear_boltzmann_transport_equation_solver {
	class memory_manager {
		memory_manager();
		memory_manager(const size_t nx, const size_t ny, const size_t nz, const size_t na, const size_t ne);
		memory_manager(const size_t n_element, const size_t n_surfaces, const size_t na, const size_t ne);
		void solver_size_set(const size_t nx, const size_t ny, const size_t nz, const size_t na, const size_t ne);
		void solver_size_set(const size_t n_element, const size_t n_surfaces, const size_t na, const size_t ne);
		void alloc_memory();
		num_t* 
	};
}