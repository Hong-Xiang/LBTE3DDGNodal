#include "memory_manager.h"
#include <xlib.h>

//linear_boltzmann_transport_equation_solver::MemoryManager::MemoryManager()
//	: solution_p_(nullptr), source_p_(nullptr), boundary_p_(nullptr), na_(0)
//{
//	alloc_memory();
//}

linear_boltzmann_transport_equation_solver::MemoryManager::~MemoryManager()
{
	clear();
}

void linear_boltzmann_transport_equation_solver::MemoryManager::alloc()
{
	solution_p_ = xlib::mkl_ext::xcalloc<num_t>(n_element_node_*na_);
	source_p_ = xlib::mkl_ext::xcalloc<num_t>(n_element_node_*na_);
	source_ext_p_ = xlib::mkl_ext::xcalloc<num_t>(n_element_node_*na_);
	boundary_p_ = xlib::mkl_ext::xcalloc<num_t>(n_surface_node_*na_);
}


#include "globals_mul.h"

void linear_boltzmann_transport_equation_solver::MemoryManager::clear()
{
	if (Global::verbose() >= 3)
	{
		std::cout << "[LBTE::MemoryManager]\tClear of memory manager." << std::endl;
	}
	
	if (solution_p_ != nullptr)
		mkl_free(solution_p_);
	if (source_p_ != nullptr)
		mkl_free(source_p_);
	if (boundary_p_ != nullptr)
		mkl_free(boundary_p_);
	if (source_ext_p_ != nullptr)
		mkl_free(source_ext_p_);
}

void linear_boltzmann_transport_equation_solver::MemoryManager::solution_pos_inc(const size_t ia, num_p& ptr, size_t& inc) const
{
	ptr = solution_p_ + ia;
	inc = na_;
}


void linear_boltzmann_transport_equation_solver::MemoryManager::source_pos_inc(const size_t ia, num_p& ptr, size_t& inc) const
{
	ptr = source_p_ + ia;
	inc = na_;
}



void linear_boltzmann_transport_equation_solver::MemoryManager::source_ext_pos_inc(const size_t ia, num_p& ptr, size_t& inc) const
{
	ptr = source_ext_p_ + ia;
	inc = na_;
}

void linear_boltzmann_transport_equation_solver::MemoryManager::boundary_pos_inc(const size_t ia, num_p& ptr, size_t& inc) const
{
	ptr = boundary_p_ + ia;
	inc = na_;
}

linear_boltzmann_transport_equation_solver::MemoryManager::MemoryManager(const size_t total_element_node, const size_t total_surface_node, const size_t na)
	:	solution_p_(nullptr), source_p_(nullptr), boundary_p_(nullptr),
		na_(na), n_element_node_(total_element_node), n_surface_node_(total_surface_node)
{
	alloc();
}



