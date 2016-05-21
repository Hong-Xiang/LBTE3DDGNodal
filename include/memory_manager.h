#pragma once
#include "types_mul.h"

namespace linear_boltzmann_transport_equation_solver {
	/*!
	* \class MemoryManager
	*
	* \ingroup Groupname
	*
	* \brief memory allocation manager, provides methods to access those pointers
	*
	* solution_pos_inc()
	* source_pos_inc()
	* source_ext_pos_inc()
	* boundary_pos_inc()
	*
	* \note
	*
	* \author HongXiang
	*
	* \version 1.0
	*
	* \date May 2016
	*
	* Contact: hx.hongxiang@gmail.com
	*
	*/
	class MemoryManager {
		//MemoryManager();
	public:
		MemoryManager(const size_t total_element_node, const size_t total_surface_node, const size_t na);
		virtual ~MemoryManager();

		//void solver_size_set(const size_t total_node, const size_t na);
		
		

		void solution_pos_inc(const size_t ia, num_p& ptr, size_t& inc) const;
		void source_pos_inc(const size_t ia, num_p& ptr, size_t& inc) const;
		void source_ext_pos_inc(const size_t ia, num_p& ptr, size_t& inc) const;
		void boundary_pos_inc(const size_t ia, num_p& ptr, size_t& inc) const;

		num_p solution_ptr() const { return solution_p_; }
		num_p source_ptr() const { return source_p_; }
		num_p source_ext_ptr() const { return source_ext_p_; }
		num_p boundary_ptr() const { return boundary_p_; }
		
		const size_t na() const { return na_; }
		const size_t total_element_basis() const { return n_element_node_; }
		const size_t total_surface_basis() const { return n_surface_node_; }

		MemoryManager(MemoryManager const&) = delete;
		MemoryManager& operator=(MemoryManager const&) = delete;
	private:
		
		void alloc();
		void clear();
	private:
		num_p solution_p_;
		num_p source_p_;
		num_p source_ext_p_;
		num_p boundary_p_;
		
		size_t na_;
		size_t n_element_node_;
		size_t n_surface_node_;
	};
}