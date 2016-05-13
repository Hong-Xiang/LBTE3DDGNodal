#pragma once
#include "Global.h"
#include "Utils.h"
#include <boost/multi_array.hpp>
#include "Elements.h"
#include "Interfaces.h"

namespace DG3DNodal {


	class Mesh {
	public:
		Mesh(size_t nx, size_t ny, size_t nz, num_t h, num_t xc, num_t yc, num_t zc);
		void generate_mesh();		
		//size_t calculate_total_element_memory();
		//size_t calculate_total_interface_memory();
		size_matrix_t get_element_interface_matrix() const;
		size_array_t get_sweep_order(Quadrant quad) const;
		size_array_t get_refine_level() const;
		std::vector<vector3> elements_center() const;
		//num_p get_memory_position(num_p data, size_t offset) const;
		size_t add_interface(mesh_refine_lvl_t refine_level);

		void print_index_of_elements(std::ostream& os);
		std::vector<bool> get_boundary_mark() const;

		size_t total_element() const { return el.size(); }
		size_t total_interface() const { return il.size(); }

	private:
		size_t nx_, ny_, nz_, h_, xc_, yc_, zc_;
		std::vector<Element3D> el;
		std::vector<Interface> il;
		boost::multi_array<size_t, 3> index3;
		std::vector<size_t> m_ix, m_iy, m_iz;
	};


	class mesh {
	public:
		mesh();
	};
}