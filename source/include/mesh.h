#pragma once
#include <boost/multi_array.hpp>
#include "types.h"
#include "system_matrices.h"

namespace dgn {
	class surface;
	class element;
	
	class mesh {
	public:
		mesh(size_t nx, size_t ny, size_t nz, num_t h, num_t xc, num_t yc, num_t zc);

		void generate();		

		void cache_sweep_order();
		std::vector<size_t> sweep_order(quadrant quad) const;

		boost::multi_array<size_t,2> element_interface_matrix() const;

		std::vector<size_t> refine_level() const;

		std::vector<vector3> center() const;
		
		std::vector<size_t> memory_element(size_t offset) const;
		std::vector<size_t> memory_surface(size_t offset) const;
		size_t memory_element_total(size_t offset) const;
		size_t memory_surface_total(size_t offset) const;

		size_t add_surface(size_t refine_level, surface_direction dir, num_t xc, num_t yc, num_t zc, num_t h);
		size_t add_element(size_t refine_level, num_t xc, num_t yc, num_t zc, num_t h);
		
		std::vector<bool> boundary_mark() const;

		size_t elements_total() const { return element_list_.size(); }
		size_t surfaces_total() const { return surface_list_.size(); }
		
		const element& element_ref(size_t id) const { return element_list_.at(id); }
		const surface& surface_ref(size_t id) const { return surface_list_.at(id); }

		
	private:
		void generate_coarse();

		//void mark_refine(num_t x0, num_t x1, num_t y0, num_t y1, num_t z0, num_t z1);
		//void generate_refine();

	private:
		size_t nx_, ny_, nz_;
		size_t element_basis_;
		size_t surface_basis_;
		num_t h_, xc_, yc_, zc_;
		std::vector<element> element_list_;
		std::vector<surface> surface_list_;		
		std::vector<std::vector<size_t>> sweep_order_cache_;
	};

	std::ostream& operator<<(std::ostream& os, const mesh& m);

	//There are two types of surfaces.
	class surface {
	public:


	public:
		surface(size_t id, size_t refine_level, surface_direction direction, num_t xc, num_t yc, num_t zc, num_t h);

		surface* handle() { return this; }
		const surface* handle() const { return this; }
		num_t xc() const { return xc_; }
		num_t yc() const { return yc_; }
		num_t zc() const { return zc_; }
		num_t h() const { return h_; }
		size_t id() const { return id_; }

		

		const element* element_pre(size_t local_id) const { 
			return element_pre_list_.at(local_id); 
		}
		const element* element_inc(size_t local_id) const { 
			return element_inc_list_.at(local_id); 
		}
		const size_t memory_position() const { 
			return memory_position_;
		}
		void set_memory_position(size_t pos) {
			memory_position_ = pos;
		}
		
		void split(std::vector<surface*> child_surfaces) {
			//TODO: 
		}
		const size_t refine_level() const {
			return refine_level_;
		}
		surface_direction direction() const {
			return direction_;
		}
		void link_element(element_direction ed, const element* e, size_t local_id = 0);

		bool is_boundary() const { return is_boundary_; }
		void mark_as_boundary(bool flag = true) { is_boundary_ = flag; }

		size_t pre_element_total() const { return element_pre_list_.size(); }
		size_t inc_element_total() const { return element_inc_list_.size(); }

	private:
		size_t id_;
		size_t refine_level_;
		size_t memory_position_;		
		std::vector<const element*> element_pre_list_;
		std::vector<const element*> element_inc_list_;
		num_t xc_, yc_, zc_, h_;
		surface_direction direction_;
		bool is_boundary_;
	};


	class element {
	public:
		element(size_t eId, size_t refine_level, num_t xc, num_t yc, num_t zc, num_t h);
	public:

		element* handle() { return this; }
		const element* handle() const { return this; }

		num_t xc() const { return xc_; }
		num_t yc() const { return yc_; }
		num_t zc() const { return zc_; }
		num_t h() const { return h_; }
		size_t id() const { return id_; }



		void split(std::vector<element*> child_elements){
			//TODO;
		}
		
		void link_surface(interface_direction direction, const surface* surface_ptr, size_t local_id = 0);

		//Get surface of one direction
		const surface* surface_adjacent(interface_direction dir) const;

		const size_t memory_position() const { 
			return memory_position_; 
		}
		void set_memory_position(size_t pos) { 
			memory_position_ = pos; 
		}

		bool is_to_refine() const {
			return refine_mark_;
		}
		void mark_to_refine(bool flag = true) {	
			refine_mark_ = flag; 
		}
		const size_t refine_level() const { 
			return refine_level_; 
		}

		bool is_active() const { return true; }
		
	public:
		size_t id_;
		size_t refine_level_;
		num_t xc_, yc_, zc_, h_;
		size_t memory_position_;

		std::vector<const surface*> surface_list_;
		std::vector<const element*> child_elements;
		std::vector<size_t> surface_local_id;
		bool refine_mark_;
		
	};
}