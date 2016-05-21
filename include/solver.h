#pragma once
#include "types.h"
#include "system_matrices.h"
#include "mesh.h"
#include "globals.h"

namespace dgn {	
	class source_generator {
	public:
		source_generator(const system_matrix_angle& s, const mesh& m);
		~source_generator();

		const vector_t x() const { return x_; }
		const vector_t y() const { return y_; }
		const vector_t z() const { return z_; }

		void calculate(vector_t& data) const;
	private:
		const system_matrix_angle* s_;
		const mesh* m_;
		void calculate_coordinates();
		vector_t x_, y_, z_;
		num_p xp_, yp_, zp_;
	};

	class boundary_generator {
	public:
		boundary_generator(const system_matrix_angle& s, const mesh& m);
		~boundary_generator();

		const vector_t x() const { return x_; }
		const vector_t y() const { return y_; }
		const vector_t z() const { return z_; }

		void calculate(vector_t data) const;
	private:
		vector_t x_, y_, z_;
		num_p xp_, yp_, zp_;
		const system_matrix_angle* s_;
		const mesh* m_;
		void calculate_coordinates();
	};

	class solver {
	public:
		static const size_t total_node_element(const mesh& m, size_t np);
		static const size_t total_node_surface(const mesh& m, size_t np);

		const size_t total_node_element() const;
		const size_t total_node_surface() const;
	public:
		solver(const system_matrix_angle& s, const mesh& m);		
		~solver();
		//void source_set(vector_t s);
		//void boudnary_set(vector_t b);
		//vector_t source_get() const;
		//vector_t boudnary_get() const;

		void solve();

		void solve_test(size_t id);
		
		void solve_analytical(size_t id);

		num_p boundary_ptr() const { return boundary_p_; }
		num_p source_ptr() const { return source_p_; }
		num_p solution_ptr() const { return solution_p_; }

		vector_t x() const { return sg_.x(); }
		vector_t y() const { return sg_.y(); }
		vector_t z() const { return sg_.z(); }
		vector_t solution() const { return solution_; }
		vector_t x_b() const { return bg_.x(); }
		vector_t y_b() const { return bg_.y(); }
		vector_t z_b() const { return bg_.z(); }
		void solution_mean(vector_t data) const;

		size_t memory_total_solution() const {
			return m_->memory_element_total(s_->basis_element());
		}
	private:
		source_generator sg_;
		boundary_generator bg_;
		num_p source_p_;
		num_p boundary_p_;
		num_p solution_p_;
		num_p boundary_mask_p_;

		vector_t source_;
		vector_t boundary_;
		vector_t solution_;
		vector_t boundary_mask_;
		const system_matrix_angle* s_;
		const mesh* m_;

	};

}