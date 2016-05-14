#pragma once
#include "types.h"
#include "system_matrices.h"
#include "mesh.h"
namespace dgn {
	class solver {
	public:
		solver(const system_matrix_angle& s, const mesh& m) {
			//TODO
		}
	};

	class source_generator {
	public:
		source_generator(const system_matrix_angle& s, const mesh& m);
		~source_generator();

		const vector_t x() const { return x_; }
		const vector_t y() const { return y_; }
		const vector_t z() const { return z_; }

		void calculate(f3d_p f, vector_t data) const;
		

		source_generator(const source_generator&) {}
		source_generator& operator=(const source_generator&)
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

		void calculate(f3d_p f, vector_t data) const;
	private:
		vector_t x_, y_, z_;
		num_p xp_, yp_, zp_;
	};
}