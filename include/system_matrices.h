#pragma once
#include "types.h"
namespace dgn {
	class invalid_dimension : public std::exception {};
	/*!
	 * \class matrices_data
	 *
	 * \brief Provide matrices related data for 1D standard element. Support up to 20 order.
	 *		  Singleton pattern, use instance() method to handle it.
	 *
	 * \author HongXiang
	 * \date May 2016
	 */
	class matrices_data {
	public:
		static const matrices_data& instance();
		const vector_t x(size_t np) const;
		const matrix_t m(size_t np) const;
		const matrix_t s(size_t np) const;		
	private:
		matrices_data();
		void read_MAT_files();
		static matrices_data* handle_;
		std::vector<num_p> x_;
		std::vector<num_p> s_;
		std::vector<num_p> m_;

	private:
		const size_t MINNP = 2;
		const size_t MAXNP = 20;
		
	};

	class matrices_base {
	public:

		matrices_base(num_t h, size_t np);
		void reset(num_t h, size_t np);
		virtual const size_t basis_total() const = 0;
		virtual const size_t basis_lower_dim_total() const = 0;

		const matrix_t lift_matrix(interface_direction idr) const { return lift_matrices_.at(static_cast<size_t>(idr)); }
		const matrix_t flux_matrix(interface_direction idr) const { return flux_matrices_.at(static_cast<size_t>(idr)); }

		void initialize();

		virtual size_t dimension() const { 
			//throw(invalid_dimension());
			return 0; 
		}
		virtual ~matrices_base();
	protected:
		virtual void calculate() = 0;
		virtual void memory_alloc();
	
	protected:
		num_t h_;
		size_t np_;

		std::vector<num_p> lift_matrices_ptr_;
		std::vector<num_p> flux_matrices_ptr_;
		std::vector<matrix_t> lift_matrices_;
		std::vector<matrix_t> flux_matrices_;

		size_t dim_;
		bool is_init_;
	};

	class matrices1d : public matrices_base {
	public:
		static size_t basis_total(size_t np) {
			return np;
		}
		static size_t basis_lower_dim_total(size_t np) {
			return 1;
		}
	public:
		
		matrices1d(num_t h, size_t np);
		const vector_t x() const { return x_; }
		const matrix_t m() const { return m_; }
		const matrix_t s() const { return s_; }
		virtual const size_t basis_total() const { return basis_total(np_); }
		virtual const size_t basis_lower_dim_total() const { return basis_lower_dim_total(np_); }
		virtual size_t dimension() const { return 1; }
		~matrices1d();
	protected:
		virtual void calculate();
		virtual void memory_alloc();

	private:
		num_p xp_;
		num_p mp_;
		num_p sp_;		


		vector_t x_;
		matrix_t m_;
		matrix_t s_;
	};


	class matrices2d : public matrices_base {
	public:
		static size_t basis_total(size_t np) {
			return np*np;
		}
		static size_t basis_lower_dim_total(size_t np) {
			return np;
		}
	public:
		matrices2d(num_t h, size_t np);
		const size_t basis_total() const { return np_*np_; }
		virtual const size_t basis_lower_dim_total() const { return np_; }
		virtual size_t dimension() const { return 2; }
		~matrices2d();
	protected:
		virtual void calculate();
		virtual void memory_alloc();
	public:
		const vector_t x() const { return x_; }
		const vector_t y() const { return y_; }
		const matrix_t m() const { return m_; }
		const matrix_t sx() const { return sx_; }
		const matrix_t sy() const { return sy_; } 		
	private:
		num_p xp_;
		num_p yp_;
		num_p mp_;
		num_p sxp_;
		num_p syp_;

		vector_t x_;
		vector_t y_;
		matrix_t m_;
		matrix_t sx_;
		matrix_t sy_;
	};

	class matrices3d : public matrices_base {
	public:
		static size_t basis_total(size_t np) {
			return np*np*np;
		}
		static size_t basis_lower_dim_total(size_t np) {
			return np*np;
		}
	public:
		matrices3d(num_t h, size_t np);
		virtual const size_t basis_total() const { return np_*np_*np_; }
		virtual const size_t basis_lower_dim_total() const { return np_*np_; }
		virtual size_t dimension() const { return 3; }
		~matrices3d();
	protected:
		virtual void calculate();
		virtual void memory_alloc();
	public:
		const vector_t x() const { return x_; }
		const vector_t y() const { return y_; }
		const vector_t z() const { return z_; }
		const matrix_t m() const { return m_; }
		const matrix_t sx() const { return sx_; }
		const matrix_t sy() const { return sy_; }
		const matrix_t sz() const { return sz_; }
	private:
		num_p xp_;
		num_p yp_;
		num_p zp_;
		num_p mp_;
		num_p sxp_;
		num_p syp_;
		num_p szp_;		


		vector_t x_;
		vector_t y_;
		vector_t z_;
		matrix_t m_;
		matrix_t sx_;
		matrix_t sy_;
		matrix_t sz_;
	};

	class system_matrix_angle {			
	public:
		static size_t basis_total(size_t np) {
			return matrices3d::basis_total(np);
		}
		static size_t basis_lower_dim_total(size_t np) {
			return matrices3d::basis_lower_dim_total(np);
		}
	public:
		system_matrix_angle(num_t h, size_t np, num_t mu, num_t xi, num_t eta, num_t sigma);
		const matrix_t system_matrix() const { return sysm_; }
		const xlib::mkl_ext::matrix_lu<num_t> system_matrixlu() const { return sysmlu_; }
		const matrix_t system_matrixiv() const { return sysmiv_; }
		const matrix_t mass_matrix() const { return massm_; }
		const matrix_t lift_matrix(interface_direction dir) const {
			return lift_matrix_.at(static_cast<size_t>(dir));
		}
		const matrix_t reorder() const {
			return reorder_;
		}
		
		const matrix_t back_order() const {
			return back_order_;
		}
		const matrix_t flux_matrix(interface_direction dir) const {
			return flux_matrix_.at(static_cast<size_t>(dir));
		}
		const vector_t x() const { return matrices_.x(); }
		const vector_t y() const { return matrices_.y(); }
		const vector_t z() const { return matrices_.z(); }

		const vector_t xs() const { return surface_matrices_.x(); }
		const vector_t ys() const { return surface_matrices_.y(); }

		const size_t np() const { return np_; }
		num_t mu() const { return mu_; }
		num_t xi() const { return xi_; }
		num_t eta() const { return eta_; }
		num_t sigma() const { return sigma_; }

		size_t basis_element() const { return matrices_.basis_total(); }
		size_t basis_surface() const { return matrices_.basis_lower_dim_total(); }
		~system_matrix_angle();
	private:
		matrices3d matrices_;
		matrices2d surface_matrices_;
		num_t h_, mu_, xi_, eta_, sigma_;
		size_t np_;
		num_p sysmp_;
		matrix_t sysm_;
		
		num_p sysmlup_;
		xlib::mkl_ext::matrix_lu<num_t> sysmlu_;
		lapack_int* ipiv;

		num_p sysmivp_;
		matrix_t sysmiv_;
		num_p massmp_;
		matrix_t massm_;
		std::vector<num_p> lift_matrix_p_;
		std::vector<num_p> flux_matrix_p_;
		std::vector<matrix_t> lift_matrix_;
		std::vector<matrix_t> flux_matrix_;		

		num_p reorder_p_;
		matrix_t reorder_;
		num_p back_order_p_;
		matrix_t back_order_;
	};
	std::ostream& operator<<(std::ostream& os , const system_matrix_angle& s);
}
