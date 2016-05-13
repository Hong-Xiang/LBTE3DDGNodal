#pragma once
#include "Global.h"
#include <boost\multi_array.hpp>
#include <string>
#include <memory>

namespace DG3DNodal {



/*!
 * \class sovler_DG3DNodal
 *
 * \ingroup solver
 *
 * \brief sovler_DG3DNodal : Core solver for problem:
 *		\Omega \cdot \nabla \psi + \sigma \psi = q
 *			problem parameters:
 *				<\mu> <\xi> <\eta> : \Omega = (\mu, \xi, \eta)
 *				\sigma
 *				xc, yc, zc, h_cell, nx, ny, nz
 *
 * TODO: use method:
 *		Initialization:
 *			A.	setup problem by:
 *					set_problem(mu, xi, eta, sigma, xc, yc, zc, h, nx, ny, nz)
 *				then calculate all data need for solver
 *					assemble_system()
 *			B.	load solver 
 *					load_solver(filename)
 *		Calculate or set source and boundary
 *		
 *
 *			
 *
 *
 *
 * \note 
 *
 * \author HongXiang
 *
 * \version 1.0
 *
 * \date May 2016
 *
 *
 */
	class mesh;
	class system_matrices;

	class sovler_DG3DNodal {
	public:
		enum statues {
			SUCCEED = 0,
			FAILED = 1,
		};

		enum phase {
			PHASE_PRE_MESH = 0,
			PHASE_PRE_CONDITION = 1,
			PHASE_PRE_SOLVE = 2
		};

		enum boundary_type {
			DIRECLET_VACUUM = 0,
			DIRECLET_NORMAL = 1
		};

	public:		
		statues mesh_set(num_t xc, num_t yc, num_t zc, num_t h, size_t nx, size_t ny, size_t nz);
		statues problem_set(num_t mu, num_t xi, num_t eta, num_t sigma);
		
		//Generate mesh, calculate and save all mesh data
		virtual statues assemble_system();
		statues mesh_calculate();
		statues mesh_load(std::string filename);
		statues mesh_save(std::string filename) const;

		statues memory_bind(num_p solution_p, num_p source_p, num_p interface_p);
						
		statues solve() const;
		
		//Helper functions
		//available on PHASE_PRE_CONDITION
		const_num_p x_e() const;
		const_num_p y_e() const;
		const_num_p z_e() const;

		const_num_p x_i() const;
		const_num_p y_i() const;
		const_num_p z_i() const;


		size_t elements_total() const;
		size_t interfaces_total() const;
		size_t basis_total() const;

		size_t source_size() const;
		size_t solution_size() const;
		size_t boundary_size() const;

		//available on PHASE_PRE_SOLVE
		num_p source() const;
		num_p solution() const;
		num_p interfaces() const;

	//bind pointers
	private:
		num_p solution_;
		num_p source_;
		num_p interfaces_;
	//flags:
	private:
		bool _is_system_assembled;
		bool _is_memory_bind;
		bool _is_boundary_source_calculated;
		
	private:
		num_p xc;
		num_p yc;
		num_p zc;

		

	private:
		boost::multi_array<size_t, 2> _element_interface_matrix;		
		boost::multi_array<size_t, 2> _sweep_order;
		std::vector<bool> _is_boundary;
	public:
		virtual ~sovler_DG3DNodal();

	private:
		std::shared_ptr<mesh> _mesh;
		std::vector<std::shared_ptr<system_matrices>> _system_matrices_list;
	};
}