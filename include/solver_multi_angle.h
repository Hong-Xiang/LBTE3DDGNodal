#pragma once
#include <string>
#include <vector>
#include "xlib.h"

#include "types_mul.h"

#include "memory_manager.h"

namespace linear_boltzmann_transport_equation_solver {


	//TODO: add more methods to input parameters, i.e. from console or file.
	class ProblemDefinition {
	public:
		ProblemDefinition() = default;
		ProblemDefinition(size_t nx, size_t ny, size_t nz, size_t np, size_t na, num_t sigma, std::string file_scattering_matrix, std::string file_angle_info);

		const size_t nx() const { return nx_; }
		const size_t ny() const { return ny_; }
		const size_t nz() const { return nz_; }
		const size_t np() const { return np_; }
		const size_t na() const { return na_; }
		const std::string file_name_scattering_matrix() const { return file_scattering_matrix_; }
		const std::string file_name_angle_info() const { return file_angle_info_; }

		const num_t h() const { return h_; }
		const num_t xc() const { return xc_; }
		const num_t yc() const { return yc_; }
		const num_t zc() const { return zc_; }
		const num_t sigma() const { return sigma_; }

		const size_t verbose() const { return verbose_; }
		
		const size_t iteration_number() const { return iteration_number_; }
	private:
		size_t nx_, ny_, nz_; //space cell of initial Mesh
		size_t np_; // basis per dimension
		size_t na_; // total number of angles
		num_t sigma_; // total cross section
		num_t h_; //cell size of initial Mesh
		num_t xc_, yc_, zc_; //center of spatial Mesh
		std::string file_scattering_matrix_;
		std::string file_angle_info_;
		size_t verbose_;
		size_t iteration_number_;
	};

/*!
 * \class SolverMultiAngle
 *
 * \ingroup linear_boltzmann_transport_equation_solver
 *
 * \brief Users' interface of LBTE solver
 *
 * User should provide problem definition only. Without any details about Mesh, angle or other things.
 * All those information are read-only by SolverMultiAngle.
 * Currently SolverMultiAngle only provide solution of mean value of voxels.
 * 
 *
 * \note 
 *
 * \author HongXiang
 *
 * \version 1.0
 *
 * \date may 2016
 *
 * Contact: hx.hongxiang@gmail.com
 *
 */	
	class SolverMultiAngle {
	public:
		SolverMultiAngle(const ProblemDefinition pd, const Mesh& Mesh);

		virtual ~SolverMultiAngle();

	public:
		std::vector<vector3> Center() const;
		
		//Vector solution_mean_voxel() const;

		//Vector solution_fix_angle(size_t ia) const;

		void SaveSolution(std::string file_name, std::string variable_name);

		void Solve();

		//Analytical tests:
		
		void SolveAnalyticalByDirect();

		void SolveAnalyticalBySolve();		
		void CalculateSourceAnalytical();
		void CalculateBoundaryAnalytical();

	private:
		void Initialization();

		void Iteration();
		void CalculateScattering();
		void SolveSingleAngle(size_t ia);
		
	private:
		void load_scattering_matrix();
		void load_angles();
		void ConstructSystemMatrices();
		void ConstructBasicDiscontinuesGalerkinNodalSolvers();
		void ConstructCoordinates();
		void BindVectors();
		void clear();
		

		
	private:
		SolverMultiAngle(SolverMultiAngle const&) = delete;
		SolverMultiAngle& operator=(SolverMultiAngle const&) = delete;
	//helper functions:
	public:
		const Matrix scattering_matrix() const { return scattering_matrix_; }
		const Vector mu() const { return mu_; }
		const Vector xi() const { return xi_; }
		const Vector eta() const { return eta_; }
	private:
		

	private:
		ProblemDefinition problem_definition_;
		const Mesh& mesh_;
		std::vector<SystemMatrixSingleAngle*> system_matrices_;
		std::vector<SolverSingleAngle*> solvers_angle_;
		Matrix scattering_matrix_;
		Vector mu_, xi_, eta_;
		MemoryManager memory_manager_;				
		bool isInitialized_;

	//phase 1 cache variables: usable after Initialization()
	private:
		std::vector<Vector> solutions_single_; //solution handle in solver of different angle
		std::vector<Vector> solutions_multiple_; //solution handle in multi angle solver
		std::vector<Vector> sources_single_;
		std::vector<Vector> sources_multiple_;
		std::vector<Vector> sources_ext_multiple_;
		std::vector<Vector> boundary_single_;
		std::vector<Vector> boundary_multiple_;
		Matrix solution_m_;
		Matrix source_m_;
		Vector solution_v_;
		Vector source_v_;
		Vector source_ext_v_;
		Vector boundary_v_;

	//phase 2 cache variables: usable after Solve()
	private:
		Vector x_solution_, y_solution_, z_solution_;
		Vector mu_solution_, xi_solution_, eta_solution_;

		Vector x_boundary_, y_boundary_, z_boundary_;
		Vector mu_boundary_, xi_boundary_, eta_boundary_;

	private:
		size_t node_per_angle_;
		size_t node_total;
	};


}