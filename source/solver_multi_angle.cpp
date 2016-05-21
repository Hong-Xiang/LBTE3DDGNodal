#include "solver_multi_angle.h"
#include "globals_mul.h"
#include "system_matrices.h"
#include "solver.h"
#include <ctime>
namespace linear_boltzmann_transport_equation_solver {

	ProblemDefinition::ProblemDefinition(size_t nx, size_t ny, size_t nz, size_t np, size_t na, num_t sigma, std::string file_scattering_matrix, std::string file_angle_info)
		: nx_(nx), ny_(ny), nz_(nz), np_(np), na_(na), sigma_(sigma),
		file_scattering_matrix_(file_scattering_matrix),
		file_angle_info_(file_angle_info)
	{
		h_ = 2.0 / nx_;
		xc_ = 0.0, yc_ = 0.0, zc_ = 0.0;
		verbose_ = 2;
		iteration_number_ = 10;
	}

	void SolverMultiAngle::Initialization()
	{		
		load_scattering_matrix();
		std::cout << "load scattering matrix finished" << std::endl;
		load_angles();
		std::cout << "load angle finished" << std::endl;
		ConstructSystemMatrices();
		std::cout << "construct system matrix finished" << std::endl;
		ConstructBasicDiscontinuesGalerkinNodalSolvers();
		std::cout << "construct solver finished" << std::endl;
		BindVectors();
		std::cout << "bind vector finished" << std::endl;
		ConstructCoordinates();		
		std::cout << "construct coordinate finished" << std::endl;
	}

	SolverMultiAngle::~SolverMultiAngle()
	{
		clear();
	}

	std::vector<dgn::vector3> SolverMultiAngle::Center() const
	{
		return mesh_.center();
	}

	void SolverMultiAngle::SaveSolution(std::string file_name, std::string variable_name)
	{
		std::vector<size_t> sz;
		sz.resize(2);
		sz.at(0) = problem_definition_.na();
		sz.at(1) = mesh_.memory_element_total(SystemMatrixSingleAngle::basis_total(problem_definition_.np()));
		xlib::imatlab::save_MAT_data(solution_v_.ptr(), sz, file_name.c_str(), variable_name.c_str());
	}

	void SolverMultiAngle::Solve()
	{
		for (size_t i = 0; i < problem_definition_.iteration_number(); i++)
		{
			time_t t0 = std::clock();
			Iteration();
			time_t t1 = std::clock();
			std::cout << "iteration " << i << " end, cost " << (double)(t1 - t0) / CLOCKS_PER_SEC << std::endl;
		}
	}

	void SolverMultiAngle::SolveAnalyticalByDirect()
	{
		for (size_t ia = 0; ia < problem_definition_.na(); ia++)
		{
			num_t mu = mu_(ia), xi = xi_(ia), eta = eta_(ia);
			for (size_t inode = 0; inode < solutions_multiple_.at(ia).size(); inode++)
			{
				solutions_multiple_.at(ia)(inode) = AnalyticalSolution::solution(
					x_solution_(inode),
					y_solution_(inode),
					z_solution_(inode),
					mu,
					xi,
					eta,
					problem_definition_.sigma()
				);
			}
		}
	}

	void SolverMultiAngle::SolveAnalyticalBySolve()
	{
		CalculateSourceAnalytical();
		CalculateBoundaryAnalytical();
		Solve();
	}

	void SolverMultiAngle::CalculateSourceAnalytical()
	{
		for (size_t ia = 0; ia < problem_definition_.na(); ia++)
		{
			num_t mu = mu_(ia), xi = xi_(ia), eta = eta_(ia);
			for (size_t inode = 0; inode < sources_ext_multiple_.at(ia).size(); inode++)
			{
				sources_ext_multiple_.at(ia)(inode) = AnalyticalSolution::source(
					x_solution_(inode),
					y_solution_(inode),
					z_solution_(inode),
					mu,
					xi,
					eta,
					problem_definition_.sigma()
				);
			}
		}
	}

	void SolverMultiAngle::CalculateBoundaryAnalytical()
	{
		for (size_t ia = 0; ia < problem_definition_.na(); ia++)
		{
			num_t mu = mu_(ia), xi = xi_(ia), eta = eta_(ia);
			for (size_t inode = 0; inode < boundary_multiple_.at(ia).size(); inode++)
			{
				boundary_multiple_.at(ia)(inode) = AnalyticalSolution::solution(
					x_boundary_(inode),
					y_boundary_(inode),
					z_boundary_(inode),
					mu,
					xi,
					eta,
					problem_definition_.sigma()
				);
			}
		}
	}

	void SolverMultiAngle::Iteration()
	{
		
		CalculateScattering();
		for (size_t i = 0; i < problem_definition_.na(); i++)
		{			
			sources_single_.at(i).copy(sources_multiple_.at(i));			
			sources_single_.at(i).add(sources_ext_multiple_.at(i));						
			boundary_single_.at(i).copy(boundary_multiple_.at(i));
			SolveSingleAngle(i);
			solutions_multiple_.at(i).copy(solutions_single_.at(i));
//#ifdef _DEBUG
//			std::cout << "ia = " << i << std::endl;
//			std::cout << "source = " << std::endl << sources_single_.at(i);
//			std::cout << "boundary = " << std::endl << boundary_single_.at(i);
//			std::cout << "solution = " << std::endl << solutions_single_.at(i);
//			std::cout << "solution multiple = " << std::endl << solutions_multiple_.at(i);
//			std::cout << "end of ia = " << i << std::endl;
//#endif
		}

	}

	void SolverMultiAngle::CalculateScattering()
	{
		
		source_m_.gemm(scattering_matrix_, solution_m_);

	}

	void SolverMultiAngle::SolveSingleAngle(size_t ia)
	{
		solvers_angle_.at(ia)->solve();
	}

	void SolverMultiAngle::load_scattering_matrix()
	{
		scattering_matrix_.alloc(problem_definition_.na(), problem_definition_.na());
		xlib::imatlab::load_MAT_data(scattering_matrix_.ptr(), problem_definition_.file_name_scattering_matrix().c_str(), "scam");
	}

	void SolverMultiAngle::load_angles()
	{
		mu_.alloc(problem_definition_.na());
		xi_.alloc(problem_definition_.na());
		eta_.alloc(problem_definition_.na());
		xlib::imatlab::load_MAT_data(mu_.ptr(), problem_definition_.file_name_angle_info().c_str(), "mu");
		xlib::imatlab::load_MAT_data(xi_.ptr(), problem_definition_.file_name_angle_info().c_str(), "xi");
		xlib::imatlab::load_MAT_data(eta_.ptr(), problem_definition_.file_name_angle_info().c_str(), "eta");
	}

	void SolverMultiAngle::ConstructSystemMatrices()
	{
		for each (auto ptr in system_matrices_)
		{
			if (ptr != nullptr)
				delete ptr;
		}

		system_matrices_.clear();
		system_matrices_.resize(problem_definition_.na());

		std::cout << "Start to construct system matrices" << std::endl;
		for (size_t i = 0; i < problem_definition_.na(); i++)
		{
			num_t mu = mu_(i), xi = xi_(i), eta = eta_(i);
			num_t h = problem_definition_.h();
			size_t np = problem_definition_.np();
			num_t sigma = problem_definition_.sigma();
			std::cout << "ia = " << i 
				<< " mu= " << mu << "xi= " << xi << "eta= " << eta
				<< "h= " << h << "np= " << np << "sigma= " << sigma << std::endl;
			system_matrices_.at(i) = new SystemMatrixSingleAngle(
				problem_definition_.h(),
				problem_definition_.np(),
				mu,
				xi,
				eta,
				problem_definition_.sigma()
			);
		}
	}

	void SolverMultiAngle::ConstructBasicDiscontinuesGalerkinNodalSolvers()
	{
		for each (auto ptr in solvers_angle_)
		{
			if (ptr != nullptr)
				delete ptr;
		}
		solvers_angle_.clear();
		solvers_angle_.resize(problem_definition_.na());
		for (size_t i = 0; i < problem_definition_.na(); i++)
		{
			num_t mu = mu_(i), xi = xi_(i), eta = eta_(i);
			solvers_angle_.at(i) = new SolverSingleAngle(
				*system_matrices_.at(i),
				mesh_
			);
		}
	}

	void SolverMultiAngle::ConstructCoordinates()
	{
		size_t total_node_element = mesh_.memory_element_total(SystemMatrixSingleAngle::basis_total(problem_definition_.np()));
		size_t total_node_surface = mesh_.memory_surface_total(SystemMatrixSingleAngle::basis_lower_dim_total(problem_definition_.np()));		
		
		x_solution_.alloc(total_node_element);
		y_solution_.alloc(total_node_element);
		z_solution_.alloc(total_node_element);
		mu_solution_.alloc(total_node_element);
		xi_solution_.alloc(total_node_element);
		eta_solution_.alloc(total_node_element);

		x_boundary_.alloc(total_node_surface);
		y_boundary_.alloc(total_node_surface);
		z_boundary_.alloc(total_node_surface);
		mu_boundary_.alloc(total_node_surface);
		xi_boundary_.alloc(total_node_surface);
		eta_boundary_.alloc(total_node_surface);





		x_solution_.copy(solvers_angle_.at(0)->x());
		y_solution_.copy(solvers_angle_.at(0)->y());
		z_solution_.copy(solvers_angle_.at(0)->z());
		x_boundary_.copy(solvers_angle_.at(0)->x_b());
		y_boundary_.copy(solvers_angle_.at(0)->y_b());
		z_boundary_.copy(solvers_angle_.at(0)->z_b());

		
		
	}

	void SolverMultiAngle::BindVectors()
	{
		size_t total_element_node_per_angle = mesh_.memory_element_total(SystemMatrixSingleAngle::basis_total(problem_definition_.np()));
		size_t total_surface_node_per_angle = mesh_.memory_surface_total(SystemMatrixSingleAngle::basis_lower_dim_total(problem_definition_.np()));
		size_t na_ = problem_definition_.na();

		solutions_single_.clear();
		solutions_single_.resize(problem_definition_.na());
		sources_single_.clear();
		sources_single_.resize(problem_definition_.na());
		boundary_single_.clear();
		boundary_single_.resize(problem_definition_.na());
		solutions_multiple_.clear();
		solutions_multiple_.resize(problem_definition_.na());
		sources_multiple_.clear();
		sources_multiple_.resize(problem_definition_.na());
		sources_ext_multiple_.clear();
		sources_ext_multiple_.resize(problem_definition_.na());
		boundary_multiple_.clear();
		boundary_multiple_.resize(problem_definition_.na());

		for (size_t ia = 0; ia < problem_definition_.na(); ia++)
		{
			solutions_single_.at(ia).bind(
				solvers_angle_.at(ia)->total_node_element(),
				solvers_angle_.at(ia)->solution_ptr()
			);

			sources_single_.at(ia).bind(
				solvers_angle_.at(ia)->total_node_element(),
				solvers_angle_.at(ia)->source_ptr()
			);

			boundary_single_.at(ia).bind(
				solvers_angle_.at(ia)->total_node_surface(),
				solvers_angle_.at(ia)->boundary_ptr()
			);

						
			num_p ptr; size_t inc;			
			memory_manager_.solution_pos_inc(ia, ptr, inc);
			solutions_multiple_.at(ia).bind(total_element_node_per_angle,ptr,inc);
			memory_manager_.source_pos_inc(ia, ptr, inc);
			sources_multiple_.at(ia).bind(total_element_node_per_angle, ptr, inc);
			memory_manager_.source_ext_pos_inc(ia, ptr, inc);
			sources_ext_multiple_.at(ia).bind(total_element_node_per_angle, ptr, inc);
			memory_manager_.boundary_pos_inc(ia, ptr, inc);
			boundary_multiple_.at(ia).bind(total_surface_node_per_angle, ptr, inc);

			solution_m_.bind(na_, total_element_node_per_angle, memory_manager_.solution_ptr());
			source_m_.bind(na_, total_element_node_per_angle, memory_manager_.source_ptr());

			solution_v_.bind(na_*total_element_node_per_angle, memory_manager_.solution_ptr());
			source_v_.bind(na_*total_element_node_per_angle, memory_manager_.source_ptr());
			source_ext_v_.bind(na_*total_element_node_per_angle, memory_manager_.source_ext_ptr());
			boundary_v_.bind(na_*total_surface_node_per_angle, memory_manager_.boundary_ptr());
		}
	}

	void SolverMultiAngle::clear()
	{
		scattering_matrix_.clear();
		mu_.clear();
		xi_.clear();
		eta_.clear();
		for each (auto ptr in system_matrices_)
		{
			if(ptr != nullptr)
				delete ptr;
		}
		for each (auto ptr in solvers_angle_)
		{
			if (ptr != nullptr)
				delete ptr;
		}
	}

	SolverMultiAngle::SolverMultiAngle(const ProblemDefinition pd, const Mesh& mesh)
		:	problem_definition_(pd), 
			mesh_(mesh),
			memory_manager_(
				mesh_.memory_element_total(SystemMatrixSingleAngle::basis_total(pd.np())),
				mesh_.memory_surface_total(SystemMatrixSingleAngle::basis_lower_dim_total(pd.np())),
				pd.na()
			)									
	{
		Initialization();
	}


}

