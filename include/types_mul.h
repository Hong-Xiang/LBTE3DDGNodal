#pragma once
#include "types.h"
#include "mesh.h"
#include "solver.h"
#include "system_matrices.h"

namespace linear_boltzmann_transport_equation_solver {

	typedef dgn::num_t num_t;
	typedef dgn::num_p num_p;
	typedef dgn::const_num_p const_num_p;
	typedef dgn::matrix_t Matrix;
	typedef dgn::vector_t Vector;

	typedef dgn::mesh Mesh;
	typedef dgn::solver SolverSingleAngle;
	typedef dgn::system_matrix_angle SystemMatrixSingleAngle;
	
}