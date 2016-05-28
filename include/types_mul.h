#pragma once
#include "types.h"
#include "Mesh.h"
#include "solver.h"
#include "system_matrices.h"

namespace linear_boltzmann_transport_equation_solver {

	typedef num_t num_t;
	typedef num_p num_p;
	typedef const_num_p const_num_p;
	typedef matrix_t Matrix;
	typedef vector_t Vector;

	typedef Mesh Mesh;
	typedef solver SolverSingleAngle;
	typedef system_matrix_angle SystemMatrixSingleAngle;
	
}