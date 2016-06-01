//////#include <mkl.h>
//////#include <iostream>
//////#include <fstream>
//////#include <random>
//////#include "system_matrices.h"
//////#include "globals.h"
//////#include "types.h"
//////#include "Utilities.h"
//////#include "xlib_mkl_ext.hpp"
//////#include "Mesh.h"
//////#include "solver.h"
//////#include "system_matrices.h"
//////#include "Mesh.h"
//////
//////#include <mat.h>
//////#include <mex.h>
//////#include <engine.h>
////
////
//////////
////////////int main_() {	
////////////	num_t h = 1, mu = 0.0, xi = 0.0, eta = 1.0, sigma = 3;
////////////	size_t Np = 1;	
////////////	size_t Nx = 3;
////////////	std::ofstream fout("result.txt");
////////////	
////////////	Mesh m(Nx, Nx, Nx, h, 0, 0, 0);
////////////	m.generate_mesh();
////////////	m.get_element_interface_matrix();
////////////
////////////	size_matrix_t mim = m.get_element_interface_matrix();
////////////	size_t total_element = m.total_element();
////////////	size_t total_interface = m.total_interface();
////////////	m.print_index_of_elements(fout);
////////////	fout << std::endl;
////////////	fout << "element_interface_matrix" << std::endl;
////////////	for (size_t ie = 0; ie < total_element; ie++)
////////////	{
////////////		fout << "id= " << ie << "\t";
////////////		for (size_t j = 0; j < 6; j++)
////////////		{
////////////			fout << mim(ie, j) << " ";
////////////		}
////////////		fout << std::endl;
////////////	}
////////////
////////////	fout << std::endl << "is boundary:" << std::endl;
////////////	std::vector<bool> bdmark = m.get_boundary_mark();
////////////	for (size_t i = 0; i < bdmark.size(); i++)
////////////	{
////////////		fout << "ii= " << i << "mark= " << bdmark.at(i) << std::endl;
////////////	}
////////////
////////////	fout << "sweep order" << std::endl;
////////////	size_array_t ans;
////////////	ans = m.get_sweep_order(Quadrant1);
////////////	for (size_t i = 0; i < total_element; i++)
////////////	{
////////////		fout << ans.at(i) << " ";
////////////	}
////////////	fout << std::endl;
////////////
////////////	ans = m.get_sweep_order(Quadrant2);
////////////	for (size_t i = 0; i < total_element; i++)
////////////	{
////////////		fout << ans.at(i) << " ";
////////////	}
////////////	fout << std::endl;
////////////
////////////	ans = m.get_sweep_order(Quadrant3);
////////////	for (size_t i = 0; i < total_element; i++)
////////////	{
////////////		fout << ans.at(i) << " ";
////////////	}
////////////	fout << std::endl;
////////////
////////////	ans = m.get_sweep_order(Quadrant4);
////////////	for (size_t i = 0; i < total_element; i++)
////////////	{
////////////		fout << ans.at(i) << " ";
////////////	}
////////////	fout << std::endl;
////////////
////////////	ans = m.get_sweep_order(Quadrant5);
////////////	for (size_t i = 0; i < total_element; i++)
////////////	{
////////////		fout << ans.at(i) << " ";
////////////	}
////////////	fout << std::endl;
////////////
////////////	ans = m.get_sweep_order(Quadrant6);
////////////	for (size_t i = 0; i < total_element; i++)
////////////	{
////////////		fout << ans.at(i) << " ";
////////////	}
////////////	fout << std::endl;
////////////
////////////	ans = m.get_sweep_order(Quadrant7);
////////////	for (size_t i = 0; i < total_element; i++)
////////////	{
////////////		fout << ans.at(i) << " ";
////////////	}
////////////	fout << std::endl;
////////////
////////////	ans = m.get_sweep_order(Quadrant8);
////////////	for (size_t i = 0; i < total_element; i++)
////////////	{
////////////		fout << ans.at(i) << " ";
////////////	}
////////////	fout << std::endl;
////////////
////////////	std::vector<vector3> ec = m.elements_center();
////////////	fout << "Element center" << std::endl;
////////////	for (size_t i = 0; i < ec.size(); i++)
////////////	{
////////////		fout << ec.at(i).x << "\t" << ec.at(i).y << "\t" << ec.at(i).z << std::endl;
////////////	}
////////////	fout.close();
////////////	system("pause");
////////////	return 0;
////////////}
////////////
////////////int main__() {
////////////	for (size_t i = 0; i < 6; i++)
////////////	{
////////////		std::cout
////////////			<< (InterfaceDirection_list.at(i) == SurfaceElementRelation::B) << " "
////////////			<< (InterfaceDirection_list.at(i) == SurfaceElementRelation::F) << " "
////////////			<< (InterfaceDirection_list.at(i) == SurfaceElementRelation::L) << " "
////////////			<< (InterfaceDirection_list.at(i) == SurfaceElementRelation::R) << " "
////////////			<< (InterfaceDirection_list.at(i) == SurfaceElementRelation::D) << " "
////////////			<< (InterfaceDirection_list.at(i) == SurfaceElementRelation::U) << " "
////////////			<< std::endl;
////////////	}
////////////	system("pause");
////////////	return 0;
////////////}
//////////
//////////
//////////int main___() {
//////////	
//////////	std::ofstream fout("result.txt");
//////////	system_matrix_angle m(2, 2, 1, 0, 0, 1);
//////////	
//////////	fout << "x:" << std::endl << m.x();
//////////	fout << "y:" << std::endl << m.y();
//////////	fout << "z:" << std::endl << m.z();
//////////	fout << "sys:" << std::endl << m.system_matrix();
//////////	fout << "sys_lu:" << std::endl << m.system_matrixlu();
//////////	fout << "sys_iv:" << std::endl << m.system_matrixiv();
//////////
//////////	for (size_t i = 0; i < kInterfaceTotal.at(3); i++)
//////////	{
//////////		fout << "lift_dir_" << i << ":" << std::endl << m.lift_matrix(InterfaceDirection_list.at(i));
//////////		fout << "flux_dir_" << i << ":" << std::endl << m.flux_matrix(InterfaceDirection_list.at(i));
//////////	}
//////////
//////////	std::ofstream fout1("result1.txt");
//////////	matrices3d m3d(2, 2);
//////////	m3d.initialize();
//////////	fout1 << "x:" << std::endl << m3d.x();
//////////	fout1 << "y:" << std::endl << m3d.y();
//////////	fout1 << "z:" << std::endl << m3d.z();
//////////	fout1 << "m:" << std::endl << m3d.m();
//////////	fout1 << "sx:" << std::endl << m3d.sx();
//////////	fout1 << "sy:" << std::endl << m3d.sy();
//////////	fout1 << "sz:" << std::endl << m3d.sz();
//////////	for (size_t i = 0; i < kInterfaceTotal.at(3); i++)
//////////	{
//////////		fout1 << "lift_dir_" << i << ":" << std::endl << m3d.lift_matrix(InterfaceDirection_list.at(i));
//////////		fout1 << "flux_dir_" << i << ":" << std::endl << m3d.flux_matrix(InterfaceDirection_list.at(i));
//////////	}
//////////	fout.close();
//////////	system("pause");
//////////	return 0;
//////////}
//////////
//////////
//////////
//////////#include <ctime>
//////////
//////////void savemat(const double* data, std::vector<size_t> sz, const char* file_name, const char* variable_name, bool clear_file = false) {
//////////	MATFile* pmat = NULL;
//////////
//////////	if (clear_file == false)
//////////	{
//////////		pmat = matOpen(file_name, "r");
//////////
//////////		if (pmat != NULL) {
//////////			std::vector<const char*> variable_name_read_list;
////////////			std::vector<mxArray*> mxarray_read_list;
////////////			variable_name_read_list.clear();
////////////			mxarray_read_list.clear();
////////////			const char* tmp_name;
////////////			mxarray_read_list.push_back(matGetNextVariable(pmat, &tmp_name));
////////////			auto itr = mxarray_read_list.end() - 1;
////////////			while (*itr != NULL) {
////////////				variable_name_read_list.push_back(tmp_name);
////////////				mxarray_read_list.push_back(matGetNextVariable(pmat, &tmp_name));
////////////				itr = mxarray_read_list.end() - 1;
////////////			}
////////////			matClose(pmat);
////////////			pmat = matOpen(file_name, "w");
////////////			for (size_t i = 0; i < mxarray_read_list.size(); i++)
////////////			{
////////////				if (mxarray_read_list.at(i) != NULL)
////////////				{
////////////					matPutVariable(pmat, variable_name_read_list.at(i), mxarray_read_list.at(i));
////////////				}
////////////			}
////////////
////////////		}
////////////		if (pmat == NULL)
////////////			pmat = matOpen(file_name, "w");
////////////	}
////////////	else {
////////////		pmat = matOpen(file_name, "w");
////////////	}
////////////	size_t* dims = new size_t[sz.size()];
////////////	for (size_t i = 0; i < sz.size(); i++)
////////////	{
////////////		dims[i] = sz.at(i);
////////////	}
////////////
////////////	mxClassID cid;
////////////	mxArray* pa = mxCreateNumericArray(sz.size(), dims, mxDOUBLE_CLASS, mxREAL);
////////////	if (dims != nullptr)
////////////		delete[] dims;
////////////	if (pa == NULL)
////////////		throw(std::exception("can not create output mxArray."));
////////////	double* par = mxGetPr(pa);
////////////	for (size_t i = 0; i < mxGetNumberOfElements(pa); i++)
////////////	{
////////////		par[i] = data[i];
////////////	}
////////////
////////////	matPutVariable(pmat, variable_name, pa);
////////////	if (pmat != NULL)
////////////		matClose(pmat);
////////////	if (pa != nullptr)
////////////		mxDestroyArray(pa);
////////////}
////////////
//////#include <ctime>
//////
//////int main() {
//////	////std::ofstream fout("Mesh.txt");
//////	//Mesh m(15, 15, 15, 2.0/3, 0, 0, 0);
//////	//m.generate();
//////	//m.cache_sweep_order();
//////	////fout << m;
//////	////fout.close();
//////
//////	////std::cout << "id\tix\tiy\tiz" << std::endl;
//////	////for (size_t i = 0; i < 3*3*3; i++)
//////	////{
//////	////	size_t ix, iy, iz;
//////	////	xlib::ind2sub(3, 3, 3, i, ix, iy, iz);
//////	////	std::cout << i << "\t" << ix << "\t" << iy << "\t" << iz << std::endl;
//////	////}
//////
//////	////std::cout << m.memory_element(8);
//////	////std::cout << m.memory_surface(4);
//////	//std::cout << m.memory_element_total(8) << std::endl;
//////	//std::cout << m.memory_surface_total(4) << std::endl;
//////
//////	//const MKL_INT N = 3;
//////	//double* a = xlib::mkl_ext::xcalloc<double>(N);
//////	//for (size_t i = 0; i < N; i++)
//////	//{
//////	//	a[i] = (double)rand() / RAND_MAX;
//////	//}
//////	//double* b = xlib::mkl_ext::xcalloc<double>(N);
//////	//double* y = xlib::mkl_ext::xcalloc<double>(N);
//////	//double scalea = 1.0, scaleb = 0.0, r , shiftb = 1.0;
//////	//std::cin >> r;
//////	//Vector va(N, a);
//////	//Vector vy(N, y);
//////	//vdLinearFrac(N, a, b, scalea, r, scaleb, shiftb, y);
//////	//std::cout << va;
//////	//std::cout << vy;
//////	//mkl_free(a);
//////	// 
//////
//////	std::cout << "Solver Start." << std::endl;
//////	size_t testid = 3;
//////	std::ofstream fout("result.txt");
//////	std::cout << "input test case id (default= " << testid << ")" << std::endl;
//////	std::cin >> testid;
//////	analytical_solution::set_test_id(testid);
//////	//size_t nx = 15, ny = 15, nz = 15;
//////	size_t nx = 2, ny = 1, nz = 1;
//////	std::cout << "intput nx, ny and nz:" << std::endl;
//////	std::cin >> nx >> ny >> nz;
//////	size_t np = 2;
//////	std::cout << "input np" << std::endl;
//////	std::cin >> np;
//////	num_t h = 2.0 / 15;
//////	//num_t h = 1.0;
//////	std::cout << "input h" << std::endl;
//////	std::cin >> h;
//////	num_t mu = 0.0, xi = 0.0, eta = 0.0, sigma = 1.0;
//////	//std::cout << "input h (default= " << h << ")" << std::endl;
//////	//std::cin >> h;
//////	std::srand(std::clock());
//////	mu = (num_t)rand() / RAND_MAX;
//////	xi = (num_t)rand() / RAND_MAX;
//////	eta = (num_t)rand() / RAND_MAX;
//////	mu = 2 * mu - 1; xi = 2 * xi - 1; eta = 2 * eta - 1;
//////	mu = -0.5; xi = 0.5; eta = 0.5;
//////	//mu = -1.0; xi = -1.0; eta = -1.0;
//////	//mu = -1.0; xi = 0; eta = 0;
//////	
//////	num_t r = mu*mu + xi*xi + eta*eta;
//////	r = std::pow(r, 0.5);
//////	//mu = mu / r; xi = xi / r; eta = eta / r;	
//////	std::cout << "input mu, xi, eta " << std::endl;
//////	std::cin >> mu >> xi >> eta;
//////	std::cout << "mu= " << mu << "\txi= " << xi << "\teta= " << eta << std::endl;
//////	
//////	
//////
//////	system_matrix_angle sma(h, np, mu, xi, eta, sigma);
//////		
//////	std::cout << "system matrix generate start." << std::endl;
//////	Mesh m(nx, ny, nz, h, 0.0, 0.0, 0.0);
//////	std::cout << "Mesh generate start." << std::endl;
//////	m.generate();
//////	std::cout << "Mesh generated." << std::endl;
//////	m.cache_sweep_order();
//////	std::cout << "Mesh sweep order finished." << std::endl;
//////	//fout << m;
//////	fout << sma;
//////	boundary_generator bdg(sma, m);
//////	fout << "boundary x" << std::endl << bdg.x();
//////	fout << "boundary y" << std::endl << bdg.y();
//////	fout << "boundary z" << std::endl << bdg.z();
//////	source_generator sg(sma, m);
//////	fout << "source x" << std::endl << sg.x();
//////	fout << "source y" << std::endl << sg.y();
//////	fout << "source z" << std::endl << sg.z();
//////	std::cout << "source boundary generated." << std::endl;
//////	//fout << sma;
//////
//////	//std::cout << "meshinfo" << m;
//////
//////	//fout << "source " << std::endl;
//////	analytical_solution::set_test_id(testid);
//////	num_p ana_test_p = xlib::mkl_ext::xcalloc<num_t>(sg.x().size());
//////	Vector v(sg.x().size(), ana_test_p);
//////	sg.calculate(v);
//////	//fout << v << std::endl;
//////
//////	solver sol(sma, m);
//////
//////	num_p datam = xlib::mkl_ext::xcalloc<num_t>(m.elements_total());
//////	Vector solm(m.elements_total(), datam);
//////	fout << "analytical solution" << std::endl;
//////	sol.solve_analytical(testid);
//////	sol.solution_mean(solm);
//////	fout << sol.solution() << std::endl;
//////	//fout << solm << std::endl;
//////
//////	std::vector<size_t> sz;
//////	sz.push_back(solm.size());
//////
//////
//////
//////	xlib::imatlab::save_MAT_data(solm.ptr(), sz, "result.mat", "ustd", true);
//////
//////
//////	fout << "numerical solution" << std::endl;
//////	std::time_t t0 = std::clock();
//////	sol.solve_test(testid);
//////	sol.solution_mean(solm);
//////	std::time_t t1 = std::clock();
//////	std::cout << "solve costs " << (double)(t1 - t0) / CLOCKS_PER_SEC;
//////	fout << sol.solution() << std::endl;
//////	//fout << solm << std::endl;
//////	xlib::imatlab::save_MAT_data(solm.ptr(), sz, "result.mat", "u");
//////	fout.close();
//////
//////	Engine* ep;
//////	ep = engOpen("\0");
//////	engEvalString(ep, "cd C:/Workspace/LBTE3DDGNodal/build");
//////	engEvalString(ep, "compare");	
//////	//engClose(ep);
//////	return 0;
//////}
//////
//////#include "memory_manager.h"
//////#include "solver_multi_angle.h"
//////#include <string>
//////
//////#include "globals_mul.h"
//////
//////int main__() {
//////	size_t nx = 45, ny = 45, nz = 45, np = 2, na = 38;
//////	linear_boltzmann_transport_equation_solver::num_t sigma = 1.0;
//////	std::string file_scattering_matrix = "scattering_matrix_c.mat";
//////	std::string file_angle_info = "angleinfo_c.mat";
//////	linear_boltzmann_transport_equation_solver::ProblemDefinition pd(
//////		nx,
//////		ny,
//////		nz,
//////		np,
//////		na,
//////		sigma,
//////		file_scattering_matrix,
//////		file_angle_info
//////	);
//////	linear_boltzmann_transport_equation_solver::AnalyticalSolution::test_id_set(2);
//////	linear_boltzmann_transport_equation_solver::Mesh Mesh(
//////		pd.nx(), pd.ny(), pd.nz(), pd.h(), pd.xc(), pd.yc(), pd.zc());
//////	std::cout << "Mesh start" << std::endl;
//////	Mesh.generate();
//////	Mesh.cache_sweep_order();
//////	std::cout << "init start" << std::endl;
//////	linear_boltzmann_transport_equation_solver::SolverMultiAngle solver(pd, Mesh);
//////	std::cout << "init finished" << std::endl;
//////	solver.SolveAnalyticalBySolve();
//////	solver.SaveSolution("result.mat","u");
//////	//solver.SolveAnalyticalByDirect();
//////	//solver.SaveSolution("result.mat", "ustd");
//////	std::cout << "solve finished" << std::endl;
//////	std::ofstream fout("result.txt");
//////	fout << "scattering matrix:" << std::endl << solver.scattering_matrix();
//////	fout << "mu:" << std::endl << solver.mu();
//////	fout << "xi:" << std::endl << solver.xi();
//////	fout << "eta:" << std::endl << solver.eta();
//////	fout.close();
//////	return 0;
//////}
//////} 
//////} 

#include <fstream>
#include <ctime>

#include "global.h"
#include "mesh.h"
#include "matrices.h"
#include "system_matrices.h"
#include "solver.h"


using namespace discontinues_galerkin_nodal_solver;

void mesh_tests() {
  std::ofstream fout("result.txt");
  size_t nx, ny, nz;
  num_t xc, yc, zc, h;
  nx = ny = nz = 15;
  xc = yc = zc = 0;
  h = 1;
  Mesh m(nx, ny, nz, h, xc, yc, zc);
  fout << m << std::endl;
  fout.close();
}

void matrices_tests() {
  //size_t i = 4;
  //Vector x{ i };
  //std::vector<size_t> sz{ i, 1 };  
  //xlib::imatlab::MatFileHelper::LoadCellDataRaw((num_p)x.p(), "matrices.mat",
  //  "x", i, sz);
  //std::cout << x << std::endl;
  std::ofstream fout("result.txt");
  num_t h = 1;
  size_t np = 2;
  const MatricesData& md = MatricesData::instance();
  Matrices1D m1d(h, np);
  m1d.Initialize();
  Matrices2D m2d(h, np);
  m2d.Initialize();
  Matrices3D m3d(h, np);
  m3d.Initialize();
  fout << md;
  fout << m1d;
  fout << m2d;
  fout << m3d;
  fout.close();
}

//void bet() {    
//  //BOOST_THROW_EXCEPTION( ExceptionOutOfRange() );
//  //throw ExceptionOutOfRange();
//  BOOST_THROW_EXCEPTION(ExceptionOutOfRange() << InfoOutOfRange(1, 2, 3)) ;
//    //InfoMinSizeT(1) << InfoMaxSizeT(2) << InfoWrongSizeT(4);
//}
//} 

void system_matrix_test() {
  std::ofstream fout("result.txt");
  size_t np = 2;
  num_t h = 1;
  num_t mu = 0.5, xi = 0.5, eta = 0.5;
  num_t sigma = 2.0;
 
  
  SystemMatrix sm(h, np, mu, xi, eta, sigma);
  fout << sm;
  fout.close();
}

#include "source.h"
void source_test() {
  std::ofstream fout("result.txt");
  size_t np = 2;
  num_t h = 1;
  num_t mu = -0.3, xi = 0.7, eta = -0.5;
  num_t sigma = 2.0;
  size_t nx, ny, nz;
  num_t xc, yc, zc;
  nx = ny = nz = 3;
  xc = yc = zc = 0;
  Mesh mesh{ nx, ny, nz, h, xc, yc, zc };
  SystemMatrix system{ h, np, mu, xi, eta, sigma };
  AnalyticalSource sor{ mesh, system, 0};
  fout << "x of sor1:\n" << sor.x();
  fout << "y of sor1:\n" << sor.y();
  fout << "z of sor1:\n" << sor.z();
  sor.CalculateValue();
  fout << "source of test id " << sor.id() << " : " << std::endl;
  fout << sor.Value();
  AnalyticalSource sor2{ mesh, system, 1, sor.x(), sor.y(), sor.z() };
  sor2.CalculateValue();
  fout << "source of test id " << sor2.id() << " : " << std::endl;
  fout << sor2.Value();
  fout.close();
}

void source_test2() {
  std::ofstream fout("result.txt");
  size_t np = 2;
  num_t h = 1;
  num_t mu = -0.3, xi = 0.7, eta = -0.5;
  num_t sigma = 2.0;
  size_t nx, ny, nz;
  num_t xc, yc, zc;
  nx = ny = nz = 3;
  xc = yc = zc = 0;
  MeshWithCoordinate mesh{ nx, ny, nz, h, xc, yc, zc,
    SystemMatrix{ h, np, mu, xi, eta, sigma } };  
  SystemMatrix system{ h, np, mu, xi, eta, sigma };
  AnalyticalSourceMeshCoordinates sor{ mesh, system, 0 };
  fout << "x of sor1:\n" << sor.x();
  fout << "y of sor1:\n" << sor.y();
  fout << "z of sor1:\n" << sor.z();
  sor.CalculateValue();
  fout << "source of test id " << sor.id() << " : " << std::endl;
  fout << sor.Value();
  AnalyticalSourceMeshCoordinates sor2{ mesh, system, 1 };
  sor2.CalculateValue();
  fout << "source of test id " << sor2.id() << " : " << std::endl;
  fout << sor2.Value();
  fout.close();
}

#include "boundary.h"
void boundary_test() {
  std::ofstream fout("result.txt");
  size_t np = 2;
  num_t h = 1;
  num_t mu = -0.3, xi = 0.7, eta = -0.5;
  num_t sigma = 2.0;
  size_t nx, ny, nz;
  num_t xc, yc, zc;
  nx = ny = nz = 3;
  xc = yc = zc = 0;
  Mesh mesh{ nx, ny, nz, h, xc, yc, zc };
  SystemMatrix system{ h, np, mu, xi, eta, sigma };
  AnalyticalBoundary bnd{ mesh, system, 0 };
  fout << "x of sor1:\n" << bnd.x();
  fout << "y of sor1:\n" << bnd.y();
  fout << "z of sor1:\n" << bnd.z();
  bnd.CalculateValue();
  fout << "source of test id " << bnd.id() << " : " << std::endl;
  fout << bnd.Value();
  AnalyticalBoundary bnd2{ mesh, system, 1, bnd.x(), bnd.y(), bnd.z(), bnd.boundary_mask() };
  bnd2.CalculateValue();
  fout << "source of test id " << bnd2.id() << " : " << std::endl;
  fout << bnd2.Value();
  fout.close();
}

void boundary_test2() {
  std::ofstream fout("result.txt");
  size_t np = 2;
  num_t h = 1;
  num_t mu = -0.3, xi = 0.7, eta = -0.5;
  num_t sigma = 2.0;
  size_t nx, ny, nz;
  num_t xc, yc, zc;
  nx = ny = nz = 3;
  xc = yc = zc = 0;
  MeshWithCoordinate mesh{ nx, ny, nz, h, xc, yc, zc,
    SystemMatrix{ h, np, mu, xi, eta, sigma } };
  SystemMatrix system{ h, np, mu, xi, eta, sigma };
  AnalyticalBoundaryMeshCoordinate bnd{ mesh, system, 0 };
  fout << "x of sor1:\n" << bnd.x();
  fout << "y of sor1:\n" << bnd.y();
  fout << "z of sor1:\n" << bnd.z();
  bnd.CalculateValue();
  fout << "source of test id " << bnd.id() << " : " << std::endl;
  fout << bnd.Value();
  AnalyticalBoundaryMeshCoordinate bnd2{ mesh, system, 1};
  bnd2.CalculateValue();
  fout << "source of test id " << bnd2.id() << " : " << std::endl;
  fout << bnd2.Value();
  fout.close();
}

void mesh_with_cood_test() {
  std::ofstream fout("result.txt");
  size_t nx, ny, nz;
  num_t xc, yc, zc;
  nx = ny = nz = 3;
  xc = yc = zc = 0;
  num_t h = 1;
  num_t mu = -0.3, xi = 0.7, eta = -0.5;
  num_t sigma = 2.0;
  size_t np = 2;
  SystemMatrix system{ h, np, mu, xi, eta, sigma };
  std::vector<Vector> xs, ys, zs;
  xs.resize(3); ys.resize(3), zs.resize(3);
  for (size_t i = 0; i < 3; ++i)
  {
    xs.at(i) = system.xs(kSurfaceDirectionList.at(i));
    ys.at(i) = system.ys(kSurfaceDirectionList.at(i));
    zs.at(i) = system.zs(kSurfaceDirectionList.at(i));
  }
  //MeshWithCoordinate mesh{ nx, ny, nz, h, xc, yc, zc,
  //    system.xe(), system.ye(), system.ze(), xs, ys, zs };
  MeshWithCoordinate mesh{ nx, ny, nz, h, xc, yc, zc, 
      SystemMatrix{h, np, mu, xi, eta, sigma} };
  fout << "xe:\n" << mesh.xe();
  fout << "ye:\n" << mesh.ye();
  fout << "ze:\n" << mesh.ze();
  fout << "xs:\n" << mesh.xs();
  fout << "ys:\n" << mesh.ys();
  fout << "zs:\n" << mesh.zs();
  fout.close();
}

void dgn_solve_test() {
  ProblemDefinitionSingleAngle pd;
  pd.Load("input.txt");
  SystemMatrix system{ pd.h(), pd.np(), pd.mu(), pd.xi(), pd.eta(), pd.sigma() };
  MeshWithCoordinate mesh{ pd.nx(), pd.ny(), pd.nz(), pd.h(), pd.xc(), pd.yc(), 
        pd.zc(), system};
  SolverSingleAngle solver(mesh);
  solver.Initialization();
  solver.Procceed();
  solver.SetCondition();
  solver.Procceed();
  solver.Solve();
  Vector x = solver.xe();
  Vector y = solver.ye();
  Vector z = solver.ze();
  Vector u = solver.solution();
  std::vector<size_t> sz{ static_cast<size_t>(x.row()), 1 };
  std::string filename = "result.mat";
  xlib::imatlab::MatFileHelper::SaveDataRaw(x.p(), filename, "x", sz);    
  xlib::imatlab::MatFileHelper::SaveDataRaw(y.p(), filename, "y", sz);
  xlib::imatlab::MatFileHelper::SaveDataRaw(z.p(), filename, "z", sz);
  xlib::imatlab::MatFileHelper::SaveDataRaw(u.p(), filename, "u", sz);

  AnalyticalSolutionMeshCoordinates sol{ mesh, system, pd.test_id() };
  sol.CalculateValue();
  Vector ustd = sol.Value();
  xlib::imatlab::MatFileHelper::SaveDataRaw(ustd.p(), filename, "ustd", sz);
}


void dgn_solver_ext_test() {
  ProblemDefinitionSingleAngle pd;
  pd.Load("input.txt");
  SystemMatrix system{ pd.h(), pd.np(), pd.mu(), pd.xi(), pd.eta(), pd.sigma() };
  MeshWithCoordinate mesh{ pd.nx(), pd.ny(), pd.nz(), pd.h(), pd.xc(), pd.yc(),
    pd.zc(), system };
  AnalyticalBoundaryMeshCoordinate bd{ mesh, system, pd.test_id() };
  AnalyticalSourceMeshCoordinates sor{ mesh, system, pd.test_id() };
  bd.CalculateValue();
  sor.CalculateValue();
  SolverSingleAngle solver(mesh, sor.Value(), bd.Value());

  time_t t1 = std::clock();

  solver.Initialization();
  solver.Procceed();

  time_t t2 = std::clock();

  solver.SetCondition();
  solver.Procceed();

  time_t t3 = std::clock();
  solver.Solve();
  time_t t4 = std::clock();
#ifdef MEATIME
  std::cout << "init cost = " << t2 - t1 << std::endl;
  std::cout << "condition cost = " << t3 - t2 << std::endl;
  std::cout << "solve cost = " << t4 - t3 << std::endl;
#endif
  Vector x = solver.xe();
  Vector y = solver.ye();
  Vector z = solver.ze();
  Vector u = solver.solution();
  std::vector<size_t> sz{ static_cast<size_t>(x.row()), 1 };
  std::string filename = "result.mat";
  xlib::imatlab::MatFileHelper::SaveDataRaw(x.p(), filename, "x", sz);
  xlib::imatlab::MatFileHelper::SaveDataRaw(y.p(), filename, "y", sz);
  xlib::imatlab::MatFileHelper::SaveDataRaw(z.p(), filename, "z", sz);
  xlib::imatlab::MatFileHelper::SaveDataRaw(u.p(), filename, "u", sz);

  AnalyticalSolutionMeshCoordinates sol{ mesh, system, pd.test_id() };
  sol.CalculateValue();
  Vector ustd = sol.Value();
  xlib::imatlab::MatFileHelper::SaveDataRaw(ustd.p(), filename, "ustd", sz);
}

void repeat_test_ker() {
  ProblemDefinitionSingleAngle pd;
  pd.Load("input.txt");
  SystemMatrix system{ pd.h(), pd.np(), pd.mu(), pd.xi(), pd.eta(), pd.sigma() };
}

void repeat_test() {
  const size_t N = 100;
  for (size_t i = 0; i < N; ++i)
    repeat_test_ker();
}
//#define EXCEDEBUG
int main() {

  //mesh_with_cood_test();
#ifdef EXCEDEBUG
  try {
   dgn_solver_ext_test();
   // repeat_test_ker();
  //system_matrix_test();
  }
  catch (boost::exception& e){
    std::cout << boost::diagnostic_information(e);
    throw;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    throw;
  }
#else
  dgn_solver_ext_test();
#endif
  return 0;
}