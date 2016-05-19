#include <mkl.h>
#include <iostream>
#include <fstream>
#include <random>

//int main_() {	
//	num_t h = 1, mu = 0.0, xi = 0.0, eta = 1.0, sigma = 3;
//	size_t Np = 1;	
//	size_t Nx = 3;
//	std::ofstream fout("result.txt");
//	
//	Mesh m(Nx, Nx, Nx, h, 0, 0, 0);
//	m.generate_mesh();
//	m.get_element_interface_matrix();
//
//	size_matrix_t mim = m.get_element_interface_matrix();
//	size_t total_element = m.total_element();
//	size_t total_interface = m.total_interface();
//	m.print_index_of_elements(fout);
//	fout << std::endl;
//	fout << "element_interface_matrix" << std::endl;
//	for (size_t ie = 0; ie < total_element; ie++)
//	{
//		fout << "id= " << ie << "\t";
//		for (size_t j = 0; j < 6; j++)
//		{
//			fout << mim(ie, j) << " ";
//		}
//		fout << std::endl;
//	}
//
//	fout << std::endl << "is boundary:" << std::endl;
//	std::vector<bool> bdmark = m.get_boundary_mark();
//	for (size_t i = 0; i < bdmark.size(); i++)
//	{
//		fout << "ii= " << i << "mark= " << bdmark.at(i) << std::endl;
//	}
//
//	fout << "sweep order" << std::endl;
//	size_array_t ans;
//	ans = m.get_sweep_order(Quadrant1);
//	for (size_t i = 0; i < total_element; i++)
//	{
//		fout << ans.at(i) << " ";
//	}
//	fout << std::endl;
//
//	ans = m.get_sweep_order(Quadrant2);
//	for (size_t i = 0; i < total_element; i++)
//	{
//		fout << ans.at(i) << " ";
//	}
//	fout << std::endl;
//
//	ans = m.get_sweep_order(Quadrant3);
//	for (size_t i = 0; i < total_element; i++)
//	{
//		fout << ans.at(i) << " ";
//	}
//	fout << std::endl;
//
//	ans = m.get_sweep_order(Quadrant4);
//	for (size_t i = 0; i < total_element; i++)
//	{
//		fout << ans.at(i) << " ";
//	}
//	fout << std::endl;
//
//	ans = m.get_sweep_order(Quadrant5);
//	for (size_t i = 0; i < total_element; i++)
//	{
//		fout << ans.at(i) << " ";
//	}
//	fout << std::endl;
//
//	ans = m.get_sweep_order(Quadrant6);
//	for (size_t i = 0; i < total_element; i++)
//	{
//		fout << ans.at(i) << " ";
//	}
//	fout << std::endl;
//
//	ans = m.get_sweep_order(Quadrant7);
//	for (size_t i = 0; i < total_element; i++)
//	{
//		fout << ans.at(i) << " ";
//	}
//	fout << std::endl;
//
//	ans = m.get_sweep_order(Quadrant8);
//	for (size_t i = 0; i < total_element; i++)
//	{
//		fout << ans.at(i) << " ";
//	}
//	fout << std::endl;
//
//	std::vector<vector3> ec = m.elements_center();
//	fout << "element center" << std::endl;
//	for (size_t i = 0; i < ec.size(); i++)
//	{
//		fout << ec.at(i).x << "\t" << ec.at(i).y << "\t" << ec.at(i).z << std::endl;
//	}
//	fout.close();
//	system("pause");
//	return 0;
//}
//
//int main__() {
//	for (size_t i = 0; i < 6; i++)
//	{
//		std::cout
//			<< (interface_direction_list.at(i) == interface_direction::B) << " "
//			<< (interface_direction_list.at(i) == interface_direction::F) << " "
//			<< (interface_direction_list.at(i) == interface_direction::L) << " "
//			<< (interface_direction_list.at(i) == interface_direction::R) << " "
//			<< (interface_direction_list.at(i) == interface_direction::D) << " "
//			<< (interface_direction_list.at(i) == interface_direction::U) << " "
//			<< std::endl;
//	}
//	system("pause");
//	return 0;
//}


#include "system_matrices.h"
#include "globals.h"
#include "types.h"
#include "utilities.h"
#include "xlib_mkl_ext.hpp"
#include "mesh.h"
int main___() {
	
	std::ofstream fout("result.txt");
	dgn::system_matrix_angle m(2, 2, 1, 0, 0, 1);
	
	fout << "x:" << std::endl << m.x();
	fout << "y:" << std::endl << m.y();
	fout << "z:" << std::endl << m.z();
	fout << "sys:" << std::endl << m.system_matrix();
	fout << "sys_lu:" << std::endl << m.system_matrixlu();
	fout << "sys_iv:" << std::endl << m.system_matrixiv();

	for (size_t i = 0; i < dgn::INTERFACE_TOTAL.at(3); i++)
	{
		fout << "lift_dir_" << i << ":" << std::endl << m.lift_matrix(dgn::interface_direction_list.at(i));
		fout << "flux_dir_" << i << ":" << std::endl << m.flux_matrix(dgn::interface_direction_list.at(i));
	}

	std::ofstream fout1("result1.txt");
	dgn::matrices3d m3d(2, 2);
	m3d.initialize();
	fout1 << "x:" << std::endl << m3d.x();
	fout1 << "y:" << std::endl << m3d.y();
	fout1 << "z:" << std::endl << m3d.z();
	fout1 << "m:" << std::endl << m3d.m();
	fout1 << "sx:" << std::endl << m3d.sx();
	fout1 << "sy:" << std::endl << m3d.sy();
	fout1 << "sz:" << std::endl << m3d.sz();
	for (size_t i = 0; i < dgn::INTERFACE_TOTAL.at(3); i++)
	{
		fout1 << "lift_dir_" << i << ":" << std::endl << m3d.lift_matrix(dgn::interface_direction_list.at(i));
		fout1 << "flux_dir_" << i << ":" << std::endl << m3d.flux_matrix(dgn::interface_direction_list.at(i));
	}
	fout.close();
	system("pause");
	return 0;
}


#include "solver.h"
#include "system_matrices.h"
#include "mesh.h"

#include <mat.h>
#include <mex.h>
#include <engine.h>

#include <ctime>

void savemat(const double* data, std::vector<size_t> sz, const char* file_name, const char* variable_name, bool clear_file = false) {
	MATFile* pmat = NULL;

	if (clear_file == false)
	{
		pmat = matOpen(file_name, "r");

		if (pmat != NULL) {
			std::vector<const char*> variable_name_read_list;
			std::vector<mxArray*> mxarray_read_list;
			variable_name_read_list.clear();
			mxarray_read_list.clear();
			const char* tmp_name;
			mxarray_read_list.push_back(matGetNextVariable(pmat, &tmp_name));
			auto itr = mxarray_read_list.end() - 1;
			while (*itr != NULL) {
				variable_name_read_list.push_back(tmp_name);
				mxarray_read_list.push_back(matGetNextVariable(pmat, &tmp_name));
				itr = mxarray_read_list.end() - 1;
			}
			matClose(pmat);
			pmat = matOpen(file_name, "w");
			for (size_t i = 0; i < mxarray_read_list.size(); i++)
			{
				if (mxarray_read_list.at(i) != NULL)
				{
					matPutVariable(pmat, variable_name_read_list.at(i), mxarray_read_list.at(i));
				}
			}

		}
		if (pmat == NULL)
			pmat = matOpen(file_name, "w");
	}
	else {
		pmat = matOpen(file_name, "w");
	}
	size_t* dims = new size_t[sz.size()];
	for (size_t i = 0; i < sz.size(); i++)
	{
		dims[i] = sz.at(i);
	}

	mxClassID cid;
	mxArray* pa = mxCreateNumericArray(sz.size(), dims, mxDOUBLE_CLASS, mxREAL);
	if (dims != nullptr)
		delete[] dims;
	if (pa == NULL)
		throw(std::exception("can not create output mxArray."));
	double* par = mxGetPr(pa);
	for (size_t i = 0; i < mxGetNumberOfElements(pa); i++)
	{
		par[i] = data[i];
	}

	matPutVariable(pmat, variable_name, pa);
	if (pmat != NULL)
		matClose(pmat);
	if (pa != nullptr)
		mxDestroyArray(pa);
}

#include <ctime>

int main() {
	////std::ofstream fout("mesh.txt");
	//dgn::mesh m(15, 15, 15, 2.0/3, 0, 0, 0);
	//m.generate();
	//m.cache_sweep_order();
	////fout << m;
	////fout.close();

	////std::cout << "id\tix\tiy\tiz" << std::endl;
	////for (size_t i = 0; i < 3*3*3; i++)
	////{
	////	size_t ix, iy, iz;
	////	xlib::ind2sub(3, 3, 3, i, ix, iy, iz);
	////	std::cout << i << "\t" << ix << "\t" << iy << "\t" << iz << std::endl;
	////}

	////std::cout << m.memory_element(8);
	////std::cout << m.memory_surface(4);
	//std::cout << m.memory_element_total(8) << std::endl;
	//std::cout << m.memory_surface_total(4) << std::endl;

	//const MKL_INT N = 3;
	//double* a = xlib::mkl_ext::xcalloc<double>(N);
	//for (size_t i = 0; i < N; i++)
	//{
	//	a[i] = (double)rand() / RAND_MAX;
	//}
	//double* b = xlib::mkl_ext::xcalloc<double>(N);
	//double* y = xlib::mkl_ext::xcalloc<double>(N);
	//double scalea = 1.0, scaleb = 0.0, r , shiftb = 1.0;
	//std::cin >> r;
	//dgn::vector_t va(N, a);
	//dgn::vector_t vy(N, y);
	//vdLinearFrac(N, a, b, scalea, r, scaleb, shiftb, y);
	//std::cout << va;
	//std::cout << vy;
	//mkl_free(a);
	// 

	std::cout << "Solver Start." << std::endl;
	
	//std::ofstream fout("result.txt");

	//size_t nx = 15, ny = 15, nz = 15;
	size_t nx = 15, ny = 15, nz = 15;
	std::cout << "intput nx, ny and nz:" << std::endl;
	std::cin >> nx >> ny >> nz;
	dgn::num_t h = 2.0/nx, mu = 0.0, xi = 0.0, eta = 0.0, sigma = 1.0;
	std::srand(std::clock());
	mu = (dgn::num_t)rand() / RAND_MAX;
	xi = (dgn::num_t)rand() / RAND_MAX;
	eta = (dgn::num_t)rand() / RAND_MAX;
	mu = 2 * mu - 1; xi = 2 * xi - 1; eta = 2 * eta - 1;
	//mu = 0.7; xi = 0.5; eta = -0.3;
	//mu = 0.0; xi = 0.0; eta = -1.0;
	dgn::num_t r = mu*mu + xi*xi + eta*eta;
	r = std::pow(r, 0.5);
	mu = mu / r; xi = xi / r; eta = eta / r;
	std::cout << "mu= " << mu << "\txi= " << xi << "\teta= " << eta << std::endl;
	size_t np = 2;

	dgn::system_matrix_angle sma(h, np, mu, xi, eta, sigma);
		
	
	dgn::mesh m(nx, ny, nz, h, 0, 0, 0);
	std::cout << "mesh generate start." << std::endl;
	m.generate();
	std::cout << "mesh generated." << std::endl;
	m.cache_sweep_order();
	std::cout << "mesh sweep order finished." << std::endl;
	//fout << m;
	dgn::boundary_generator bdg(sma, m);
	//fout << "boundary x" << std::endl << bdg.x();
	//fout << "boundary y" << std::endl << bdg.y();
	//fout << "boundary z" << std::endl << bdg.z();
	dgn::source_generator sg(sma, m);
	//fout << "source x" << std::endl << sg.x();
	//fout << "source y" << std::endl << sg.y();
	//fout << "source z" << std::endl << sg.z();
	std::cout << "source boundary generated." << std::endl;
	//fout << sma;

	size_t testid = 4;

	//fout << "source " << std::endl;
	dgn::analytical_solution::set_test_id(testid);
	dgn::num_p ana_test_p = xlib::mkl_ext::xcalloc<dgn::num_t>(sg.x().size());
	dgn::vector_t v(sg.x().size(), ana_test_p);
	sg.calculate(v);
	//fout << v << std::endl;

	dgn::solver sol(sma, m);

	dgn::num_p datam = xlib::mkl_ext::xcalloc<dgn::num_t>(m.elements_total());
	dgn::vector_t solm(m.elements_total(), datam);
	//fout << "analytical solution" << std::endl;
	sol.solve_analytical(testid);
	sol.solution_mean(solm);
	//fout << sol.solution() << std::endl;
	//fout << solm << std::endl;

	std::vector<size_t> sz;
	sz.push_back(solm.size());



	xlib::imatlab::save_MAT_data(solm.ptr(), sz, "result.mat", "ustd", true);


	//fout << "numerical solution" << std::endl;
	std::time_t t0 = std::clock();
	sol.solve_test(testid);
	sol.solution_mean(solm);
	std::time_t t1 = std::clock();
	std::cout << "solve costs " << (double)(t1 - t0) / CLOCKS_PER_SEC;
	//fout << sol.solution() << std::endl;
	//fout << solm << std::endl;
	xlib::imatlab::save_MAT_data(solm.ptr(), sz, "result.mat", "u");
	//fout.close();

	Engine* ep;
	ep = engOpen("\0");
	engEvalString(ep, "cd C:/Workspace/LBTE3DDGNodal/build");
	engEvalString(ep, "compare");
	engClose(ep);
	return 0;
}

int main_____() {
	dgn::num_t h = 1;
	size_t np = 2;
	std::ofstream fout("result.txt");
	dgn::num_t mu = 0.0, xi = 0.0, eta = -1.0, sigma = 1.0;
	dgn::system_matrix_angle sm1(h, np, mu, xi, eta, sigma);
	fout << sm1;

	mu = 0.0, xi = 0.0, eta = 1.0, sigma = 1.0;
	dgn::system_matrix_angle sm2(h, np, mu, xi, eta, sigma);
	fout << sm2;

	mu = 0.0, xi = 1.0, eta = 0.0, sigma = 1.0;
	dgn::system_matrix_angle sm3(h, np, mu, xi, eta, sigma);
	fout << sm3;

	mu = 0.0, xi = -1.0, eta = 0.0, sigma = 1.0;
	dgn::system_matrix_angle sm4(h, np, mu, xi, eta, sigma);
	fout << sm4;

	fout.close();
	return 0;
}