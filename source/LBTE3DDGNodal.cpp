#include "Utils.h"
#include "Matrices.h"
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

	const MKL_INT N = 3;
	double* a = xlib::mkl_ext::xcalloc<double>(N);
	for (size_t i = 0; i < N; i++)
	{
		a[i] = (double)rand() / RAND_MAX;
	}
	double* b = xlib::mkl_ext::xcalloc<double>(N);
	double* y = xlib::mkl_ext::xcalloc<double>(N);
	double scalea = 1.0, scaleb = 0.0, r , shiftb = 1.0;
	std::cin >> r;
	dgn::vector_t va(N, a);
	dgn::vector_t vy(N, y);
	vdLinearFrac(N, a, b, scalea, r, scaleb, shiftb, y);
	std::cout << va;
	std::cout << vy;
	mkl_free(a);
	return 0;
}