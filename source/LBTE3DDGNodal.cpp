#include "Utils.h"
#include "Matrices.h"
#include <mkl.h>
#include <iostream>
#include <fstream>
#include <IMATLABC.h>
#include "Global.h"
#include <random>
#include "mkl_vsl.h"
#include <iomanip>
#include "SystemMatrices.h"
using namespace DG3DNodal;

#include "Meshes.h"

int main_() {	
	num_t h = 1, mu = 0.0, xi = 0.0, eta = 1.0, sigma = 3;
	size_t Np = 1;	
	size_t Nx = 3;
	std::ofstream fout("result.txt");
	
	Mesh m(Nx, Nx, Nx, h, 0, 0, 0);
	m.generate_mesh();
	m.get_element_interface_matrix();

	size_matrix_t mim = m.get_element_interface_matrix();
	size_t total_element = m.total_element();
	size_t total_interface = m.total_interface();
	m.print_index_of_elements(fout);
	fout << std::endl;
	fout << "element_interface_matrix" << std::endl;
	for (size_t ie = 0; ie < total_element; ie++)
	{
		fout << "id= " << ie << "\t";
		for (size_t j = 0; j < 6; j++)
		{
			fout << mim(ie, j) << " ";
		}
		fout << std::endl;
	}

	fout << std::endl << "is boundary:" << std::endl;
	std::vector<bool> bdmark = m.get_boundary_mark();
	for (size_t i = 0; i < bdmark.size(); i++)
	{
		fout << "ii= " << i << "mark= " << bdmark.at(i) << std::endl;
	}

	fout << "sweep order" << std::endl;
	size_array_t ans;
	ans = m.get_sweep_order(Quadrant1);
	for (size_t i = 0; i < total_element; i++)
	{
		fout << ans.at(i) << " ";
	}
	fout << std::endl;

	ans = m.get_sweep_order(Quadrant2);
	for (size_t i = 0; i < total_element; i++)
	{
		fout << ans.at(i) << " ";
	}
	fout << std::endl;

	ans = m.get_sweep_order(Quadrant3);
	for (size_t i = 0; i < total_element; i++)
	{
		fout << ans.at(i) << " ";
	}
	fout << std::endl;

	ans = m.get_sweep_order(Quadrant4);
	for (size_t i = 0; i < total_element; i++)
	{
		fout << ans.at(i) << " ";
	}
	fout << std::endl;

	ans = m.get_sweep_order(Quadrant5);
	for (size_t i = 0; i < total_element; i++)
	{
		fout << ans.at(i) << " ";
	}
	fout << std::endl;

	ans = m.get_sweep_order(Quadrant6);
	for (size_t i = 0; i < total_element; i++)
	{
		fout << ans.at(i) << " ";
	}
	fout << std::endl;

	ans = m.get_sweep_order(Quadrant7);
	for (size_t i = 0; i < total_element; i++)
	{
		fout << ans.at(i) << " ";
	}
	fout << std::endl;

	ans = m.get_sweep_order(Quadrant8);
	for (size_t i = 0; i < total_element; i++)
	{
		fout << ans.at(i) << " ";
	}
	fout << std::endl;

	std::vector<vector3> ec = m.elements_center();
	fout << "element center" << std::endl;
	for (size_t i = 0; i < ec.size(); i++)
	{
		fout << ec.at(i).x << "\t" << ec.at(i).y << "\t" << ec.at(i).z << std::endl;
	}
	fout.close();
	system("pause");
	return 0;
}

int main() {
	for (size_t i = 0; i < 6; i++)
	{
		std::cout
			<< (interface_direction_list.at(i) == interface_direction::B) << " "
			<< (interface_direction_list.at(i) == interface_direction::F) << " "
			<< (interface_direction_list.at(i) == interface_direction::L) << " "
			<< (interface_direction_list.at(i) == interface_direction::R) << " "
			<< (interface_direction_list.at(i) == interface_direction::D) << " "
			<< (interface_direction_list.at(i) == interface_direction::U) << " "
			<< std::endl;
	}
	system("pause");
	return 0;
}