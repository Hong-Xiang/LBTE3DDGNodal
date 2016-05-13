#pragma once
#include <string>
#include <iostream>
#include <vector>
#include "Global.h"

namespace DG3DNodal {
	void print_vector(std::ostream& os, const double* data, size_t total_element, size_t offset = 1);
	void print_matrix(std::ostream& os, const double* data, size_t total_row, size_t total_column, size_t offset_row = 1, size_t offset_column = 1);

	//Convert integer to string
	const std::string i2s(int i);

	inline size_t column_major_index(size_t irow, size_t icolumn, size_t nrow) {
		return irow + icolumn*nrow;
	}

	void ind2sub2(const size_t Np, const size_t ind, size_t& ix, size_t& iy);
	void ind2sub3(const size_t Np, const size_t ind, size_t& ix, size_t& iy, size_t& iz);
	void sub2ind2(const size_t Np, size_t& ind, const size_t ix, const size_t iy);
	void sub2ind3(const size_t Np, size_t& ind, const size_t ix, const size_t iy, const size_t iz);

	void ind2sub2(const size_t Np, const size_t N, const size_t* ind, size_t* ix, size_t* iy);
	void ind2sub3(const size_t Np, const size_t N, const size_t* ind, size_t* ix, size_t* iy, size_t* iz);
	void sub2ind2(const size_t Np, const size_t N, const size_t* ix, const size_t* iy, size_t* ind);
	void sub2ind3(const size_t Np, const size_t N, const size_t* ix, const size_t* iy, const size_t* iz, size_t* ind);

	void ind2sub2(const size_t Np, const std::vector<size_t>& ind, std::vector<size_t>& ix, std::vector<size_t>& iy);
	void ind2sub3(const size_t Np, const std::vector<size_t>& ind, std::vector<size_t>& ix, std::vector<size_t>& iy, std::vector<size_t>& iz);
	void sub2ind2(const size_t Np, const std::vector<size_t>& ix, const std::vector<size_t>& iy, std::vector<size_t>& ind);
	void sub2ind3(const size_t Np, const std::vector<size_t>& ix, const std::vector<size_t>& iy, const std::vector<size_t>& iz, std::vector<size_t> ind);

	class size_matrix_t : public std::vector < std::vector <size_t> > {
	public:
		size_t& operator()(size_t ix, size_t iy) {
			return at(ix).at(iy);
		}
		size_t operator()(size_t ix, size_t iy) const {
			return at(ix).at(iy);
		}
	};

	class size_array_t : public  std::vector < size_t > {		
		size_t& operator()(size_t ix) {
			return at(ix);
		}
		size_t operator()(size_t ix) const {
			return at(ix);
		}
	};


	struct vector3
	{
		double x, y, z;
	};
	typedef int mesh_refine_lvl_t;

	quadrant quadrant_of_angle(num_t mu, num_t xi, num_t eta);
}