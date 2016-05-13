#include "Utils.h"

#include <sstream>
void DG3DNodal::print_vector(std::ostream & os, const double * data, size_t total_element, size_t offset)
{
	for (size_t i = 0; i < total_element; i++)
	{
		os << data[i*offset];
	}
	os << std::endl;
}

void DG3DNodal::print_matrix(std::ostream & os, const double * data, size_t total_row, size_t total_column, size_t offset_row, size_t offset_column)
{
	for (size_t i = 0; i < total_row; i++)
	{
		for (size_t j = 0; j < total_column; j++)
		{
			os << data[column_major_index(i*offset_row, j*offset_column, total_row*offset_row)] << " ";
		}
		os << std::endl;
	}
}

const std::string DG3DNodal::i2s(int i)
{
	std::stringstream ss;
	ss << i;
	return ss.str();
}

void DG3DNodal::ind2sub2(const size_t Np, const size_t ind, size_t & ix, size_t & iy)
{
	ix = ind % Np;
	iy = ind / Np;
}

void DG3DNodal::ind2sub3(const size_t Np, const size_t ind, size_t & ix, size_t & iy, size_t & iz)
{
	ix = ind % Np;
	iy = ind / Np % Np;
	iz = ind / Np / Np;
}

void DG3DNodal::sub2ind2(const size_t Np, size_t& ind, const size_t ix, const size_t iy)
{
	ind = ix + iy * Np;
}

void DG3DNodal::sub2ind3(const size_t Np, size_t & ind, const size_t ix, const size_t iy, const size_t iz)
{
	ind = ix + iy * Np + iz * Np * Np;
}





void DG3DNodal::ind2sub2(const size_t Np, const size_t N, const size_t * ind, size_t * ix, size_t * iy)
{
	for (size_t i = 0; i < N; i++)
	{
		ind2sub2(Np, ind[i], ix[i], iy[i]);
	}
}

void DG3DNodal::ind2sub3(const size_t Np, const size_t N, const size_t * ind, size_t * ix, size_t * iy, size_t * iz)
{
	for (size_t i = 0; i < N; i++)
	{
		ind2sub3(Np, ind[i], ix[i], iy[i], iz[i]);
	}
}

void DG3DNodal::sub2ind2(const size_t Np, const size_t N, const size_t * ix, const size_t * iy, size_t * ind)
{
	for (size_t i = 0; i < N; i++)
	{
		sub2ind2(Np, ind[i], ix[i], iy[i]);
	}
}

void DG3DNodal::sub2ind3(const size_t Np, const size_t N, const size_t * ix, const size_t * iy, const size_t * iz, size_t * ind)
{
	for (size_t i = 0; i < N; i++)
	{
		sub2ind3(Np, ind[i], ix[i], iy[i], iz[i]);
	}
}

void DG3DNodal::ind2sub2(const size_t Np, const std::vector<size_t>& ind, std::vector<size_t>& ix, std::vector<size_t>& iy)
{
	for (size_t i = 0; i < Np; i++)
	{
		ind2sub2(Np, ind.at(i), ix.at(i), iy.at(i));
	}
}

void DG3DNodal::ind2sub3(const size_t Np, const std::vector<size_t>& ind, std::vector<size_t>& ix, std::vector<size_t>& iy, std::vector<size_t>& iz)
{
	for (size_t i = 0; i < Np; i++)
	{
		ind2sub3(Np, ind.at(i), ix.at(i), iy.at(i), iz.at(i));
	}
}

void DG3DNodal::sub2ind2(const size_t Np, const std::vector<size_t>& ix, const std::vector<size_t>& iy, std::vector<size_t>& ind)
{
	for (size_t i = 0; i < Np; i++)
	{
		sub2ind2(Np, ind.at(i), ix.at(i), iy.at(i));
	}
}

void DG3DNodal::sub2ind3(const size_t Np, const std::vector<size_t>& ix, const std::vector<size_t>& iy, const std::vector<size_t>& iz, std::vector<size_t> ind)
{
	for (size_t i = 0; i < Np; i++)
	{
		sub2ind3(Np, ind.at(i), ix.at(i), iy.at(i), iz.at(i));
	}
}

DG3DNodal::quadrant DG3DNodal::quadrant_of_angle(num_t mu, num_t xi, num_t eta)
{
	if (mu >= 0.0 && xi >= 0.0 && eta >= 0.0)
		return quadrant::quadrant1;
	if (mu < 0.0 && xi >= 0.0 && eta >= 0.0)
		return quadrant::quadrant2;
	if (mu < 0.0 && xi < 0.0 && eta >= 0.0)
		return quadrant::quadrant3;
	if (mu >= 0.0 && xi < 0.0 && eta >= 0.0)
		return quadrant::quadrant4;
	if (mu >= 0.0 && xi >= 0.0 && eta < 0.0)
		return quadrant::quadrant5;
	if (mu < 0.0 && xi >= 0.0 && eta < 0.0)
		return quadrant::quadrant6;
	if (mu < 0.0 && xi < 0.0 && eta < 0.0)
		return quadrant::quadrant7;
	if (mu >= 0.0 && xi < 0.0 && eta < 0.0)
		return quadrant::quadrant8;
}

