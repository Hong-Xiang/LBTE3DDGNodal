#include "Matrices.h"
using namespace DG3DNodal;

SystemMatrixStandard1D* SystemMatrixStandard1D::handle_ = nullptr;
SystemMatrixStandard2D* SystemMatrixStandard2D::handle_ = nullptr;
SystemMatrixStandard3D* SystemMatrixStandard3D::handle_ = nullptr;

SystemMatrixStandard1D & DG3DNodal::SystemMatrixStandard1D::getHandle()
{
	if (handle_ == nullptr)
		handle_ = new SystemMatrixStandard1D();
	return *handle_;
}

const vector_t DG3DNodal::SystemMatrixStandard1D::getX(size_t Np) const
{
	if (Np < MINNP || Np > MAXNP)
		throw SystemMatrixStandardOutOfOrder();
	return vector_t(Np, x.at(Np - MINNP));
}

const matrix_t DG3DNodal::SystemMatrixStandard1D::getM(size_t Np) const
{
	if (Np < MINNP || Np > MAXNP)
		throw SystemMatrixStandardOutOfOrder();
	return matrix_t(Np, Np, m.at(Np - MINNP));
}

const matrix_t DG3DNodal::SystemMatrixStandard1D::getS(size_t Np) const
{
	if (Np < MINNP || Np > MAXNP)
		throw SystemMatrixStandardOutOfOrder();
	return matrix_t(Np, Np, s.at(Np - MINNP));
}

DG3DNodal::SystemMatrixStandard1D::SystemMatrixStandard1D()
{
	x.resize(MAXNP - MINNP + 1);
	m.resize(MAXNP - MINNP + 1);
	s.resize(MAXNP - MINNP + 1);
	for (size_t i = MINNP; i <= MAXNP; i++)
	{
		x.at(i - MINNP) = new double[i];
		m.at(i - MINNP) = new double[i*i];
		s.at(i - MINNP) = new double[i*i];
	}
	ReadMatFiles();
}

#include <IMATLABC.h>
#include <string>
#include <iostream>
#include <mkl.h>
#include "Utils.h"
void DG3DNodal::SystemMatrixStandard1D::ReadMatFiles()
{
	for (size_t i = MINNP; i <= MAXNP; i++)
	{
		std::cout << "Reading Mat-File for Np = " << i << std::endl;
		std::string filename_prefix = "matrix1D_";
		std::string filename_Np = i2s(i);
		std::string filename_suffix = ".mat";
		std::string filename = filename_prefix + filename_Np + filename_suffix;
		IMATLABC::iMC_readMATFile(x[i - MINNP], filename.c_str(), "x");
		IMATLABC::iMC_readMATFile(m[i - MINNP], filename.c_str(), "M");
		IMATLABC::iMC_readMATFile(s[i - MINNP], filename.c_str(), "S");
	}
}

SystemMatrixStandard2D & DG3DNodal::SystemMatrixStandard2D::getHandle()
{
	if (handle_ == nullptr)
		handle_ = new SystemMatrixStandard2D();
	return *handle_;
}

const vector_t DG3DNodal::SystemMatrixStandard2D::getX(size_t Np) const
{
	return vector_t(getTotalBasis(Np), x.at(Np - MINNP));
}

const vector_t DG3DNodal::SystemMatrixStandard2D::getY(size_t Np) const
{
	return vector_t(getTotalBasis(Np), y.at(Np - MINNP));
}

const matrix_t DG3DNodal::SystemMatrixStandard2D::getM(size_t Np) const
{
	return matrix_t(getTotalBasis(Np), getTotalBasis(Np), m.at(Np - MINNP));
}

const matrix_t DG3DNodal::SystemMatrixStandard2D::getSx(size_t Np) const
{
	return matrix_t(getTotalBasis(Np), getTotalBasis(Np), sx.at(Np - MINNP));
}

const matrix_t DG3DNodal::SystemMatrixStandard2D::getSy(size_t Np) const
{
	return matrix_t(getTotalBasis(Np), getTotalBasis(Np), sy.at(Np - MINNP));
}

DG3DNodal::SystemMatrixStandard2D::SystemMatrixStandard2D()
{
	size_t total_order = MAXNP - MINNP + 1;
	x.resize(total_order);
	y.resize(total_order);
	m.resize(total_order);
	sx.resize(total_order);
	sy.resize(total_order);
	for (size_t i = MINNP; i <= MAXNP; i++)
	{
		size_t total_basis = getTotalBasis(i);
		size_t index_order = i - MINNP;
		x.at(index_order) = new double[total_basis];
		y.at(index_order) = new double[total_basis];
		m.at(index_order) = new double[total_basis*total_basis];
		sx.at(index_order) = new double[total_basis*total_basis];
		sy.at(index_order) = new double[total_basis*total_basis];
	}
	calculateMatrices();
}

#include <vector>

void DG3DNodal::SystemMatrixStandard2D::calculateMatrices()
{
	const DG3DNodal::SystemMatrixStandard1D& sm1d = DG3DNodal::SystemMatrixStandard1D::getHandle();
	for (size_t Np = MINNP; Np < MAXNP; Np++)
	{
		const vector_t x1d = sm1d.getX(Np);
		const matrix_t m1d = sm1d.getM(Np);
		const matrix_t s1d = sm1d.getS(Np);

		size_t total_basis = getTotalBasis(Np);
		size_t index_order = Np - MINNP;


		size_t* index_basis = new size_t[total_basis];
		for (size_t i = 0; i < total_basis; i++)
		{
			index_basis[i] = i;
		}
		size_t* ix = new size_t[total_basis];
		size_t* iy = new size_t[total_basis];
		DG3DNodal::ind2sub2(Np, total_basis, index_basis, ix, iy);

		vector_t x_ref(total_basis, x.at(index_order));
		vector_t y_ref(total_basis, y.at(index_order));
		matrix_t m_ref(total_basis, total_basis, m.at(index_order));
		matrix_t sx_ref(total_basis, total_basis, sx.at(index_order));
		matrix_t sy_ref(total_basis, total_basis, sy.at(index_order));
		for (size_t i = 0; i < total_basis; i++)
		{
			size_t ix_c = ix[i];
			size_t iy_c = iy[i];
			x_ref(i) = x1d(ix_c);
			y_ref(i) = x1d(iy_c);
			for (size_t j = 0; j < total_basis; j++)
			{
				size_t ix_cj = ix[j];
				size_t iy_cj = iy[j];
				m_ref(i, j) = m1d(ix_c, ix_cj) * m1d(iy_c, iy_cj);
				sx_ref(i, j) = s1d(ix_c, ix_cj) * m1d(iy_c, iy_cj);
				sy_ref(i, j) = m1d(ix_c, ix_cj) * s1d(iy_c, iy_cj);
			}
		}

		delete[] index_basis;
		delete[] ix;
		delete[] iy;
	}
}

SystemMatrixStandard3D & DG3DNodal::SystemMatrixStandard3D::getHandle()
{
	if (handle_ == nullptr)
		handle_ = new SystemMatrixStandard3D();
	return *handle_;
}

const vector_t DG3DNodal::SystemMatrixStandard3D::getX(size_t Np) const
{
	return vector_t(getTotalBasis(Np), x.at(Np - MINNP));
}

const vector_t DG3DNodal::SystemMatrixStandard3D::getY(size_t Np) const
{
	return vector_t(getTotalBasis(Np), y.at(Np - MINNP));
}

const vector_t DG3DNodal::SystemMatrixStandard3D::getZ(size_t Np) const
{
	return vector_t(getTotalBasis(Np), z.at(Np - MINNP));
}

const matrix_t DG3DNodal::SystemMatrixStandard3D::getM(size_t Np) const
{
	return matrix_t(getTotalBasis(Np), getTotalBasis(Np), m.at(Np - MINNP));
}

const matrix_t DG3DNodal::SystemMatrixStandard3D::getSx(size_t Np) const
{
	return matrix_t(getTotalBasis(Np), getTotalBasis(Np), sx.at(Np - MINNP));
}

const matrix_t DG3DNodal::SystemMatrixStandard3D::getSy(size_t Np) const
{
	return matrix_t(getTotalBasis(Np), getTotalBasis(Np), sy.at(Np - MINNP));
}

const matrix_t DG3DNodal::SystemMatrixStandard3D::getSz(size_t Np) const
{
	return matrix_t(getTotalBasis(Np), getTotalBasis(Np), sz.at(Np - MINNP));
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrixStandard3D::getLb(size_t Np) const
{
	return matrix_t(getTotalBasis(Np), getTotalBasis(Np), lb.at(Np - MINNP));
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrixStandard3D::getLf(size_t Np) const
{
	return matrix_t(getTotalBasis(Np), getTotalBasis(Np), lf.at(Np - MINNP));
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrixStandard3D::getLl(size_t Np) const
{
	return matrix_t(getTotalBasis(Np), getTotalBasis(Np), ll.at(Np - MINNP));
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrixStandard3D::getLr(size_t Np) const
{
	return matrix_t(getTotalBasis(Np), getTotalBasis(Np), lr.at(Np - MINNP));
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrixStandard3D::getLu(size_t Np) const
{
	return matrix_t(getTotalBasis(Np), getTotalBasis(Np), lu.at(Np - MINNP));
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrixStandard3D::getLt(size_t Np) const
{
	return matrix_t(getTotalBasis(Np), getTotalBasis(Np), lt.at(Np - MINNP));
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrixStandard3D::getBb(size_t Np) const
{
	return matrix_t(getTotalBasis2d(Np), getTotalBasis(Np), bb.at(Np - MINNP));
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrixStandard3D::getBf(size_t Np) const
{
	return matrix_t(getTotalBasis2d(Np), getTotalBasis(Np), bf.at(Np - MINNP));
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrixStandard3D::getBl(size_t Np) const
{
	return matrix_t(getTotalBasis2d(Np), getTotalBasis(Np), bl.at(Np - MINNP));
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrixStandard3D::getBr(size_t Np) const
{
	return matrix_t(getTotalBasis2d(Np), getTotalBasis(Np), br.at(Np - MINNP));
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrixStandard3D::getBu(size_t Np) const
{
	return matrix_t(getTotalBasis2d(Np), getTotalBasis(Np), bu.at(Np - MINNP));
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrixStandard3D::getBt(size_t Np) const
{
	return matrix_t(getTotalBasis2d(Np), getTotalBasis(Np), bt.at(Np - MINNP));
}

DG3DNodal::SystemMatrixStandard3D::SystemMatrixStandard3D()
{
	size_t total_order = MAXNP - MINNP + 1;
	x.resize(total_order);
	y.resize(total_order);
	z.resize(total_order);
	m.resize(total_order);
	sx.resize(total_order);
	sy.resize(total_order);
	sz.resize(total_order);

	lb.resize(total_order);
	lf.resize(total_order);
	ll.resize(total_order);
	lr.resize(total_order);
	lu.resize(total_order);
	lt.resize(total_order);

	bb.resize(total_order);
	bf.resize(total_order);
	bl.resize(total_order);
	br.resize(total_order);
	bu.resize(total_order);	
	bt.resize(total_order);
	for (size_t i = MINNP; i <= MAXNP; i++)
	{
		size_t total_basis = getTotalBasis(i);
		size_t index_order = i - MINNP;
		x.at(index_order) = new double[total_basis]();
		y.at(index_order) = new double[total_basis]();
		z.at(index_order) = new double[total_basis]();
		m.at(index_order) = new double[total_basis*total_basis]();
		sx.at(index_order) = new double[total_basis*total_basis]();
		sy.at(index_order) = new double[total_basis*total_basis]();
		sz.at(index_order) = new double[total_basis*total_basis]();
		lb.at(index_order) = new double[total_basis*total_basis]();
		lf.at(index_order) = new double[total_basis*total_basis]();
		ll.at(index_order) = new double[total_basis*total_basis]();
		lr.at(index_order) = new double[total_basis*total_basis]();
		lu.at(index_order) = new double[total_basis*total_basis]();
		lt.at(index_order) = new double[total_basis*total_basis]();
		bb.at(index_order) = new double[total_basis*total_basis]();
		bf.at(index_order) = new double[total_basis*total_basis]();
		bl.at(index_order) = new double[total_basis*total_basis]();
		br.at(index_order) = new double[total_basis*total_basis]();
		bu.at(index_order) = new double[total_basis*total_basis]();
		bt.at(index_order) = new double[total_basis*total_basis]();
		
	}
	calculateMatrices();
}

void DG3DNodal::SystemMatrixStandard3D::calculateMatrices()
{
	const DG3DNodal::SystemMatrixStandard1D& sm1d = DG3DNodal::SystemMatrixStandard1D::getHandle();
	const DG3DNodal::SystemMatrixStandard2D& sm2d = DG3DNodal::SystemMatrixStandard2D::getHandle();

	for (size_t Np = MINNP; Np < MAXNP; Np++)
	{
		const vector_t x1d = sm1d.getX(Np);
		const matrix_t m1d = sm1d.getM(Np);
		const matrix_t s1d = sm1d.getS(Np);

		const matrix_t m2d = sm2d.getM(Np);
		
		size_t total_basis = getTotalBasis(Np);
		size_t total_basis2 = sm2d.getTotalBasis(Np);
		size_t index_order = Np - MINNP;

		 
		size_t* index_basis = new size_t[total_basis];
		for (size_t i = 0; i < total_basis; i++)
		{
			index_basis[i] = i;
		}
		size_t* ix = new size_t[total_basis];
		size_t* iy = new size_t[total_basis];
		size_t* iz = new size_t[total_basis];
		DG3DNodal::ind2sub3(Np, total_basis, index_basis, ix, iy, iz);

		vector_t x_ref(total_basis, x.at(index_order));
		vector_t y_ref(total_basis, y.at(index_order));
		vector_t z_ref(total_basis, z.at(index_order));
		matrix_t m_ref(total_basis, total_basis, m.at(index_order));
		matrix_t sx_ref(total_basis, total_basis, sx.at(index_order));
		matrix_t sy_ref(total_basis, total_basis, sy.at(index_order));
		matrix_t sz_ref(total_basis, total_basis, sz.at(index_order));

		matrix_t lb_ref(total_basis, total_basis, lb.at(index_order));
		matrix_t lf_ref(total_basis, total_basis, lf.at(index_order));
		matrix_t ll_ref(total_basis, total_basis, ll.at(index_order));
		matrix_t lr_ref(total_basis, total_basis, lr .at(index_order));
		matrix_t lu_ref(total_basis, total_basis, lu.at(index_order));
		matrix_t lt_ref(total_basis, total_basis, lt.at(index_order));

		matrix_t bb_ref(total_basis2, total_basis, bb.at(index_order));
		matrix_t bf_ref(total_basis2, total_basis, bf.at(index_order));
		matrix_t bl_ref(total_basis2, total_basis, bl.at(index_order));
		matrix_t br_ref(total_basis2, total_basis, br.at(index_order));
		matrix_t bu_ref(total_basis2, total_basis, bu.at(index_order));
		matrix_t bt_ref(total_basis2, total_basis, bt.at(index_order));

		std::vector<size_t> ibb(total_basis);
		std::vector<size_t> ibf(total_basis);
		std::vector<size_t> ibl(total_basis);
		std::vector<size_t> ibr(total_basis);
		std::vector<size_t> ibu(total_basis);
		std::vector<size_t> ibt(total_basis);
		size_t ibbc = 1, ibfc = 1, iblc = 1, ibrc = 1, ibuc = 1, ibtc = 1;
		for (size_t i = 0; i < total_basis; ++i)
		{
			size_t ix_c = ix[i];
			size_t iy_c = iy[i];
			size_t iz_c = iz[i];
			if (ix_c == 0)
			{
				ibb.at(i) = ibbc;
				ibbc++;
			}
			if (ix_c == Np-1)
			{
				ibf.at(i) = ibfc;
				ibfc++;
			}
			if (iy_c == 0)
			{
				ibl.at(i) = iblc;
				iblc++;
			}
			if (iy_c == Np - 1)
			{
				ibr.at(i) = ibrc;
				ibrc++;
			}
			if (iz_c == 0)
			{
				ibu.at(i) = ibuc;
				ibuc++;
			}
			if (iz_c == Np - 1)
			{
				ibt.at(i) = ibtc;
				ibtc++;
			}
		}

		for (size_t i = 0; i < total_basis; ++i)
		{
			size_t ix_c = ix[i];
			size_t iy_c = iy[i];
			size_t iz_c = iz[i];
			x_ref(i) = x1d(ix_c);
			y_ref(i) = x1d(iy_c);
			z_ref(i) = x1d(iz_c);
			if (ibb.at(i) > 0)
				bb_ref(ibb.at(i) - 1, i) = 1.0;
			if (ibf.at(i) > 0)
				bf_ref(ibf.at(i) - 1, i) = 1.0;
			if (ibl.at(i) > 0)
				bl_ref(ibl.at(i) - 1, i) = 1.0;
			if (ibr.at(i) > 0)
				br_ref(ibr.at(i) - 1, i) = 1.0;
			if (ibu.at(i) > 0)
				bu_ref(ibu.at(i) - 1, i) = 1.0;
			if (ibt.at(i) > 0)
				bt_ref(ibt.at(i) - 1, i) = 1.0;
			for (size_t j = 0; j < total_basis; j++)
			{
				size_t ix_cj = ix[j];
				size_t iy_cj = iy[j];
				size_t iz_cj = iz[j];
				m_ref(i, j) = m1d(ix_c, ix_cj) * m1d(iy_c, iy_cj) * m1d(iz_c, iz_cj);
				sx_ref(i, j) = s1d(ix_c, ix_cj) * m1d(iy_c, iy_cj) * m1d(iz_c, iz_cj);
				sy_ref(i, j) = m1d(ix_c, ix_cj) * s1d(iy_c, iy_cj) * m1d(iz_c, iz_cj);
				sz_ref(i, j) = m1d(ix_c, ix_cj) * m1d(iy_c, iy_cj) * s1d(iz_c, iz_cj);
				if (ibb.at(i) > 0 && ibb.at(j) > 0)		
					lb_ref(i, j) = m2d(ibb.at(i) - 1, ibb.at(j) - 1);
				if (ibf.at(i) > 0 && ibf.at(j) > 0)
					lf_ref(i, j) = m2d(ibf.at(i) - 1, ibf.at(j) - 1);
				if (ibl.at(i) > 0 && ibl.at(j) > 0)
					ll_ref(i, j) = m2d(ibl.at(i) - 1, ibl.at(j) - 1);
				if (ibr.at(i) > 0 && ibr.at(j) > 0)
					lr_ref(i, j) = m2d(ibr.at(i) - 1, ibr.at(j) - 1);
				if (ibu.at(i) > 0 && ibu.at(j) > 0)
					lu_ref(i, j) = m2d(ibu.at(i) - 1, ibu.at(j) - 1);
				if (ibt.at(i) > 0 && ibt.at(j) > 0)
					lt_ref(i, j) = m2d(ibt.at(i) - 1, ibt.at(j) - 1);								
			}
		}

		delete[] index_basis;
		delete[] ix;
		delete[] iy;
		delete[] iz;
		
	}
}

//const BoundaryMask1D & DG3DNodal::BoundaryMask1D::getHandle()
//{
//	if (handle_ == nullptr)
//		handle_ = new BoundaryMask1D();
//	return *handle_;
//}
//
//const std::vector<bool> DG3DNodal::BoundaryMask1D::getB(size_t Np) const
//{
//	return is_boundary_b.at(Np - MINNP);
//}
//
//const std::vector<bool> DG3DNodal::BoundaryMask1D::getF(size_t Np) const
//{
//	return is_boundary_f.at(Np - MINNP);
//}
//
//DG3DNodal::BoundaryMask1D::BoundaryMask1D()
//{
//	is_boundary_b.resize(MAXNP - MINNP + 1);
//	is_boundary_f.resize(MAXNP - MINNP + 1);
//	for (size_t Np = MINNP; Np < MAXNP; Np++)
//	{
//		is_boundary_b.at(Np - MINNP).resize(getTotalBasis(Np));	
//		is_boundary_f.at(Np - MINNP).resize(getTotalBasis(Np));
//	}
//	calculate();
//}
//
//void DG3DNodal::BoundaryMask1D::calculate()
//{
//	for (size_t Np = MINNP; Np < MAXNP; Np++)
//	{
//		for (size_t i = 0; i < getTotalBasis(Np); ++i)
//		{
//			if (i == 0)
//				is_boundary_b.at(Np).at(i) = true;
//			else
//				is_boundary_b.at(Np).at(i) = false;
//			if (i == Np - 1)
//				is_boundary_f.at(Np).at(i) = true;
//			else
//				is_boundary_f.at(Np).at(i) = false;			
//		}		
//	}
//}
//
//const BoundaryMask2D & DG3DNodal::BoundaryMask2D::getHandle()
//{
//	if (handle_ == nullptr)
//		handle_ = new BoundaryMask2D();
//	return *handle_;
//}
//
//const std::vector<bool> DG3DNodal::BoundaryMask2D::getB(size_t Np) const
//{
//	return is_boundary_b.at(Np - MINNP);
//}
//
//const std::vector<bool> DG3DNodal::BoundaryMask2D::getF(size_t Np) const
//{
//	return is_boundary_f.at(Np - MINNP);
//}
//
//const std::vector<bool> DG3DNodal::BoundaryMask2D::getL(size_t Np) const
//{
//	return is_boundary_l.at(Np - MINNP);
//}
//
//const std::vector<bool> DG3DNodal::BoundaryMask2D::getR(size_t Np) const
//{
//	return is_boundary_r.at(Np - MINNP);
//}
//
//DG3DNodal::BoundaryMask2D::BoundaryMask2D()
//{
//	is_boundary_b.resize(MAXNP - MINNP + 1);
//	is_boundary_f.resize(MAXNP - MINNP + 1);
//	is_boundary_l.resize(MAXNP - MINNP + 1);
//	is_boundary_r.resize(MAXNP - MINNP + 1);
//	for (size_t Np = MINNP; Np < MAXNP; Np++)
//	{
//		is_boundary_b.at(Np - MINNP).resize(getTotalBasis(Np));
//		is_boundary_f.at(Np - MINNP).resize(getTotalBasis(Np));
//		is_boundary_l.at(Np - MINNP).resize(getTotalBasis(Np));
//		is_boundary_r.at(Np - MINNP).resize(getTotalBasis(Np));
//	}
//	calculate();
//}
//
//void DG3DNodal::BoundaryMask2D::calculate()
//{
//	for (size_t Np = MINNP; Np < MAXNP; Np++)
//	{
//		for (size_t i = 0; i < getTotalBasis(Np); ++i)
//		{
//			size_t ix, iy;
//			ind2sub2(Np, i, ix, iy);
//			if (ix == 0)
//				is_boundary_b.at(Np).at(i) = true;
//			else
//				is_boundary_b.at(Np).at(i) = false;
//			if (ix == Np-1)
//				is_boundary_f.at(Np).at(i) = true;
//			else
//				is_boundary_f.at(Np).at(i) = false;
//			if (iy == 0)
//				is_boundary_l.at(Np).at(i) = true;
//			else
//				is_boundary_l.at(Np).at(i) = false;
//			if (iy == Np - 1)
//				is_boundary_r.at(Np).at(i) = true;
//			else
//				is_boundary_r.at(Np).at(i) = false;
//		}
//	}
//}
//
//const BoundaryMask3D & DG3DNodal::BoundaryMask3D::getHandle()
//{
//	if (handle_ == nullptr)
//		handle_ = new BoundaryMask3D();
//	return *handle_;
//}
//
////const std::vector<bool> DG3DNodal::BoundaryMask3D::getB(size_t Np) const
////{
////	return is_boundary_b.at(Np - MINNP);
////}
////
////const std::vector<bool> DG3DNodal::BoundaryMask3D::getF(size_t Np) const
////{
////	return is_boundary_f.at(Np - MINNP);
////}
////
////const std::vector<bool> DG3DNodal::BoundaryMask3D::getL(size_t Np) const
////{
////	return is_boundary_l.at(Np - MINNP);
////}
////
////const std::vector<bool> DG3DNodal::BoundaryMask3D::getR(size_t Np) const
////{
////	return is_boundary_r.at(Np - MINNP);
////}
////
////const std::vector<bool> DG3DNodal::BoundaryMask3D::getU(size_t Np) const
////{
////	return is_boundary_u.at(Np - MINNP);
////}
////
////const std::vector<bool> DG3DNodal::BoundaryMask3D::getT(size_t Np) const
////{
////	return is_boundary_t.at(Np - MINNP);
////}
//
//DG3DNodal::BoundaryMask3D::BoundaryMask3D()
//{
//	b.resize(MAXNP - MINNP + 1);
//	f.resize(MAXNP - MINNP + 1);
//	l.resize(MAXNP - MINNP + 1);
//	r.resize(MAXNP - MINNP + 1);
//	u.resize(MAXNP - MINNP + 1);
//	t.resize(MAXNP - MINNP + 1);
//	for (size_t Np = MINNP; Np < MAXNP; Np++)
//	{
//		b.at(Np - MINNP)= new size_t[getTotalBasis(Np)];
//		f.at(Np - MINNP)= new size_t[getTotalBasis(Np)];
//		l.at(Np - MINNP)= new size_t[getTotalBasis(Np)];
//		r.at(Np - MINNP)= new size_t[getTotalBasis(Np)];
//		u.at(Np - MINNP)= new size_t[getTotalBasis(Np)];
//		t.at(Np - MINNP)= new size_t[getTotalBasis(Np)];
//	}
//	calculate();
//}
//
//void DG3DNodal::BoundaryMask3D::calculate()
//{
//	for (size_t Np = MINNP; Np < MAXNP; Np++)
//	{
//		for (size_t i = 0; i < getTotalBasis(Np); ++i)
//		{
//			//size_t ix, iy, iz;
//			//ind2sub3(Np, i, ix, iy, iz);
//			//if (ix == 0)
//			//	b.at(Np).at(i) = true;
//			//else
//			//	b.at(Np).at(i) = false;
//			//if (ix == Np - 1)
//			//	f.at(Np).at(i) = true;
//			//else
//			//	f.at(Np).at(i) = false;
//			//if (iy == 0)
//			//	l.at(Np).at(i) = true;
//			//else
//			//	l.at(Np).at(i) = false;
//			//if (iy == Np - 1)
//			//	r.at(Np).at(i) = true;
//			//else
//			//	r.at(Np).at(i) = false;
//			//if (iz == 0)
//			//	u.at(Np).at(i) = true;
//			//else
//			//	u.at(Np).at(i) = false;
//			//if (iz == Np - 1)
//			//	t.at(Np).at(i) = true;
//			//else
//			//	t.at(Np).at(i) = false;
//		}
//	}
//}
//
//
