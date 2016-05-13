#include "SystemMatrices.h"
using namespace DG3DNodal;

#include "Matrices.h"

DG3DNodal::SystemMatrix3D::SystemMatrix3D(num_t h, num_t mu, num_t xi, num_t eta, num_t sigma, size_t Np)
	: h_(h), mu_(mu), xi_(xi), eta_(eta), sigma_(sigma), np_(Np)
{
	const SystemMatrixStandard3D& sm3 = SystemMatrixStandard3D::getHandle();
	size_t total_basis = sm3.getTotalBasis(np_);
	size_t total_basis2 = sm3.getTotalBasis2d(np_);

	xp = calloc(total_basis);
	yp = calloc(total_basis);
	zp = calloc(total_basis);
	sysmp = calloc(total_basis*total_basis);
	mxp = calloc(total_basis*total_basis2);
	myp = calloc(total_basis*total_basis2);
	mzp = calloc(total_basis*total_basis2);
	fxp = calloc(total_basis2*total_basis);
	fyp = calloc(total_basis2*total_basis);
	fzp = calloc(total_basis2*total_basis);

	sysmlp = calloc(total_basis*total_basis);
	sysmup = calloc(total_basis*total_basis);
	sysmip = calloc(total_basis*total_basis);

	x.bind(total_basis, xp);
	y.bind(total_basis, yp);
	z.bind(total_basis, zp);
	sysm.bind(total_basis, total_basis, sysmp);
	mx.bind(total_basis, total_basis2, mxp);
	my.bind(total_basis, total_basis2, myp);
	mz.bind(total_basis, total_basis2, mzp);
	fx.bind(total_basis2, total_basis, fxp);
	fy.bind(total_basis2, total_basis, fyp);
	fz.bind(total_basis2, total_basis, fzp);

	sysml.bind(total_basis, total_basis, sysmlp);
	sysmu.bind(total_basis, total_basis, sysmup);
	sysmi.bind(total_basis, total_basis, sysmip);

	ipiv = (lapack_int*)mkl_calloc(total_basis, sizeof(int), 64);	
	
	calculate();
}

DG3DNodal::SystemMatrix3D::~SystemMatrix3D()
{
	mkl_free(xp);
	mkl_free(yp);
	mkl_free(zp);
	mkl_free(sysmp);
	mkl_free(mxp);
	mkl_free(myp);
	mkl_free(mzp);
	mkl_free(fxp);
	mkl_free(fyp);
	mkl_free(fzp);

	mkl_free(sysmlp);
	mkl_free(sysmup);
	mkl_free(sysmip);
}

const DG3DNodal::SystemMatrix3D& DG3DNodal::SystemMatrix3D::getHandle() const
{
	return *this;
}

const DG3DNodal::vector_t DG3DNodal::SystemMatrix3D::getX() const
{
	return vector_t(x);
}

const DG3DNodal::vector_t DG3DNodal::SystemMatrix3D::getY() const
{
	return vector_t(y);
}

const DG3DNodal::vector_t DG3DNodal::SystemMatrix3D::getZ() const
{
	return vector_t(z);
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrix3D::getSysM() const
{
	return matrix_t(sysm);
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrix3D::getMx() const
{
	return matrix_t(mx);
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrix3D::getMy() const
{
	return matrix_t(my);
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrix3D::getMz() const
{
	return matrix_t(mz);
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrix3D::getFx() const
{
	return matrix_t(fx);
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrix3D::getFy() const
{
	return matrix_t(fy);
}

const DG3DNodal::matrix_t DG3DNodal::SystemMatrix3D::getFz() const
{
	return matrix_t(fz);
}

size_t DG3DNodal::SystemMatrix3D::getRefineLevel() const
{
	return size_t(refine_level);
}

void DG3DNodal::SystemMatrix3D::calculate()
{
	const DG3DNodal::SystemMatrixStandard3D& sm3 = SystemMatrixStandard3D::getHandle();
	const_num_p m_s3d = sm3.getM(np_).getPtr();
	const_num_p sx_s3d = sm3.getSx(np_).getPtr();
	const_num_p sy_s3d = sm3.getSy(np_).getPtr();
	const_num_p sz_s3d = sm3.getSz(np_).getPtr();
	const_num_p lb_s3d, lf_s3d, ll_s3d, lr_s3d, lu_s3d, lt_s3d;
	const_num_p bb_s3d, bf_s3d, bl_s3d, br_s3d, bu_s3d, bt_s3d;

	lb_s3d = sm3.getLb(np_).getPtr();
	lf_s3d = sm3.getLf(np_).getPtr();
	ll_s3d = sm3.getLl(np_).getPtr();
	lr_s3d = sm3.getLr(np_).getPtr();
	lu_s3d = sm3.getLu(np_).getPtr();
	lt_s3d = sm3.getLt(np_).getPtr();

	bb_s3d = sm3.getBb(np_).getPtr();
	bf_s3d = sm3.getBf(np_).getPtr();
	bl_s3d = sm3.getBl(np_).getPtr();
	br_s3d = sm3.getBr(np_).getPtr();
	bu_s3d = sm3.getBu(np_).getPtr();
	bt_s3d = sm3.getBt(np_).getPtr();	

	size_t nb = getTotalBasis();
	size_t ni = sm3.getTotalBasis2d(np_);
	num_t c3 = h_*h_*h_ / 8.0;
	num_t c2 = h_*h_ / 4.0;
	//sysm = sigma * c3 * M 
	mkl_domatadd('C', 'N', 'N', nb, nb, 0.0, sysmp, nb, sigma_*c3,	m_s3d, nb, sysmp, nb);

	//sysm = sysm + c2*mu*Sx
	mkl_domatadd('C', 'N', 'N', nb, nb, 1.0, sysmp, nb, mu_*c2,	sx_s3d, nb, sysmp, nb);

	//sysm = sysm + c2*xi*Sy
	mkl_domatadd('C', 'N', 'N', nb, nb, 1.0, sysmp, nb, xi_*c2,	sy_s3d, nb, sysmp, nb);

	//sysm = sysm + c2*eta*Sz
	mkl_domatadd('C', 'N', 'N', nb, nb, 1.0, sysmp, nb, eta_*c2,	sz_s3d, nb, sysmp, nb);

	/** mx, my, mz: nb-by-ni matrix
			calculate income flux by m(xyz) * data_interface(xyz)
			mx = mu * lb * bb' /  mu * lf * bf'
	*/
	
	const_num_p lx, ly, lz;
	if (mu_ >= 0.0)
	{
		//lx = bb
		lx = lb_s3d;
		//fxp = c2 * bf
		mkl_domatadd('C', 'N', 'N', ni, nb, 0.0, fxp, ni, 1.0, bf_s3d, ni, fxp, ni);

		//mxp =  mu * c2 * lx * bb
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, nb, ni, nb, mu_*c2, lx, nb, bb_s3d, ni, 0.0, mxp, nb);
	}
	else
	{
		lx = lf_s3d;
		mkl_domatadd('C', 'N', 'N', ni, nb, 0.0, fxp, ni, 1.0, bb_s3d, ni, fxp, ni);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, nb, ni, nb, mu_*c2, lx, nb, bf_s3d, ni, 0.0, mxp, nb);
	}

	if (xi_ >= 0.0)
	{
		ly = ll_s3d;
		mkl_domatadd('C', 'N', 'N', ni, nb, 0.0, fyp, ni, 1.0, bl_s3d, ni, fyp, ni);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, nb, ni, nb, xi_*c2, ly, nb, bl_s3d, ni, 0.0, myp, nb);
	}
	else
	{
		ly = lr_s3d;
		mkl_domatadd('C', 'N', 'N', ni, nb, 0.0, fyp, ni, 1.0, bl_s3d, ni, fyp, ni);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, nb, ni, nb, xi_*c2, ly, nb, br_s3d, ni, 0.0, myp, nb);
	}

	if (eta_ >= 0.0)
	{		
		lz = lu_s3d;		
		mkl_domatadd('C', 'N', 'N', ni, nb, 0.0, fzp, ni, 1.0, bt_s3d, ni, fzp, ni);		
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, nb, ni, nb, eta_*c2, lz, nb, bu_s3d, ni, 0.0, mzp, nb);
	}
	else
	{
		lz = lt_s3d;
		mkl_domatadd('C', 'N', 'N', ni, nb, 0.0, fzp, ni, 1.0, bu_s3d, ni, fzp, ni);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, nb, ni, nb, eta_*c2, lz, nb, bt_s3d, ni, 0.0, mzp, nb);
	}

	mkl_domatadd('C', 'N', 'N', nb, nb, 1.0, sysmp, nb, c2*mu_, lx, nb, sysmp, nb);
	mkl_domatadd('C', 'N', 'N', nb, nb, 1.0, sysmp, nb, c2*xi_, ly, nb, sysmp, nb);
	mkl_domatadd('C', 'N', 'N', nb, nb, 1.0, sysmp, nb, c2*eta_, lz, nb, sysmp, nb);
	
	cblas_dcopy(nb*nb, sysm.getPtr(), 1, sysml.getPtr(), 1);
	LAPACKE_dgetrf(LAPACK_COL_MAJOR, nb, nb, sysml.getPtr(), nb, ipiv);
	cblas_dcopy(nb*nb, sysml.getPtr(), 1, sysmi.getPtr(), 1);
	LAPACKE_dgetri(LAPACK_COL_MAJOR, nb, sysmi.getPtr(), nb, ipiv);
}
