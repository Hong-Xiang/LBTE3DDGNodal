#include <stdexcept>
#include "system_matrices.h"
#include "utilities.h"

using namespace dgn;
//========	Implement of matrices_data	========================================

dgn::matrices_data* dgn::matrices_data::handle_ = nullptr;

const dgn::matrices_data & dgn::matrices_data::instance()
{
	if (handle_ == nullptr)
		handle_ = new dgn::matrices_data();
	return *handle_;
}

const vector_t dgn::matrices_data::x(size_t Np) const
{
	if (Np < MINNP || Np > MAXNP)
		throw(std::out_of_range("Np = " + xlib::i2s(Np)));
	return vector_t(Np, x_.at(Np - MINNP));
}

const matrix_t dgn::matrices_data::m(size_t Np) const
{
	if (Np < MINNP || Np > MAXNP)
		throw(std::out_of_range("Np = " + xlib::i2s(Np)));
	return matrix_t(Np, Np, m_.at(Np - MINNP));
}

const matrix_t dgn::matrices_data::s(size_t Np) const
{
	if (Np < MINNP || Np > MAXNP)
		throw(std::out_of_range("Np = " + xlib::i2s(Np)));
	return matrix_t(Np, Np, s_.at(Np - MINNP));
}

dgn::matrices_data::matrices_data()
{
	x_.resize(MAXNP - MINNP + 1);
	m_.resize(MAXNP - MINNP + 1);
	s_.resize(MAXNP - MINNP + 1);
	for (size_t i = MINNP; i <= MAXNP; i++)
	{
#ifdef _DEBUG
		std::cout << "Reading matrices files " << i << std::endl;
#endif
		x_.at(i - MINNP) = xlib::mkl_ext::xcalloc<num_t>(i);
		m_.at(i - MINNP) = xlib::mkl_ext::xcalloc<num_t>(i*i);
		s_.at(i - MINNP) = xlib::mkl_ext::xcalloc<num_t>(i*i);
	}
	read_MAT_files();
}

void dgn::matrices_data::read_MAT_files()
{
	for (size_t i = MINNP; i <= MAXNP; i++)
	{
//		std::cout << "Reading Mat-File for Np = " << i << std::endl;
		const char* filename = "matrices.mat";
		double* xtmp = new double[i];
		double* mtmp = new double[i*i];
		double* stmp = new double[i*i];
		xlib::imatlab::load_MAT_cell_data(xtmp, filename, "x", i-1);
		xlib::imatlab::load_MAT_cell_data(mtmp, filename, "M", i-1);
		xlib::imatlab::load_MAT_cell_data(stmp, filename, "S", i-1);
		for (size_t j = 0; j < i; j++)
			x_.at(i - MINNP)[j] = xtmp[j];
		for (size_t j = 0; j < i*i; j++)
		{
			m_.at(i - MINNP)[j] = mtmp[j];
			s_.at(i - MINNP)[j] = stmp[j];
		}
		delete[] xtmp;
		delete[] mtmp;
		delete[] stmp;
	}
}

//==============================================================================




//========	Implement of matrices_base	========================================

dgn::matrices_base::matrices_base(num_t h, size_t np)
	: h_(h), np_(np), is_init_(false), dim_(dimension())
{
}

void dgn::matrices_base::reset(num_t h, size_t np)
{
	h_ = h;
	np_ = np;
	initialize();
}

void dgn::matrices_base::initialize()
{
	dim_ = dimension();
	memory_alloc();
	calculate();
	is_init_ = true;
}

dgn::matrices_base::~matrices_base()
{
	for (size_t i = 0; i < INTERFACE_TOTAL.at(dim_); i++)
	{
		mkl_free(lift_matrices_ptr_.at(i));
		mkl_free(flux_matrices_ptr_.at(i));		
		lift_matrices_.at(i).unbind();
		flux_matrices_.at(i).unbind();
	}
}

void dgn::matrices_base::memory_alloc()
{
	size_t ift = INTERFACE_TOTAL.at(dim_);
	lift_matrices_ptr_.resize(ift);
	flux_matrices_ptr_.resize(ift);
	lift_matrices_.resize(ift);
	flux_matrices_.resize(ift);
	for (size_t i = 0; i < ift; i++)
	{
		if (lift_matrices_ptr_.at(i) != nullptr)
			mkl_free(lift_matrices_ptr_.at(i));
		lift_matrices_ptr_.at(i) = xlib::mkl_ext::xcalloc<num_t>(basis_total()*basis_total());
		lift_matrices_.at(i).bind(basis_total(), basis_total(), lift_matrices_ptr_.at(i));
		if (flux_matrices_ptr_.at(i) != nullptr)
			mkl_free(flux_matrices_ptr_.at(i));
		flux_matrices_ptr_.at(i) = xlib::mkl_ext::xcalloc<num_t>(basis_lower_dim_total()*basis_total());
		flux_matrices_.at(i).bind(basis_lower_dim_total(), basis_total(), flux_matrices_ptr_.at(i));
	}
}

//==============================================================================

//========	Implement of matrices1d	============================================
dgn::matrices1d::matrices1d(num_t h, size_t np)
	: matrices_base(h, np), xp_(nullptr), mp_(nullptr), sp_(nullptr)
{
}

dgn::matrices1d::~matrices1d()
{
	if(xp_ != nullptr)
		mkl_free(xp_);
	if(mp_ != nullptr)
		mkl_free(mp_);
	if(sp_ != nullptr)
		mkl_free(sp_);
}

void dgn::matrices1d::calculate()
{
	const matrices_data& md = matrices_data::instance();
	vector_t xmd = md.x(np_);
	matrix_t mmd = md.m(np_);
	matrix_t smd = md.s(np_);

	for (size_t i = 0; i < basis_total(); i++)
	{
		x_(i) = xmd(i) * h_ / 2.0;
		for (size_t j = 0; j < basis_total(); j++)
		{
			m_(i, j) = mmd(i, j)*h_ / 2.0;
			s_(i, j) = smd(i, j);
		}
	}

	
	lift_matrices_.at(static_cast<size_t>(interface_direction::B))(1, 1) = 1.0;
	lift_matrices_.at(static_cast<size_t>(interface_direction::F))(np_, np_) = 1.0;

}

void dgn::matrices1d::memory_alloc()
{
	matrices_base::memory_alloc();
	if (xp_ != nullptr)
		mkl_free(xp_);
	xp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total());
	if (mp_ != nullptr)
		mkl_free(mp_);
	mp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total()*basis_total());
	if (sp_ != nullptr)
		mkl_free(sp_);
	sp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total()*basis_total());
	x_.bind(basis_total(), xp_);
	m_.bind(basis_total(), basis_total(), mp_);
	s_.bind(basis_total(), basis_total(), sp_);
}

//==============================================================================







//========	Implement of matrices2d	============================================

dgn::matrices2d::matrices2d(num_t h, size_t np)
	: matrices_base(h, np), 
	  xp_(nullptr), yp_(nullptr), 
	  mp_(nullptr), sxp_(nullptr), syp_(nullptr)
{
}

dgn::matrices2d::~matrices2d()
{
	mkl_free(xp_);
	mkl_free(yp_);
	mkl_free(mp_);
	mkl_free(sxp_);
	mkl_free(syp_);
}

void dgn::matrices2d::calculate()
{
	const matrices_data& md = matrices_data::instance();
	const vector_t xmd = md.x(np_);
	const matrix_t mmd = md.m(np_);
	const matrix_t smd = md.s(np_);

	for (size_t i = 0; i < basis_total(); i++)
	{
		size_t ix, iy;
		xlib::ind2sub(np_, np_, i, ix, iy);
		x_(i) = xmd(ix) * h_ / 2.0;
		y_(i) = xmd(iy) * h_ / 2.0;
		for (size_t j = 0; j < basis_total(); j++)
		{
			size_t jx, jy;
			xlib::ind2sub(np_, np_, j, jx, jy);
			m_(i, j) = mmd(ix, jx)*mmd(iy, jy) * h_ * h_ / 4.0;
			sx_(i, j) = smd(ix, jx)*mmd(iy, jy) * h_ / 2.0;
			sy_(i, j) = mmd(ix, jx)*smd(iy, jy) * h_ / 2.0;
		}
	}

	matrices1d m1d(h_, np_);
	m1d.initialize();
	for (size_t i = 0; i < INTERFACE_TOTAL.at(dimension()); i++)
	{
		interface_direction idir = interface_direction_list.at(i);
		std::vector<size_t> idx;
		std::vector<bool> bdmask;
		size_t cid = 0;
		idx.resize(basis_total());
		bdmask.resize(basis_total());
		for (size_t j = 0; j < basis_total(); j++)
		{
			size_t ix, iy;
			xlib::ind2sub(np_, np_, j, ix, iy);
			bdmask.at(j) = dgn::utilities::is_boundary(idir, np_, ix, iy);
			if (bdmask.at(j))
			{
				idx.at(j) = cid;
				cid++;
			}
		}
		for (size_t j = 0; j < basis_total(); j++)
		{
			if (bdmask.at(j))
				flux_matrices_.at(i)(idx.at(j), j) = 1.0;
			for (size_t k = 0; k < basis_total(); k++)
				if (bdmask.at(j) && bdmask.at(k))
					lift_matrices_.at(i)(j, k) = m1d.m()(idx.at(j), idx.at(k));
		}
	}
}

void dgn::matrices2d::memory_alloc()
{
	matrices_base::memory_alloc();
	if (xp_ != nullptr)
		mkl_free(xp_);
	xp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total());
	if (yp_ != nullptr)
		mkl_free(yp_);
	yp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total());
	if (mp_ != nullptr)
		mkl_free(mp_);
	mp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total()*basis_total());
	if (sxp_ != nullptr)
		mkl_free(sxp_);
	sxp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total()*basis_total());
	if (syp_ != nullptr)
		mkl_free(syp_);
	syp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total()*basis_total());

	x_.bind(basis_total(), xp_);
	y_.bind(basis_total(), yp_);
	m_.bind(basis_total(), basis_total(), mp_);
	sx_.bind(basis_total(), basis_total(), sxp_);
	sy_.bind(basis_total(), basis_total(), syp_);
}

//==============================================================================


//========	Implement of matrices3d	============================================

dgn::matrices3d::matrices3d(num_t h, size_t np)
	: matrices_base(h, np),
	  xp_(nullptr), yp_(nullptr), zp_(nullptr),
	  mp_(nullptr), sxp_(nullptr), syp_(nullptr), szp_(nullptr)
{
}

dgn::matrices3d::~matrices3d()
{
	mkl_free(xp_);
	mkl_free(yp_);
	mkl_free(zp_);
	mkl_free(mp_);
	mkl_free(sxp_);
	mkl_free(syp_);
	mkl_free(szp_);
}


void dgn::matrices3d::calculate()
{
	const matrices_data& md = matrices_data::instance();
	const vector_t xmd = md.x(np_);
	const matrix_t mmd = md.m(np_);
	const matrix_t smd = md.s(np_);

	for (size_t i = 0; i < basis_total(); i++)
	{
		size_t ix, iy, iz;
		xlib::ind2sub(np_, np_, np_, i, ix, iy, iz);
		x_(i) = xmd(ix) * h_ / 2.0;
		y_(i) = xmd(iy) * h_ / 2.0;
		z_(i) = xmd(iz) * h_ / 2.0;
		for (size_t j = 0; j < basis_total(); j++)
		{
			size_t jx, jy, jz;
			xlib::ind2sub(np_, np_, np_, j, jx, jy, jz);
			m_(i, j) = mmd(ix, jx) * mmd(iy, jy) * mmd(iz, jz) * h_ * h_ * h_ / 8.0;
			sx_(i, j) = smd(ix, jx) * mmd(iy, jy) * mmd(iz, jz) * h_ * h_ / 4.0;
			sy_(i, j) = mmd(ix, jx) * smd(iy, jy) * mmd(iz, jz) * h_ * h_ / 4.0;
			sz_(i, j) = mmd(ix, jx) * mmd(iy, jy) * smd(iz, jz) * h_ * h_ / 4.0;
		}
	}

	matrices2d m2d(h_, np_);
	m2d.initialize();
	for (size_t i = 0; i < INTERFACE_TOTAL.at(dimension()); i++)
	{
		interface_direction idir = interface_direction_list.at(i);
		std::vector<size_t> idx;
		std::vector<bool> bdmask;
		size_t cid = 0;
		idx.resize(basis_total());
		bdmask.resize(basis_total());
		for (size_t j = 0; j < basis_total(); j++)
		{
			size_t ix, iy, iz;
			xlib::ind2sub(np_, np_, np_, j, ix, iy, iz);
			bdmask.at(j) = dgn::utilities::is_boundary(idir, np_, ix, iy, iz);
			if (bdmask.at(j))
			{
				idx.at(j) = cid;
				cid++;
			}
		}
		for (size_t j = 0; j < basis_total(); j++)
		{
			if (bdmask.at(j))
				flux_matrices_.at(i)(idx.at(j), j) = 1.0;
			for (size_t k = 0; k < basis_total(); k++)			
				if (bdmask.at(j) && bdmask.at(k))
					lift_matrices_.at(i)(j, k) = m2d.m()(idx.at(j), idx.at(k));			
		}
	}
}

void dgn::matrices3d::memory_alloc()
{
	matrices_base::memory_alloc();
	if (xp_ != nullptr)
		mkl_free(xp_);
	xp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total());
	if (yp_ != nullptr)
		mkl_free(yp_);
	yp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total());
	if (zp_ != nullptr)
		mkl_free(zp_);
	zp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total());
	if (mp_ != nullptr)
		mkl_free(mp_);
	mp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total()*basis_total());
	if (sxp_ != nullptr)
		mkl_free(sxp_);
	sxp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total()*basis_total());
	if (syp_ != nullptr)
		mkl_free(syp_);
	syp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total()*basis_total());
	if (szp_ != nullptr)
		mkl_free(szp_);
	szp_ = xlib::mkl_ext::xcalloc<num_t>(basis_total()*basis_total());
	x_.bind(basis_total(), xp_);
	y_.bind(basis_total(), yp_);
	z_.bind(basis_total(), zp_);
	m_.bind(basis_total(), basis_total(), mp_);
	sx_.bind(basis_total(), basis_total(), sxp_);
	sy_.bind(basis_total(), basis_total(), syp_);
	sz_.bind(basis_total(), basis_total(), szp_);
}

//==============================================================================





//========	Implement of system_matrix_angle	================================

/*
dgn::system_matrix_angle::system_matrix_angle(num_t h, size_t np, num_t mu, num_t xi, num_t eta, num_t sigma)
	: matrices_(h, np), surface_matrices_(h, np), h_(h), np_(np), mu_(mu), xi_(xi), eta_(eta), sigma_(sigma)
{
	matrices_.initialize();
	surface_matrices_.initialize();
	size_t tbs = matrices_.basis_total();

	sysmp_ = xlib::mkl_ext::xcalloc<num_t>(tbs*tbs);
	sysm_.bind(tbs, tbs, sysmp_);
	sysmlup_ = xlib::mkl_ext::xcalloc<num_t>(tbs*tbs);
	sysmlu_.bind(tbs, tbs, sysmlup_);
	sysmivp_ = xlib::mkl_ext::xcalloc<num_t>(tbs*tbs);
	sysmiv_.bind(tbs, tbs, sysmivp_);
	
	massmp_ = xlib::mkl_ext::xcalloc<num_t>(tbs*tbs);
	massm_.bind(tbs, tbs, massmp_);

	massm_.omatcopy(matrices_.m());

	size_t tis = INTERFACE_TOTAL.at(3);
	lift_matrix_p_.resize(tis);
	flux_matrix_p_.resize(tis);
	lift_matrix_.resize(tis);
	flux_matrix_.resize(tis);

	size_t tbls = matrices_.basis_lower_dim_total();
	for (size_t i = 0; i < tis; i++)
	{
		interface_direction dir = interface_direction_list.at(i);
		lift_matrix_p_.at(i) = xlib::mkl_ext::xcalloc<num_t>(tbs*tbls);
		flux_matrix_p_.at(i) = xlib::mkl_ext::xcalloc<num_t>(tbls*tbs);
		lift_matrix_.at(i).bind(tbs, tbls, lift_matrix_p_.at(i));
		lift_matrix_.at(i).gemm(matrices_.lift_matrix(dir), matrices_.flux_matrix(dir), 1.0, 0.0, CblasNoTrans, CblasTrans);
		
		flux_matrix_.at(i).bind(tbls, tbs, flux_matrix_p_.at(i));
		flux_matrix_.at(i).omatcopy(matrices_.flux_matrix(dir));
		flux_matrix_.at(i).imatcopy('T');
		
		if (dir == interface_direction::B || dir == interface_direction::F) {
			lift_matrix_.at(i).imatcopy('N', mu_);
			flux_matrix_.at(i).imatcopy('N', 1.0);
		}
		if (dir == interface_direction::L || dir == interface_direction::R) {
			lift_matrix_.at(i).imatcopy('N', xi_);
			flux_matrix_.at(i).imatcopy('N', 1.0);
		}
		if (dir == interface_direction::D || dir == interface_direction::U) {
			lift_matrix_.at(i).imatcopy('N', eta_);
			flux_matrix_.at(i).imatcopy('N', 1.0);
		}
		
	}
	
	sysm_.omatadd(matrices_.m(), sysm_, sigma_, 0.0);
	
	sysm_.omatadd(matrices_.sx(), sysm_, mu_, 1.0);
	
	sysm_.omatadd(matrices_.sy(), sysm_, xi_, 1.0);
	
	sysm_.omatadd(matrices_.sz(), sysm_, eta_, 1.0);
	

	num_p tmpp = xlib::mkl_ext::xcalloc<num_t>(tbs*tbs);
	matrix_t tmpm(tbs, tbs, tmpp);

	size_t idir;
	if (mu_ >= 0)
	{
		interface_direction dir = interface_direction::B;
		idir = static_cast<size_t>(dir);	
		tmpm.gemm(lift_matrix_.at(idir), matrices_.flux_matrix(dir));		
		sysm_.omatadd(tmpm, sysm_, 1.0, 1.0);
	}
	else
	{
		interface_direction dir = interface_direction::F;
		idir = static_cast<size_t>(dir);
		tmpm.gemm(lift_matrix_.at(idir), matrices_.flux_matrix(dir));
		sysm_.omatadd(tmpm, sysm_, 1.0, 1.0);
	}
	
	if (xi_ >= 0)
	{
		interface_direction dir = interface_direction::L;
		idir = static_cast<size_t>(dir);
		tmpm.gemm(lift_matrix_.at(idir), matrices_.flux_matrix(dir));
		sysm_.omatadd(tmpm, sysm_, 1.0, 1.0);
	}
	else
	{
		interface_direction dir = interface_direction::R;
		idir = static_cast<size_t>(dir);
		tmpm.gemm(lift_matrix_.at(idir), matrices_.flux_matrix(dir));
		sysm_.omatadd(tmpm, sysm_, 1.0, 1.0);
	}
	
	if (eta_ >= 0)
	{
		interface_direction dir = interface_direction::D;
		idir = static_cast<size_t>(dir);
		tmpm.gemm(lift_matrix_.at(idir), matrices_.flux_matrix(dir));
		sysm_.omatadd(tmpm, sysm_, 1.0, 1.0);
	}
	else
	{
		interface_direction dir = interface_direction::U;
		idir = static_cast<size_t>(dir);
		tmpm.gemm(lift_matrix_.at(idir), matrices_.flux_matrix(dir));
		sysm_.omatadd(tmpm, sysm_, 1.0, 1.0);
	}

	mkl_free(tmpp);
	
	ipiv = new MKL_INT[tbs];
	sysmlu_.lu_factorization(sysm_, ipiv);

	sysmiv_.getri(sysm_);
	

	

}
*/

dgn::system_matrix_angle::system_matrix_angle(num_t h, size_t np, num_t mu, num_t xi, num_t eta, num_t sigma)
	: matrices_(h, np), surface_matrices_(h, np), h_(h), np_(np), mu_(mu), xi_(xi), eta_(eta), sigma_(sigma)
{
	matrices_.initialize();
	surface_matrices_.initialize();
	size_t tbs = matrices_.basis_total();

	sysmp_ = xlib::mkl_ext::xcalloc<num_t>(tbs*tbs);
	sysm_.bind(tbs, tbs, sysmp_);
	sysmlup_ = xlib::mkl_ext::xcalloc<num_t>(tbs*tbs);
	sysmlu_.bind(tbs, tbs, sysmlup_);
	sysmivp_ = xlib::mkl_ext::xcalloc<num_t>(tbs*tbs);
	sysmiv_.bind(tbs, tbs, sysmivp_);

	massmp_ = xlib::mkl_ext::xcalloc<num_t>(tbs*tbs);
	massm_.bind(tbs, tbs, massmp_);

	massm_.omatcopy(matrices_.m());

	reorder_p_ = xlib::mkl_ext::xcalloc<num_t>(tbs*tbs);
	back_order_p_ = xlib::mkl_ext::xcalloc<num_t>(tbs*tbs);

	reorder_.bind(tbs, tbs, reorder_p_);
	back_order_.bind(tbs, tbs, back_order_p_);


	std::vector<size_t> pos;
	pos.clear();
	pos.resize(tbs);
	for (size_t i = 0; i < tbs; i++)
	{
		size_t ix, iy, iz;
		xlib::ind2sub(np, np, np, i, ix, iy, iz);
		if (mu < 0)
			ix = np - 1 - ix;
		if (xi < 0)
			iy = np - 1 - iy;
		if (eta < 0)
			iz = np - 1 - iz;
		size_t nid;
		xlib::sub2ind(np, np, np, nid, ix, iy, iz);
		pos.at(i) = nid;
	}

	for (size_t i = 0; i < tbs; i++)
	{
		reorder_(i, pos.at(i)) = 1.0;
		back_order_(pos.at(i), i) = 1.0;
	}

	size_t tis = INTERFACE_TOTAL.at(3);
	lift_matrix_p_.resize(tis);
	flux_matrix_p_.resize(tis);
	lift_matrix_.resize(tis);
	flux_matrix_.resize(tis);

	size_t tbls = matrices_.basis_lower_dim_total();
	for (size_t i = 0; i < tis; i++)
	{
		interface_direction dir = interface_direction_list.at(i);
		lift_matrix_p_.at(i) = xlib::mkl_ext::xcalloc<num_t>(tbs*tbls);
		flux_matrix_p_.at(i) = xlib::mkl_ext::xcalloc<num_t>(tbls*tbs);
		lift_matrix_.at(i).bind(tbs, tbls, lift_matrix_p_.at(i));
		lift_matrix_.at(i).gemm(matrices_.lift_matrix(dir), matrices_.flux_matrix(dir), 1.0, 0.0, CblasNoTrans, CblasTrans);

		flux_matrix_.at(i).bind(tbls, tbs, flux_matrix_p_.at(i));
		flux_matrix_.at(i).omatcopy(matrices_.flux_matrix(dir));
		flux_matrix_.at(i).imatcopy('T');

		if (dir == interface_direction::B || dir == interface_direction::F) {
			lift_matrix_.at(i).imatcopy('N', abs(mu_));
			flux_matrix_.at(i).imatcopy('N', 1.0);
		}
		if (dir == interface_direction::L || dir == interface_direction::R) {
			lift_matrix_.at(i).imatcopy('N', abs(xi_));
			flux_matrix_.at(i).imatcopy('N', 1.0);
		}
		if (dir == interface_direction::D || dir == interface_direction::U) {
			lift_matrix_.at(i).imatcopy('N', abs(eta_));
			flux_matrix_.at(i).imatcopy('N', 1.0);
		}

	}

	sysm_.omatadd(matrices_.m(), sysm_, sigma_, 0.0);

	sysm_.omatadd(matrices_.sx(), sysm_, abs(mu_), 1.0);

	sysm_.omatadd(matrices_.sy(), sysm_, abs(xi_), 1.0);

	sysm_.omatadd(matrices_.sz(), sysm_, abs(eta_), 1.0);


	num_p tmpp = xlib::mkl_ext::xcalloc<num_t>(tbs*tbs);
	matrix_t tmpm(tbs, tbs, tmpp);

	size_t idir;

	{
		interface_direction dir = interface_direction::B;
		idir = static_cast<size_t>(dir);
		tmpm.gemm(lift_matrix_.at(idir), matrices_.flux_matrix(dir));
		sysm_.omatadd(tmpm, sysm_, 1.0, 1.0);
	}



	{
		interface_direction dir = interface_direction::L;
		idir = static_cast<size_t>(dir);
		tmpm.gemm(lift_matrix_.at(idir), matrices_.flux_matrix(dir));
		sysm_.omatadd(tmpm, sysm_, 1.0, 1.0);
	}




	{
		interface_direction dir = interface_direction::D;
		idir = static_cast<size_t>(dir);
		tmpm.gemm(lift_matrix_.at(idir), matrices_.flux_matrix(dir));
		sysm_.omatadd(tmpm, sysm_, 1.0, 1.0);
	}


	mkl_free(tmpp);

	ipiv = new MKL_INT[tbs];
	sysmlu_.lu_factorization(sysm_, ipiv);

	sysmiv_.getri(sysm_);




}

dgn::system_matrix_angle::~system_matrix_angle()
{
	if (ipiv != nullptr)
		delete[] ipiv;
	if (sysmp_ != nullptr)
		mkl_free(sysmp_);
	if (sysmlup_ != nullptr)
		mkl_free(sysmlup_);
	if (sysmivp_ != nullptr)
		mkl_free(sysmivp_);
	if (massmp_ != nullptr)
		mkl_free(massmp_);
	for (size_t i = 0; i < lift_matrix_p_.size(); i++)	
		if (lift_matrix_p_.at(i) != nullptr)	
			mkl_free(lift_matrix_p_.at(i));

	for (size_t i = 0; i < flux_matrix_p_.size(); i++)	
		if (flux_matrix_p_.at(i) != nullptr)
			mkl_free(flux_matrix_p_.at(i));	

	if (reorder_p_ != nullptr)
		mkl_free(reorder_p_);
	if (back_order_p_ != nullptr)
		mkl_free(back_order_p_);
}

//==============================================================================

#include <iomanip>

std::ostream & dgn::operator<<(std::ostream & os, const system_matrix_angle & s)
{
	os << std::setprecision(5);
	os << "information of system matrices" << std::endl;
	os << "total basis of element = " << s.basis_element() << std::endl;
	os << "total basis of surface = " << s.basis_surface() << std::endl;

	os << "========   problem info:  ================" << std::endl;
	os << "mu\txi\teta\tsigma" << std::endl;
	os << s.mu() << "\t" << s.xi() << "\t" << s.eta() << "\t" << s.sigma() << std::endl;
	os << "==========================================" << std::endl;

	os << "========   node info:  ===================" << std::endl;
	os << "x" << std::endl << s.x();
	os << "y" << std::endl << s.y();
	os << "z" << std::endl << s.z();
	os << "xs" << std::endl << s.xs();
	os << "ys" << std::endl << s.ys();	
	os << "==========================================" << std::endl;

	os << "========   system matrix:  ===============" << std::endl;
	os << "system matrix:" << std::endl << s.system_matrix();
	os << "inverse of system matrix:" << std::endl << s.system_matrixiv();
	os << "==========================================" << std::endl;

	os << "========   mass matrix: ==================" << std::endl;
	os << "mass matrix:" << std::endl << s.mass_matrix();	
	os << "==========================================" << std::endl;

	os << "========  lift matries:  =================" << std::endl;
	for (size_t i = 0; i < INTERFACE_TOTAL.at(3); i++)
	{
		os << "direction " << interface_direction_name_list.at(i).c_str() << ":" << std::endl;
		os << s.lift_matrix(interface_direction_list.at(i));
	}	
	os << "==========================================" << std::endl;
	
	os << "========  flux matries:  =================" << std::endl;
	for (size_t i = 0; i < INTERFACE_TOTAL.at(3); i++)
	{
		os << "direction " << interface_direction_name_list.at(i).c_str() << ":" << std::endl;
		os << s.flux_matrix(interface_direction_list.at(i));
	}
	os << "==========================================" << std::endl;

	os << "========  reorder matries:  =================" << std::endl;
	os << "reorder " << std::endl << s.reorder();
	os << "back order" << std::endl << s.back_order();
	os << "==========================================" << std::endl;
	return os;
}
