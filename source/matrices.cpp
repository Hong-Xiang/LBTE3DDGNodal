#include "matrices.h"

namespace discontinues_galerkin_nodal_solver {
  /************************************************************************/
  /*          MatricesData                                                */
  /************************************************************************/
  MatricesData* MatricesData::handle_ = nullptr;

  const MatricesData & MatricesData::instance()
  {
    if (handle_ == nullptr)
      handle_ = new MatricesData();
    return *handle_;
  }

  const Vector MatricesData::x(size_t Np) const
  {
    if (Np < kMinNp || Np > kMaxNp)
      BOOST_THROW_EXCEPTION(ExceptionOutOfRange() << 
          InfoOutOfRange(Np, kMinNp, kMaxNp));
    return x_.at(Np - kMinNp);
  }

  const Matrix MatricesData::m(size_t Np) const
  {
    if (Np < kMinNp || Np > kMaxNp)
      BOOST_THROW_EXCEPTION(ExceptionOutOfRange() << 
          InfoOutOfRange(Np, kMinNp, kMaxNp));
    return m_.at(Np - kMinNp);
  }

  const Matrix MatricesData::s(size_t Np) const
  {
    if (Np < kMinNp || Np > kMaxNp)
      BOOST_THROW_EXCEPTION(ExceptionOutOfRange() << 
          InfoOutOfRange(Np, kMinNp, kMaxNp));
    return s_.at(Np - kMinNp);
  }

  MatricesData::MatricesData()
  {
    x_.resize(kMaxNp - kMinNp + 1);
    m_.resize(kMaxNp - kMinNp + 1);
    s_.resize(kMaxNp - kMinNp + 1);
    for (size_t i = kMinNp; i <= kMaxNp; i++)
    {
      x_.at(i - kMinNp).New(i);
      m_.at(i - kMinNp).New(i, i);
      s_.at(i - kMinNp).New(i, i);
    }
    read_MAT_files();
  }

  void MatricesData::read_MAT_files()
  {
    for (size_t i = kMinNp; i <= kMaxNp; i++)
    {
#ifdef DGNSDEBUG
      std::cout << "Reading matrices files " << i << std::endl;
#endif
      const char* filename = "matrices.mat";
      size_t pos = i - kMinNp;
      xlib::imatlab::MatFileHelper::LoadCellDataRaw(x_.at(pos).p(), filename,
          "x", i-1, std::vector<size_t>{i, 1});
      xlib::imatlab::MatFileHelper::LoadCellDataRaw(m_.at(pos).p(), filename,
          "M", i-1, std::vector<size_t>{i, i});
      xlib::imatlab::MatFileHelper::LoadCellDataRaw(s_.at(pos).p(), filename,
          "S", i-1, std::vector<size_t>{i, i});
    }
  }

  std::ostream& operator<<(std::ostream& os, const MatricesData& m)
  {
    os << "Info of MatricesData:" << std::endl;
    for (size_t i = m.min_np(); i < m.max_np(); ++i)
    {
      os << "Matrices data of order " << i << "." << std::endl;
      os << "x:\n" << m.x(i);
      os << "m:\n" << m.m(i);
      os << "s:\n" << m.s(i);
    }
    return os;
  }

  //========================================================================




 /************************************************************************/
 /*                   MatricesBase                                       */
 /************************************************************************/

  MatricesBase::MatricesBase(num_t h, size_t np)
    : h_(h), np_(np), is_init_(false), dim_(0)
  {    
  }

  void MatricesBase::reset(num_t h, size_t np)
  {
    h_ = h;
    np_ = np;
    Initialize();
  }

  void MatricesBase::Initialize()
  {
    dim_ = dimension();
    MemoryAlloc();
    Calculate();
    is_init_ = true;
  }

  MatricesBase::~MatricesBase()
  {
#ifdef DGNSDEBUG
    std::cout << "Deconstruction of Matrices base called." << std::endl;
#endif
  }

  void MatricesBase::MemoryAlloc()
  {
    size_t ift = kSurfaceTotal.at(dim_);    
    lift_matrices_.resize(ift);
    flux_matrices_.resize(ift);
    for (size_t i = 0; i < ift; i++)
    {     
      lift_matrices_.at(i).New(n_basis(), n_basis());      
      flux_matrices_.at(i).New(n_basis_lower(), n_basis());
    }
  }

  //========================================================================

  /************************************************************************/
  /*                        Matrices1D                                    */
  /************************************************************************/
  Matrices1D::Matrices1D(num_t h, size_t np)
    : MatricesBase(h, np)
  {
    dim_ = dimension();
  }

  Matrices1D::~Matrices1D()
  {
#ifdef DGNSDEBUG
    std::cout << "Deconstruction of Matrices1d base called." << std::endl;
#endif
  }

  void Matrices1D::Calculate()
  {
    const MatricesData& md = MatricesData::instance();
    size_t nb = n_basis();
    x_.AddVector(h_ / 2.0, md.x(np()));
    m_.ScalarMulMatrix(h_ / 2.0, md.m(np()), xlib::mkl::Transport::N);
    sx_.ScalarMulMatrix(1.0, md.s(np()), xlib::mkl::Transport::N);    
    lift_matrices_.at(Utilities::id(SurfaceElementRelation::B))(1, 1) = 1.0;
    lift_matrices_.at(Utilities::id(SurfaceElementRelation::F))(np(), np()) = 1.0;
  }

  void Matrices1D::MemoryAlloc()
  {
    MatricesBase::MemoryAlloc();
    size_t nb = n_basis();
    x_.New(nb);
    m_.New(nb, nb);
    sx_.New(nb, nb);
  }

  std::ostream& operator<<(std::ostream& os, const Matrices1D& m) {
    os << "Info of Matrices1D: h= " << m.h() << "; np= " << m.np() << std::endl;
    os << "x:\n" << m.x();
    os << "m:\n" << m.m();
    os << "sx:\n" << m.sx();
    for (size_t i = 0; i < kSurfaceTotal.at(m.dimension()); ++i)
    {
      SurfaceElementRelation ser = kSurfaceElementRelationList.at(i);
      os << "lift matrix of " << Utilities::toString(ser) << "\n"
         << m.lift_matrix(ser);  
      os << "flux matrix of " << Utilities::toString(ser) << "\n"
         << m.flux_matrix(ser);
    }
    return os;
  }

  //==============================================================================







  //========	Implement of matrices2d	============================================

  Matrices2D::Matrices2D(num_t h, size_t np)
    : MatricesBase(h, np)
  {
    dim_ = dimension();
  }

  Matrices2D::~Matrices2D()
  {
#ifdef DGNSDEBUG
    std::cout << "Deconstruction of Matrices2D base called." << std::endl;
#endif
  }

  void Matrices2D::Calculate()
  {
    const MatricesData& md = MatricesData::instance();
    const Vector xmd = md.x(np_);
    const Matrix mmd = md.m(np_);
    const Matrix smd = md.s(np_);

    for (size_t i = 0; i < n_basis(); i++)
    {
      size_t ix, iy;
      xlib::Utilities::ind2sub(np_, np_, i, ix, iy);
      x_(i) = xmd(ix) * h_ / 2.0;
      y_(i) = xmd(iy) * h_ / 2.0;
      for (size_t j = 0; j < n_basis(); j++)
      {
        size_t jx, jy;
        xlib::Utilities::ind2sub(np_, np_, j, jx, jy);
        m_(i, j) = mmd(ix, jx)*mmd(iy, jy) * h_ * h_ / 4.0;
        sx_(i, j) = smd(ix, jx)*mmd(iy, jy) * h_ / 2.0;
        sy_(i, j) = mmd(ix, jx)*smd(iy, jy) * h_ / 2.0;
      }
    }

    Matrices1D m1d(h_, np_);
    m1d.Initialize();
    for (size_t i = 0; i < kSurfaceTotal.at(dimension()); i++)
    {
      SurfaceElementRelation idir = kSurfaceElementRelationList.at(i);
      std::vector<size_t> idx;
      std::vector<bool> bdmask;
      size_t cid = 0;
      idx.resize(n_basis());
      bdmask.resize(n_basis());
      for (size_t j = 0; j < n_basis(); j++)
      {
        size_t ix, iy;
        xlib::Utilities::ind2sub(np_, np_, j, ix, iy);
        bdmask.at(j) = Utilities::is_boundary(idir, np_, ix, iy);
        if (bdmask.at(j))
        {
          idx.at(j) = cid;
          cid++;
        }
      }
      for (size_t j = 0; j < n_basis(); j++)
      {
        if (bdmask.at(j))
          flux_matrices_.at(i)(idx.at(j), j) = 1.0;
        for (size_t k = 0; k < n_basis(); k++)
          if (bdmask.at(j) && bdmask.at(k))
            lift_matrices_.at(i)(j, k) = m1d.m()(idx.at(j), idx.at(k));
      }
    }
  }

  void Matrices2D::MemoryAlloc()
  {
    MatricesBase::MemoryAlloc();
    size_t nb = n_basis();
    x_.New(nb);
    y_.New(nb);
    m_.New(nb, nb);
    sx_.New(nb, nb);
    sy_.New(nb, nb);
  }

  std::ostream& operator<<(std::ostream& os, const Matrices2D& m) {
    os << "Info of Matrices2D: h= " << m.h() << "; np= " << m.np() << std::endl;
    os << "x:\n" << m.x();
    os << "y:\n" << m.y();
    os << "m:\n" << m.m();
    os << "sx:\n" << m.sx();
    os << "sy:\n" << m.sy();
    for (size_t i = 0; i < kSurfaceTotal.at(m.dimension()); ++i)
    {
      SurfaceElementRelation ser = kSurfaceElementRelationList.at(i);
      os << "lift matrix of " << Utilities::toString(ser) << "\n"
        << m.lift_matrix(ser);
      os << "flux matrix of " << Utilities::toString(ser) << "\n"
        << m.flux_matrix(ser);
    }
    return os;
  }

  //==============================================================================


  //========	Implement of matrices3d	============================================

  Matrices3D::Matrices3D(num_t h, size_t np)
    : MatricesBase(h, np)
  {
    dim_ = dimension();
  }

  Matrices3D::~Matrices3D()
  {
#ifdef DGNSDEBUG
    std::cout << "Deconstruction of Matrices3D base called." << std::endl;
#endif
  }


  void Matrices3D::Calculate()
  {
    const MatricesData& md = MatricesData::instance();
    const Vector xmd = md.x(np_);
    const Matrix mmd = md.m(np_);
    const Matrix smd = md.s(np_);

    for (size_t i = 0; i < n_basis(); i++)
    {
      size_t ix, iy, iz;
      xlib::Utilities::ind2sub(np_, np_, np_, i, ix, iy, iz);
      x_(i) = xmd(ix) * h_ / 2.0;
      y_(i) = xmd(iy) * h_ / 2.0;
      z_(i) = xmd(iz) * h_ / 2.0;
      for (size_t j = 0; j < n_basis(); j++)
      {
        size_t jx, jy, jz;
        xlib::Utilities::ind2sub(np_, np_, np_, j, jx, jy, jz);
        m_(i, j) = mmd(ix, jx) * mmd(iy, jy) * mmd(iz, jz) * h_ * h_ * h_ / 8.0;
        sx_(i, j) = smd(ix, jx) * mmd(iy, jy) * mmd(iz, jz) * h_ * h_ / 4.0;
        sy_(i, j) = mmd(ix, jx) * smd(iy, jy) * mmd(iz, jz) * h_ * h_ / 4.0;
        sz_(i, j) = mmd(ix, jx) * mmd(iy, jy) * smd(iz, jz) * h_ * h_ / 4.0;
      }
    }

    Matrices2D m2d(h_, np_);
    m2d.Initialize();
    for (size_t i = 0; i < kSurfaceTotal.at(dimension()); i++)
    {
      SurfaceElementRelation idir = kSurfaceElementRelationList.at(i);
      std::vector<size_t> idx;
      std::vector<bool> bdmask;
      size_t cid = 0;
      idx.resize(n_basis());
      bdmask.resize(n_basis());
      for (size_t j = 0; j < n_basis(); j++)
      {
        size_t ix, iy, iz;
        xlib::Utilities::ind2sub(np_, np_, np_, j, ix, iy, iz);
        bdmask.at(j) = Utilities::is_boundary(idir, np_, ix, iy, iz);
        if (bdmask.at(j))
        {
          idx.at(j) = cid;
          cid++;
        }
      }
      for (size_t j = 0; j < n_basis(); j++)
      {
        if (bdmask.at(j))
          flux_matrices_.at(i)(idx.at(j), j) = 1.0;
        for (size_t k = 0; k < n_basis(); k++)
          if (bdmask.at(j) && bdmask.at(k))
            lift_matrices_.at(i)(j, k) = m2d.m()(idx.at(j), idx.at(k));
      }
    }
  }

  void Matrices3D::MemoryAlloc()
  {
    MatricesBase::MemoryAlloc();
    size_t nb = n_basis();
    x_.New(nb);
    y_.New(nb);
    z_.New(nb);
    m_.New(nb, nb);
    sx_.New(nb, nb);
    sy_.New(nb, nb);
    sz_.New(nb, nb);
  }

  std::ostream& operator<<(std::ostream& os, const Matrices3D& m) {
    os << "Info of Matrices1D of: h= " << m.h() << "; np= "<< m.np() << std::endl;
    os << "x:\n" << m.x();
    os << "y:\n" << m.y();
    os << "z:\n" << m.z();
    os << "m:\n" << m.m();
    os << "sx:\n" << m.sx();
    os << "sy:\n" << m.sy();
    os << "sz:\n" << m.sz();
    for (size_t i = 0; i < kSurfaceTotal.at(m.dimension()); ++i)
    {
      SurfaceElementRelation ser = kSurfaceElementRelationList.at(i);
      os << "lift matrix of " << Utilities::toString(ser) << "\n"
        << m.lift_matrix(ser);
      os << "flux matrix of " << Utilities::toString(ser) << "\n"
        << m.flux_matrix(ser);
    }
    return os;
  }
  //========================================================================
}