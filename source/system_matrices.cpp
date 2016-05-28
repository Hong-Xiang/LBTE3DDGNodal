#include <stdexcept>
#include "system_matrices.h"
#include "utilities.h"

namespace discontinues_galerkin_nodal_solver {

  SystemMatrix::SystemMatrix(num_t h, size_t np, num_t mu, num_t xi, num_t eta, num_t sigma)
    : h_(h), np_(np), mu_(mu), xi_(xi), eta_(eta), sigma_(sigma)
  {
    size_t nbe = Matrices3D::n_basis(np);
    size_t nbs = Matrices2D::n_basis(np);
    n_basis_element_ = nbe;
    n_basis_surface_ = nbs;
    x_e_.New(nbe);
    y_e_.New(nbe);
    z_e_.New(nbe);

    x_s_.resize(kSurfaceTotal.at(3));
    y_s_.resize(kSurfaceTotal.at(3));
    z_s_.resize(kSurfaceTotal.at(3));
    for (auto& v : x_s_) v.New(nbs);
    for (auto& v : y_s_) v.New(nbs);
    for (auto& v : z_s_) v.New(nbs);
    x_s3_.resize(3);
    y_s3_.resize(3);
    z_s3_.resize(3);
    for (auto& v : x_s3_) v.New(nbs);
    for (auto& v : y_s3_) v.New(nbs);
    for (auto& v : z_s3_) v.New(nbs);

    sysm_.New(nbe, nbe);
    sysmlu_.New(nbe, nbe);
    sysmiv_.New(nbe, nbe);
    massm_.New(nbe, nbe);

    lift_matrix_.resize(kSurfaceTotal.at(3));
    flux_matrix_.resize(kSurfaceTotal.at(3));
    for (auto& lm : lift_matrix_)
      lm.New(nbe, nbs);
    for (auto& fm : flux_matrix_)
      fm.New(nbs, nbe);

    CalculateMatrix();
    CalculateCoordinates();
  }

  SystemMatrix::~SystemMatrix()
  {
  }

  void discontinues_galerkin_nodal_solver::SystemMatrix::CalculateMatrix()
  {

    Matrices3D element_matrices(h(), np());
    element_matrices.Initialize();
    Matrices2D surface_matrices(h(), np());
    surface_matrices.Initialize();

    for (auto ser : kSurfaceElementRelationList)
    {
      size_t i = Utilities::id(ser);
      num_t absangle;
      if (ser == SurfaceElementRelation::B || ser == SurfaceElementRelation::F) {
        absangle = std::abs(mu());
      }
      if (ser == SurfaceElementRelation::L || ser == SurfaceElementRelation::R) {
        absangle = std::abs(xi());
      }
      if (ser == SurfaceElementRelation::D || ser == SurfaceElementRelation::U) {
        absangle = std::abs(eta());
      }


      lift_matrix_.at(i).AddMatrixMulMatrix(absangle,
        element_matrices.lift_matrix(ser), xlib::mkl::Transport::N,
        element_matrices.flux_matrix(ser), xlib::mkl::Transport::T,
        0.0);

      flux_matrix_.at(i).ScalarMulMatrix(1.0,
        element_matrices.flux_matrix(ser), xlib::mkl::Transport::N);
    }

    massm_.ScalarMulMatrix(sigma(), element_matrices.m(), 
        xlib::mkl::Transport::N);
    sysm_.ScalarMulMatrixAddScalarMulMatrix(sigma(),
      element_matrices.m(), xlib::mkl::Transport::N,
      0.0, sysm_, xlib::mkl::Transport::N);
    sysm_.ScalarMulMatrixAddScalarMulMatrix(mu(),
      element_matrices.sx(), xlib::mkl::Transport::N,
      1.0, sysm_, xlib::mkl::Transport::N);
    sysm_.ScalarMulMatrixAddScalarMulMatrix(xi(),
      element_matrices.sy(), xlib::mkl::Transport::N,
      1.0, sysm_, xlib::mkl::Transport::N);
    sysm_.ScalarMulMatrixAddScalarMulMatrix(eta(),
      element_matrices.sz(), xlib::mkl::Transport::N,
      1.0, sysm_, xlib::mkl::Transport::N);
    Matrix tmpm{ static_cast<MKL_INT>(n_basis_element()),
        static_cast<MKL_INT>(n_basis_element()) };

    size_t iser;
    auto quad = Utilities::quadrant_of_angle(mu(), xi(), eta());
    auto ser_x = Utilities::upwind_source_direction_x(quad);
    auto ser_y = Utilities::upwind_source_direction_y(quad);
    auto ser_z = Utilities::upwind_source_direction_z(quad);
    
    iser = Utilities::id(ser_x);
    tmpm.AddMatrixMulMatrix(1.0,
      lift_matrix_.at(iser), xlib::mkl::Transport::N,
      element_matrices.flux_matrix(ser_x), xlib::mkl::Transport::N,
      0.0);    
    sysm_.ScalarMulMatrixAddScalarMulMatrix(
      1.0, tmpm, xlib::mkl::Transport::N,
      1.0, sysm_, xlib::mkl::Transport::N);
    iser = Utilities::id(ser_y);
    tmpm.AddMatrixMulMatrix(1.0,
      lift_matrix_.at(iser), xlib::mkl::Transport::N,
      element_matrices.flux_matrix(ser_y), xlib::mkl::Transport::N,
      0.0);
    sysm_.ScalarMulMatrixAddScalarMulMatrix(
      1.0, tmpm, xlib::mkl::Transport::N,
      1.0, sysm_, xlib::mkl::Transport::N);
    iser = Utilities::id(ser_z);
    tmpm.AddMatrixMulMatrix(1.0,
      lift_matrix_.at(iser), xlib::mkl::Transport::N,
      element_matrices.flux_matrix(ser_z), xlib::mkl::Transport::N,
      0.0);
    sysm_.ScalarMulMatrixAddScalarMulMatrix(
      1.0, tmpm, xlib::mkl::Transport::N,
      1.0, sysm_, xlib::mkl::Transport::N);

    sysmlu_.LuFactorization(system_matrix());
    sysmiv_.Inverse(system_matrix());

   
  }

  void SystemMatrix::CalculateCoordinates()
  {
    Matrices3D element_matrices(h(), np());
    element_matrices.Initialize();
    Matrices2D surface_matrices(h(), np());
    surface_matrices.Initialize();

    x_e_.Copy(element_matrices.x());
    y_e_.Copy(element_matrices.y());
    z_e_.Copy(element_matrices.z());
    for (auto ser : kSurfaceElementRelationList)
    {
      size_t id = Utilities::id(ser);
      num_t xc = 0.0, yc = 0.0, zc = 0.0;

      if (ser == SurfaceElementRelation::B) {
        xc -= h() / 2;
        x_s_.at(id).Ones();
        x_s_.at(id).MulScalar(xc);
        y_s_.at(id).Copy(surface_matrices.x());
        z_s_.at(id).Copy(surface_matrices.y());
      }

      if (ser == SurfaceElementRelation::F)
      {
        xc += h() / 2;
        x_s_.at(id).Ones();
        x_s_.at(id).MulScalar(xc);
        y_s_.at(id).Copy(surface_matrices.x());
        z_s_.at(id).Copy(surface_matrices.y());
      }

      if (ser == SurfaceElementRelation::L) {
        yc -= h() / 2;
        y_s_.at(id).Ones();
        y_s_.at(id).MulScalar(yc);
        x_s_.at(id).Copy(surface_matrices.x());
        z_s_.at(id).Copy(surface_matrices.y());
      }

      if (ser == SurfaceElementRelation::R)
      {
        yc += h() / 2;
        y_s_.at(id).Ones();
        y_s_.at(id).MulScalar(yc);
        x_s_.at(id).Copy(surface_matrices.x());
        z_s_.at(id).Copy(surface_matrices.y());
      }

      if (ser == SurfaceElementRelation::D) {
        zc -= h() / 2;
        z_s_.at(id).Ones();
        z_s_.at(id).MulScalar(zc);
        x_s_.at(id).Copy(surface_matrices.x());
        y_s_.at(id).Copy(surface_matrices.y());
      }

      if (ser == SurfaceElementRelation::U)
      {
        zc += h() / 2;
        z_s_.at(id).Ones();
        z_s_.at(id).MulScalar(zc);
        x_s_.at(id).Copy(surface_matrices.x());
        y_s_.at(id).Copy(surface_matrices.y());
      }
    }

    for (auto sd : kSurfaceDirectionList)
    {
      if (sd == SurfaceDirection::X)
      {
        size_t id = Utilities::id(sd);
        x_s3_.at(id).Zeros();
        y_s3_.at(id).Copy(surface_matrices.x());
        z_s3_.at(id).Copy(surface_matrices.y());
      }
      if (sd == SurfaceDirection::Y)
      {
        size_t id = Utilities::id(sd);
        y_s3_.at(id).Zeros();
        x_s3_.at(id).Copy(surface_matrices.x());
        z_s3_.at(id).Copy(surface_matrices.y());
      }
      if (sd == SurfaceDirection::Z)
      {
        size_t id = Utilities::id(sd);
        z_s3_.at(id).Zeros();
        x_s3_.at(id).Copy(surface_matrices.x());
        y_s3_.at(id).Copy(surface_matrices.y());
      }
    }

  }

  //==============================================================================

#include <iomanip>

  std::ostream & operator<<(std::ostream & os, const SystemMatrix & s)
  {
    os << std::setprecision(5);
    os << "information of system matrices" << std::endl;
    os << "total basis of Element = " << s.n_basis_element() << std::endl;
    os << "total basis of Surface = " << s.n_basis_surface() << std::endl;

    os << "========   problem info:  ================" << std::endl;
    os << "mu\txi\teta\tsigma" << std::endl;
    os << s.mu() << "\t" << s.xi() << "\t" << s.eta() << "\t" << s.sigma() << std::endl;
    os << "==========================================" << std::endl;

    os << "========   node info:  ===================" << std::endl;
    os << "x" << std::endl << s.xe();
    os << "y" << std::endl << s.ye();
    os << "z" << std::endl << s.ze();
    for (auto ser : kSurfaceElementRelationList)
    {
      size_t id = Utilities::id(ser);
      os << "coordinate of surface " 
                << Utilities::toString(ser) << std::endl;      
      os << "xs" << std::endl << s.xs(ser);
      os << "ys" << std::endl << s.ys(ser);
      os << "zs" << std::endl << s.zs(ser);
    }
    
    os << "==========================================" << std::endl;

    os << "========   system matrix:  ===============" << std::endl;
    os << "system matrix:" << std::endl << s.system_matrix();
    os << "inverse of system matrix:" << std::endl << s.system_matrixiv();
    os << "==========================================" << std::endl;

    os << "========   mass matrix: ==================" << std::endl;
    os << "mass matrix:" << std::endl << s.mass_matrix();
    os << "==========================================" << std::endl;

    os << "========  lift matries:  =================" << std::endl;
    for (auto ser : kSurfaceElementRelationList)
    {      
      os << "surface " << Utilities::toString(ser) << ": " << std::endl;
      os << s.lift_matrix(ser);
    }
    os << "==========================================" << std::endl;

    os << "========  flux matries:  =================" << std::endl;
    for (auto ser : kSurfaceElementRelationList)
    {
      os << "surface " << Utilities::toString(ser) << ": " << std::endl;
      os << s.flux_matrix(ser);
    }
    os << "==========================================" << std::endl;

    return os;
  }
}