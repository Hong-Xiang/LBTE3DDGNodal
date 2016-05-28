#pragma once
#include "types.h"
#include "matrices.h"

namespace discontinues_galerkin_nodal_solver {
    
  class SystemMatrix {
  private:
    static size_t n_basis(size_t np) {
      return Matrices3D::n_basis(np);
    }
    static size_t n_basis_low(size_t np) {
      return Matrices3D::n_basis_lower(np);
    }
  public:
    SystemMatrix(num_t h, size_t np, num_t mu, num_t xi, num_t eta, num_t sigma);
    const Matrix system_matrix() const { return sysm_; }
    const MatrixLu system_matrixlu() const { return sysmlu_; }
    const Matrix system_matrixiv() const { return sysmiv_; }
    const Matrix mass_matrix() const { return massm_; }
    const Matrix lift_matrix(SurfaceElementRelation dir) const {
      return lift_matrix_.at(Utilities::id(dir));
    }
    const Matrix flux_matrix(SurfaceElementRelation dir) const {
      return flux_matrix_.at(Utilities::id(dir));
    }
    const Vector& xe() const { return x_e_; }
    const Vector& ye() const { return y_e_; }
    const Vector& ze() const { return z_e_; }

    const Vector& xs(SurfaceElementRelation ser) const {
      return x_s_.at(Utilities::id(ser));
    }
    const Vector& ys(SurfaceElementRelation ser) const {
      return y_s_.at(Utilities::id(ser));
    }
    const Vector& zs(SurfaceElementRelation ser) const {
      return z_s_.at(Utilities::id(ser));
    }

    const Vector& xs(SurfaceDirection sd) const {
      return x_s3_.at(Utilities::id(sd));
    }

    const Vector& ys(SurfaceDirection sd) const {
      return y_s3_.at(Utilities::id(sd));
    }

    const Vector& zs(SurfaceDirection sd) const {
      return z_s3_.at(Utilities::id(sd));
    }

    const size_t np() const { return np_; }
    num_t mu() const { return mu_; }
    num_t xi() const { return xi_; }
    num_t eta() const { return eta_; }
    num_t sigma() const { return sigma_; }
    num_t h() const { return h_; }

    size_t n_basis_element() const { return n_basis_element_; }
    size_t n_basis_surface() const { return n_basis_surface_; }
    ~SystemMatrix();
  private:
    void CalculateMatrix();
    void CalculateCoordinates();
  private:
    num_t h_, mu_, xi_, eta_, sigma_;
    size_t np_;

    //Matrices3D matrices_;
    //Matrices2D surface_matrices_;

    Matrix sysm_;
    MatrixLu sysmlu_;
    Matrix sysmiv_;
    Matrix massm_;
    std::vector<Matrix> lift_matrix_;
    std::vector<Matrix> flux_matrix_;

    size_t n_basis_element_;
    size_t n_basis_surface_;

    Vector x_e_, y_e_, z_e_;
    std::vector<Vector> x_s_, y_s_, z_s_;
    std::vector<Vector> x_s3_, y_s3_, z_s3_;

  };
  std::ostream& operator<<(std::ostream& os, const SystemMatrix& s);
}
