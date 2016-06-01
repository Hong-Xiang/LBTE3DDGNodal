#pragma once
#include "types.h"

namespace discontinues_galerkin_nodal_solver {
  class ProblemDefinitionSingleAngle {
  public:
    void Load(std::string filename);
    num_t xc() const { return xc_; }
    num_t yc() const { return yc_; }
    num_t zc() const { return zc_; }
    num_t mu() const { return mu_; }
    num_t xi() const { return xi_; }
    num_t eta() const { return eta_; }
    num_t sigma() const { return sigma_; }
    num_t h() const { return h_; }
    size_t nx() const { return nx_; }
    size_t ny() const { return ny_; }
    size_t nz() const { return nz_; }
    size_t np() const { return np_; }
    size_t test_id() const { return test_id_; }
  private:
    num_t xc_, yc_, zc_;
    num_t h_;
    size_t nx_, ny_, nz_;
    num_t mu_, xi_, eta_;
    num_t sigma_;
    size_t np_;
    size_t test_id_;
  };
}