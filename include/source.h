#pragma once
#include "global.h"
#include "mesh.h"
#include "system_matrices.h"

namespace discontinues_galerkin_nodal_solver {
  class Source {
  public:
    Source(const Mesh& mesh, const SystemMatrix& system);
    Vector Value() const { return data_; }
    void Copy(const Vector& v);
    size_t n_node() const { return n_node_; }
    const Mesh& mesh() const { return mesh_; }
    const SystemMatrix& system_matrices() const { return system_matrices_; }
  protected:
    Vector data_;
  private:
    size_t n_node_;         
    const Mesh& mesh_;
    const SystemMatrix& system_matrices_;
  };

  class AnalyticalSource : public Source{
  public:
    AnalyticalSource(const Mesh& mesh, const SystemMatrix& system, size_t id);
    AnalyticalSource(const Mesh& mesh, const SystemMatrix& system, size_t id,
        const Vector& x, const Vector& y, const Vector& z);
    void CalculateValue();
    Vector x() const { return x_; }
    Vector y() const { return y_; }
    Vector z() const { return z_; }
    num_t mu() const { return mu_; }
    num_t xi() const { return xi_; }
    num_t eta() const { return eta_; }
    num_t sigma() const { return sigma_; }
    size_t id() const { return id_; }
  private:
    void CalculateCoordinates();

  private:
    size_t id_;
    Vector x_;
    Vector y_;
    Vector z_;
    num_t mu_, xi_, eta_, sigma_;
  };


  class AnalyticalSourceMeshCoordinates : public Source {
  public:
    AnalyticalSourceMeshCoordinates(const MeshWithCoordinate& mesh,
        const SystemMatrix& system, size_t id);
    void CalculateValue();
    Vector x() const { return x_; }
    Vector y() const { return y_; }
    Vector z() const { return z_; }
    num_t mu() const { return mu_; }
    num_t xi() const { return xi_; }
    num_t eta() const { return eta_; }
    num_t sigma() const { return sigma_; }
    size_t id() const { return id_; }
  private:
    size_t id_;
    Vector x_;
    Vector y_;
    Vector z_;
    num_t mu_, xi_, eta_, sigma_;
  };
}