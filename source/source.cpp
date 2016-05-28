#include "source.h"

namespace discontinues_galerkin_nodal_solver {
  Source::Source(const Mesh & mesh, const SystemMatrix & system)
    : mesh_(mesh), system_matrices_(system)
  {
    n_node_ = mesh.n_elements()*system_matrices_.n_basis_element();
    data_.New( n_node() );
  }

  void Source::Copy(const Vector & v)
  {
    data_.Copy(v);
  }

  AnalyticalSource::AnalyticalSource(const Mesh & mesh,
      const SystemMatrix & system, size_t id)
    : Source(mesh, system)
  {
    x_.New(n_node());
    y_.New(n_node());
    z_.New(n_node());
    mu_ = system_matrices().mu();
    xi_ = system_matrices().xi();
    eta_ = system_matrices().eta();
    sigma_ = system_matrices().sigma();
    id_ = id;
    CalculateCoordinates();
  }

  AnalyticalSource::AnalyticalSource(const Mesh & mesh, 
      const SystemMatrix & system, size_t id, 
      const Vector & x, const Vector & y, const Vector & z)
    : Source(mesh, system)
  {
    x_ = x;
    y_ = y;
    z_ = z;
    mu_ = system_matrices().mu();
    xi_ = system_matrices().xi();
    eta_ = system_matrices().eta();
    sigma_ = system_matrices().sigma();
    id_ = id;
  }

  void AnalyticalSource::CalculateValue()
  {
    for (size_t i = 0; i < n_node(); i++)
    {
      data_(i) = AnalyticalSolution::source(x_(i), y_(i), z_(i), 
          mu(), xi(), eta(), sigma(), id());
    }
  }

  void AnalyticalSource::CalculateCoordinates()
  {
    size_t n_basis = system_matrices().n_basis_element();
    std::vector<Vector3> ec = mesh().CenterElements();
    for (size_t i = 0; i < mesh().n_elements(); ++i)
    {
      Vector xt, yt, zt;
      xt.Bind(n_basis, x().p() + i*n_basis);
      yt.Bind(n_basis, y().p() + i*n_basis);
      zt.Bind(n_basis, z().p() + i*n_basis);
      xt.ScalarMulVectorAddScalar(1.0, system_matrices().xe(), ec.at(i).x());
      yt.ScalarMulVectorAddScalar(1.0, system_matrices().ye(), ec.at(i).y());
      zt.ScalarMulVectorAddScalar(1.0, system_matrices().ze(), ec.at(i).z());
    }
  }

  AnalyticalSourceMeshCoordinates::AnalyticalSourceMeshCoordinates(const MeshWithCoordinate & mesh, const SystemMatrix & system, size_t id)
    : Source(mesh, system)
  {
    x_ = mesh.xe();
    y_ = mesh.ye();
    z_ = mesh.ze();
    mu_ = system_matrices().mu();
    xi_ = system_matrices().xi();
    eta_ = system_matrices().eta();
    sigma_ = system_matrices().sigma();
    id_ = id;
  }

  void AnalyticalSourceMeshCoordinates::CalculateValue()
  {
    for (size_t i = 0; i < n_node(); i++)
    {
      data_(i) = AnalyticalSolution::source(x_(i), y_(i), z_(i),
        mu(), xi(), eta(), sigma(), id());
    }
  }


}