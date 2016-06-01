#include "boundary.h"

namespace discontinues_galerkin_nodal_solver {
  Boundary::Boundary(const Mesh & mesh, const SystemMatrix & system)
    : mesh_(mesh), system_matrices_(system)
  {
    n_node_ = mesh.n_surfaces()*system_matrices_.n_basis_surface();
    data_.New(n_node());    
  }

  void Boundary::Copy(const Vector & v)
  {
    data_.Copy(v);
  }

  AnalyticalBoundary::AnalyticalBoundary(const Mesh & mesh,
    const SystemMatrix & system, size_t id)
    : Boundary(mesh, system)
  {
    x_.New(n_node());
    y_.New(n_node());
    z_.New(n_node());
    boundary_mask_.New(n_node());
    mu_ = system_matrices().mu();
    xi_ = system_matrices().xi();
    eta_ = system_matrices().eta();
    sigma_ = system_matrices().sigma();
    id_ = id;
    CalculateCoordinates();
  }

  AnalyticalBoundary::AnalyticalBoundary(const Mesh & mesh,
    const SystemMatrix & system, size_t id,
    const Vector & x, const Vector & y, const Vector & z, const Vector& bdmask)
    : Boundary(mesh, system)
  {
    x_ = x;
    y_ = y;
    z_ = z;
    boundary_mask_ = bdmask;
    mu_ = system_matrices().mu();
    xi_ = system_matrices().xi();
    eta_ = system_matrices().eta();
    sigma_ = system_matrices().sigma();
    id_ = id;
  }

  void AnalyticalBoundary::CalculateValue()
  {
    for (size_t i = 0; i < n_node(); i++)
    {
      data_(i) = AnalyticalSolution::solution(x_(i), y_(i), z_(i),
        mu(), xi(), eta(), sigma(), id());
    }
    data_.PointMulVector(data_, boundary_mask_);
  }

  void AnalyticalBoundary::CalculateCoordinates()
  {
    size_t n_basis = system_matrices().n_basis_surface();
    std::vector<Vector3> ec = mesh().CenterSurfaces();
    std::vector<bool> bdf = mesh().BoundaryMask();
    std::vector<SurfaceDirection> sd = mesh().SurfaceType();
    for (size_t i = 0; i < mesh().n_surfaces(); ++i)
    {
      Vector xt, yt, zt, bdm;
      xt.Bind(n_basis, x().p() + i*n_basis);
      yt.Bind(n_basis, y().p() + i*n_basis);
      zt.Bind(n_basis, z().p() + i*n_basis);
      bdm.Bind(n_basis, boundary_mask_.p() + i*n_basis);
      SurfaceDirection sdv = sd.at(i);
      xt.ScalarMulVectorAddScalar(1.0, system_matrices().xs(sdv), ec.at(i).x());
      yt.ScalarMulVectorAddScalar(1.0, system_matrices().ys(sdv), ec.at(i).y());
      zt.ScalarMulVectorAddScalar(1.0, system_matrices().zs(sdv), ec.at(i).z());
      if (bdf.at(i))
        bdm.Ones();
      else
        bdm.Zeros();
    }
  }
  AnalyticalBoundaryMeshCoordinate::AnalyticalBoundaryMeshCoordinate(const MeshWithCoordinate & mesh, const SystemMatrix & system, size_t id)
    : Boundary(mesh, system)
  {
    x_ = mesh.xs();
    y_ = mesh.ys();
    z_ = mesh.zs();
    boundary_mask_ = mesh.BoundaryMaskNodal();
    mu_ = system_matrices().mu();
    xi_ = system_matrices().xi();
    eta_ = system_matrices().eta();
    sigma_ = system_matrices().sigma();
    id_ = id;
  }
  void AnalyticalBoundaryMeshCoordinate::CalculateValue()
  {
    for (size_t i = 0; i < n_node(); i++)
    {
      data_(i) = AnalyticalSolution::solution(x_(i), y_(i), z_(i),
        mu(), xi(), eta(), sigma(), id());
    }
    data_.PointMulVector(data_, boundary_mask_);
  }
}