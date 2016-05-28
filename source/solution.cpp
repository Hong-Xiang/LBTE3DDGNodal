#include "solution.h"

namespace discontinues_galerkin_nodal_solver {
  Solution::Solution(const Mesh & mesh, const SystemMatrix & system)
    : mesh_(mesh), system_matrices_(system)
  {
    n_node_ = mesh.n_elements()*system_matrices_.n_basis_element();
    data_.New(n_node());
  }

  void Solution::Copy(const Vector & v)
  {
    data_.Copy(v);
  }
  AnalyticalSolutionMeshCoordinates::AnalyticalSolutionMeshCoordinates(const MeshWithCoordinate & mesh, const SystemMatrix & system, size_t id)
    : Solution(mesh, system)
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
  void AnalyticalSolutionMeshCoordinates::CalculateValue()
  {
    for (size_t i = 0; i < n_node(); i++)
    {
      data_(i) = AnalyticalSolution::solution(x_(i), y_(i), z_(i),
        mu(), xi(), eta(), sigma(), id());
    }
  }
}