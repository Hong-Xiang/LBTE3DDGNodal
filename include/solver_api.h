#pragma once
#include "mesh.h"
namespace discontinues_galerkin_nodal_solver {
  class SingleAngleSolver {
  public:
    SingleAngleSolver(const MeshWithCoordinate& m);
    void Initialization();
    void Solve();
    Vector solution();
    Vector source();
    Vector boundary();
  };
}