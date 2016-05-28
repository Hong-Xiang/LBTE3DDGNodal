#include "problem_definition.h"
#include <fstream>

namespace discontinues_galerkin_nodal_solver {
  void ProblemDefinitionSingleAngle::Load(std::string filename)
  {
    std::ifstream fin(filename.c_str());
    fin >> xc_ >> yc_ >> zc_;
    fin >> h_;
    fin >> np_;
    fin >> nx_ >> ny_ >> nz_;
    fin >> mu_ >> xi_ >> eta_;
    fin >> sigma_;
    fin >> test_id_;
    fin.close();
  }

}