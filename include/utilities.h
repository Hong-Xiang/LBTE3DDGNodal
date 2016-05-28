#pragma once
#include "types.h"
#include <string>
namespace discontinues_galerkin_nodal_solver{

  class Utilities {
  public:
    static Quadrant quadrant_of_angle(num_t mu, num_t xi, num_t eta);
    static bool quadrant_x_positive_flag(Quadrant quad);
    static bool quadrant_y_positive_flag(Quadrant quad);
    static bool quadrant_z_positive_flag(Quadrant quad);
    static bool is_boundary(SurfaceElementRelation idir, size_t np, size_t ix, size_t iy, size_t iz);
    static bool is_boundary(SurfaceElementRelation idir, size_t np, size_t ix, size_t iy);
    static bool is_boundary(SurfaceElementRelation idir, size_t np, size_t ix);

    static SurfaceElementRelation upwind_source_direction_x(Quadrant quad);
    static SurfaceElementRelation upwind_source_direction_y(Quadrant quad);
    static SurfaceElementRelation upwind_source_direction_z(Quadrant quad);
    static SurfaceElementRelation upwind_target_direction_x(Quadrant quad);
    static SurfaceElementRelation upwind_target_direction_y(Quadrant quad);
    static SurfaceElementRelation upwind_target_direction_z(Quadrant quad);

    static std::string toString(Quadrant quad);
    static std::string toString(SurfaceElementRelation direction);
    static std::string toString(SurfaceDirection sd);
    static size_t id(SurfaceElementRelation r);
    static size_t id(SurfaceDirection r);
  };
  

}


