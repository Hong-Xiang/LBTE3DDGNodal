#include "utilities.h"
namespace discontinues_galerkin_nodal_solver {
  Quadrant Utilities::quadrant_of_angle(num_t mu, num_t xi, num_t eta)
  {
    if (mu >= 0.0 && xi >= 0.0 && eta >= 0.0)
      return Quadrant::Q1;
    if (mu < 0.0 && xi >= 0.0 && eta >= 0.0)
      return Quadrant::Q2;
    if (mu < 0.0 && xi < 0.0 && eta >= 0.0)
      return Quadrant::Q3;
    if (mu >= 0.0 && xi < 0.0 && eta >= 0.0)
      return Quadrant::Q4;
    if (mu >= 0.0 && xi >= 0.0 && eta < 0.0)
      return Quadrant::Q5;
    if (mu < 0.0 && xi >= 0.0 && eta < 0.0)
      return Quadrant::Q6;
    if (mu < 0.0 && xi < 0.0 && eta < 0.0)
      return Quadrant::Q7;
    if (mu >= 0.0 && xi < 0.0 && eta < 0.0)
      return Quadrant::Q8;
    return Quadrant::Q8;
  }

  bool Utilities::quadrant_x_positive_flag(Quadrant quad)
  {
    return (quad == Quadrant::Q1 ||
      quad == Quadrant::Q4 ||
      quad == Quadrant::Q5 ||
      quad == Quadrant::Q8);
  }
  bool Utilities::quadrant_y_positive_flag(Quadrant quad)
  {
    return (quad == Quadrant::Q1 ||
      quad == Quadrant::Q2 ||
      quad == Quadrant::Q5 ||
      quad == Quadrant::Q6);
  }
  bool Utilities::quadrant_z_positive_flag(Quadrant quad)
  {
    return (quad == Quadrant::Q1 ||
      quad == Quadrant::Q2 ||
      quad == Quadrant::Q3 ||
      quad == Quadrant::Q4);
  }

  bool Utilities::is_boundary(SurfaceElementRelation idir, size_t np, size_t ix, size_t iy, size_t iz)
  {
    bool ans = false;
    if (idir == SurfaceElementRelation::B && ix == 0)
      ans = true;
    if (idir == SurfaceElementRelation::F && ix == np - 1)
      ans = true;
    if (idir == SurfaceElementRelation::L && iy == 0)
      ans = true;
    if (idir == SurfaceElementRelation::R && iy == np - 1)
      ans = true;
    if (idir == SurfaceElementRelation::D && iz == 0)
      ans = true;
    if (idir == SurfaceElementRelation::U && iz == np - 1)
      ans = true;
    return ans;
  }

  bool Utilities::is_boundary(SurfaceElementRelation idir, size_t np, size_t ix, size_t iy)
  {
    bool ans = false;
    if (idir == SurfaceElementRelation::B && ix == 0)
      ans = true;
    if (idir == SurfaceElementRelation::F && ix == np - 1)
      ans = true;
    if (idir == SurfaceElementRelation::L && iy == 0)
      ans = true;
    if (idir == SurfaceElementRelation::R && iy == np - 1)
      ans = true;
    return ans;
  }

  bool Utilities::is_boundary(SurfaceElementRelation idir, size_t np, size_t ix)
  {
    bool ans = false;
    if (idir == SurfaceElementRelation::B && ix == 0)
      ans = true;
    if (idir == SurfaceElementRelation::F && ix == np - 1)
      ans = true;
    return ans;
  }



  SurfaceElementRelation Utilities::upwind_source_direction_x(Quadrant quad)
  {
    if (quad == Quadrant::Q1 || quad == Quadrant::Q4 || 
        quad == Quadrant::Q5 || quad == Quadrant::Q8)
      return SurfaceElementRelation::B;
    else
      return SurfaceElementRelation::F;     
  }



  SurfaceElementRelation Utilities::upwind_source_direction_y(Quadrant quad)
  {
    if (quad == Quadrant::Q1 || quad == Quadrant::Q2 || 
        quad == Quadrant::Q5 || quad == Quadrant::Q6)
      return SurfaceElementRelation::L;
    else
      return SurfaceElementRelation::R;     
  }

  SurfaceElementRelation Utilities::upwind_source_direction_z(Quadrant quad)
  {
    if (quad == Quadrant::Q1 || quad == Quadrant::Q2 || 
        quad == Quadrant::Q3 || quad == Quadrant::Q4)
      return SurfaceElementRelation::D;
    else
      return SurfaceElementRelation::U;    
  }

  SurfaceElementRelation Utilities::upwind_target_direction_x(Quadrant quad)
  {
    if (quad == Quadrant::Q1 || quad == Quadrant::Q4 ||
        quad == Quadrant::Q5 || quad == Quadrant::Q8)
      return SurfaceElementRelation::F;
    else
      return SurfaceElementRelation::B;
  }

  SurfaceElementRelation Utilities::upwind_target_direction_y(Quadrant quad)
  {
    if (quad == Quadrant::Q1 || quad == Quadrant::Q2 ||
        quad == Quadrant::Q5 || quad == Quadrant::Q6)
      return SurfaceElementRelation::R;
    else
      return SurfaceElementRelation::L;
  }

  SurfaceElementRelation Utilities::upwind_target_direction_z(Quadrant quad)
  {
    if (quad == Quadrant::Q1 || quad == Quadrant::Q2 ||
        quad == Quadrant::Q3 || quad == Quadrant::Q4)
      return SurfaceElementRelation::U;
    else
      return SurfaceElementRelation::D;
  }

  std::string Utilities::toString(Quadrant quad)
  {
    return kQuadrantNameList.at(static_cast<size_t>(quad));
  }

  std::string Utilities::toString(SurfaceElementRelation direction)
  {
    return kSurfaceElementRelationNameList.at(static_cast<size_t>(direction));
  }

  std::string Utilities::toString(SurfaceDirection sd)
  {
    std::string ans;
    switch (sd)
    {
    case discontinues_galerkin_nodal_solver::SurfaceDirection::X:
      ans = "X";
      break;
    case discontinues_galerkin_nodal_solver::SurfaceDirection::Y:
      ans = "Y";
      break;
    case discontinues_galerkin_nodal_solver::SurfaceDirection::Z:
      ans = "Z";
      break;
    default:
      break;
    }
    return ans;
  }



  size_t Utilities::id(SurfaceElementRelation r)
  {
    return static_cast<size_t>(r);
  }

  size_t Utilities::id(SurfaceDirection r)
  {
    return static_cast<size_t>(r);
  }

}