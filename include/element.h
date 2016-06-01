#pragma once
#include "global.h"

namespace discontinues_galerkin_nodal_solver {
  class Surface;
  class Element {
  public:
    Element(size_t eId, size_t refine_level, num_t xc, num_t yc, num_t zc, num_t h);

    num_t xc() const { return xc_; }
    num_t yc() const { return yc_; }
    num_t zc() const { return zc_; }
    num_t h() const { return h_; }
    size_t id() const { return id_; }

    void LinkSurface(SurfaceElementRelation direction, size_t sid, size_t local_id = 0);

    //Get Surface of one direction
    const size_t surface_adjacent(SurfaceElementRelation dir) const;

  public:
    bool is_to_refine() const {
      return refine_mark_;
    }
    void mark_to_refine(bool flag = true) {
      refine_mark_ = flag;
    }
    const size_t refine_level() const {
      return refine_level_;
    }
    bool is_active() const { return true; }

  public:
    size_t id_;
    size_t refine_level_;
    num_t xc_, yc_, zc_, h_;   

    std::vector<size_t> surface_list_;
    std::vector<size_t> child_elements_;
    
    bool refine_mark_;
    bool is_active_;
  };
}