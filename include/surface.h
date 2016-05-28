#pragma once
#include "global.h"

namespace discontinues_galerkin_nodal_solver {

  class Element;

  class Surface {
  public:
    Surface(size_t id, size_t refine_level, SurfaceDirection direction, 
        num_t xc, num_t yc, num_t zc, num_t h);

    num_t xc() const { return xc_; }
    num_t yc() const { return yc_; }
    num_t zc() const { return zc_; }
    num_t h() const { return h_; }
    size_t id() const { return id_; }

    const size_t element_pre(size_t local_id = 0) const {
      return element_pre_list_.at(local_id);
    }
    const size_t element_inc(size_t local_id = 0) const {
      return element_inc_list_.at(local_id);
    }

    SurfaceDirection direction() const {
      return direction_;
    }

    void LinkElement(const ElementDirection ed, const size_t eid, 
        const size_t local_id = 0);

    bool is_boundary() const { return is_boundary_; }
    void mark_as_boundary(bool flag = true) { is_boundary_ = flag; }

    size_t pre_element_total() const { return element_pre_list_.size(); }
    size_t inc_element_total() const { return element_inc_list_.size(); }

  private:
    size_t id_;
    
    
    std::vector<size_t> element_pre_list_;
    std::vector<size_t> element_inc_list_;
    num_t xc_, yc_, zc_, h_;
    SurfaceDirection direction_;
    bool is_boundary_;  

  //Refine related:
  public:
    const size_t refine_level() const {
      return refine_level_;
    }
    size_t refine_level_;
  
  };
}