#include "surface.h"

namespace discontinues_galerkin_nodal_solver {

  Surface::Surface(size_t id, size_t refine_level, SurfaceDirection direction,
    num_t xc, num_t yc, num_t zc, num_t h)
    : id_(id), refine_level_(refine_level), direction_(direction),
      xc_(xc), yc_(yc), zc_(zc), h_(h), is_boundary_(false)
  {}

  void Surface::LinkElement(const ElementDirection ed, const size_t eid, 
      const size_t local_id)
  {
    if (ed == ElementDirection::pre)
      element_pre_list_.push_back(eid);
    else
      element_inc_list_.push_back(eid);
  }

}