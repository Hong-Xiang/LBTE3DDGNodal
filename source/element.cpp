#include "element.h"

namespace discontinues_galerkin_nodal_solver {
  Element::Element(size_t eId, size_t refine_level, num_t xc, num_t yc, num_t zc, num_t h)
    : id_(eId), refine_level_(refine_level), xc_(xc), yc_(yc), zc_(zc), h_(h)
  {
    surface_list_.resize(kSurfaceTotal.at(3));
    for each (auto sid in surface_list_)
      sid = kInfSize;

    child_elements_.clear();   
    refine_mark_ = false;
  }


  void Element::LinkSurface(const SurfaceElementRelation direction, const size_t sid, const size_t local_id)
  {
    surface_list_.at(static_cast<size_t>(direction)) = sid;
  }

  const size_t Element::surface_adjacent(SurfaceElementRelation dir) const
  {
    return surface_list_.at(static_cast<size_t>(dir));
  }
}