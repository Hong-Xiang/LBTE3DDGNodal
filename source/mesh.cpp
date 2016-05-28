#include <stack>
#include <iomanip>
#include "mesh.h"
#include "system_matrices.h"

namespace discontinues_galerkin_nodal_solver {

  Mesh::Mesh(size_t nx, size_t ny, size_t nz, num_t h, num_t xc, num_t yc, num_t zc)
    : nx_(nx), ny_(ny), nz_(nz), h_(h), xc_(xc), yc_(yc), zc_(zc)
  {
    element_list_.clear();
    surface_list_.clear();
    sweep_orders_.clear();
    sweep_orders_.resize(kQuadrantTotal.at(3));
    Generate();
    Calculate();
  }

  void Mesh::Generate()
  {
    GenerateCoarse();
  }


  void Mesh::Calculate()
  {
    CacheCenterElement();
    CacheCenterSurfaces();
    CacheSweepOrder();
    CacheBoundaryMark();
    CacheElementSurfaceMatrix();
  }



  void Mesh::CacheElementSurfaceMatrix()
  {
    auto sz = boost::extents[element_list_.size()][kSurfaceTotal.at(3)];
    element_surface_adjacent_matrix_.resize(sz);
    boost::multi_array<size_t, 2>& esam = element_surface_adjacent_matrix_;
    for (auto el : element_list_)
      for (auto r : kSurfaceElementRelationList)
        esam[el.id()][Utilities::id(r)] = el.surface_adjacent(r);
  }

  std::vector<size_t> Mesh::RefineLevel() const
  {
    std::vector<size_t> ans;
    ans.clear();
    ans.resize(n_elements());
    for (auto el : element_list_)
      ans.at(el.id()) = el.refine_level();
    return ans;
  }


  size_t Mesh::AddSurface(size_t refine_level, SurfaceDirection dir, num_t xc, num_t yc, num_t zc, num_t h)
  {
    size_t id = surface_list_.size();
    surface_list_.push_back(Surface(id, refine_level, dir, xc, yc, zc, h));
    return id;
  }

  size_t Mesh::AddElement(size_t refine_level, num_t xc, num_t yc, num_t zc, num_t h)
  {
    size_t id = element_list_.size();
    element_list_.push_back(Element(id, refine_level, xc, yc, zc, h));
    return id;
  }

  void Mesh::GenerateCoarse()
  {
    size_t total_elements = nx_*ny_*nz_;
    size_t total_interfaces = total_elements * 3 + nx_*ny_ + ny_*nz_ + nx_*nz_;
    surface_list_.reserve(total_interfaces);
    boost::multi_array<size_t, 3> index3;
    index3.resize(boost::extents[nx_][ny_][nz_]);

    //Generate Elements	

    for (size_t i = 0; i < total_elements; i++)
    {
      size_t ix, iy, iz;
      xlib::Utilities::ind2sub(nx_, ny_, nz_, i, ix, iy, iz);
      num_t xct, yct, zct;
      xct = ix*h_ - (nx_ - 1) * h_ / 2.0 + xc_;
      yct = iy*h_ - (ny_ - 1) * h_ / 2.0 + yc_;
      zct = iz*h_ - (nz_ - 1) * h_ / 2.0 + zc_;
      index3[ix][iy][iz] = AddElement(0, xct, yct, zct, h_);
    }

    //Generate surfaces

    size_t ids;
    for (size_t iz = 0; iz < nz_; iz++)
      for (size_t iy = 0; iy < ny_; iy++)
        for (size_t ix = 0; ix < nx_; ix++)
        {
          size_t ide = index3[ix][iy][iz];
          num_t xce = element_list_.at(ide).xc();
          num_t yce = element_list_.at(ide).yc();
          num_t zce = element_list_.at(ide).zc();
          num_t he = element_list_.at(ide).h();
          ids = AddSurface(0, SurfaceDirection::X, xce - he / 2.0, yce, zce, he);
          element_list_.at(ide).LinkSurface(SurfaceElementRelation::B, surface_list_.at(ids).id());
          surface_list_.at(ids).LinkElement(ElementDirection::inc, element_list_.at(ide).id());
          ids = AddSurface(0, SurfaceDirection::Y, xce, yce - he / 2.0, zce, he);
          element_list_.at(ide).LinkSurface(SurfaceElementRelation::L, surface_list_.at(ids).id());
          surface_list_.at(ids).LinkElement(ElementDirection::inc, element_list_.at(ide).id());
          ids = AddSurface(0, SurfaceDirection::Z, xce, yce, zce - he / 2.0, he);
          element_list_.at(ide).LinkSurface(SurfaceElementRelation::D, surface_list_.at(ids).id());
          surface_list_.at(ids).LinkElement(ElementDirection::inc, element_list_.at(ide).id());
        }

    size_t ide = 0;
    size_t idf = 0;
    for (size_t ix = 0; ix < nx_; ix++)
      for (size_t iy = 0; iy < ny_; iy++)
        for (size_t iz = 0; iz < nz_; iz++)
        {
          ide = index3[ix][iy][iz];
          if (ix == nx_ - 1)
          {
            num_t xce = element_list_.at(ide).xc();
            num_t yce = element_list_.at(ide).yc();
            num_t zce = element_list_.at(ide).zc();
            num_t he = element_list_.at(ide).h();
            size_t tid = AddSurface(0, SurfaceDirection::X, xce + he / 2.0, yce, zce, he);
            element_list_.at(ide).LinkSurface(SurfaceElementRelation::F, surface_list_.at(tid).id());
            surface_list_.at(tid).LinkElement(ElementDirection::pre, element_list_.at(ide).id());
          }
          else
          {
            idf = index3[ix + 1][iy][iz];
            element_list_.at(ide).LinkSurface(SurfaceElementRelation::F, element_list_.at(idf).surface_adjacent(SurfaceElementRelation::B));
            surface_list_.at(element_list_.at(idf).surface_adjacent(SurfaceElementRelation::B)).LinkElement(ElementDirection::pre, element_list_.at(ide).id());
          }
          if (iy == ny_ - 1)
          {
            num_t xce = element_list_.at(ide).xc();
            num_t yce = element_list_.at(ide).yc();
            num_t zce = element_list_.at(ide).zc();
            num_t he = element_list_.at(ide).h();
            size_t tid = AddSurface(0, SurfaceDirection::Y, xce, yce + he / 2.0, zce, he);
            element_list_.at(ide).LinkSurface(SurfaceElementRelation::R, surface_list_.at(tid).id());
            surface_list_.at(tid).LinkElement(ElementDirection::pre, element_list_.at(ide).id());
          }
          else
          {
            idf = index3[ix][iy + 1][iz];
            element_list_.at(ide).LinkSurface(SurfaceElementRelation::R, element_list_.at(idf).surface_adjacent(SurfaceElementRelation::L));
            surface_list_.at(element_list_.at(idf).surface_adjacent(SurfaceElementRelation::L)).LinkElement(ElementDirection::pre, element_list_.at(ide).id());
          }
          if (iz == (nz_ - 1))
          {
            num_t xce = element_list_.at(ide).xc();
            num_t yce = element_list_.at(ide).yc();
            num_t zce = element_list_.at(ide).zc();
            num_t he = element_list_.at(ide).h();
            size_t tid = AddSurface(0, SurfaceDirection::Z, xce, yce, zce + he / 2.0, he);
            element_list_.at(ide).LinkSurface(SurfaceElementRelation::U, surface_list_.at(tid).id());
            surface_list_.at(tid).LinkElement(ElementDirection::pre, element_list_.at(ide).id());
          }
          else
          {
            idf = index3[ix][iy][iz + 1];
            element_list_.at(ide).LinkSurface(SurfaceElementRelation::U, element_list_.at(idf).surface_adjacent(SurfaceElementRelation::D));
            surface_list_.at(element_list_.at(idf).surface_adjacent(SurfaceElementRelation::D)).LinkElement(ElementDirection::pre, element_list_.at(ide).id());
          }
        }

    //Mark boundaries

    for (size_t ix = 0; ix < nx_; ix++)
    {
      for (size_t iy = 0; iy < ny_; iy++)
      {
        surface_list_.at(element_list_.at(index3[ix][iy][0]).surface_adjacent(SurfaceElementRelation::D)).mark_as_boundary();
        surface_list_.at(element_list_.at(index3[ix][iy][nz_ - 1]).surface_adjacent(SurfaceElementRelation::U)).mark_as_boundary();
      }
    }
    for (size_t ix = 0; ix < nx_; ix++)
    {
      for (size_t iz = 0; iz < nz_; iz++)
      {
        surface_list_.at(element_list_.at(index3[ix][0][iz]).surface_adjacent(SurfaceElementRelation::L)).mark_as_boundary();
        surface_list_.at(element_list_.at(index3[ix][ny_ - 1][iz]).surface_adjacent(SurfaceElementRelation::R)).mark_as_boundary();
      }
    }
    for (size_t iz = 0; iz < nz_; iz++)
    {
      for (size_t iy = 0; iy < ny_; iy++)
      {
        surface_list_.at(element_list_.at(index3[0][iy][iz]).surface_adjacent(SurfaceElementRelation::B)).mark_as_boundary();
        surface_list_.at(element_list_.at(index3[nx_ - 1][iy][iz]).surface_adjacent(SurfaceElementRelation::F)).mark_as_boundary();
      }
    }

  }

  void Mesh::CacheCenterElement() {
    center_elements_.clear();
    center_elements_.resize(n_elements());
    for (auto el : element_list_)
      center_elements_.at(el.id()) = Vector3{ el.xc(), el.yc(), el.zc() };
  }

  void Mesh::CacheCenterSurfaces() {
    center_surfaces_.clear();
    center_surfaces_.resize(n_surfaces());
    surface_types_.clear();
    surface_types_.resize(n_surfaces());
    for (auto sf : surface_list_)
    {
      center_surfaces_.at(sf.id()) = Vector3{ sf.xc(), sf.yc(), sf.zc() };
      surface_types_.at(sf.id()) = sf.direction();
    }
  }

  void Mesh::CacheBoundaryMark()
  {
    boundary_marks_.clear();
    boundary_marks_.resize(n_surfaces());    
    for (auto sf : surface_list_)
      boundary_marks_.at(sf.id()) = sf.is_boundary();
  }

  void Mesh::sweep_order_finder(size_t cid, std::vector<size_t>& sweep_order,
    std::vector<bool>& added,
    const std::vector<SurfaceElementRelation>& search_direction,
    const std::vector<bool>& revert_sweep) const
  {
    const Element& e_now = element_list_.at(cid);
    for (size_t i = 0; i < 3; ++i)
    {
      size_t surface_id;
      surface_id = e_now.surface_adjacent(search_direction.at(i));
      const Surface& sf_now = surface_list_.at(surface_id);
      if (sf_now.is_boundary())
        continue;
      size_t element_next_id = revert_sweep.at(i) ? sf_now.element_inc()
        : sf_now.element_pre();
      if (added.at(element_next_id))
        continue;
      else
        sweep_order_finder(element_next_id, sweep_order,
          added, search_direction, revert_sweep);
    }
    if (!added.at(cid))
    {
      sweep_order.push_back(cid);
      added.at(cid) = true;
    }
  }

  void Mesh::CacheSweepOrder()
  {

    for (size_t i = 0; i < kQuadrantTotal.at(3); i++)
    {
      Quadrant quad = kQuadrantList.at(i);
      size_t iquad = static_cast<size_t>(quad);

      std::vector<SurfaceElementRelation> search_directions;
      search_directions.clear();
      search_directions.resize(kSurfaceTotal.at(3) / 2);
      std::vector<bool> revert_sweep;
      revert_sweep.clear();
      revert_sweep.resize(kSurfaceTotal.at(3) / 2);

      search_directions.at(0) = Utilities::upwind_source_direction_x(quad);
      revert_sweep.at(0) = !Utilities::quadrant_x_positive_flag(quad);
      search_directions.at(1) = Utilities::upwind_source_direction_y(quad);
      revert_sweep.at(1) = !Utilities::quadrant_y_positive_flag(quad);
      search_directions.at(2) = Utilities::upwind_source_direction_z(quad);
      revert_sweep.at(2) = !Utilities::quadrant_z_positive_flag(quad);

      std::vector<size_t> ans;
      ans.clear();

      std::vector<bool> added;
      added.resize(n_elements());

      for (auto flag : added) flag = false;

      for (auto el : element_list_)
        sweep_order_finder(el.id(), ans, added, search_directions, revert_sweep);

      //Using stack to find all elements
      //////std::stack<size_t> stk;
      //////while (!stk.empty())
      //////  stk.pop();
      //
      //////for (size_t i = 0; i < n_elements(); i++)
      //////{
      //////  if (!added.at(i))
      //////    stk.push(i);
      //
      //////  while (!stk.empty())
      //////  {
      //////    size_t e_now = stk.top();
      //////    bool is_bd_el = true; //is boundary Element
      //////    for (size_t j = 0; j < 3; j++)
      //////    {
      //////      size_t surface_search = element_list_.at(e_now).surface_adjacent(search_directions.at(j));
      //////      const Surface& s_ref = surface_list_.at(surface_search);
      //////      if (s_ref.is_boundary())
      //////        continue;
      //////      if (!revert_sweep.at(j)) {
      //////        for (size_t k = 0; k < s_ref.pre_element_total(); k++)
      //////        {
      //////          size_t element_search = s_ref.element_pre(k);
      //////          if (added.at(element_search))
      //////            continue;
      //////          is_bd_el = false;
      //////          stk.push(element_search);
      //////        }
      //////      }
      //////      else {
      //////        for (size_t k = 0; k < s_ref.inc_element_total(); k++)
      //////        {
      //////          size_t element_search = s_ref.element_inc(k);
      //////          if (added.at(element_search))
      //////            continue;
      //////          is_bd_el = false;
      //////          stk.push(element_search);
      //////        }
      //////      }
      //////    }
      //
      //////    if (is_bd_el)
      //////    {
      //////      ans.push_back(e_now);
      //////      added.at(e_now) = true;
      //////      stk.pop();
      //////    }
      //
      //////  }
      //////}

      sweep_orders_.at(iquad) = ans;
    }
  }


  std::ostream & operator<<(std::ostream & os, const Mesh & m)
  {
    os << std::setprecision(5);
    os << "information of Mesh" << std::endl;
    os << "total elements = " << m.n_elements() << std::endl;
    os << "total interface = " << m.n_surfaces() << std::endl;
    os << "========   elements info:  ================" << std::endl;
    os << "id\txc\tyc\tzc" << std::endl;
    std::vector<Vector3> ec = m.CenterElements();
    for (size_t i = 0; i < m.n_elements(); i++)
    {
      os << i << "\t" << ec.at(i).x() << "\t" << ec.at(i).y() << "\t" << ec.at(i).z() << std::endl;
    }
    os << "===========================================" << std::endl;

    std::vector<Vector3> sc = m.CenterSurfaces();
    os << "========   Surface info:   ================" << std::endl;
    os << "id\tbd_flag\tdiret\txc\tyc\tzc" << std::endl;
    std::vector<bool> bdf = m.BoundaryMask();
    for (size_t i = 0; i < m.n_surfaces(); i++)
    {
      os << i << "\t" << bdf.at(i) << "\t" << Utilities::toString(m.surface_ref(i).direction()) << "\t" << sc.at(i).x() << "\t" << sc.at(i).y() << "\t" << sc.at(i).z() << std::endl;
    }
    os << "===========================================" << std::endl;

    os << "========   Element Surface link:   ========" << std::endl;
    boost::multi_array<size_t, 2> esam = m.ElementSurfaceMatrix();
    os << "id\tB\tF\tL\tR\tD\tU" << std::endl;
    for (size_t i = 0; i < m.n_elements(); i++)
    {
      os << i;
      for (size_t j = 0; j < kSurfaceTotal.at(3); j++)
      {
        os << "\t" << esam[i][j];
      }
      os << std::endl;
    }
    os << "===========================================" << std::endl;

    os << "========== sweep order ====================" << std::endl;
    for (auto quad : kQuadrantList)
    {
      std::vector<size_t> order = m.SweepOrder(quad);
      os << Utilities::toString(quad) << std::endl;
      for (size_t i = 0; i < order.size(); i++)
      {
        os << order.at(i) << "\t";
      }
      os << std::endl;
    }
    os << "===========================================" << std::endl;

    return os;
  }



  MeshWithCoordinate::MeshWithCoordinate(size_t nx, size_t ny, size_t nz, num_t h,
    num_t xc, num_t yc, num_t zc,
    const Vector& xe, const Vector& ye, const Vector& ze,
    const std::vector<Vector>& xs, const std::vector<Vector>& ys,
    const std::vector<Vector>& zs)
    : Mesh(nx, ny, nz, h, xc, yc, zc),
    n_basis_element_(xe.row()), n_basis_surface_(xs.at(0).row())
  {
    n_node_element_ = n_elements()*n_basis_element();
    n_node_surface_ = n_surfaces()*n_basis_surface();
    xe_.New(n_node_element());
    ye_.New(n_node_element());
    ze_.New(n_node_element());
    xs_.New(n_node_surface());
    ys_.New(n_node_surface());
    zs_.New(n_node_surface());
    boundary_mask_.New(n_node_surface());
    boundary_mask_iv_.New(n_node_surface());
    CalculateElementCoordinates(xe, ye, ze);
    CalculateSurfaceCoordinates(xs, ys, zs);
  }

  MeshWithCoordinate::MeshWithCoordinate(size_t nx, size_t ny, size_t nz, num_t h, num_t xc, num_t yc, num_t zc, const SystemMatrix & system)
    : Mesh{ nx, ny, nz, h, xc, yc, zc },
    n_basis_element_(system.n_basis_element()), n_basis_surface_(system.n_basis_surface())
  {
    n_node_element_ = n_elements()*n_basis_element();
    n_node_surface_ = n_surfaces()*n_basis_surface();
    xe_.New(n_node_element());
    ye_.New(n_node_element());
    ze_.New(n_node_element());
    xs_.New(n_node_surface());
    ys_.New(n_node_surface());
    zs_.New(n_node_surface());
    boundary_mask_.New(n_node_surface());
    boundary_mask_iv_.New(n_node_surface());
    std::vector<Vector> xs, ys, zs;
    xs.resize(3); ys.resize(3), zs.resize(3);
    for (size_t i = 0; i < 3; ++i)
    {
      xs.at(i) = system.xs(kSurfaceDirectionList.at(i));
      ys.at(i) = system.ys(kSurfaceDirectionList.at(i));
      zs.at(i) = system.zs(kSurfaceDirectionList.at(i));
    }
    CalculateElementCoordinates(system.xe(), system.ye(), system.ze());
    CalculateSurfaceCoordinates(xs, ys, zs);
  }

  void MeshWithCoordinate::CalculateElementCoordinates(const Vector& xei,
    const Vector& yei, const Vector& zei)
  {    
    std::vector<Vector3> ec =CenterElements();
    for (size_t i = 0; i < n_elements(); ++i)
    {
      Vector xt, yt, zt;
      xt.Bind(n_basis_element_, xe().p() + i*n_basis_element_);
      yt.Bind(n_basis_element_, ye().p() + i*n_basis_element_);
      zt.Bind(n_basis_element_, ze().p() + i*n_basis_element_);
      xt.ScalarMulVectorAddScalar(1.0, xei, ec.at(i).x());
      yt.ScalarMulVectorAddScalar(1.0, yei, ec.at(i).y());
      zt.ScalarMulVectorAddScalar(1.0, zei, ec.at(i).z());
    }
  }

  void MeshWithCoordinate::CalculateSurfaceCoordinates(const std::vector<Vector>& xsi, 
    const std::vector<Vector>& ysi, const std::vector<Vector>& zsi)
  {
    size_t n_basis = n_basis_surface();
    std::vector<Vector3> ec = CenterSurfaces();
    std::vector<bool> bdf = BoundaryMask();
    std::vector<SurfaceDirection> sd = SurfaceType();
    for (size_t i = 0; i < n_surfaces(); ++i)
    {
      Vector xt, yt, zt, bdm, bdmiv;
      xt.Bind(n_basis, xs().p() + i*n_basis);
      yt.Bind(n_basis, ys().p() + i*n_basis);
      zt.Bind(n_basis, zs().p() + i*n_basis);
      bdm.Bind(n_basis, boundary_mask_.p() + i*n_basis);
      bdmiv.Bind(n_basis, boundary_mask_iv_.p() + i*n_basis);
      SurfaceDirection sdv = sd.at(i);
      xt.ScalarMulVectorAddScalar(1.0, xsi.at(Utilities::id(sdv)), ec.at(i).x());
      yt.ScalarMulVectorAddScalar(1.0, ysi.at(Utilities::id(sdv)), ec.at(i).y());
      zt.ScalarMulVectorAddScalar(1.0, zsi.at(Utilities::id(sdv)), ec.at(i).z());
      if (bdf.at(i))
      {
        bdm.Ones();
        bdmiv.Zeros();
      }
      else
      {
        bdm.Zeros();
        bdmiv.Ones();
      }
    }
  }

}