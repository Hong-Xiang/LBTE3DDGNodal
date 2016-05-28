#pragma once
#include <vector>
#include <boost/multi_array.hpp>
#include "global.h"
#include "element.h"
#include "surface.h"
#include <memory>


namespace discontinues_galerkin_nodal_solver {

  /*!
  * \file Mesh.h
  * \date 2016/05/17 1:53
  *
  * \author HongXiang
  * Contact: hx.hongxiang@gmail.com
  *
  * \brief Mesh information of 3d Cartesian mesh.
  *		  Include: center of elements: xc(), yc() zc()
  *
  * Mesh is designed to include only static information to a fixed mesh, e.g. no element basis specific information is included.
  * Mesh include information to apply upwind sweep. Thus the element_interface matrix.
  * It can be combined with different elements.
  *
  *
  *
  * \note
  */
  class Mesh {
  public:
    Mesh(size_t nx, size_t ny, size_t nz, num_t h, num_t xc, num_t yc, num_t zc);
    void Generate();
    void Calculate();

  protected:
    size_t AddSurface(size_t refine_level, SurfaceDirection dir, num_t xc, num_t yc, num_t zc, num_t h);
    size_t AddElement(size_t refine_level, num_t xc, num_t yc, num_t zc, num_t h);

    void CacheSweepOrder();


  public:
    //Report functions:

    //For SystemMatrix
    num_t h_element(size_t id) const
    {
      return element_list_.at(id).h();
    }
    num_t h_surface(size_t id) const
    {
      return surface_list_.at(id).h();
    }

    //Static, cache in Mesh.
    size_t n_elements() const {
      return element_list_.size();
    }
    size_t n_surfaces() const {
      return surface_list_.size();
    }

    const std::vector<Vector3>& CenterElements() const {
      return center_elements_;
    }
    const std::vector<Vector3>& CenterSurfaces() const {
      return center_surfaces_;
    }

    const boost::multi_array<size_t, 2>& ElementSurfaceMatrix() const {
      return element_surface_adjacent_matrix_;
    }
    const std::vector<size_t>& SweepOrder(Quadrant quad) const
    {
      return sweep_orders_.at(static_cast<size_t>(quad));
    }
    const std::vector<bool>& BoundaryMask() const { return boundary_marks_; }
    const std::vector<SurfaceDirection>& SurfaceType() const { return surface_types_; }

    
    //Debug functions
    const Element& element_ref(size_t id) const { return element_list_.at(id); }
    const Surface& surface_ref(size_t id) const { return surface_list_.at(id); }
    //SurfaceDirection surface_direction(size_t id) const;

  private:
    void GenerateCoarse();
    void CacheCenterElement();
    void CacheCenterSurfaces();
    void CacheBoundaryMark();
    void CacheElementSurfaceMatrix();    
    //void mark_refine(num_t x0, num_t x1, num_t y0, num_t y1, num_t z0, num_t z1);
    //void generate_refine();

  private:
    size_t nx_, ny_, nz_;
    num_t h_, xc_, yc_, zc_;
    std::vector<Element> element_list_;
    std::vector<Surface> surface_list_;
    std::vector<std::vector<size_t>> sweep_orders_;
    boost::multi_array<size_t, 2> element_surface_adjacent_matrix_;
    std::vector<bool> boundary_marks_;
    std::vector<Vector3> center_elements_;
    std::vector<Vector3> center_surfaces_;
    std::vector<SurfaceDirection> surface_types_;
    //Refine related:
  public:
    std::vector<size_t> RefineLevel() const;

  private:
    //helper functions
    void sweep_order_finder(size_t cid, std::vector<size_t>& sweep_order,
      std::vector<bool>& added, 
      const std::vector<SurfaceElementRelation>& search_direction, 
      const std::vector<bool>& revert_sweep) const;

    
  };


  class SystemMatrix;
  class MeshWithCoordinate : public Mesh {
  public: 
    MeshWithCoordinate(size_t nx, size_t ny, size_t nz, num_t h, 
        num_t xc, num_t yc, num_t zc, 
        const Vector& xe, const Vector& ye, const Vector& ze,
        const std::vector<Vector>& xs, const std::vector<Vector>& ys, 
        const std::vector<Vector>& zs );
    MeshWithCoordinate(size_t nx, size_t ny, size_t nz, num_t h,
      num_t xc, num_t yc, num_t zc, const SystemMatrix& system);
    Vector xe() const { return xe_; }
    Vector ye() const { return ye_; }
    Vector ze() const { return ze_; }
    Vector xs() const { return xs_; }
    Vector ys() const { return ys_; }
    Vector zs() const { return zs_; }
    Vector BoundaryMaskNodal() const { return boundary_mask_; }
    Vector BoundaryMaskIvNodal() const { return boundary_mask_iv_; }
    size_t n_node_element() const { return n_node_element_; }
    size_t n_node_surface() const { return n_node_surface_; }
    size_t n_basis_element() const { return n_basis_element_; }
    size_t n_basis_surface() const { return n_basis_surface_; }
  private:
    void CalculateElementCoordinates(const Vector& xe, 
        const Vector& ye, const Vector& ze);
    void CalculateSurfaceCoordinates(const std::vector<Vector>& xs, 
        const std::vector<Vector>& ys, const std::vector<Vector>& zs);
  private:
    Vector xe_, ye_, ze_, xs_, ys_, zs_;
    Vector boundary_mask_;
    Vector boundary_mask_iv_;
    size_t n_node_element_;
    size_t n_node_surface_;
    size_t n_basis_element_;
    size_t n_basis_surface_;
  };
  std::ostream& operator<<(std::ostream& os, const Mesh& m);
}