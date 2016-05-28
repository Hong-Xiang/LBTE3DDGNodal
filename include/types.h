#pragma once
#include <array>
#include <xlib.h>

namespace discontinues_galerkin_nodal_solver {

  enum class DgnsStatues {
    kDgnsSucceed = 0,
    kDgnsFailed = 1
  };

#ifdef USE_DOUBLE_PRECISION
  typedef double num_t;
#else
  typedef float num_t;
#endif

  typedef num_t* num_p;
  typedef const num_t* const_num_p;
  typedef size_t* size_p;

  typedef xlib::mkl::Vector<num_t> Vector;
  typedef xlib::mkl::Matrix<num_t> Matrix;
  typedef xlib::mkl::MatrixLu<num_t> MatrixLu;

  using xlib::mkl::operator<<;
  using xlib::operator<<;

  class Vector3 {
  public:
    Vector3() {}
    Vector3(num_t x, num_t y, num_t z)
      : x_(x), y_(y), z_(z)
    {}
    num_t x() const { return x_; }
    num_t y() const { return y_; }
    num_t z() const { return z_; }
  private:
    num_t x_, y_, z_;
  };

  //quadrants and list for it
  enum class Quadrant {
    Q1 = 0,
    Q2 = 1,
    Q3 = 2,
    Q4 = 3,
    Q5 = 4,
    Q6 = 5,
    Q7 = 6,
    Q8 = 7,  
  };
  const std::array<Quadrant, 8> kQuadrantList = {
    Quadrant::Q1,
    Quadrant::Q2,
    Quadrant::Q3,
    Quadrant::Q4,
    Quadrant::Q5,
    Quadrant::Q6,
    Quadrant::Q7,
    Quadrant::Q8
  };
  const std::array<std::string, 8> kQuadrantNameList = {
    "quadrant1",
    "quadrant2",
    "quadrant3",
    "quadrant4",
    "quadrant5",
    "quadrant6",
    "quadrant7",
    "quadrant8"
  };

  //interface directions and list for it
  enum class SurfaceElementRelation {
    B = 0,
    F = 1,
    L = 2,
    R = 3,
    D = 4,
    U = 5
  };

  const std::array<SurfaceElementRelation, 6> kSurfaceElementRelationList = {
    SurfaceElementRelation::B,
    SurfaceElementRelation::F,
    SurfaceElementRelation::L,
    SurfaceElementRelation::R,
    SurfaceElementRelation::D,
    SurfaceElementRelation::U
  };

  const std::array<std::string, 6> kSurfaceElementRelationNameList = {
    "Back",
    "Front",
    "Left",
    "Right",
    "Down",
    "Up"
  };


  static const std::vector<size_t> kSurfaceTotal = { 1, 2, 4, 6 };
  static const std::vector<size_t> kQuadrantTotal = { 1, 2, 4, 8 };

  //normal directions 
  enum class SurfaceDirection {
    X = 0,
    Y = 1,
    Z = 2
  };

  const std::array<SurfaceDirection, 3> kSurfaceDirectionList = {
    SurfaceDirection::X,
    SurfaceDirection::Y,
    SurfaceDirection::Z
  };

  enum class ElementDirection {
    pre = 0,
    inc = 1
  };

}