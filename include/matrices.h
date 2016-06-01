#pragma once
#include "global.h"

namespace discontinues_galerkin_nodal_solver {
  class MatricesData {
  public:
    static const MatricesData& instance();
    const Vector x(size_t np) const;
    const Matrix m(size_t np) const;
    const Matrix s(size_t np) const;
    const size_t min_np() const { return kMinNp; }
    const size_t max_np() const { return kMaxNp; }
  private:
    MatricesData();
    void read_MAT_files();
    static MatricesData* handle_;
    std::vector<Vector> x_;
    std::vector<Matrix> s_;
    std::vector<Matrix> m_;

  private:
    const size_t kMinNp = 2;
    const size_t kMaxNp = 20;
  };

  std::ostream& operator<<(std::ostream& os, const MatricesData& m);
    
  
  class MatricesBase {
  public:

    MatricesBase(num_t h, size_t np);
    void reset(num_t h, size_t np);
    virtual const size_t n_basis() const = 0;
    virtual const size_t n_basis_lower() const = 0;
    const num_t h() const { return h_; }
    const size_t np() const { return np_; }
    const Matrix lift_matrix(SurfaceElementRelation idr) const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization());
      return lift_matrices_.at(static_cast<size_t>(idr)); 
    }
    const Matrix flux_matrix(SurfaceElementRelation idr) const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization()); 
      return flux_matrices_.at(static_cast<size_t>(idr)); 
    }    
    void Initialize();

    virtual size_t dimension() const {
      BOOST_THROW_EXCEPTION(ExceptionOutOfRange() 
          << InfoOutOfRange(0, 1, kInfSize));
      return 0;
    }
    virtual ~MatricesBase();
  protected:
    virtual void Calculate() = 0;
    virtual void MemoryAlloc();

  protected:
    num_t h_;
    size_t np_;

    std::vector<Matrix> lift_matrices_;
    std::vector<Matrix> flux_matrices_;

    size_t dim_;
    bool is_init_;
  };

  class Matrices1D : public MatricesBase {
  public:
    static size_t n_basis(size_t np) {
      return np;
    }
    static size_t n_basis_lower(size_t np) {
      return 1;
    }
  public:

    Matrices1D(num_t h, size_t np);
    const Vector x() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization());
      return x_; 
    }
    const Matrix m() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization());
      return m_; 
    }
    const Matrix sx() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization());
      return sx_; 
    }
    virtual const size_t n_basis() const { return n_basis(np_); }
    virtual const size_t n_basis_lower() const { return n_basis_lower(np_); }
    virtual size_t dimension() const { return 1; }
    virtual ~Matrices1D();
  protected:
    virtual void Calculate();
    virtual void MemoryAlloc();

  private:
    Vector x_;
    Matrix m_;
    Matrix sx_;
  };

  std::ostream& operator<<(std::ostream& os, const Matrices1D& m);

  class Matrices2D : public MatricesBase {
  public:
    static size_t n_basis(size_t np) {
      return np*np;
    }
    static size_t n_basis_lower(size_t np) {
      return np;
    }
  public:
    Matrices2D(num_t h, size_t np);
    const size_t n_basis() const { return np_*np_; }
    virtual const size_t n_basis_lower() const { return np_; }
    virtual size_t dimension() const { return 2; }
    virtual ~Matrices2D();
  protected:
    virtual void Calculate();
    virtual void MemoryAlloc();
  public:
    const Vector x() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization()); 
      return x_; 
    }
    const Vector y() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization()); 
      return y_; 
    }
    const Matrix m() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization()); 
      return m_; 
    }
    const Matrix sx() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization()); 
      return sx_; 
    }
    const Matrix sy() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization()); 
      return sy_; 
    }
  private:

    Vector x_;
    Vector y_;
    Matrix m_;
    Matrix sx_;
    Matrix sy_;
  };

  std::ostream& operator<<(std::ostream& os, const Matrices2D& m);

  class Matrices3D : public MatricesBase {
  public:
    static size_t n_basis(size_t np) {
      return np*np*np;
    }
    static size_t n_basis_lower(size_t np) {
      return np*np;
    }
  public:
    Matrices3D(num_t h, size_t np);
    virtual const size_t n_basis() const { return np_*np_*np_; }
    virtual const size_t n_basis_lower() const { return np_*np_; }
    virtual size_t dimension() const { return 3; }
    virtual ~Matrices3D();
  protected:
    virtual void Calculate();
    virtual void MemoryAlloc();
  public:
    const Vector x() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization()); 
      return x_; }
    const Vector y() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization()); 
      return y_; }
    const Vector z() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization()); 
      return z_; }
    const Matrix m() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization()); 
      return m_; }
    const Matrix sx() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization()); 
      return sx_; }
    const Matrix sy() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization()); 
      return sy_; }
    const Matrix sz() const { 
      if (!is_init_) BOOST_THROW_EXCEPTION(ExceptionRequireInitialization()); 
      return sz_; }
  private:
    Vector x_;
    Vector y_;
    Vector z_;
    Matrix m_;
    Matrix sx_;
    Matrix sy_;
    Matrix sz_;
  };

  std::ostream& operator<<(std::ostream& os, const Matrices3D& m);
}