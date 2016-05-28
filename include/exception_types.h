#pragma once
#include <boost/exception/all.hpp>
#include <exception>
namespace discontinues_galerkin_nodal_solver {

  typedef boost::error_info<struct tag_size_min, size_t> InfoMinSizeT;
  typedef boost::error_info<struct tag_size_max, size_t> InfoMaxSizeT;
  typedef boost::error_info<struct tag_size_wrong, size_t> InfoWrongSizeT;
  typedef boost::error_info<struct tag_size_valid, size_t> InfoValidSizeT;
  typedef boost::tuple<InfoWrongSizeT, InfoMinSizeT, InfoMaxSizeT> InfoOutOfRange;
  typedef boost::tuple<InfoWrongSizeT, InfoValidSizeT> InfoDimensionMismatch;
  
  class ExceptionBase : virtual public boost::exception,
    virtual public std::exception {};

  class ExceptionOutOfRange : virtual public ExceptionBase {};

  class ExceptionRequireInitialization : virtual public ExceptionBase {};

  class ExceptionUnknownTestId : virtual public ExceptionBase {};

  class ExceptionInvalidSolverStatus : virtual public ExceptionBase {};
}