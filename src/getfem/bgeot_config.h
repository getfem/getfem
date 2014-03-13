/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 1999-2012 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/

/**@file bgeot_config.h
   @brief defines and typedefs for namespace bgeot
   @date December 20, 1999
   @author Yves Renard
*/
#ifndef BGEOT_CONFIG_H__
#define BGEOT_CONFIG_H__

#include "getfem/getfem_arch_config.h"

#ifdef GETFEM_HAVE_FEENABLEEXCEPT
# include <fenv.h>
# define FE_ENABLE_EXCEPT { feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW); }
#else
# define FE_ENABLE_EXCEPT {}
#endif

#include "dal_config.h"
#include "gmm/gmm_kernel.h"
#include "gmm/gmm_dense_lu.h"

#ifdef GETFEM_HAVE_BOOST
# include <boost/noncopyable.hpp>
#else
# include <getfem_boost/noncopyable.hpp> 
#endif

#ifdef GETFEM_HAVE_QDLIB
// #  define NO_INLINE
#  ifdef GETFEM_QDLIB_USE_QUAD
#    include <qd/qd_real.h>
#  else
#    include <qd/dd_real.h>
#  endif
#  include <qd/fpu.h>
#endif

/// Basic Geometric Tools
namespace bgeot {

  using std::endl; using std::cout; using std::cerr;
  using std::ends; using std::cin;
  

  static const size_t ST_NIL = size_t(-1);
  typedef gmm::uint16_type dim_type;
  typedef gmm::uint16_type short_type;
  typedef size_t size_type;
  typedef double scalar_type;
  typedef std::complex<double> complex_type;
  inline double to_double(double &a) { return a; }
  inline scalar_type to_scalar(const scalar_type &a) { return a; }

#ifndef GETFEM_HAVE_QDLIB
  typedef double long_scalar_type;
  typedef double opt_long_scalar_type;
# define LONG_SCALAR_ATOF(st) (atof(st))
# define LONG_SCALAR_EPS 1E-16
# define LONG_SCAL(xx) long_scalar_type(xx)
#else
#  ifdef GETFEM_QDLIB_USE_QUAD
  typedef qd_real long_scalar_type;
  typedef qd_real opt_long_scalar_type;
  inline scalar_type to_scalar(const qd_real &a) { return to_double(a); }
# define LONG_SCALAR_ATOF(st) (long_scalar_type(st))
# define LONG_SCALAR_EPS 1E-64
#  else
  typedef dd_real long_scalar_type;
  typedef dd_real opt_long_scalar_type;
  inline scalar_type to_scalar(const dd_real &a) { return to_double(a); }
# define LONG_SCALAR_ATOF(st) (long_scalar_type(st))
# define LONG_SCALAR_EPS 1E-32
#  endif
#  define LONG_SCAL(xx) long_scalar_type(#xx) /* string assignment to preserve the precision */
#endif



  // For compatibility with Getfem 2.0

  using gmm::dimension_error;
  using gmm::file_not_found_error;
  using gmm::internal_error;
  using gmm::to_be_done_error;
  using gmm::failure_error;

}  /* end of namespace bgeot.                                             */


#endif /* BGEOT_CONFIG_H__                                                */
