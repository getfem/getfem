// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 1999-2007 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

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
#include <getfem_boost/noncopyable.hpp>

#ifdef GETFEM_HAVE_QDLIB
// #  define NO_INLINE
#  ifdef GETFEM_QDLIB_USE_QUAD
#    include <qd/qd.h>
#  else
#    include <qd/dd.h>
#  endif
#  include <qd/fpu.h>
#endif

/// Basic Geometric Tools
namespace bgeot {

  static const size_t ST_NIL = size_t(-1);
  /// Dimension type (<255)
  typedef gmm::uint8_type dim_type;
  typedef gmm::uint16_type short_type;
  typedef size_t size_type;
  typedef double scalar_type;
  typedef std::complex<double> complex_type;
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
# define LONG_SCALAR_ATOF(st) (long_scalar_type(st))
# define LONG_SCALAR_EPS 1E-64
#  else
  typedef dd_real long_scalar_type;
  typedef dd_real opt_long_scalar_type;
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
