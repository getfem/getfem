// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2007-2007 Yves Renard
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

/**@file dal_config.h
   @brief defines and typedefs for namespace dal
   @date January 30, 2007
   @author Yves Renard
*/
#ifndef DAL_CONFIG_H__
#define DAL_CONFIG_H__

#include "gmm/gmm_except.h"
#include "gmm/gmm_algobase.h"

namespace dal {

  // For compatibility with Getfem 2.0

  using gmm::dimension_error;
  using gmm::file_not_found_error;
  using gmm::internal_error;
  using gmm::to_be_done_error;
  using gmm::failure_error;

  using gmm::int8_type;
  using gmm::uint8_type;
  using gmm::int16_type;
  using gmm::uint16_type;
  using gmm::int32_type;
  using gmm::uint32_type;
  using gmm::int64_type;
  using gmm::uint64_type;

  using gmm::lexicographical_less;
  using gmm::approx_less;
  using gmm::uclock_sec;

  
  // inline void DAL_THROW() IS_DEPRECATED;
  inline void DAL_THROW() {}
  inline void DAL_INTERNAL_ERROR() IS_DEPRECATED;
  inline void DAL_INTERNAL_ERROR() {}


#define DAL_STANDARD_CATCH_ERROR GMM_STANDARD_CATCH_ERROR
#define DAL_THROW(a, b) { GMM_THROW(a, b); dal::DAL_THROW(); }
#define DAL_INTERNAL_ERROR(a)				\
  { GMM_INTERNAL_ERROR(a); dal::DAL_INTERNAL_ERROR(); }
#define DAL_WARNING0 GMM_WARNING0
#define DAL_WARNING1 GMM_WARNING1
#define DAL_WARNING2 GMM_WARNING2
#define DAL_WARNING3 GMM_WARNING3
#define DAL_WARNING4 GMM_WARNING4
#define DAL_PRETTY_FUNCTION GMM_PRETTY_FUNCTION
#define DAL_SET_EXCEPTION_DEBUG GMM_SET_EXCEPTION_DEBUG

}  /* end of namespace bgeot.                                             */


#endif /* DAL_CONFIG_H__                                                  */
