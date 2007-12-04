// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2007-2008 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.
//
//===========================================================================

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

  
  inline void DAL_THROW() IS_DEPRECATED;
  inline void DAL_THROW() {}

#define DAL_STANDARD_CATCH_ERROR GMM_STANDARD_CATCH_ERROR
#define DAL_THROW(a, b) { GMM_THROW(a, b); dal::DAL_THROW(); }
#define DAL_WARNING0 GMM_WARNING0
#define DAL_WARNING1 GMM_WARNING1
#define DAL_WARNING2 GMM_WARNING2
#define DAL_WARNING3 GMM_WARNING3
#define DAL_WARNING4 GMM_WARNING4
#define DAL_PRETTY_FUNCTION GMM_PRETTY_FUNCTION
#define DAL_SET_EXCEPTION_DEBUG GMM_SET_EXCEPTION_DEBUG

}  /* end of namespace bgeot.                                             */


#endif /* DAL_CONFIG_H__                                                  */
