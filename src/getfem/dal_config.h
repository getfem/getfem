/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2007-2015 Yves Renard

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

  using std::endl; using std::cout; using std::cerr;
  using std::ends; using std::cin;

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

}  /* end of namespace bgeot.                                             */


#endif /* DAL_CONFIG_H__                                                  */
