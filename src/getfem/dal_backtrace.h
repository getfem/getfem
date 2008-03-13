// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 1995-2008 Yves Renard
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

/**@file dal_backtrace.h
   @author Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date June 01, 2003.
   @brief Get debug information.
*/
#ifndef DAL_BACKTRACE
#define DAL_BACKTRACE

#include "getfem/getfem_arch_config.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>


namespace dal {
  std::string demangle(const char *mangled_name);
  inline void dump_backtrace() {
#ifdef GETFEM_HAVE_BACKTRACE
    void dump_glibc_backtrace();
    dump_glibc_backtrace();
#endif
  }
}

#endif
