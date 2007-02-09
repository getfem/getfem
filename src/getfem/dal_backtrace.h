// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 1995-2007 Yves Renard
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

/**@file dal_backtrace.h
   @author Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date June 01, 2003.
   @brief Get debug information.
*/
#ifndef DAL_BACKTRACE
#define DAL_BACKTRACE

#include "getfem/getfem_arch_config.h"
#include <stdio.h>
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
