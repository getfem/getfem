// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Dynamic Array Library (dal)
// File    : dal_backtrace.h : For debugging information
//           
// Date    : June 01, 2003.
// Authors : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1995-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================

/**@file dal_backtrace.h
   @brief Get debug information.
*/
#ifndef DAL_BACKTRACE
#define DAL_BACKTRACE

#include <getfem_arch_config.h>
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
