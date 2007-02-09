// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2006-2006 Yves Renard, Julien Pommier.
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

#ifndef GETFEMINT_GSPARSE_MISC_H
#define GETFEMINT_GSPARSE_MISC_H

#include <getfemint_gsparse.h>

namespace getfemint {
  void spmat_set_diag(gsparse &gsp, mexargs_in& in, bool create_matrix);
  void spmat_load(mexargs_in& in, mexargs_out& out, mexarg_out::output_sparse_fmt fmt);
}

#endif /* GETFEMINT_GSPARSE_MISC_H */
