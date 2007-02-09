// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2001-2006 Yves Renard, Julien Pommier.
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


#ifndef GETFEMINT_PGT_H__
#define GETFEMINT_PGT_H__

#include <getfemint_std.h>
#include <getfem/bgeot_geometric_trans.h>

namespace getfemint
{
  id_type ind_pgt(bgeot::pgeometric_trans);
  bgeot::pgeometric_trans addr_pgt(id_type);
  bool exists_pgt(id_type);
}  /* end of namespace getfemint.                                          */

#endif /* GETFEMINT_PGT_H__                                               */
