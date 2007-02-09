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


#ifndef GETFEMINT_MATELEM_H__
#define GETFEMINT_MATELEM_H__

#include <getfemint_std.h>
#include <getfem/getfem_mat_elem.h>


namespace getfemint
{
  id_type ind_matelem(getfem::pmat_elem_computation);
  getfem::pmat_elem_computation addr_matelem(id_type);
  bool exists_matelem(id_type);
  
}  /* end of namespace getfemint.                                          */

#endif /* GETFEMINT_MATELEM_H__                                            */
