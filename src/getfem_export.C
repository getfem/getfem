/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_export.C : export and import data.                    */
/*     									   */
/* Date : October 15, 2001.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001-2002  Yves Renard.                                   */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */


#include <getfem_export.h>

namespace getfem
{
  void classical_mesh_fem(mesh_fem& mf, short_type K)
  {
    dal::bit_vector bv = mf.linked_mesh().convex_index();
    for (size_type cv = bv.take_first(); cv != size_type(-1); cv << bv) {
      bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
      mf.set_finite_element(cv, classical_fem(pgt,K), exact_classical_im(pgt));
    }
  }


}  /* end of namespace getfem.                                             */


