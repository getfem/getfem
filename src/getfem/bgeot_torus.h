/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2014-2017 Liang Jin Lim

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
/**
@file bgeot_torus.h
@brief Provides mesh of torus.
@date May 2014
@author Liang Jin Lim
*/

#pragma once

#ifndef BGEOT_TORUS_H__
#define BGEOT_TORUS_H__

#include "getfem/bgeot_geometric_trans.h"

namespace bgeot{

/**An adaptor that adapts a two dimensional geometric_trans to include radial dimension.*/
struct torus_geom_trans : public geometric_trans{

  virtual void poly_vector_val(const base_node &, bgeot::base_vector &) const;
  virtual void poly_vector_val(const base_node &, const bgeot::convex_ind_ct &,
    bgeot::base_vector &) const;
  virtual void poly_vector_grad(const base_node &, bgeot::base_matrix &) const;
  inline virtual void poly_vector_grad(const base_node &,
    const bgeot::convex_ind_ct &, bgeot::base_matrix &) const;
  inline virtual void compute_K_matrix
    (const bgeot::base_matrix &, const bgeot::base_matrix &, bgeot::base_matrix &) const;

  virtual void poly_vector_hess(const base_node &, bgeot::base_matrix &) const;
  virtual void project_into_reference_convex(base_node &) const;

  torus_geom_trans(bgeot::pgeometric_trans poriginal_trans);

  pgeometric_trans get_original_transformation() const;

private:
  pgeometric_trans poriginal_trans_;
};

pconvex_structure torus_structure_descriptor(pconvex_structure ori_structure);

bool is_torus_structure(pconvex_structure cvs);

pgeometric_trans torus_geom_trans_descriptor(pgeometric_trans poriginal_trans);

bool is_torus_geom_trans(pgeometric_trans pgt);

}

#endif /* BGEOT_TORUS_H__  */
