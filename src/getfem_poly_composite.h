// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_poly_composite.h : polynomials by parts
//           
// Date    : August 26, 2002.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2002-2005 Yves Renard
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

/**@file getfem_poly_composite.h
   @brief Handle composite polynomials.

   Composite polynomials are used in hierarchical FEM.
*/

#ifndef GETFEM_POLY_COMPOSITE_H__
#define GETFEM_POLY_COMPOSITE_H__

#include <bgeot_poly.h>
#include <bgeot_imbricated_box.h>
#include <getfem_mesh.h>

namespace getfem
{

  struct mesh_precomposite {

    typedef dal::dynamic_tree_sorted<base_node,
      bgeot::imbricated_box_less> PTAB;

    const mesh *msh;
    PTAB vertexes;
    std::vector<base_matrix> gtrans;
    std::vector<scalar_type> det;
    std::vector<base_node> orgs;
    mutable std::vector<bool> elt;
    
    const mesh &linked_mesh(void) const { return *msh; }
    size_type nb_convex(void) const { return gtrans.size(); }
    dim_type dim(void) const { return msh->dim(); }
    bgeot::pgeometric_trans trans_of_convex(size_type ic) const
    { return msh->trans_of_convex(ic); }
    
    mesh_precomposite(const mesh &m);
  };

  typedef const mesh_precomposite *pmesh_precomposite;

  class polynomial_composite {

  protected :
    const mesh_precomposite *mp;
    std::vector<bgeot::base_poly> polytab;

  public :
    
    template <class ITER> scalar_type eval(const ITER &it) const;
    scalar_type eval(const base_node &pt) const;
    void derivative(short_type k);
    base_poly &poly_of_subelt(size_type l) { return polytab[l]; }
    const base_poly &poly_of_subelt(size_type l) const { return polytab[l]; }
    

    polynomial_composite(void) {}
    polynomial_composite(const mesh_precomposite &m);

  };

  template <class ITER>
  scalar_type polynomial_composite::eval(const ITER &it) const {
    base_node pt(mp->dim());
    std::copy(it, it + mp->dim(), pt.begin());
    return eval(pt);
  }

  void structured_mesh_for_convex(bgeot::pconvex_ref cvr, short_type k,
				  pmesh &pm, pmesh_precomposite &pmp, bool force_simplexification=false);

  /** simplexify a convex_ref.
      @param cvr the convex_ref.
      @param k the refinement level.
      @return a pointer to a statically allocated mesh. Do no free it!
  */
  const mesh *
  refined_simplex_mesh_for_convex(bgeot::pconvex_ref cvr, short_type k);

  /** simplexify the faces of a convex_ref

      @param cvr the convex_ref.

      @param k the refinement level.

      @return vector of pointers to a statically allocated
      mesh_structure objects. Do no free them! The point numbers in
      the mesh_structure refer to the points of the mesh given by
      refined_simplex_mesh_for_convex.
  */      
  const std::vector<bgeot::mesh_structure*>&
  refined_simplex_mesh_for_convex_faces(bgeot::pconvex_ref cvr, short_type k);
}  /* end of namespace getfem.                                            */


#endif
