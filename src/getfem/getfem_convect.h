// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2009-2009 Yves Renard
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
// As a special exception, you  may use  this file  as it is a part of a free
// software  library  without  restriction.  Specifically,  if   other  files
// instantiate  templates  or  use macros or inline functions from this file,
// or  you compile this  file  and  link  it  with other files  to produce an
// executable, this file  does  not  by itself cause the resulting executable
// to be covered  by the GNU Lesser General Public License.  This   exception
// does not  however  invalidate  any  other  reasons why the executable file
// might be covered by the GNU Lesser General Public License.
//
//===========================================================================

/**@file getfem_convect.h
   @author Yves Renard <Yves.Renard@insa-lyon.fr>,
   @date October 27, 2009.
   @brief Compute the convection of a quantity with respect to a vector field.
*/
#ifndef GETFEM_CONVECT_H__
#define GETFEM_CONVECT_H__

#include "getfem_mesh_fem.h"
#include "getfem_interpolation.h"
#include "gmm/gmm_dense_qr.h"

namespace getfem {

  // Specific interpolation
  template<typename VECTU, typename VECTV>
  void interpolation_convect(const mesh_fem &mf_source,
			     mesh_trans_inv &mti, const VECTU &UU, VECTV &V) {

    typedef typename gmm::linalg_traits<VECTU>::value_type T;
    const mesh &msh(mf_source.linked_mesh());
    dim_type qdim_s = mf_source.get_qdim();
    dim_type qqdim = dim_type(gmm::vect_size(UU)/mf_source.nb_dof());

    std::vector<T> U(mf_source.nb_basic_dof()*qqdim);
    mf_source.extend_vector(UU, U);

    mti.distribute(2);
    std::vector<size_type> itab;    
    base_matrix G;

    /* interpolation */
    dal::bit_vector dof_done; dof_done.add(0, mti.nb_points());
    std::vector<T> val(qdim_s);
    std::vector<std::vector<T> > coeff;
    base_tensor Z;
    std::vector<size_type> dof_source;

    for (dal::bv_visitor cv(mf_source.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt = msh.trans_of_convex(cv);
      mti.points_on_convex(cv, itab);
      if (itab.size() == 0) continue;

      pfem pf_s = mf_source.fem_of_element(cv);
      if (pf_s->need_G()) 
	bgeot::vectors_to_base_matrix(G, msh.points_of_convex(cv));

      fem_interpolation_context ctx(pgt, pf_s, base_node(), G, cv,
				    size_type(-1));
      coeff.resize(qqdim);
      for (size_type qq=0; qq < qqdim; ++qq) {
	coeff[qq].resize(mf_source.nb_basic_dof_of_element(cv));
	gmm::sub_index SUBI(mf_source.ind_basic_dof_of_element(cv));
	gmm::copy(gmm::sub_vector(U, SUBI), coeff[qq]);
      }
      for (size_type i = 0; i < itab.size(); ++i) {
	size_type dof_t = itab[i];
	if (dof_done.is_in(dof_t)) {
	  dof_done.sup(dof_t);
	  ctx.set_xref(mti.reference_coords()[dof_t]);
	  size_type pos = dof_t * qdim_s;
	  for (size_type qq=0; qq < qqdim; ++qq) {           
	    pf_s->interpolation(ctx, coeff[qq], val, qdim_s);
	    for (size_type k=0; k < qdim_s; ++k)
	      V[(pos + k)*qqdim+qq] = val[k];
	  }
	}
      }
    }
    if (dof_done.card() != 0)
      GMM_WARNING2("WARNING : in interpolation (different meshes),"
		   << dof_done.card() << " dof of target mesh_fem have "
		   << " been missed\nmissing dofs : " << dof_done);
    
  }




  /** Compute the convection of a quantity on a getfem::mesh_fem with respect
      to a velocity field
      @param mf the source mesh_fem. Should be of Lagrange type.
      @param U the source field.
      @param mf_v the mesh_fem on which the vector field is described
      @param V contains the vector field described on mf_v.
      @param nt number of time integration step.
  */
  template<class VECT1, class VECT2>
  void convect(const mesh_fem &mf, VECT1 &U, const mesh_fem &mf_v,
	       const VECT2 &V, scalar_type dt, size_type nt) {
    // Should be robustified on the boundaries -> control of the nodes going
    //   out the mesh.
    // Should control that the point do not move to fast ( < h/2 for instance).
    // Could be extended to non-lagragian fem with a projection (moving the 
    //   gauss points).
    // Can take into account a source term by integration on the
    //   caracteristics.

    typedef typename gmm::linalg_traits<VECT1>::value_type T;
    
    if (nt == 0) return;

    GMM_ASSERT1(!(mf.is_reduced()),
		"This convection algorithm work only on pure Lagrange fems");
    /* test if mf is really of Lagrange type.         */
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished();++cv) {
      pfem pf_t = mf.fem_of_element(cv);
      GMM_ASSERT1(pf_t->target_dim() == 1 && pf_t->is_lagrange(),
		  "This convection algorithm work only on pure Lagrange fems");
    }

    // Get the nodes of mf
    const mesh &msh(mf.linked_mesh());
    size_type N = msh.dim();
    getfem::mesh_trans_inv mti(msh, 1E-10);
    size_type qdim = mf.get_qdim();
    size_type nbpts = mf.nb_basic_dof() / qdim;
    std::vector<base_node> nodes(nbpts);
    for (size_type i = 0; i < nbpts; ++i)
      nodes[i] = mf.point_of_basic_dof(i * qdim);

    // Convect the nodes with respect to v
    size_type qqdimt = (gmm::vect_size(V) / mf_v.nb_dof()) * mf_v.get_qdim();
    GMM_ASSERT1(qqdimt == N, "The velocity field should be a vector field "
	       "of the smae dimension as the mesh");
    std::vector<T> VI(nbpts*N);

    scalar_type ddt = dt / scalar_type(nt);
    for (size_type i = 0; i < nt; ++i) {
      mti.clear();
      mti.add_points(nodes);
      interpolation_convect(mf_v, mti, V, VI);

      for (size_type j = 0; j < nbpts; ++j) {
	gmm::add(gmm::scaled(gmm::sub_vector(VI, gmm::sub_interval(N*j, N)),
			     -ddt), nodes[j]);
      }
    }

    // 3 interpolation finale
    std::vector<T> UI(nbpts*qdim);
    mti.clear();
    mti.add_points(nodes);
    interpolation_convect(mf, mti, U, UI);
    gmm::copy(UI, U);
  }


}


#endif
