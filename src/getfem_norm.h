/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_norm.h: computes various norms on pde solutions.      */
/*     									   */
/* Date : November 17, 2000.                                               */
/* Authors : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                     */
/*           Julien Pommier, pommier@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2000-2002  Yves Renard.                                   */
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
 

#ifndef GETFEM_NORM_H__
#define GETFEM_NORM_H__

#include <getfem_mesh_fem.h>
#include <getfem_mesh_slicers.h>
#include <bgeot_geotrans_inv.h>
#include <getfem_export.h>
#include <bgeot_rtree.h>

namespace getfem
{
  enum { L2_NORM=1, H1_SEMI_NORM=2, H1_NORM=3, LINF_NORM=4 };

  template<typename VECT1, typename VECT2>
  class slicer_compute_norm_diff : public slicer_mesh_with_mesh {
    const mesh_fem &mf1, &mf2;
    const VECT1& U1;
    const VECT2& U2;
    papprox_integration im;
    int what; /* combiantion of L2_NORM, H1_SEMI_NORM .. */    
    bgeot::pgeotrans_precomp gp;
    bgeot::geotrans_inv_convex gti;
    scalar_type maxd;
  public:
    scalar_type l2_norm_sqr, h1_semi_norm_sqr, linf_norm;
    slicer_compute_norm_diff(const mesh_fem& mf1_, const VECT1 &U1_, 
			     const mesh_fem& mf2_, const VECT2 &U2_, 
			     papprox_integration im_, int what_=L2_NORM) 
      : slicer_mesh_with_mesh(mf2_.linked_mesh()), mf1(mf1_), mf2(mf2_), U1(U1_), U2(U2_), 
	im(im_), what(what_), gp(0) {
	l2_norm_sqr = 0.; h1_semi_norm_sqr = 0.; linf_norm = 0.; maxd=0.;
      }
    scalar_type norm(int w) { 
      scalar_type s = 0;
      if (w & L2_NORM) s += sqrt(l2_norm_sqr);
      if (w & H1_SEMI_NORM) s += sqrt(h1_semi_norm_sqr);
      if (w & LINF_NORM) s+= linf_norm;
      return s;
    }
    /*virtual*/ void intersection_callback(mesh_slicer &ms, size_type slmcv) {
      if (ms.splx_in.card() == 0) return;
      const getfem_mesh &m1 = mf1.linked_mesh();
      const getfem_mesh &m2 = mf2.linked_mesh();
      bgeot::pgeometric_trans pgt1 = m1.trans_of_convex(ms.cv);
      bgeot::pgeometric_trans pgt2 = m2.trans_of_convex(slmcv);
      pfem pf1 = mf1.fem_of_element(ms.cv); 
      pfem pf2 = mf2.fem_of_element(slmcv);
      base_matrix G1, G2;
      size_type qdim = mf1.get_qdim(), mdim = m1.dim();
      base_vector val1(qdim), val2(qdim);
      base_matrix gval1(qdim,mdim), gval2(qdim,mdim);

      vectors_to_base_matrix(G1,m1.points_of_convex(ms.cv));
      vectors_to_base_matrix(G2,m2.points_of_convex(slmcv));
      fem_interpolation_context ctx1(pgt1,pf1,base_node(),G1,ms.cv);
      fem_interpolation_context ctx2(pgt2,pf2,base_node(),G2,slmcv);

      /* coordinates on the ref convex of each mesh */
      std::vector<base_node> nodes1(im->nb_points_on_convex()), 
	nodes2(im->nb_points_on_convex());
      gti.init(m2.convex(slmcv), pgt2);
      //cout << "hello, doing cv " << ms.cv << " <-> " << slmcv << " [" << ms.splx_in.card() << " simplexes]\n";
      for (dal::bv_visitor is(ms.splx_in); !is.finished(); ++is) {
	const slice_simplex &s = ms.simplexes[is];
	if (gp == 0 || s.dim() != im->dim()) 
	  gp = bgeot::geotrans_precomp(bgeot::simplex_geotrans(s.dim(),1), &im->integration_points());
	//DAL_THROW(dal::dimension_error, "incompatible dimension of the slice");
	base_matrix M(s.dim(),s.dim());

	for (size_type i=0; i < s.dim(); ++i) 
	  for (size_type j=0; j < s.dim(); ++j)
	    M(i,j) = ms.nodes[s.inodes[i+1]].pt[j] - ms.nodes[s.inodes[0]].pt[j];
	
	scalar_type J = dal::abs(gmm::lu_det(M));
	
	/* build nodes1 and nodes2 */
	std::vector<base_node> ptref(s.dim()+1);

	for (size_type k=0; k <= s.dim(); ++k) 
	  gti.invert(ms.nodes[s.inodes[k]].pt, ptref[k]);
	base_matrix G; vectors_to_base_matrix(G,ptref);
	for (size_type k=0; k < nodes1.size(); ++k) 
	  nodes2[k] = gp->transform(k,G);

	for (size_type k=0; k <= s.dim(); ++k) 
	  ptref[k] = ms.nodes[s.inodes[k]].pt_ref;	
	vectors_to_base_matrix(G,ptref);
	for (size_type k=0; k < nodes1.size(); ++k) 
	  nodes1[k] = gp->transform(k,G);
	
	base_vector coeff1(mf1.nb_dof_of_element(ms.cv)), coeff2(mf2.nb_dof_of_element(slmcv));
	gmm::copy(gmm::sub_vector(U1,gmm::sub_index(mf1.ind_dof_of_element(ms.cv))),coeff1);
	gmm::copy(gmm::sub_vector(U2,gmm::sub_index(mf2.ind_dof_of_element(slmcv))),coeff2);

	for (size_type i=0; i < im->nb_points_on_convex(); ++i) {
	  ctx1.set_xref(nodes1[i]); ctx2.set_xref(nodes2[i]);
	  if (what & (L2_NORM | LINF_NORM)) {
	    pf1->interpolation(ctx1, coeff1, val1, qdim);
	    pf2->interpolation(ctx2, coeff2, val2, qdim);
	    for (size_type q=0; q < qdim; ++q) {
	      l2_norm_sqr += dal::sqr(val1[q]-val2[q])*J*im->coeff(i);
	      linf_norm = std::max(linf_norm, dal::abs(val1[q]-val2[q]));
	    }
	  }
	  if (what & H1_SEMI_NORM) {
	    pf1->interpolation_grad(ctx1, coeff1, gval1, qdim);
	    pf2->interpolation_grad(ctx2, coeff2, gval2, qdim);
	    for (size_type q=0; q < qdim*mdim; ++q) {
	      scalar_type v = dal::sqr(gval1[q]-gval2[q])*J*im->coeff(i);
	      if (v > maxd) { cout << "new maxd: " << v << ", J=" << J << " at " << ctx1.xreal() << ", " << ctx2.xreal() << ", v1=" << gval1[q] << ", v2=" << gval2[q] << "\n"; maxd = v; }
	      h1_semi_norm_sqr += dal::sqr(gval1[q]-gval2[q])*J*im->coeff(i);
	    }
	  }
	}
      }
    }
  };

  template<typename VECT1, typename VECT2>
  void solutions_distance(const mesh_fem& mf1, const VECT1& U1, 
			  const mesh_fem& mf2, const VECT2& U2,
			  papprox_integration im, scalar_type *pl2=0, scalar_type *psh1=0) {
    mesh_slicer slicer(mf1.linked_mesh()); 
    slicer_compute_norm_diff<VECT1,VECT2> 
      cn(mf1,U1,mf2,U2,im, (pl2 ? L2_NORM : 0) +
	 (psh1 ? H1_SEMI_NORM : 0));
    slicer.push_back_action(cn);
    slicer.exec();
    if (pl2) *pl2 = cn.norm(L2_NORM);
    if (psh1) *psh1 = cn.norm(H1_SEMI_NORM);
  }
}
#endif 
