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
#include <getfem_mesh_slice.h>
#include <bgeot_geotrans_inv.h>
#include <getfem_export.h>
namespace getfem
{
  enum { L2_NORM=1, H1_SEMI_NORM=2 };

  template<typename VECT1, typename VECT2>
  scalar_type L2_or_H1_dist(const mesh_fem& mf1, const VECT1& U1, 
		      const mesh_fem& mf2, const VECT2& U2,
		      papprox_integration im, int nrefine=1) {
    int what = L2_NORM;
    const getfem_mesh &m1 = mf1.linked_mesh(), &m2 = mf2.linked_mesh();
    size_type mdim = m1.dim();
    size_type qdim = mf1.get_qdim();
    if (mdim != m2.dim()) DAL_THROW(dal::dimension_error,"");
    if (qdim != mf2.get_qdim()) DAL_THROW(dal::dimension_error,"different values of Qdim");
    
    mesh_slice sl(m1);
    cout << "L2_or_H1_dist building slice\n";
    if (&m1 != &m2) {
      slicer_mesh s(m2); sl.build(&s,nrefine);
    } else {
      slicer_none s; sl.build(&s,nrefine);
    }
    cout << "slice done\n";
    exit(1);
    std::vector<scalar_type> V1, V2, J;
    V1.resize(qdim * sl.nb_simplexes(im->dim()) * im->nb_points_on_convex());
    V2.resize(qdim * sl.nb_simplexes(im->dim()) * im->nb_points_on_convex());
    J.reserve(sl.nb_simplexes(im->dim()));
    bgeot::geotrans_inv gti;
    bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(im->dim(),1);
    for (dal::bv_visitor cv(m1.convex_index()); !cv.finished(); ++cv) {
      const mesh_slice::cs_nodes_ct& nodes = sl.nodes(cv);
      const mesh_slice::cs_simplexes_ct& splxs = sl.simplexes(cv);
      for (size_type scnt=0; scnt < splxs.size(); ++scnt) {
	const slice_simplex& s = splxs[scnt];
	if (s.dim() != im->dim()) continue;	
	base_matrix M(s.dim(),s.dim());
	for (size_type i=0; i < s.dim(); ++i) 
	  for (size_type j=0; j < s.dim(); ++j)
	    M(i,j) = nodes[s.inodes[i+1]].pt[j] - nodes[s.inodes[0]].pt[j];
	J.push_back(dal::abs(gmm::lu_det(M)));
	std::vector<base_node> pts(s.dim()+1);	
	for (size_type i=0; i < s.dim()+1; ++i) {
	  pts[i] = nodes[s.inodes[i]].pt;
	  //cout << "pts[" << i << "]=" << pts[i] << "\n";
	}
	for (size_type i=0; i < im->nb_points_on_convex(); ++i) {
	  base_node n = pgt->transform(im->point(i),pts.begin());
	  //cout << "n=" << n << " " << im->point(i) << "\n";
	  gti.add_point(n);
	}
      }
    }
    //cout << "U1.size()" <<U1.size() << "gti.nb_points()"<<gti.nb_points()<<"COEFF.size()"<<COEFF.size()<<"\n";
    scalar_type res = 0;
    assert(V1.size() == gti.nb_points()*qdim);
    cout << "L2_or_H1_dist interpolating\n";
    if (what == H1_SEMI_NORM) {
      std::vector<scalar_type> gV1(V1.size()*mdim*qdim);
      std::vector<scalar_type> gV2(V2.size()*mdim*qdim);
      interpolation_solution(mf1,gti,U1,V1,&gV1);
      interpolation_solution(mf2,gti,U2,V2,&gV2);
      for (size_type i=0, pos=0; i < J.size(); ++i) {
	for (size_type j=0; j < im->nb_points_on_convex(); ++j) {
	  for (size_type q=0; q < qdim*mdim; ++q, ++pos) {
	    res += dal::sqr(gV1[pos]-gV2[pos])*J[i]*im->coeff(j);
	  }
	}
      }
    /*      gmm::add(gmm::scaled(gV2,-1.),gV1);
      for (size_type i=0; i < COEFF.size(); ++i) {
	res += COEFF[i]*gmm::vect_norm2_sqr(
	       gmm::sub_vector(gV1,gmm::sub_interval(i*mdim,mdim))); 
      }
    */
    } else {
      interpolation_solution(mf1,gti,U1,V1);
      interpolation_solution(mf2,gti,U2,V2);
    }
    cout << "L2_or_H1_dist done\n";
    for (size_type i=0, pos=0; i < J.size(); ++i) {
      //cout << "J[" << i << "]=" << J[i] << ", res=" << res << ", V1[pos]=" << V1[pos] << ", V2[pos]=" << V2[pos] << "\n";
      for (size_type j=0; j < im->nb_points_on_convex(); ++j) {
	for (size_type q=0; q < qdim; ++q, ++pos) {
	  res += dal::sqr(V1[pos]-V2[pos])*J[i]*im->coeff(j);
	}
      }
    }
    return sqrt(res);
  }
}
#endif 
