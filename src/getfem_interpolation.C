/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (Getfem)             */
/* File    :  getfem_interpolation.C                                       */
/*     									   */
/* Date : June 1, 2004.                                                    */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*          Julien Pommier, pommier@gmm.insa-tlse.fr                       */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001-2004  Yves Renard, Julien Pommier.                   */
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


#include <getfem_interpolation.h>

namespace getfem {

  void mesh_trans_inv::points_on_convex(size_type i,
					std::vector<size_type> &itab) const {
    itab.resize(pts_cvx[i].size()); size_type j = 0;
    for (map_iterator it = pts_cvx[i].begin(); it != pts_cvx[i].end(); ++it)
      itab[j++] = it->first;
  }
  
  void mesh_trans_inv::ventilate(bool projection) {
    size_type nbpts = nb_points();
    dal::bit_vector nn = mesh.convex_index(), npt;
    size_type nbcvx = nn.last_true() + 1;
    ref_coords.resize(nbpts); dist.resize(nbpts); cvx_pts.resize(nbpts);
    pts_cvx.resize(nbcvx);
    base_node min, max, pt_ref; /* bound of the box enclosing the convex */
    bgeot::kdtree_tab_type boxpts;


    for (size_type j = nn.take_first(); j != size_type(-1); j << nn) {
      pts_cvx[j].clear();
      bgeot::pgeometric_trans pgt = mesh.trans_of_convex(j);
      bounding_box(min, max, mesh.points_of_convex(j), pgt);
      for (size_type k=0; k < min.size(); ++k) { min[k]-=EPS; max[k]+=EPS; }
      gic.init(mesh.convex(j), pgt);
      points_in_box(boxpts, min, max);
      for (size_type l = 0; l < boxpts.size(); ++l) {
	bool gicisin = gic.invert(boxpts[l].n, pt_ref, EPS);
	bool toadd = projection || gicisin;
	double isin = pgt->convex_ref()->is_in(pt_ref);
	size_type ind = boxpts[l].i;
	if (toadd && npt[ind]) {
	  if (isin < dist[ind]) pts_cvx[cvx_pts[ind]].erase(ind);
	  else toadd = false;
	}
	if (toadd) {
	  ref_coords[ind] = pt_ref;
	  dist[ind] = isin;
	  cvx_pts[ind] = j;
	  pts_cvx[j][ind] = void_type();
	  npt[ind] = true;
	}
      }
    }
  }
  
}  /* end of namespace getfem.                                             */

