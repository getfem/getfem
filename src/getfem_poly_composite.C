/* *********************************************************************** */
/*                                                                         */
/* Library : GEneric Tool for Finite Element Methods (getfem)              */
/* File    : getfem_poly_composite.C : polynomials by parts                */
/*                                                                         */
/* Date : August 26, 2002.                                                 */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; version 2 of the License.                 */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software Foundation, */
/* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.         */
/*                                                                         */
/* *********************************************************************** */


#include <getfem_poly_composite.h>

namespace getfem
{ 
  mesh_precomposite::mesh_precomposite(const getfem_mesh &m) {
    mesh = &m;
    elt.resize(m.nb_convex());
    orgs.resize(m.nb_convex());
    gtrans.resize(m.nb_convex());
    dal::bit_vector nn = m.points().index();
    size_type i;
    for (i = 0; i <= nn.last_true(); ++i) {
      vertexes.add(m.points()[i]);
    }
    nn = m.convex_index();
    for (i << nn; i != size_type(-1); i << nn) {
      
      bgeot::pgeometric_trans pgt = m.trans_of_convex(i);
      size_type N = pgt->structure()->dim();
      size_type P = m.dim();
      if (!(pgt->is_linear()) || N != P) 
	DAL_THROW(internal_error, "Bad geometric transformations.");
      
      base_poly PO;
      base_matrix a(P, pgt->nb_points());
      base_matrix pc(pgt->nb_points() , N);
      base_matrix grad(P, N), TMP1(N,N), B0(N, P);
      
      for (size_type j = 0; j < pgt->nb_points(); ++j)
	for (size_type k = 0; k < P; ++k)
	  a(i,j) = (m.points_of_convex(i)[j])[i];
      
      for (size_type k = 0; k < pgt->nb_points(); ++k)
	for (dim_type n = 0; n < N; ++n)
	  { PO = pgt->poly_vector()[k]; PO.derivative(n); pc(k,n) = PO[0]; }
      bgeot::mat_product(a, pc, grad);
      bgeot::mat_gauss_inverse(grad, TMP1); B0 = grad;
      gtrans[i] = B0;
      orgs[i] = m.points_of_convex(i)[0];
    }
  }

  void polynomial_composite::derivative(short_type k) {
    std::vector<bgeot::base_poly>::iterator it = polytab.begin();
    std::vector<bgeot::base_poly>::iterator ite = polytab.end();
    for ( ; it != ite; ++it) it->derivative(k);
  }
  
}  /* end of namespace getfem.                                            */
