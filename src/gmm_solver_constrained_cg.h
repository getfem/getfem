/* -*- c++ -*- (enables emacs c++ mode)                                    */
//=======================================================================
// Copyright (C) 1997-2001
// Authors: Andrew Lumsdaine <lums@osl.iu.edu> 
//          Lie-Quan Lee     <llee@osl.iu.edu>
//
// This file is part of the Iterative Template Library
//
// You should have received a copy of the License Agreement for the
// Iterative Template Library along with the software;  see the
// file LICENSE.  
//
// Permission to modify the code and to distribute modified code is
// granted, provided the text of this NOTICE is retained, a notice that
// the code was modified is included with the above COPYRIGHT NOTICE and
// with the COPYRIGHT NOTICE in the LICENSE file, and that the LICENSE
// file is distributed with the modified code.
//
// LICENSOR MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
// By way of example, but not limitation, Licensor MAKES NO
// REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY
// PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS
// OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS
// OR OTHER RIGHTS.
//=======================================================================
/* *********************************************************************** */
/*                                                                         */
/* Library : Generic Matrix Methods  (gmm)                                 */
/* File    : gmm_solver_constrained_cg.h : constrained conjugate, gradient */
/*           Modified version of I.T.L. conjugated gradient.               */
/*     									   */
/* Date : October 13, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
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

// a unifier avec cg

#ifndef __GMM_SOLVER_CCG_H
#define __GMM_SOLVER_CCG_H

namespace gmm {

template <class CMatrix, class CINVMatrix, class VectorX>
void pseudo_inverse(const CMatrix C, CINVMatrix CINV, VectorX&)
{
  // compute the pseudo inverse of the non-square matrix C such
  // CINV = inv(C * trans(C)) * C.
  // based on a conjugate gradient method.

  // optimisable : copie de la ligne, precalcul de C * trans(C).

  typedef VectorX TmpVec;
  typedef size_t size_type;
  typedef typename linalg_traits<VectorX>::value_type value_type;

  TmpVec d(mat_nrows(C)), e(mat_nrows(C)), l(mat_ncols(C));
  TmpVec p(mat_nrows(C)), q(mat_nrows(C)), r(mat_nrows(C));
  value_type rho, rho_1, alpha;
  clear(C);

  for (size_type i = 0; i < mat_nrows(C); ++i)
  {
    d[i] = 1.0; rho = 1.0;
    std::fill(e.begin(), e.end(), 0.0);
    copy(d, r); gmm::copy(d, p);
    
    while (rho >= 1E-38) /* conjugate gradient to compute e               */
    {                    /* which is the i nd row of inv(C * trans(C))    */
      mult(gmm::transposed(C), p, l);
      mult(C, l, q);	  
      alpha = rho / vect_sp(p, q);
      add(e, scaled(p, alpha), e);  
      add(r, scaled(q, -alpha), r); 
      rho_1 = rho;
      rho = vect_sp(r, r);
      add(r, scaled(p, rho / rho_1), p);
    }
    
    mult(transposed(C), e, l); /* l is the i nd row of CINV     */
    for (size_type j = 0; j < mat_ncols(C); ++j) /* copy of the row       */
      if (dal::abs(l[j]) >  1E-15) {
	CINV(i, j) = l[j];
	// std::cout << "i = " << i << " j = " << j << " : " << l[j] << endl;
      }

   d[i] = 0.0;
  }

/*   for (size_type i = 0; i < C.nrows(); i++) */
/*   {  */
/*     std::cout << "ligne " << i << " [ "; */
/*     for (size_type j = 0; j < C.nrows(); j++) */
/*     { */
/*       double al = mtl::dot(C[j], CINV[i]); */
/*       if (al != 0.0) */
/* 	std::cout << "(" << j << "," << al << ")  "; */
/*     } */
/*     std::cout << "]" << endl; */
/*   } */

}

template < class Matrix, class CMatrix, class VectorX, class VectorB, 
           class Preconditioner >
int constrained_cg(const Matrix& A, const CMatrix& C, VectorX& x,
		   const VectorB& b, const Preconditioner& M,
		   int itemax, double residu, int noisy)
{
  typedef typename temporary_plain_vector<VectorX>::vector_type TmpVec;
  typedef typename temporary_plain_vector<typename
    linalg_traits<CMatrix>::sub_row_type>::vector_type TmpCVec;
  typedef row_matrix<TmpCVec> TmpCmat;
  
  typedef size_t size_type;
  typedef typename linalg_traits<VectorX>::value_type value_type;
  value_type rho = 1.0, rho_1, lambda, gamma, norm_b = vect_norm2(b);
  TmpVec p(vect_size(x)), q(vect_size(x)), q2(vect_size(x)),
    r(vect_size(x)), old_z(vect_size(x)), z(vect_size(x)), memox(vect_size(x));
  std::vector<bool> satured(mat_nrows(C));
  clear(p);
  int iter = 0;

  TmpCmat CINV(mat_nrows(C), mat_ncols(C));
  clear(CINV);
  pseudo_inverse(C, CINV, x);

  while(true)
  {
    // computation of residu
    copy(z, old_z);
    copy(x, memox);
    mult(A, scaled(x, -1.0), b, r);
    mult(M, r, z); // ...
    bool transition = false;
    for (size_type i = 0; i < mat_nrows(C); ++i) {
      value_type al = vect_sp(mat_row(C, i), x);
      if (al >= -1.0E-15)
      {
	if (!satured[i]) { satured[i] = true; transition = true; }
	value_type bb = vect_sp(mat_row(CINV, i), z);
	if (bb > 0.0)
	  add(scaled(mat_row(C, i), -bb), z);
 
/* 	bb = itl::dot(mtl::rows(CINV)[i], r); */
/* 	if (bb > 0.0) // itl::add(r, itl::scaled(mtl::rows(C)[i], -bb), r); */
/* 	  add_vect_sparse__(r, mtl::rows(C)[i].begin(), */
/* 			       mtl::rows(C)[i].end(), -bb); */
      }
      else
	satured[i] = false;
    }
    
    // descent direction
    // rho_1 = rho; rho = itl::dot(r, r);
    rho_1 = rho; rho = vect_sp(z, z); // ...
    // std::cout << "norm of residu : " << rho << endl; getchar();


    if (sqrt(dal::abs(rho)) < residu * norm_b) break;

    if (noisy && transition) std::cout << "transition\n";
    if (transition || (iter == 0)) gamma = 0.0;
    else
      // gamma = std::max(0.0, (rho - itl::dot(old_r, r) ) / rho_1);
      // gamma = rho / rho_1;
      gamma = std::max(0.0, (rho - vect_sp(old_z, z) ) / rho_1); // ...
    // std::cout << "gamma = " << gamma << endl;
    // itl::add(r, itl::scaled(p, gamma), p);
    add(z, scaled(p, gamma), p); // ...
 
    if (++iter >= itemax) return -1;
    // one dimensionnal optimization
    mult(A, p, q2);
    mult(M, q2, q);
    lambda = rho / vect_sp(q,p);
    for (size_type i = 0; i < mat_nrows(C); ++i)
      if (!satured[i])
      {
	value_type bb = vect_sp(mat_row(C, i), p);
	if (bb > 0.0)
	  lambda = std::min(lambda, -vect_sp(mat_row(C, i), x) / bb);
      }
    add(x, scaled(p, lambda), x);
    add(memox, scaled(x, -1.0), memox);
    if (noisy > 0)  cout << "iter " << iter << " residu "
			 << sqrt(dal::abs(rho)) / norm_b << endl;
    
  }
  return iter;
}

}





#endif //  __GMM_SOLVER_CCG_H
