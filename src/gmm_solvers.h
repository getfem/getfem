/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_solvers.h : generic solvers.                             */
/*     									   */
/* Date : October 13, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001  Yves Renard.                                        */
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


#ifndef __GMM_SOLVERS_H
#define __GMM_SOLVERS_H

#include <gmm_dense_lu.h>

#include <gmm_iter.h>


namespace gmm {

  // Needed by ilut and choleskyt
  template<class T> struct _elt_rsvector_value_less {
    inline bool operator()(const _elt_rsvector<T>& a, 
			   const _elt_rsvector<T>& b) const
    { return (dal::abs(a.e) > dal::abs(b.e)); }
  };
}

#include <gmm_precond_diagonal.h>
#include <gmm_precond_cholesky.h>
#include <gmm_precond_choleskyt.h>
#include <gmm_precond_mr_approx_inverse.h>
#include <gmm_precond_ilu.h>
#include <gmm_precond_ilut.h>



#include <gmm_solver_cg.h>
#include <gmm_solver_bicgstab.h>
#include <gmm_solver_qmr.h>
#include <gmm_solver_cheby.h>
#include <gmm_solver_constrained_cg.h>
#include <gmm_solver_Schwarz_additive.h>
#include <gmm_given_rotation.h>
#include <gmm_modified_gram_schmidt.h>
#include <gmm_tri_solve.h>
#include <gmm_solver_gmres.h>



#endif //  __GMM_SOLVERS_H
