/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_solvers_Schwarz_additive.h : generic solver.             */
/*     									   */
/* Date : October 13, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*          Michel Fournie, fournie@mip.ups-tlse.fr                        */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard, Michel Fournie.                        */
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


#ifndef __GMM_SOLVERS_SCHWARZ_ADDITIVE_H
#define __GMM_SOLVERS_SCHWARZ_ADDITIVE_H

namespace gmm {
      
  /* ******************************************************************** */
  /*		Schwartz Additive method                                  */
  /* ******************************************************************** */

  template <class Matrix1, class Matrix2, class Matrix3, class SUBI>
  struct schwarz_additif_matrix {
    typedef typename linalg_traits<Matrix2>::value_type value_type;
    typedef typename plain_vector_type<value_type>::vector_type vector_type; 
    const Matrix1 *A;
    const std::vector<Matrix2> *ml1;
    const std::vector<Matrix3> *ml2;
    const std::vector<SUBI> *cor;
    int itemax, noisy;
    mutable int itebilan;
    mutable double residu_act;
    double residu;
    std::vector<vector_type> *gi;
    std::vector<vector_type> *fi;
  };

  template <class Matrix1, class Matrix2, class Matrix3, class Matrix4,
	    class Matrix5, class Matrix6, class SUBI, class Vector2, class Vector3>
  int schwarz_additif(const Matrix1 &A,
		      Vector3 &u, const Matrix4 &CO,
		      const std::vector<Matrix2> &ml1,
		      const std::vector<Matrix6> &mco1, 
		      const std::vector<Matrix3> &ml2,
		      const std::vector<Matrix5> &mco2, 
		      const std::vector<SUBI> &cor,
		      const Vector2 &f,
		      int itemax,  double residu, int noisy = 1) {

    typedef typename linalg_traits<Matrix2>::value_type value_type;
    typedef typename plain_vector_type<value_type>::vector_type vector_type;

    size_type nb_sub = ml1.size() + ml2.size();
    int itebilan = 0;
    std::vector<vector_type> gi(nb_sub);
    std::vector<vector_type> fi(nb_sub);

    residu /= nb_sub; // utile ?

    size_type ms = ml1.size();

    for (size_type i = 0; i < nb_sub; ++i) {
      size_type k = i < ms ? mat_nrows(ml1[i]) : mat_nrows(ml2[i-ms]);
      fi[i] = gi[i] = vector_type(k);
      clear(gi[i]);
    }

    size_type nb_dof = f.size();
    global_to_local(f, fi, cor);

    for (size_type i = 0; i < ms; ++i) {
      // cout << "mco1[" << i << "] = " << mco1[i] << endl;
      itebilan = std::max(itebilan,
		constrained_cg(ml1[i], mco1[i], gi[i], fi[i],
			       identity_matrix(),
			       identity_matrix(), itemax, residu, noisy - 1));
      // cout << "gi[" << i << "] = " << gi[i] << endl;
    }
    for (size_type i = 0; i < ml2.size(); ++i) {
      itebilan = std::max(itebilan,
		constrained_cg(ml2[i], mco2[i], gi[i+ms], fi[i+ms],
			       identity_matrix(),
			       identity_matrix(), itemax, residu, noisy - 1));
    }

    vector_type g(nb_dof);
    local_to_global(gi, g, cor);
    
    schwarz_additif_matrix<Matrix1, Matrix2, Matrix3, SUBI> SAM;
    SAM.A = &A; SAM.ml1 = &ml1; SAM.ml2 = &ml2; SAM.cor = &cor;
    SAM.itemax = itemax;
    // SAM.residu_act = 1E-2;
    SAM.residu_act = residu;
    SAM.residu = residu;
    SAM.noisy = noisy; SAM.gi = &gi; SAM.fi = &fi; SAM.itebilan = itebilan;
   
    constrained_cg(SAM, CO, u, g, A, identity_matrix(), itemax, residu, noisy);

    return SAM.itebilan;
  }
  
  template <class Matrix1, class Matrix2, class Matrix3, class SUBI,
    class Vector2, class Vector3>
  void mult(const schwarz_additif_matrix<Matrix1, Matrix2, Matrix3, SUBI> &M,
	    const Vector2 &p, Vector3 &q) {

    int itebilan = 0;
    size_type ms = (M.ml1)->size();
    mult(*(M.A), p, q);
    global_to_local(q, *(M.fi), *(M.cor));
    for (size_type i = 0; i < (M.ml1)->size(); ++i)
      itebilan = std::max(itebilan, 
			  cg((*(M.ml1))[i], (*(M.gi))[i], (*(M.fi))[i],
		      identity_matrix(), M.itemax, M.residu_act, M.noisy-1));

    for (size_type i = 0; i < (M.ml2)->size(); ++i)
      itebilan = std::max(itebilan, 
			  cg((*(M.ml2))[i],(*(M.gi))[i+ms],
			     (*(M.fi))[i+ms], identity_matrix(),
			     M.itemax, M.residu_act, M.noisy-1));

    local_to_global(*(M.gi), q, *(M.cor));
    cout << "itebloc = " << itebilan << endl;
    M.itebilan += itebilan;
    M.residu_act += (M.residu - M.residu_act) * 0.5;
  }

  template <class Matrix1, class Matrix2, class Matrix3, class SUBI,
    class Vector2, class Vector3, class Vector4>
  void mult(const schwarz_additif_matrix<Matrix1, Matrix2, Matrix3, SUBI> &M,
	    const Vector2 &p, const Vector3 &p2, Vector4 &q)
  { mult(M, p, q); add(p2, q); }

  template <class SUBI, class Vector2, class Vector3>
  void global_to_local(const Vector2 &f, std::vector<Vector3> &fi,
		       const std::vector<SUBI> &cor) {
    for (size_type i = 0; i < fi.size(); ++i) {
      typename linalg_traits<Vector3>::iterator it2 = fi[i].begin();
      for (size_type j = 0, l = cor[i].size(); j < l; ++j , ++it2)
        *it2 = f[cor[i].index(j)];
    }
  }

  template <class SUBI, class Vector2, class Vector3>
  void local_to_global(const std::vector<Vector3> &fi, Vector2 &f, 
		       const std::vector<SUBI> &cor) {
    clear(f);
    for (size_type i = 0; i < fi.size(); ++i) {
      typename linalg_traits<Vector3>::const_iterator it2=fi[i].begin();
      for (size_type j = 0, l = cor[i].size(); j < l; ++j, ++it2) {
	f[cor[i].index(j)] += *it2;
      }
    }
  }
  
}


#endif //  __GMM_SOLVERS_SCHWARZ_ADDITIVE_H
