/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_solvers_Schwarz_additive.h : generic solver.             */
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


#ifndef __GMM_SOLVERS_SCHWARZ_ADDITIVE_H
#define __GMM_SOLVERS_SCHWARZ_ADDITIVE_H

namespace gmm {
      
  /* ******************************************************************** */
  /*		Schwartz Additive method                                  */
  /* ******************************************************************** */

  template <class Matrix1, class Matrix2, class Matrix3, class Vector1>
  struct schwarz_additif_matrix {
    typedef typename linalg_traits<Matrix2>::value_type value_type;
    const Matrix1 *A;
    const std::vector<Matrix2> *ml1;
    const std::vector<Matrix3> *ml2;
    const std::vector<Vector1> *cor;
    int itemax, noisy;
    double residu;
    std::vector< std::vector<value_type> > *gi;
    std::vector< std::vector<value_type> > *fi;
  };

  template <class Matrix1, class Matrix2, class Matrix3,
    class Vector1, class Vector2, class Vector3>
  void schwarz_additif(const Matrix1 &A,
		       Vector3 &u,
		       const std::vector<Matrix2> &ml1,
		       const std::vector<Matrix3> &ml2,
		       const std::vector<Vector1> &cor,
		       const Vector2 &f,
		       int itemax,  double residu, int noisy = 1) {
    typedef typename linalg_traits<Matrix2>::value_type value_type;
    
    size_type nb_sub = ml1.size() + ml2.size();
    std::vector< std::vector<value_type> > gi(nb_sub);
    std::vector< std::vector<value_type> > fi(nb_sub);

    cout << "nb sub domains : " << nb_sub << endl;

    size_type ms = ml1.size();

    for (int i = 0; i < nb_sub; ++i) {
      size_type k = i < ms ? mat_nrows(ml1[i]) : mat_nrows(ml2[i-ms]);
      cout << "Taille du sous système " << i << " = " << k << endl;
      fi[i] = gi[i] = std::vector<value_type>(k);
    }

    size_type nb_dof = f.size();
    global_to_local(f, fi, cor);

    for (int i = 0; i < nb_sub; ++i)
      cg(i < ms ? ml1[i] : ml2[i-ms], gi[i], fi[i], itemax, residu, noisy - 1);

    std::vector<value_type> g(nb_dof);
    local_to_global(g, gi, cor);
    schwarz_additif_matrix<Matrix1, Matrix2, Matrix3, Vector1> SAM;
    SAM.A = &A; SAM.ml1 = &ml1; SAM.ml2 = &ml2; SAM.cor = &cor;
    SAM.itemax = itemax; SAM.residu = residu; SAM.noisy = noisy;
    SAM.gi = &gi; SAM.fi = &fi;
    cg(SAM, u, g, A, identity_matrix(), itemax, residu, noisy - 1);
  }

  
  template <class Matrix1, class Matrix2, class Matrix3, class Vector1,
    class Vector2, class Vector3>
  void mult(const schwarz_additif_matrix<Matrix1, Matrix2, Matrix3,Vector1> &M,
	    const Vector2 &p, Vector3 &q) {

    size_type ms = (M.ml1)->size();
    mult(*(M.A), p, q);
    global_to_local(q, *(M.fi), *(M.cor));
    for (int i = 0; i < (M.ml1)->size(); ++i)
      cg((*(M.ml1))[i], (*(M.gi))[i], (*(M.fi))[i], M.itemax, M.residu,
	 M.noisy-1);

    for (int i = 0; i < (M.ml2)->size(); ++i)
      cg((*(M.ml2))[i],(*(M.gi))[i+ms],(*(M.fi))[i+ms],M.itemax, M.residu,
	 M.noisy-1);

    local_to_global(q, *(M.gi), *(M.cor));
  }

  template <class Matrix1, class Matrix2, class Matrix3, class Vector1,
    class Vector2, class Vector3, class Vector4>
  void mult(const schwarz_additif_matrix<Matrix1, Matrix2, Matrix3,Vector1> &M,
	    const Vector2 &p, const Vector3 &p2, Vector4 &q) {
    
    size_type ms = (M.ml1)->size();
    mult(*(M.A), p, q);
    global_to_local(q, *(M.fi), *(M.cor));
    for (int i = 0; i < (M.ml1)->size(); ++i)
      cg((*(M.ml1))[i], (*(M.gi))[i], (*(M.fi))[i], M.itemax, M.residu,
	 M.noisy-1);

    for (int i = 0; i < (M.ml2)->size(); ++i)
      cg((*(M.ml2))[i],(*(M.gi))[i+ms],(*(M.fi))[i+ms],M.itemax,M.residu,
	 M.noisy-1);

    local_to_global(q, *(M.gi), *(M.cor));
    add(p2, q);
  }


  template <class Vector1, class Vector2, class T>
  void global_to_local(const Vector2 &f, std::vector<std::vector<T> > &fi,
		  const std::vector<Vector1> &cor) {
    for (int i = 0; i < fi.size(); ++i) {
      typename Vector1::const_iterator it = cor[i].begin(), ite = cor[i].end();
      typename std::vector<T>::iterator it2 = fi[i].begin(),
	ite2 = fi[i].end();
      for (; it != ite; ++it, ++it2) *it2 = f[*it]; 
    }
  }

  template <class Vector1, class Vector2, class T>
  void local_to_global(Vector2 &f, const std::vector<std::vector<T> > &fi,
		  const std::vector<Vector1> &cor) {
    clear(f);
    for (int i = 0; i < fi.size(); ++i) {
      typename Vector1::const_iterator it = cor[i].begin(), ite = cor[i].end();
      typename std::vector<T>::const_iterator it2=fi[i].begin(),
	ite2=fi[i].end();
      for (; it != ite; ++it, ++it2) f[*it] += *it2; 
    }
  }
  
}


#endif //  __GMM_SOLVERS_SCHWARZ_ADDITIVE_H
