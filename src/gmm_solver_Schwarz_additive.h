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
/* This file is a part of GMM++                                            */
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


  #define PRECOND choleskyt_precond

  template <class Matrix1, class Matrix2, class Matrix3, class SUBI>
  struct schwarz_additif_matrix {
    typedef typename linalg_traits<Matrix2>::value_type value_type;
    typedef typename dense_vector_type<value_type>::vector_type vector_type; 
    const Matrix1 *A;
    const std::vector<Matrix2> *ml1;
    const std::vector<Matrix3> *ml2;
    const std::vector<SUBI> *cor;
    mutable iteration iter;
    double residu;
    mutable size_t itebilan;
    std::vector<vector_type> *gi;
    std::vector<vector_type> *fi;
    
    std::vector<PRECOND<Matrix2> > *precond1;
    std::vector<PRECOND<Matrix3> > *precond2;

  };

  template <class Matrix1, class Matrix2, class Matrix3, class SUBI,
	    class Vector2, class Vector3>
  int schwarz_additif(const Matrix1 &A, Vector3 &u,
		      const std::vector<Matrix2> &ml1,
		      const std::vector<Matrix3> &ml2,
		      const std::vector<SUBI> &cor,
		      const Vector2 &f,
		      iteration &iter) {

    typedef typename linalg_traits<Matrix2>::value_type value_type;
    typedef typename dense_vector_type<value_type>::vector_type vector_type;

    size_type nb_sub = ml1.size() + ml2.size();
    size_t itebilan = 0;
    std::vector<vector_type> gi(nb_sub);
    std::vector<vector_type> fi(nb_sub);
    
    std::vector<PRECOND<Matrix2> > precond1(ml1.size());
    std::vector<PRECOND<Matrix3> > precond2(ml2.size());

    for (size_type i = 0; i < ml1.size(); ++i)
      precond1[i] = PRECOND<Matrix2>(ml1[i], 10, 1E-7);
    for (size_type i = 0; i < ml2.size(); ++i)
      precond2[i] = PRECOND<Matrix2>(ml2[i], 10, 1E-7);

    iter.set_rhsnorm(vect_norm2(f));
    if (iter.get_rhsnorm() == 0.0) { clear(u); return 0; }

    size_type ms = ml1.size();

    for (size_type i = 0; i < nb_sub; ++i) {
      size_type k = i < ms ? mat_nrows(ml1[i]) : mat_nrows(ml2[i-ms]);
      fi[i] = gi[i] = vector_type(k);
      clear(gi[i]);
    }

    size_type nb_dof = f.size();
    global_to_local(f, fi, cor);

    iteration iter2 = iter;
    iter2.reduce_noisy();

    for (size_type i = 0; i < ms; ++i) {
      iter2.init();
      cg(ml1[i], gi[i], fi[i], identity_matrix(), precond1[i], iter2);
      itebilan = std::max(itebilan, iter2.get_iteration());
    }
    for (size_type i = 0; i < ml2.size(); ++i) {
      iter2.init();
      cg(ml2[i], gi[i+ms], fi[i+ms],
	 identity_matrix(), precond2[i], iter2);
      itebilan = std::max(itebilan, iter2.get_iteration());
    }

    vector_type g(nb_dof);
    local_to_global(gi, g, cor);
    
    schwarz_additif_matrix<Matrix1, Matrix2, Matrix3, SUBI> SAM;
    SAM.A = &A; SAM.ml1 = &ml1; SAM.ml2 = &ml2; SAM.cor = &cor;
    SAM.precond1 = &precond1; SAM.precond2 = &precond2; 
    iter2.init();
    SAM.iter = iter2;
    SAM.residu = iter.get_resmax();
    // SAM.residu_act = 1E-2;
    SAM.gi = &gi; SAM.fi = &fi; SAM.itebilan = itebilan;
   
    cg(SAM, u, g, A, identity_matrix(), iter);

    return SAM.itebilan;
  }
  
  template <class Matrix1, class Matrix2, class Matrix3, class SUBI,
    class Vector2, class Vector3>
  void mult(const schwarz_additif_matrix<Matrix1, Matrix2, Matrix3, SUBI> &M,
	    const Vector2 &p, Vector3 &q) {

    size_type itebilan = 0;
    size_type ms = (M.ml1)->size();
    mult(*(M.A), p, q);
    global_to_local(q, *(M.fi), *(M.cor));
    for (size_type i = 0; i < (M.ml1)->size(); ++i) {
      M.iter.init();
      cg((*(M.ml1))[i], (*(M.gi))[i], (*(M.fi))[i], (*(M.precond1))[i], M.iter);
      itebilan = std::max(itebilan, M.iter.get_iteration());
    }

    for (size_type i = 0; i < (M.ml2)->size(); ++i) {
      M.iter.init();
      cg((*(M.ml2))[i],(*(M.gi))[i+ms], (*(M.fi))[i+ms], (*(M.precond2))[i],
	 M.iter);
      itebilan = std::max(itebilan, M.iter.get_iteration());
    }

    local_to_global(*(M.gi), q, *(M.cor));
    cout << "itebloc = " << itebilan << endl;
    M.itebilan += itebilan;
    M.iter.set_resmax((M.iter.get_resmax() + M.residu) * 0.5);
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

  template <class SUBI, class Vector2, class Vector3>
  void small_local_to_global(const std::vector<Vector3> &fi, Vector2 &f, 
		       const std::vector<SUBI> &cor, size_type i) {
    clear(f);
    typename linalg_traits<Vector3>::const_iterator it2=fi[i].begin();
    for (size_type j = 0, l = cor[i].size(); j < l; ++j, ++it2) {
      f[cor[i].index(j)] = *it2;
    }
  }
  
  // CO global constraint matrix (CO * U <= cof)
  // f  RHS

  template <class Matrix1, class Matrix2, class Matrix3, class Matrix4,
	    class Matrix5, class Matrix6, class SUBI, class Vector2,
	    class Vector3, class Vector4, class Vector5>
  int schwarz_with_constraints(const Matrix1 &A,
			       Vector3 &u, const Matrix4 &CO,
			       const std::vector<Matrix2> &ml1,
			       const std::vector<Matrix6> &mco1, 
			       const std::vector<Matrix3> &ml2,
			       const std::vector<Matrix5> &mco2, 
			       const std::vector<SUBI> &cor,
			       const Vector2 &f,
			       const Vector4 &cof,
			       const std::vector<Vector5> &cofi,
			       iteration &iter) {

    typedef typename linalg_traits<Matrix2>::value_type value_type;
    typedef typename dense_vector_type<value_type>::vector_type vector_type;
    
    size_type nb_sub = ml1.size() + ml2.size();
    size_t itebilan = 0;
    std::vector<vector_type> gi(nb_sub);
    std::vector<vector_type> fi(nb_sub);
    std::vector<vector_type> ui(nb_sub);
    std::vector<vector_type> wi(nb_sub);
    iter.set_rhsnorm(vect_norm2(f));

    size_type ms = ml1.size();

    for (size_type i = 0; i < nb_sub; ++i) {
      size_type k = i < ms ? mat_nrows(ml1[i]) : mat_nrows(ml2[i-ms]);
      size_type l = i < ms ? mat_nrows(mco1[i]) : mat_nrows(mco2[i]); 
      ui[i] = vector_type(k);   fi[i] = vector_type(k);
      gi[i] = vector_type(k);   wi[i] = vector_type(k);
    }

    vector_type w(vect_size(u));
    global_to_local(f, fi, cor);
    // global_to_local(cof, cofi, cor); // pas bon, il faudrait un cor pour les contraintes ...

    iteration iter2 = iter;
    iter2.reduce_noisy();

    for (;;) {

      // Step 1
      gmm::mult(A, u, w);
      global_to_local(w, wi, cor);
      global_to_local(u, ui, cor);
      
      for (size_type i = 0; i < nb_sub; ++i) {
	gmm::add(fi[i], gmm::scaled(wi[i], -1.0), wi[i]);
	clear(gi[i]);
      }
      
      for (size_type i = 0; i < ms; ++i) {
	iter2.init();
	vector_type cofloc(mat_nrows(mco1[i]));
	if (mat_nrows(mco1[i]) > 0) 
	  gmm::mult(mco1[i], gmm::scaled(ui[i], -1.0), cofi[i], cofloc);
	constrained_cg(ml1[i], mco1[i], gi[i], wi[i], cofloc,
		       identity_matrix(), identity_matrix(), iter2);
	itebilan = std::max(itebilan, iter2.get_iteration());
      }
      
      for (size_type i = 0; i < ml2.size(); ++i) {
	iter2.init();
	vector_type cofloc(mat_nrows(mco2[i]));
	gmm::mult(mco2[i], gmm::scaled(ui[i+ms], -1.0), cofi[i+ms], cofloc); 
	constrained_cg(ml2[i], mco2[i], gi[i+ms], wi[i+ms], cofloc,
		       identity_matrix(), identity_matrix(), iter2);
	itebilan = std::max(itebilan, iter2.get_iteration());
      }
      
      // Step 2

      gmm::row_matrix<std::vector<value_type> > global_sm(nb_sub, nb_sub);
      gmm::col_matrix<std::vector<value_type> > 
	global_CO1(mat_nrows(CO), nb_sub);
      gmm::row_matrix<std::vector<value_type> > 
	global_CO2(mat_nrows(CO), nb_sub);
      std::vector<value_type> global_cof(mat_nrows(CO));
      std::vector<value_type> global_f(nb_sub), alpha(nb_sub);
      std::vector<gmm::wsvector<value_type> > Gi(nb_sub);
      gmm::wsvector<value_type> W(vect_size(u));
      for (size_type i = 0; i < nb_sub; ++i) {
	Gi[i] = gmm::wsvector<value_type>(vect_size(u));
	small_local_to_global(gi, Gi[i], cor, i);
      } // to be optimized (passer à des produits locaux)
      
      for (size_type i = 0; i < nb_sub; ++i) {
	if (mat_nrows(CO) > 0) {
	  gmm::mult(CO, Gi[i], gmm::mat_col(global_CO1, i));
	  gmm::mult(CO, u, cof, global_cof);
	}
	gmm::mult(A, Gi[i], W);
	global_f[i] = gmm::vect_sp(f, Gi[i]) - gmm::vect_sp(w, Gi[i]) ;
	for (size_type j = 0; j <= i; ++j)
	  global_sm(i,j) = global_sm(j,i) = gmm::vect_sp(W, Gi[j]);
      } // to be optimized (symmetrie et produits locaux)

      cout << "global_sm = " << global_sm << endl;

      if (mat_nrows(CO) > 0) gmm::copy(global_CO1, global_CO2);
      
      
      size_type nbconst = 0;
      for (size_type i = 0; i < mat_nrows(CO); ++i) {
	if (gmm::vect_norm2(mat_row(global_CO2, i)) > 1E-10)
	  nbconst++;
      }
      gmm::row_matrix<std::vector<value_type> > 
	global_CO3(nbconst, nb_sub);
      std::vector<value_type> global_cof2(nbconst);
      for (size_type i = 0, k = 0; i < mat_nrows(CO); ++i) {
	if (gmm::vect_norm2(mat_row(global_CO2, i)) > 1E-10) {
	  global_cof2[k] = global_cof[i];
	  copy(mat_row(global_CO2, i), mat_row(global_CO3, ++k));
	}
      }
     
      iteration iter3 = iter;
      iter3.reduce_noisy();
      iter3.init();
      gmm::clear(alpha);

//        global_CO3(nbconst, 0) = -1.0; global_cof2[nbconst] = -10.0;
//        alpha[0] = 10;

      cout << "global_CO3 = " << global_CO3 << endl;

      constrained_cg(global_sm, global_CO3, alpha, global_f, global_cof2,
 		     identity_matrix(), identity_matrix(), iter3);
      value_type res = 0, sum_alphai = 0;
      for (size_type i = 0; i < nb_sub; ++i) {
	cout << "alpha[" << i << "] = " << alpha[i] << endl;
	// cout << "u[" << i << "] = " << Gi[i] << endl;
	if (alpha[i] < 0)
	  cout << "WARNING : alpha[" << i << "] = " << alpha[i] << endl;
	res += alpha[i] * vect_norm2(Gi[i]);
	sum_alphai += alpha[i];
	gmm::add(u, gmm::scaled(Gi[i], alpha[i]), u);
      }
      cout << "sum alpha_i = " << sum_alphai << endl;

      ++iter;
      if (iter.finished(res)) break;
    }

    return itebilan;
  }








  
}


#endif //  __GMM_SOLVERS_SCHWARZ_ADDITIVE_H
