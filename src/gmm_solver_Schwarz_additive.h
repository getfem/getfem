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
/* Copyright (C) 2002-2004  Yves Renard, Michel Fournié.                   */
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


#ifndef GMM_SOLVERS_SCHWARZ_ADDITIVE_H__
#define GMM_SOLVERS_SCHWARZ_ADDITIVE_H__

#include <gmm_kernel.h>
#include <gmm_precond_diagonal.h>
#include <gmm_superlu_interface.h>

namespace gmm {
      
  /* ******************************************************************** */
  /*		Sequential Linear Additive Schwarz method                 */
  /* ******************************************************************** */
  /* ref : Domain decomposition algorithms for the p-version finite       */
  /*       element method for elliptic problems, Luca F. Pavarino,        */
  /*       PhD thesis, Courant Institute of Mathematical Sciences, 1992.  */
  /* ******************************************************************** */

  /* ******************************************************************** */
  /*		Additive Schwarz interfaced local solvers                 */
  /* ******************************************************************** */

  struct using_cg {};
  struct using_gmres {};
  struct using_bicgstab {};

  template <typename P, typename local_solver, typename Matrix>
  struct actual_precond {
    typedef P APrecond;
    static APrecond &transform(const P &PP) { return PP; }
  };

  template <typename Matrix1, typename Precond, typename Vector> 
  void SA_local_solve(using_cg, const Matrix1 &A, Vector &x, const Vector &b,
		 const Precond &P, iteration &iter)
  { cg(A, x, b, P, iter); }

  template <typename Matrix1, typename Precond, typename Vector> 
  void SA_local_solve(using_gmres, const Matrix1 &A, Vector &x,
		      const Vector &b, const Precond &P, iteration &iter)
  { gmres(A, x, b, P, 100, iter); }
  
  template <typename Matrix1, typename Precond, typename Vector> 
  void SA_local_solve(using_bicgstab, const Matrix1 &A, Vector &x,
		      const Vector &b, const Precond &P, iteration &iter)
  { bicgstab(A, x, b, P, iter); }

#ifdef GMM_USES_SUPERLU
  struct using_superlu {};

  template <typename P, typename Matrix>
  struct actual_precond<P, using_superlu, Matrix> {
    typedef typename linalg_traits<Matrix>::value_type value_type;
    typedef SuperLU_factor<value_type> APrecond;
    template <typename PR>
    static APrecond transform(const PR &) { return APrecond(); }
    static APrecond &transform(const APrecond &PP) { return PP; }
  };

  template <typename Matrix1, typename Precond, typename Vector> 
  void SA_local_solve(using_superlu, const Matrix1 &, Vector &x,
		      const Vector &b, const Precond &P, iteration &iter)
  { P.solve(x, b); iter.set_iteration(1); }
#endif

  /* ******************************************************************** */
  /*		Additive Schwarz Linear system                            */
  /* ******************************************************************** */

  template <typename Matrix1, typename Matrix2, typename Precond,
	    typename local_solver>
  struct schwadd_mat{
    typedef typename linalg_traits<Matrix1>::value_type value_type;
    typedef typename std::vector<value_type> vector_type;
    typedef typename actual_precond<Precond, local_solver, Matrix1>::APrecond
            aprecond;

    const Matrix1 *A;
    const std::vector<Matrix2> *vB;
    std::vector<Matrix2> vAloc;
    mutable iteration iter;
    double residu;
    mutable size_type itebilan;
    mutable std::vector<vector_type> gi, fi;
    std::vector<aprecond> precond1;

    void init(const Matrix1 &A_, const std::vector<Matrix2> &vB_,
	      iteration iter_, const Precond &P, double residu_);

    schwadd_mat(void) {}
    schwadd_mat(const Matrix1 &A_, const std::vector<Matrix2> &vB_,
		iteration iter_, const Precond &P, double residu_)
    { init(A_, vB_, iter_, P, residu_); }
  };

  template <typename Matrix1, typename Matrix2, typename Precond,
	    typename local_solver>
  void schwadd_mat<Matrix1, Matrix2, Precond, local_solver>::init(
       const Matrix1 &A_, const std::vector<Matrix2> &vB_,
       iteration iter_, const Precond &P, double residu_) {

    vB = &vB_; A = &A_; iter = iter_;
    residu = residu_;
    
    size_type nb_sub = vB->size();
    vAloc.resize(nb_sub);
    gi.resize(nb_sub); fi.resize(nb_sub);
    precond1.resize(nb_sub);
    std::fill(precond1.begin(), precond1.end(),
	      actual_precond<Precond, local_solver, Matrix1>::transform(P));
    itebilan = 0;
    
    if (iter.get_noisy()) cout << "Init pour sub dom ";
    for (size_type i = 0; i < nb_sub; ++i) {
      if (iter.get_noisy()) cout << i << " " << std::flush;
      Matrix2 Maux(mat_nrows((*vB)[i]), mat_ncols((*vB)[i])),
	BT(mat_ncols((*vB)[i]), mat_nrows((*vB)[i]));
      
      gmm::copy(gmm::transposed((*vB)[i]), BT);
      gmm::resize(vAloc[i], mat_nrows((*vB)[i]), mat_nrows((*vB)[i]));      
      gmm::mult((*vB)[i], *A, Maux);
      gmm::mult(Maux, BT, vAloc[i]);
      precond1[i].build_with(vAloc[i]);
      gmm::resize(fi[i], mat_nrows((*vB)[i]));
      gmm::resize(gi[i], mat_nrows((*vB)[i]));
    }
    if (iter.get_noisy()) cout << "\n";
  }


  /* ******************************************************************** */
  /*		Additive Schwarz interfaced global solvers                */
  /* ******************************************************************** */

  template <typename Matrix1, typename Matrix2, typename Precond,
	    typename local_solver, typename Vect>
  void SA_global_solve(using_cg, const schwadd_mat<Matrix1, Matrix2, Precond,
		  local_solver> &ASM, Vect &x, const Vect &b, iteration &iter)
  { cg(ASM, x, b, *(ASM.A),  identity_matrix(), iter); }

  template <typename Matrix1, typename Matrix2, typename Precond,
	    typename local_solver, typename Vect>
  void SA_global_solve(using_gmres, const schwadd_mat<Matrix1, Matrix2,
	Precond, local_solver> &ASM, Vect &x, const Vect &b, iteration &iter)
  { gmres(ASM, x, b, identity_matrix(), 100, iter); }

  template <typename Matrix1, typename Matrix2, typename Precond,
	    typename local_solver, typename Vect>
  void SA_global_solve(using_bicgstab, const schwadd_mat<Matrix1, Matrix2,
	 Precond, local_solver> &ASM, Vect &x, const Vect &b, iteration &iter)
  { bicgstab(ASM, x, b, identity_matrix(), 100, iter); }

#ifdef GMM_USES_SUPERLU
  template <typename Matrix1, typename Matrix2, typename Precond,
	    typename local_solver, typename Vect>
  void SA_global_solve(using_superlu, const schwadd_mat<Matrix1, Matrix2,
	    Precond, local_solver> &, Vect &, const Vect &, iteration &) {
    DAL_THROW(failure_error,
     "You cannot use SuperLU as global solver in additive Schwarz meethod");
  }
#endif
  
  /* ******************************************************************** */
  /*		Additive Schwarz algorithm                                */
  /* ******************************************************************** */

  /** Function to call if the ASM matrix is precomputed for successive solve
   * with the same system.
   */
  template <typename Matrix1, typename Matrix2,
	    typename Vector2, typename Vector3, typename Precond,
	    typename local_solver, typename global_solver>
  int sequential_additive_schwarz(
    schwadd_mat<Matrix1, Matrix2, Precond, local_solver> &ASM, Vector3 &u,
    const Vector2 &f, iteration &iter, const global_solver&) {

    typedef typename linalg_traits<Matrix1>::value_type value_type;

    size_type nb_sub = ASM.vB->size(), nb_dof = gmm::vect_size(f);
    ASM.itebilan = 0;
    std::vector<value_type> g(nb_dof);

    for (size_type i = 0; i < nb_sub; ++i) {
      gmm::mult((*(ASM.vB))[i], f, ASM.fi[i]);
      ASM.iter.init();
      SA_local_solve(local_solver(), ASM.vAloc[i], ASM.gi[i], ASM.fi[i],
		     ASM.precond1[i], ASM.iter);
      // local_solver::solve(ASM.vAloc[i], ASM.gi[i], ASM.fi[i],
      //		  ASM.precond1[i], ASM.iter);
      ASM.itebilan = std::max(ASM.itebilan, ASM.iter.get_iteration());
      gmm::mult(gmm::transposed((*(ASM.vB))[i]), ASM.gi[i], g, g);
    }

    double utime = uclock_sec();
    SA_global_solve(global_solver(), ASM, u, g, iter);
    if (iter.get_noisy() > 0)
      cout << "Time for global resolution : "
	   << uclock_sec() - utime << " seconds" << endl;
    return ASM.itebilan;
  }

  /** Global function. Compute the ASM matrix and call the previous function.
      The ASM matrix represent the preconditionned linear system.
   */
  template <typename Matrix1, typename Matrix2,
	    typename Vector2, typename Vector3, typename Precond,
	    typename local_solver, typename global_solver>
  int sequential_additive_schwarz(const Matrix1 &A, Vector3 &u,
				  const Vector2 &f, const Precond &P,
				  const std::vector<Matrix2> &vB,
				  iteration &iter, const local_solver&,
				  const global_solver &GS) {
    iter.set_rhsnorm(vect_norm2(f));
    if (iter.get_rhsnorm() == 0.0) { gmm::clear(u); return 0; }
    iteration iter2 = iter; iter2.reduce_noisy();
    iter2.set_maxiter(size_type(-1));
    // iter2.set_resmax(iter.get_resmax() / 100.0);

    schwadd_mat<Matrix1, Matrix2, Precond, local_solver>
      ASM(A, vB, iter2, P, iter.get_resmax());

    return sequential_additive_schwarz(ASM, u, f, iter, GS);
  }
  
  template <typename Matrix1, typename Matrix2, typename Precond,
	    typename Vector2, typename Vector3, typename local_solver>
  void mult(const schwadd_mat<Matrix1, Matrix2, Precond, local_solver> &M,
	    const Vector2 &p, Vector3 &q) {
    size_type itebilan = 0;
    mult(*(M.A), p, q);
    globaltolocal(q, M.fi, *(M.vB));
    for (size_type i = 0; i < M.fi.size(); ++i) {
      M.iter.init();
      SA_local_solve(local_solver(), (M.vAloc)[i], (M.gi)[i],
		     (M.fi)[i],(M.precond1)[i],M.iter);
      // local_solver::solve((M.vAloc)[i], (M.gi)[i],
      //		       (M.fi)[i],(M.precond1)[i],M.iter);
      itebilan = std::max(itebilan, M.iter.get_iteration());
    }
    localtoglobal(M.gi, q, *(M.vB));
    if (M.iter.get_noisy() > 0) cout << "itebloc = " << itebilan << endl;
    M.itebilan += itebilan;
    M.iter.set_resmax((M.iter.get_resmax() + M.residu) * 0.5);
  }

  template <typename Matrix1, typename Matrix2, typename Precond,
	    typename Vector2, typename Vector3, typename Vector4,
	    typename local_solver>
  void mult(const schwadd_mat<Matrix1, Matrix2, Precond, local_solver> &M,
	    const Vector2 &p, const Vector3 &p2, Vector4 &q)
  { mult(M, p, q); add(p2, q); }

  template <typename Matrix2, typename Vector2, typename Vector3>
  void globaltolocal(const Vector2 &f, std::vector<Vector3> &fi,
		       const std::vector<Matrix2> &vB) {
    for (size_type i = 0; i < fi.size(); ++i) gmm::mult(vB[i], f, fi[i]);
  }

  template <typename Matrix2, typename Vector2, typename Vector3>
  void localtoglobal(const std::vector<Vector3> &fi, Vector2 &f, 
		     const std::vector<Matrix2> &vB) {
    gmm::clear(f);
    for (size_type i = 0; i < fi.size(); ++i)
      gmm::mult(gmm::transposed(vB[i]), fi[i], f, f);
  }




  /* ******************************************************************** */
  /*		Sequential Non-Linear Additive Schwarz method             */
  /* ******************************************************************** */
  /* ref : Nonlinearly Preconditionned Inexact Newton Algorithms,         */
  /*       Xiao-Chuan Cai, David E. Keyes,                                */
  /*       SIAM J. Sci. Comp. 24: p183-200. .                             */
  /* ******************************************************************** */

#ifdef sdlkfjfdslfkj

  template <typename Matrixt, typename MatrixBi> 
  class NewtonAS_struct {
    
    typedef Matrixt tangent_matrix_type;
    typedef typename linalg_traits<Matrix1>::value_type value_type;
    
  public :
    size_type nb_subdomains(void);
    size_type size(void);
    const std::vector<MatrixBi> &get_vB();
    
    void compute_tangent_matrix(Matrixt &M, std::vector<value_type> x);
    
    // compute Bi^T grad(F(X)) Bi
    void compute_sub_tangent_matrix(Matrixt &Mloc, 
				    std::vector<value_type> x, size_type i);
    
    // compute Bi^T F(X)
    void compute_sub_F(std::vector<value_type> fi,
		       std::vector<value_type> x, size_type i);
    
  };
  
  
  template <typename NewtonAS_struct_, typename Vector, typename Precond
	    typename local_solver, typename global_solver>
  Newton_additive_Schwarz(NewtonAS_struct_ &NS, const Vector &u_,
			  iteration &iter, const Precond &P,
			  local_solver &LS, global_solver &GS) {
    Vector &u = const_cast<Vector &>(u_);
    typedef typename linalg_traits<Vector>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type m_type;
    typedef typename NewtonAS_struct_::tangent_matrix_type Matrixt;
    typedef typename actual_precond<Precond, local_solver, Matrix1>::APrecond
      aprecond;

    aprecond PP;
    PP = actual_precond<Precond, local_solver, Matrix1>::transform(P);

    iter.set_rhsnorm(m_type(1));
    iteration iter2(iter);
    iter2.reduce_noisy(); iter2.set_maxiter(size_type(-1));
    iteration iter3(iter); iter3.reduce_noisy();
    
    std::vector<value_type> rhs(NS.size()), x(NS.size()), d(NS.size());
    std::vector<value_type> xi, xii, fi, di;
    Matrixt Mloc, M(NS.size(), NS.size());
    
    for(;;) { // Newton global iteration
      
      gmm::clear(rhs); 
      for (size_type isd = 0; isd < NS.nb_subdomains(); ++isd) {
	// solving non-linear local problem isd
	typename NewtonAS_struct_::Matrixloc &Bi = (NS.get_vB())[isd];
	size_type si = mat_ncols(Bi);
	gmm::resize(Mloc, si, si);
	xi.resize(si); xii.resize(si); fi.resize(si); di.resize(si);
	
	iter2.init(); 
	gmm::clear(xi);
	gmm::copy(u, x);
	NS.compute_sub_F(fi, x, isd);
	m_type r = gmm::norm2(fi), r_t;
	
	while(!iter2.finished_vect(r)) {
	  
	  NS.compute_sub_tangent_matrix(Mloc, x, isd);
	  P.init_with(Mloc);
	  iter3.init();
	  SA_local_solve(local_solver(), Mloc, di, fi, P, iter3);
	  // local_solver::solve(Mloc, di, fi, P, iter3);
	  
	  for (m_type alpha(1); alpha>=m_type(1)/m_type(8); alpha/=m_type(2)) {
	    gmm::add(xi, gmm::scaled(di, alpha), xii);
	    gmm::mult(Bi, xii, u, x);
	    NS.compute_sub_F(fi, x, isd);
	    if ((r_t = gmm::norm2(fi)) <= r) break;
	  }
	  r = r_t;
	  gmm::copy(xii, xi);
	}
	gmm::mult(Bi, gmm::scaled(xii, value_type(-1)), rhs, rhs);
      }
      
      if (iter2.finished_vect(rhs)) break;
      
      // solving linear system for the global Newton method
      NS.compute_tangent_matrix(M, u);

      schwadd_mat<Matrix1, Matrix2, Precond, local_solver>
	ASM(M, NS.get_vB(), iter2, P, iter.get_resmax());

      SA_global_solve(global_solver(), ASM, d, rhs, iter);
      gmm::add(d, u);
    }
  }



#endif














  /* ******************************************************************** */
  /*		Old version, obsolete.                                    */
  /* ******************************************************************** */

  #define PRECOND ildltt_precond

  template <typename Matrix1, typename Matrix2, typename Matrix3,
	    typename SUBI>
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

  template <typename Matrix1, typename Matrix2, typename Matrix3,
	    typename SUBI,    typename Vector2, typename Vector3>
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
      precond2[i] = PRECOND<Matrix3>(ml2[i], 10, 1E-7);

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
  
  template <typename Matrix1, typename Matrix2, typename Matrix3,
	    typename SUBI,    typename Vector2, typename Vector3>
  void mult(const schwarz_additif_matrix<Matrix1, Matrix2, Matrix3, SUBI> &M,
	    const Vector2 &p, Vector3 &q) {

    size_type itebilan = 0;
    size_type ms = (M.ml1)->size();
    mult(*(M.A), p, q);
    global_to_local(q, *(M.fi), *(M.cor));
    for (size_type i = 0; i < (M.ml1)->size(); ++i) {
      M.iter.init();
      cg((*(M.ml1))[i], (*(M.gi))[i], (*(M.fi))[i],(*(M.precond1))[i], M.iter);
      itebilan = std::max(itebilan, M.iter.get_iteration());
    }

    for (size_type i = 0; i < (M.ml2)->size(); ++i) {
      M.iter.init();
      cg((*(M.ml2))[i],(*(M.gi))[i+ms], (*(M.fi))[i+ms], (*(M.precond2))[i],
	 M.iter);
      itebilan = std::max(itebilan, M.iter.get_iteration());
    }

    local_to_global(*(M.gi), q, *(M.cor));
    if (M.iter.get_noisy()) cout << "itebloc = " << itebilan << endl;
    M.itebilan += itebilan;
    M.iter.set_resmax((M.iter.get_resmax() + M.residu) * 0.5);
  }

  template <typename Matrix1, typename Matrix2, typename Matrix3,
	    typename SUBI,    typename Vector2, typename Vector3,
	    typename Vector4>
  void mult(const schwarz_additif_matrix<Matrix1, Matrix2, Matrix3, SUBI> &M,
	    const Vector2 &p, const Vector3 &p2, Vector4 &q)
  { mult(M, p, q); add(p2, q); }

  template <typename SUBI, typename Vector2, typename Vector3>
  void global_to_local(const Vector2 &f, std::vector<Vector3> &fi,
		       const std::vector<SUBI> &cor) {
    for (size_type i = 0; i < fi.size(); ++i) {
      typename linalg_traits<Vector3>::iterator it2 = fi[i].begin();
      for (size_type j = 0, l = cor[i].size(); j < l; ++j , ++it2)
        *it2 = f[cor[i].index(j)];
    }
  }

  template <typename SUBI, typename Vector2, typename Vector3>
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


#endif //  GMM_SOLVERS_SCHWARZ_ADDITIVE_H__
