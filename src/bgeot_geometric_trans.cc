// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_geometric_trans.cc : geometric transformations on
//           convexes.
// Date    : December 20, 2000.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2000-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================


#include <dal_singleton.h>
#include <dal_tree_sorted.h>
#include <ftool_naming.h>
#include <bgeot_geometric_trans.h>
#include <bgeot_precomp.h>

namespace bgeot
{
  /* a quick hack for deallocation of allocated fem on program termination
     (avoid false leak alerts from valgrind) */
  struct cleanup_allocated_geotrans : public dal::ptr_collection<geometric_trans> {};
  static pgeometric_trans remember_for_cleanup(geometric_trans *p) {
    dal::singleton<cleanup_allocated_geotrans>::instance().push_back(p);
    return p;
  }

  const base_node& geotrans_interpolation_context::xref() const { 
    if (!have_xref()) 
      if (have_pgp()) xref_ = pgp_->get_point_tab()[ii_];
      else DAL_THROW(dal::failure_error, "missing xref");
    return xref_; 
  }

  const base_node& geotrans_interpolation_context::xreal() const {
    if (!have_xreal()) {
      if (have_pgp()) {
	xreal_ = pgp_->transform(ii_, G());
      } else xreal_ = pgt()->transform(xref(),G());
    }
    return xreal_;
  }

  const base_matrix& geotrans_interpolation_context::B() const {
    if (!have_B()) {
      if (!have_G() || !have_pgt()) {
	DAL_THROW(dal::failure_error, "unable to compute B\n");
      } else {
	size_type P = pgt_->structure()->dim();
	B_.resize(N(), P);
	base_matrix K(N(),P), CS(P,P);
	if (have_pgp()) {
	  //const base_matrix& gr = pgp_->grad(ii_);
	  //gmm::mult(gmm::transposed(pgp_->grad(ii_)), gmm::transposed(G()), K);/*O*/
	  gmm::mult(G(), pgp_->grad(ii_), K);
	} else {
	  base_matrix pc(pgt()->nb_points(), P); 
	  pgt()->gradient(xref(), pc);
	  //gmm::mult(gmm::transposed(pc), gmm::transposed(G()), K);/*O*/
	  gmm::mult(G(),pc,K);
	}
	if (P != N()) {
	  gmm::mult(gmm::transposed(K), K, CS);/*O*/
	  J_ = ::sqrt(gmm::lu_inverse(CS));
	  gmm::mult(K, CS, B_);/*O*/
	} else {
	  //J_ = gmm::abs(gmm::lu_inverse(K)); B_.swap(K); /*O*/
	 J_ = gmm::abs(gmm::lu_inverse(K)); gmm::copy(gmm::transposed(K), B_);
	}
      }
    }
    return B_;
  }

  const base_matrix& geotrans_interpolation_context::B3() const {
    if (!have_B3()) {
      const base_matrix &BB = B(); 
      size_type P=gmm::mat_ncols(BB), N_=gmm::mat_nrows(BB);
      B3_.resize(N_*N_, P*P);
      for (short_type i = 0; i < P; ++i)
	for (short_type j = 0; j < P; ++j)
	  for (short_type k = 0; k < N_; ++k)
	    for (short_type l = 0; l < N_; ++l)
	      B3_(k + N_*l, i + P*j) = BB(k, i) * BB(l, j);
    }
    return B3_;
  }

  const base_matrix& geotrans_interpolation_context::B32() const {
    if (!have_B32()) {
      const base_matrix &BB = B(); 
      size_type P=gmm::mat_ncols(BB), N_=gmm::mat_nrows(BB);	
      B32_.resize(N_*N_, P);
      if (!pgt()->is_linear()) {
	base_matrix B2(P*P, P), Htau(N_, P*P);
	const base_matrix& BB3 = B3();
	if (have_pgp()) {
	  gmm::mult(G(), pgp_->hessian(ii_), Htau);
	} else {
	  DAL_THROW(dal::to_be_done_error,"to be done..");
	}
	for (short_type i = 0; i < P; ++i)
	  for (short_type j = 0; j < P; ++j)
	    for (short_type k = 0; k < P; ++k)
	      for (short_type l = 0; l < N_; ++l)
		B2(i + P*j, k) += Htau(l, i + P*j) * BB(l,k);
	gmm::mult(BB3, B2, B32_);
      } else gmm::clear(B32_);
    }
    return B32_;
  }
  
 void geotrans_interpolation_context::set_ii(size_type ii__) { 
    if (ii_ == ii__) return;
    if (have_B() && !pgt()->is_linear()) { B_.resize(0,0); }
    if (have_B3() && !pgt()->is_linear()) { 
      B3_.resize(0,0); B32_.resize(0,0); }
    xref_.clear(); xreal_.clear();
    ii_=ii__; 
  }

  void geotrans_interpolation_context::set_xref(const base_node& P) {
    xref_ = P;
    if (have_B() && !pgt()->is_linear()) { B_.resize(0,0); }
    if (have_B3() && !pgt()->is_linear()) { 
      B3_.resize(0,0); B32_.resize(0,0); }
    xreal_.clear(); ii_ = size_type(-1);
  }

  geotrans_interpolation_context::geotrans_interpolation_context() :
    G_(0), pgt_(0), pgp_(0), ii_(size_type(-1)) {}
  geotrans_interpolation_context::geotrans_interpolation_context
  (bgeot::pgeotrans_precomp pgp__, size_type ii__, const base_matrix& G__) :
    G_(&G__), pgt_(pgp__->get_trans()), pgp_(pgp__), ii_(ii__) {}
  geotrans_interpolation_context::geotrans_interpolation_context
  (bgeot::pgeometric_trans pgt__, const base_node& xref__,const base_matrix& G__) :
    xref_(xref__), G_(&G__), pgt_(pgt__), pgp_(0), ii_(size_type(-1)) {}
 

  typedef ftool::naming_system<geometric_trans>::param_list gt_param_list;

  base_node geometric_trans::transform(const base_node &pt, 
				       const base_matrix &G) const {
    size_type N = G.nrows(), k = nb_points();
    base_node P(N); P.fill(0.0);
    base_matrix::const_iterator git = G.begin();
    for (size_type l = 0; l < k; ++l) {
      scalar_type a = poly_vector()[l].eval(pt.begin());
      base_node::iterator pit = P.begin(), pite = P.end();
      for (; pit != pite; ++git, ++pit) *pit += a * (*git);
    }
    return P;
  }

  void geometric_trans::gradient(const base_node& x, base_matrix& pc) const {
    base_poly PP;
    pc.resize(nb_points(),dim());
    for (size_type i = 0; i < nb_points(); ++i)
      for (dim_type n = 0; n < dim(); ++n) {
	PP = poly_vector()[i];
	//if (!is_linear()) {
	  PP.derivative(n);
	  pc(i, n) = PP.eval(x.begin());
	  /*} else {
	  pc(i, n) = PP[n+1];
	  }*/
      }
  }

  /* ******************************************************************** */
  /* transformation on simplex.                                           */
  /* ******************************************************************** */

  struct simplex_trans_ : public geometric_trans
  {
    void calc_base_func(base_poly &p, size_type i, short_type K) const
    {
      dim_type N = dim();
      base_poly l0(N, 0), l1(N, 0);
      power_index w(N+1);
      l0.one(); l1.one(); p = l0;
      for (int nn = 0; nn < N; ++nn) l0 -= base_poly(N, 1, nn);
      
      w[0] = K;
      for (int nn = 1; nn <= N; ++nn) {
	w[nn]=int(floor(0.5+((cvr->points()[i])[nn-1]*double(K))));
	w[0]-=w[nn];
      }
      
      for (int nn = 0; nn <= N; ++nn)
	for (int j = 0; j < w[nn]; ++j)
	  if (nn == 0)
	    p *= (l0 * (scalar_type(K) / scalar_type(j+1))) 
	       - (l1 * (scalar_type(j) / scalar_type(j+1)));
	  else
	    p *= (base_poly(N, 1, nn-1) * (scalar_type(K) / scalar_type(j+1))) 
	       - (l1 * (scalar_type(j) / scalar_type(j+1)));
    }

    simplex_trans_(dim_type nc, short_type k)
    {
      cvr = simplex_of_reference(nc, k);
      size_type R = cvr->structure()->nb_points();
      is_lin = (k == 1);
      trans.resize(R);
      for (size_type r = 0; r < R; ++r) calc_base_func(trans[r], r, k);
    }
  };

  static pgeometric_trans PK_gt(gt_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    if (n < 0 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");
    return remember_for_cleanup(new simplex_trans_(n, k));
  }

  /* ******************************************************************** */
  /* direct product transformation                                        */
  /* ******************************************************************** */

  struct cv_pr_t_ : public geometric_trans
  {
    cv_pr_t_(pgeometric_trans a, pgeometric_trans b)
    {
      cvr = convex_ref_product(a->convex_ref(), b->convex_ref());
      is_lin = false;

      size_type n1 = a->nb_points(), n2 = b->nb_points();
      trans.resize(n1 * n2);
      for (size_type i1 = 0; i1 < n1; ++i1)
	for (size_type i2 = 0; i2 < n2; ++i2)
	{
	  trans[i1 + i2 * n1] = a->poly_vector()[i1];
	  trans[i1 + i2 * n1].direct_product(b->poly_vector()[i2]);
	}
    }
  };

  static pgeometric_trans product_gt(gt_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	  "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 1 || params[1].type() != 1)
      DAL_THROW(failure_error, "Bad type of parameters");
    pgeometric_trans a = params[0].method();
    pgeometric_trans b = params[1].method();
    return remember_for_cleanup(new cv_pr_t_(a, b));
  }

  /* ******************************************************************** */
  /* linear direct product transformation.                                */
  /* ******************************************************************** */

  struct cv_pr_tl_ : public geometric_trans
  {
    cv_pr_tl_(pgeometric_trans a, pgeometric_trans b)
    {
      if (!(a->is_linear() && b->is_linear()))
	DAL_THROW(not_linear_error, 
		  "linear product of non-linear transformations");
      cvr = convex_ref_product(a->convex_ref(), b->convex_ref());
      is_lin = true;

      trans.resize(a->nb_points() * b->nb_points());
      std::fill(trans.begin(), trans.end(), null_poly(dim()));

      std::stringstream name;
      name << "GT_PK(" << int(dim()) << ",1)";
      pgeometric_trans pgt = geometric_trans_descriptor(name.str());

      for (size_type i = 0; i <= dim(); ++i)
	trans[cvr->structure()->ind_dir_points()[i]] 
	  = pgt->poly_vector()[i];
    }
  };

  static pgeometric_trans linear_product_gt(gt_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	  "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 1 || params[1].type() != 1)
      DAL_THROW(failure_error, "Bad type of parameters");
    pgeometric_trans a = params[0].method();
    pgeometric_trans b = params[1].method();
    return remember_for_cleanup(new cv_pr_tl_(a, b));
  }

  /* ******************************************************************** */
  /* parallelepiped transformation.                                       */
  /* ******************************************************************** */

  static pgeometric_trans QK_gt(gt_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    if (n <= 0 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");

    std::stringstream name;
    if (n == 1)
      name << "GT_PK(1," << k << ")";
    else 
      name << "GT_PRODUCT(GT_QK(" << n-1 << "," << k << "),GT_PK(1,"
	   << k << "))";
    return geometric_trans_descriptor(name.str());
  }

  static pgeometric_trans prism_gt(gt_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    if (n <= 1 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");

    std::stringstream name;
    name << "GT_PRODUCT(GT_PK(" << n-1 << "," << k << "),GT_PK(1,"
	 << k << "))";
    return geometric_trans_descriptor(name.str());
  }

  static pgeometric_trans linear_qk(gt_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	  "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    return parallelepiped_linear_geotrans(n);
  }

//   pgeometric_trans associated_trans(pconvex_structure cvs)
//   {
//     DAL::THROW(internal_error, "Obsolete function");
//     size_type n = cvs->dim(), nbp = cvs->nb_points();
//     if (nbp == n+1 && cvs == bgeot::simplex_structure(n))
//       return simplex_trans(n, 1);

//     if (nbp == (size_type(1) << n) && cvs==bgeot::parallelepiped_structure(n))
//       return parallelepiped_trans(n, 1);

//     if (nbp == 2 * n && cvs == bgeot::prism_structure(n))
// 	return prism_trans(n, 1);
    
//     // To be completed
    
//     DAL_THROW(to_be_done_error, 
// 	      "This element is not taken into account. Contact us");   
//     return NULL;
//   }

  /* norm of returned vector is the ratio between the face surface on
     the reel element and the face surface on the reference element 
     IT IS NOT UNITARY
     
     pt is the position of the evaluation point on the reference element
  */
  base_small_vector compute_normal(const geotrans_interpolation_context& c,
				   size_type face) {
    if (c.G().ncols() != c.pgt()->nb_points())
      DAL_THROW(dimension_error, "dimensions mismatch");
    base_small_vector un(c.N());
    gmm::mult(c.B(), c.pgt()->normals()[face], un);
    return un;
  }

  /*
    return the local basis (i.e. the norm in the first column, and the
    tangent vectors in the other columns 
  */
  base_matrix 
  compute_local_basis(const geotrans_interpolation_context& c,
		      size_type face) {
    if (c.G().ncols() != c.pgt()->nb_points()) DAL_THROW(dimension_error, "dimensions mismatch");
    base_small_vector up = c.pgt()->normals()[face];
    base_small_vector un(c.N());
    size_type P = c.pgt()->structure()->dim();

    base_matrix baseP(P, P);
    gmm::copy(gmm::identity_matrix(), baseP);
    size_type i0 = 0;
    for (size_type i=1; i < P; ++i) 
      if (gmm::abs(up[i])>gmm::abs(up[i0])) i0=i;
    if (i0) gmm::copy(gmm::mat_col(baseP, 0), gmm::mat_col(baseP, i0));
    gmm::copy(up, gmm::mat_col(baseP, 0));

    base_matrix baseN(c.N(), P);
    gmm::mult(c.B(), baseP, baseN);

    /* modified gram-schmidt */
    for (size_type k=0; k < P; ++k) {
      for (size_type l=0; l < k; ++l) {
	gmm::add(gmm::scaled(gmm::mat_col(baseN,l), 
			     -gmm::vect_sp(gmm::mat_col(baseN,l),gmm::mat_col(baseN,k))), 
		 gmm::mat_col(baseN,k));
      }
      gmm::scale(gmm::mat_col(baseN,k), 1./gmm::vect_norm2(gmm::mat_col(baseN,k)));
    }
    /* TODO: for cases where P < N,
       complete the basis */
    /* ensure that the baseN is direct */
    if (c.N() == P && c.N()>1 && gmm::lu_det(baseN) < 0) {
      gmm::scale(gmm::mat_col(baseN,1),-1.);
    }
    return baseN;
  }


  /* ******************************************************************** */
  /*    Naming system                                                     */
  /* ******************************************************************** */

  struct geometric_trans_naming_system : public ftool::naming_system<geometric_trans> {
    geometric_trans_naming_system() : 
      ftool::naming_system<geometric_trans>("GT") {
      add_suffix("PK", PK_gt);
      add_suffix("QK", QK_gt);
      add_suffix("PRISM", prism_gt);
      add_suffix("PRODUCT", product_gt);
      add_suffix("LINEAR_PRODUCT", linear_product_gt);
      add_suffix("LINEAR_QK", linear_qk);
    }
  };
  
  pgeometric_trans geometric_trans_descriptor(std::string name) {
    size_type i=0;
    return dal::singleton<geometric_trans_naming_system>::instance().method(name, i);
  }

  std::string name_of_geometric_trans(pgeometric_trans p) {
    return dal::singleton<geometric_trans_naming_system>::instance().shorter_name_of_method(p);
  }

  /* Fonctions pour la ref. directe.                                     */
  
  pgeometric_trans simplex_geotrans(size_type n, short_type k) {
    static pgeometric_trans pgt = 0;
    static size_type d = size_type(-2);
    static short_type r = short_type(-2);
    if (d != n || r != k) {
      std::stringstream name;
      name << "GT_PK(" << n << "," << k << ")";
      pgt = geometric_trans_descriptor(name.str());
      d = n; r = k;
    }
    return pgt;
  }

  pgeometric_trans parallelepiped_geotrans(size_type n, short_type k) {
    static pgeometric_trans pgt = 0;
    static size_type d = size_type(-2);
    static short_type r = short_type(-2);
    if (d != n || r != k) {
      std::stringstream name;
      name << "GT_QK(" << n << "," << k << ")";
      pgt = geometric_trans_descriptor(name.str());
      d = n; r = k;
    }
    return pgt;
  }

  static std::string name_of_linear_qk_trans(int dim) {
    switch (dim) {
    case 1: return "GT_PK(1,1)";
    default: return std::string("GT_LINEAR_PRODUCT(")
			   + name_of_linear_qk_trans(dim-1) 
			   + std::string(",GT_PK(1,1))");
    }
  }

  pgeometric_trans parallelepiped_linear_geotrans(size_type n) {
    static pgeometric_trans pgt = 0;
    static size_type d = size_type(-2);
    if (d != n) {
      std::stringstream name(name_of_linear_qk_trans(n));
      pgt = geometric_trans_descriptor(name.str());
      d = n;
    }
    return pgt;
  }

  pgeometric_trans prism_linear_geotrans(size_type n) {
    static pgeometric_trans pgt = 0;
    static size_type d = size_type(-2);
    if (d != n) {
      std::stringstream name;
      name << "GT_LINEAR_PRODUCT(GT_PK(" << (n-1) << ", 1), GT_PK(1,1))";
      pgt = geometric_trans_descriptor(name.str());
      d = n;
    }
    return pgt;
  }

  pgeometric_trans linear_product_geotrans(pgeometric_trans pg1,
					   pgeometric_trans pg2) {
    std::stringstream name;
    name << "GT_LINEAR_PRODUCT(" << name_of_geometric_trans(pg1) << "," 
	 << name_of_geometric_trans(pg2) << ")";
    return geometric_trans_descriptor(name.str());
  }

  pgeometric_trans prism_geotrans(size_type n, short_type k) {
    static pgeometric_trans pgt = 0;
    static size_type d = size_type(-2);
    static short_type r = short_type(-2);
    if (d != n || r != k) {
      std::stringstream name;
      name << "GT_PRISM(" << n << "," << k << ")";
      pgt = geometric_trans_descriptor(name.str());
      d = n; r = k;
    }
    return pgt;
  }

  pgeometric_trans product_geotrans(pgeometric_trans pg1,
				    pgeometric_trans pg2) {
    static pgeometric_trans pgt = 0;
    static pgeometric_trans pg1_ = 0;
    static pgeometric_trans pg2_ = 0;
    if (pg1 != pg1_ || pg2 != pg2_) {
      std::stringstream name;
      name << "GT_PRODUCT(" << name_of_geometric_trans(pg1) << "," 
	   << name_of_geometric_trans(pg2) << ")";
      pgt = geometric_trans_descriptor(name.str());
      pg1_ = pg1; pg2_ = pg2;
    }
    return pgt;    
  }

}  /* end of namespace bgeot.                                            */

