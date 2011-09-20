// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2000-2008 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================


#include "getfem/dal_singleton.h"
#include "getfem/dal_tree_sorted.h"
#include "getfem/dal_naming_system.h"
#include "getfem/bgeot_geometric_trans.h"
#include "getfem/bgeot_poly_composite.h"

namespace bgeot {

  const base_node& geotrans_interpolation_context::xref() const {
    if (!have_xref()) {
      if (pspt_) xref_ = (*pspt_)[ii_];
      else GMM_ASSERT1(false, "missing xref");
    }
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

  void geotrans_interpolation_context::compute_J(void) const {
    GMM_ASSERT1(have_G() && have_pgt(), "unable to compute B\n");
    size_type P = pgt_->structure()->dim();
    base_matrix CS(P,P);
    if (P != N()) {
      gmm::mult(gmm::transposed(K()), K(), CS);
      // gmm::abs below because on flat convexes determinant could be -1e-27.
      J_ = ::sqrt(gmm::abs(gmm::lu_det(CS)));
    }
    else
      J_ = gmm::abs(gmm::lu_det(K()));
  }

  const base_matrix& geotrans_interpolation_context::K() const {
    if (!have_K()) {
      GMM_ASSERT1(have_G() && have_pgt(), "unable to compute K\n");
      size_type P = pgt_->structure()->dim();
      K_.resize(N(), P);
      if (have_pgp()) {

        if (&pgp_->grad(ii_) == 0) { cerr << "OULA!! " << ii_ << "\n"; }
        else if (pgp_->grad(ii_).size() == 0) { cerr << "OUCH\n"; }

        assert(ii_ < pgp_->get_point_tab().size());

        gmm::mult(G(), pgp_->grad(ii_), K_);
      } else {
        base_matrix pc(pgt()->nb_points(), P);
        pgt()->poly_vector_grad(xref(), pc);
        gmm::mult(G(),pc,K_);
      }
    }
    return K_;
  }

  const base_matrix& geotrans_interpolation_context::B() const {
    if (!have_B()) {
      GMM_ASSERT1(have_G() && have_pgt(), "unable to compute K\n");
      size_type P = pgt_->structure()->dim();
      B_.resize(N(), P);
      if (P != N()) {
        base_matrix CS(P,P);
        gmm::mult(gmm::transposed(K()), K(), CS);
        // gmm::abs below because on flat convexes determinant could be -1e-27.
        J_ = ::sqrt(gmm::abs(gmm::lu_inverse(CS)));
        gmm::mult(K(), CS, B_);
      } else {
        gmm::copy(gmm::transposed(K()), B_);
        J_ = gmm::abs(gmm::lu_inverse(B_));
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
        if (have_pgp()) {
          gmm::mult(G(), pgp_->hessian(ii_), Htau);
        } else {
          /* very inefficient of course... */
          base_matrix pc(pgt()->nb_points(), P*P);
          pgt()->poly_vector_hess(xref(), pc);
          gmm::mult(G(), pc, Htau);
        }
        for (short_type i = 0; i < P; ++i)
          for (short_type j = 0; j < P; ++j)
            for (short_type k = 0; k < P; ++k)
              for (short_type l = 0; l < N_; ++l)
                B2(i + P*j, k) += Htau(l, i + P*j) * BB(l,k);
        gmm::mult(B3(), B2, B32_);
      } else gmm::clear(B32_);
    }
    return B32_;
  }

 void geotrans_interpolation_context::set_ii(size_type ii__) {
    if (ii_ == ii__) return;
    if (have_K() && !pgt()->is_linear()) { K_.resize(0,0); }
    if (have_B() && !pgt()->is_linear()) { B_.resize(0,0); }
    if (have_B3() && !pgt()->is_linear()) {
      B3_.resize(0,0); B32_.resize(0,0); }
    xref_.clear(); xreal_.clear();
    ii_=ii__;  J_ = scalar_type(-1);
  }

  void geotrans_interpolation_context::set_xref(const base_node& P) {
    xref_ = P;
    if (have_K() && !pgt()->is_linear()) { K_.resize(0,0); }
    if (have_B() && !pgt()->is_linear()) { B_.resize(0,0); }
    if (have_B3() && !pgt()->is_linear()) {
      B3_.resize(0,0); B32_.resize(0,0); }
    xreal_.clear(); ii_ = size_type(-1); J_ = scalar_type(-1); pspt_ = 0;
  }

  geotrans_interpolation_context::geotrans_interpolation_context() :
    G_(0), pgt_(0), pgp_(0), pspt_(0), ii_(size_type(-1)), J_(-1) {}
  geotrans_interpolation_context::geotrans_interpolation_context
  (bgeot::pgeotrans_precomp pgp__, size_type ii__, const base_matrix& G__) :
    G_(&G__), pgt_(pgp__->get_trans()), pgp_(pgp__),
    pspt_(&pgp__->get_point_tab()), ii_(ii__), J_(-1) {}
  geotrans_interpolation_context::geotrans_interpolation_context
  (bgeot::pgeometric_trans pgt__, bgeot::pstored_point_tab pspt__,
   size_type ii__,  const base_matrix& G__) :
    G_(&G__), pgt_(pgt__), pgp_(0),
    pspt_(pspt__), ii_(ii__), J_(-1) {}
  geotrans_interpolation_context::geotrans_interpolation_context
  (bgeot::pgeometric_trans pgt__, const base_node& xref__,
   const base_matrix& G__) :
    xref_(xref__), G_(&G__), pgt_(pgt__), pgp_(0), pspt_(0),
    ii_(size_type(-1)), J_(-1) {}


  typedef dal::naming_system<geometric_trans>::param_list gt_param_list;

  base_node geometric_trans::transform(const base_node &pt,
                                       const base_matrix &G) const {
    size_type N = G.nrows(), k = nb_points();
    base_node P(N); base_vector val(k);
    poly_vector_val(pt, val);
    base_matrix::const_iterator git = G.begin();
    for (size_type l = 0; l < k; ++l) {
      scalar_type a = val[l];
      base_node::iterator pit = P.begin(), pite = P.end();
      for (; pit != pite; ++git, ++pit) *pit += a * (*git);
    }
    return P;
  }

  void geometric_trans::fill_standard_vertices(void) {
    vertices_.resize(0);
    for (size_type ip = 0; ip < nb_points(); ++ip) {
      bool vertex = true;
      for (size_type i = 0; i < cvr->points()[ip].size(); ++i)
        if (gmm::abs(cvr->points()[ip][i]) > 1e-10
            && gmm::abs(cvr->points()[ip][i]-1.0) > 1e-10)
          { vertex = false; break; }
      if (vertex) vertices_.push_back(ip);
    }
    assert(vertices_.size() >= dim());
  }

  /* ******************************************************************** */
  /* Instantied geometric transformations.                                */
  /* ******************************************************************** */

  template <class FUNC>
  struct igeometric_trans : public geometric_trans {

    std::vector<FUNC> trans;

    virtual void poly_vector_val(const base_node &pt, base_vector &val) const {
      val.resize(nb_points());
      for (size_type k = 0; k < nb_points(); ++k)
        val[k] = trans[k].eval(pt.begin());
    }

    virtual void poly_vector_val(const base_node &pt, const convex_ind_ct &ind_ct,
                                 base_vector &val) const {
      size_type nb_funcs=ind_ct.size();
      val.resize(nb_funcs);
      for (size_type k = 0; k < nb_funcs; ++k)
        val[k] = trans[ind_ct[k]].eval(pt.begin());
    }

    virtual void poly_vector_grad(const base_node &pt, base_matrix &pc) const {
      FUNC PP;
      pc.resize(nb_points(),dim());
      for (size_type i = 0; i < nb_points(); ++i)
        for (dim_type n = 0; n < dim(); ++n) {
          PP = trans[i];
          PP.derivative(n);
          pc(i, n) = PP.eval(pt.begin());
        }
    }

    virtual void poly_vector_grad(const base_node &pt, const convex_ind_ct &ind_ct,
                                  base_matrix &pc) const {
      FUNC PP;
      size_type nb_funcs=ind_ct.size();
      pc.resize(nb_funcs,dim());
      for (size_type i = 0; i < nb_funcs; ++i)
        for (dim_type n = 0; n < dim(); ++n) {
          PP = trans[ind_ct[i]];
          PP.derivative(n);
          pc(i, n) = PP.eval(pt.begin());
        }
    }

    virtual void poly_vector_hess(const base_node &pt, base_matrix &pc) const {
      FUNC PP, QP;
      pc.resize(nb_points(),dim()*dim());
      for (size_type i = 0; i < nb_points(); ++i)
        for (dim_type n = 0; n < dim(); ++n) {
          QP = trans[i]; QP.derivative(n);
          for (dim_type m = 0; m <= n; ++m) {
            PP = QP; PP.derivative(m);
            pc(i, n*dim()+m) = pc(i, m*dim()+n) = PP.eval(pt.begin());
          }
        }
    }

  };

  typedef igeometric_trans<base_poly> poly_geometric_trans;
  typedef igeometric_trans<polynomial_composite> comppoly_geometric_trans;

  /* ******************************************************************** */
  /* transformation on simplex.                                           */
  /* ******************************************************************** */

  struct simplex_trans_ : public poly_geometric_trans {
    void calc_base_func(base_poly &p, size_type i, short_type K) const {
      dim_type N = dim();
      base_poly l0(N, 0), l1(N, 0);
      power_index w(short_type(N+1));
      l0.one(); l1.one(); p = l0;
      for (short_type nn = 0; nn < N; ++nn) l0 -= base_poly(N, 1, nn);

      w[0] = K;
      for (int nn = 1; nn <= N; ++nn) {
        w[nn]=short_type(floor(0.5+(((cvr->points())[i])[nn-1]*double(K))));
        w[0]=short_type(w[0]-w[nn]);
      }

      for (short_type nn = 0; nn <= N; ++nn)
        for (short_type j = 0; j < w[nn]; ++j)
          if (nn == 0)
            p *= (l0 * (scalar_type(K) / scalar_type(j+1)))
               - (l1 * (scalar_type(j) / scalar_type(j+1)));
          else
            p *= (base_poly(N, 1, short_type(nn-1)) * (scalar_type(K) / scalar_type(j+1)))
               - (l1 * (scalar_type(j) / scalar_type(j+1)));
    }

    simplex_trans_(dim_type nc, short_type k) {
      cvr = simplex_of_reference(nc, k);
      size_type R = cvr->structure()->nb_points();
      is_lin = (k == 1);
      complexity_ = k;
      trans.resize(R);
      for (size_type r = 0; r < R; ++r) calc_base_func(trans[r], r, k);
      fill_standard_vertices();
    }
  };

  static pgeometric_trans
  PK_gt(gt_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0,
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    GMM_ASSERT1(n >= 0 && n < 100 && k >= 0 && k <= 150 &&
                double(n) == params[0].num() && double(k) == params[1].num(),
                "Bad parameters");
    dependencies.push_back(simplex_of_reference(dim_type(n), dim_type(k)));
    return new simplex_trans_(dim_type(n), dim_type(k));
  }

  /* ******************************************************************** */
  /* direct product transformation                                        */
  /* ******************************************************************** */

  struct cv_pr_t_ : public poly_geometric_trans {
    cv_pr_t_(const poly_geometric_trans *a, const poly_geometric_trans *b) {
      cvr = convex_ref_product(a->convex_ref(), b->convex_ref());
      is_lin = false;
      complexity_ = a->complexity() * b->complexity();

      size_type n1 = a->nb_points(), n2 = b->nb_points();
      trans.resize(n1 * n2);
      for (size_type i1 = 0; i1 < n1; ++i1)
        for (size_type i2 = 0; i2 < n2; ++i2) {
          trans[i1 + i2 * n1] = a->trans[i1];
          trans[i1 + i2 * n1].direct_product(b->trans[i2]);
        }
      for (size_type i2 = 0; i2 < b->nb_vertices(); ++i2)
        for (size_type i1 = 0; i1 < a->nb_vertices(); ++i1)
          vertices_.push_back(a->vertices()[i1] + b->vertices()[i2] * n1);
    }
  };

  static pgeometric_trans product_gt(gt_param_list &params,
                  std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 1 && params[1].type() == 1,
                "Bad type of parameters");
    pgeometric_trans a = params[0].method();
    pgeometric_trans b = params[1].method();
    dependencies.push_back(a); dependencies.push_back(b);
    dependencies.push_back(convex_ref_product(a->convex_ref(),
                                              b->convex_ref()));
    const poly_geometric_trans *aa
      = dynamic_cast<const poly_geometric_trans *>(a.get());
    const poly_geometric_trans *bb
      = dynamic_cast<const poly_geometric_trans *>(b.get());
    GMM_ASSERT1(aa && bb, "The product of geometric transformations "
                "is only defined for polynomial ones");
    return new cv_pr_t_(aa, bb);
  }

  /* ******************************************************************** */
  /* linear direct product transformation.                                */
  /* ******************************************************************** */

  struct cv_pr_tl_ : public poly_geometric_trans {
    cv_pr_tl_(const poly_geometric_trans *a, const poly_geometric_trans *b) {
      GMM_ASSERT1(a->is_linear() && b->is_linear(),
                  "linear product of non-linear transformations");
      cvr = convex_ref_product(a->convex_ref(), b->convex_ref());
      is_lin = true;
      complexity_ = std::max(a->complexity(), b->complexity());

      trans.resize(a->nb_points() * b->nb_points());
      std::fill(trans.begin(), trans.end(), null_poly(dim()));

      std::stringstream name;
      name << "GT_PK(" << int(dim()) << ",1)";
      pgeometric_trans pgt_ = geometric_trans_descriptor(name.str());
      const poly_geometric_trans *pgt
      = dynamic_cast<const poly_geometric_trans *>(pgt_.get());

      for (size_type i = 0; i <= dim(); ++i)
        trans[cvr->structure()->ind_dir_points()[i]]
          = pgt->trans[i];
      for (size_type i2 = 0; i2 < b->nb_vertices(); ++i2)
        for (size_type i1 = 0; i1 < a->nb_vertices(); ++i1)
          vertices_.push_back(a->vertices()[i1]
                              + b->vertices()[i2] * a->nb_points());
    }
  };

  static pgeometric_trans linear_product_gt(gt_param_list &params,
        std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 1 && params[1].type() == 1,
                "Bad type of parameters");
    pgeometric_trans a = params[0].method();
    pgeometric_trans b = params[1].method();
    dependencies.push_back(a); dependencies.push_back(b);
    dependencies.push_back(convex_ref_product(a->convex_ref(),
                                              b->convex_ref()));
    const poly_geometric_trans *aa
      = dynamic_cast<const poly_geometric_trans *>(a.get());
    const poly_geometric_trans *bb
      = dynamic_cast<const poly_geometric_trans *>(b.get());
    GMM_ASSERT1(aa && bb, "The product of geometric transformations "
                "is only defined for polynomial ones");
    return new cv_pr_tl_(aa, bb);
  }

  /* ******************************************************************** */
  /* parallelepiped transformation.                                       */
  /* ******************************************************************** */

  static pgeometric_trans QK_gt(gt_param_list &params,
        std::vector<dal::pstatic_stored_object> &) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0,
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    GMM_ASSERT1(n > 0 && n < 100 && k >= 0 && k <= 150 &&
                double(n) == params[0].num() && double(k) == params[1].num(),
                "Bad parameters");
    std::stringstream name;
    if (n == 1)
      name << "GT_PK(1," << k << ")";
    else
      name << "GT_PRODUCT(GT_QK(" << n-1 << "," << k << "),GT_PK(1,"
           << k << "))";
    return geometric_trans_descriptor(name.str());
  }

  static pgeometric_trans prism_gt(gt_param_list &params,
        std::vector<dal::pstatic_stored_object> &) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
                << params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0,
                "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    GMM_ASSERT1(n > 0 && n < 100 && k >= 0 && k <= 150 &&
                double(n) == params[0].num() && double(k) == params[1].num(),
                "Bad parameters");
    std::stringstream name;
    name << "GT_PRODUCT(GT_PK(" << n-1 << "," << k << "),GT_PK(1,"
         << k << "))";
    return geometric_trans_descriptor(name.str());
  }

  static pgeometric_trans linear_qk(gt_param_list &params,
        std::vector<dal::pstatic_stored_object> &) {
    GMM_ASSERT1(params.size() == 1, "Bad number of parameters : "
                << params.size() << " should be 1.");
    GMM_ASSERT1(params[0].type() == 0, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    return parallelepiped_linear_geotrans(n);
  }

  /* norm of returned vector is the ratio between the face surface on
     the reel element and the face surface on the reference element
     IT IS NOT UNITARY

     pt is the position of the evaluation point on the reference element
  */
  base_small_vector compute_normal(const geotrans_interpolation_context& c,
                                   size_type face) {
    GMM_ASSERT1(c.G().ncols() == c.pgt()->nb_points(), "dimensions mismatch");
    base_small_vector un(c.N());
    gmm::mult(c.B(), c.pgt()->normals()[face], un);
    return un;
  }

  /*
    return the local basis (i.e. the normal in the first column, and the
    tangent vectors in the other columns
  */
  base_matrix
  compute_local_basis(const geotrans_interpolation_context& c,
                      size_type face) {
    GMM_ASSERT1(c.G().ncols() == c.pgt()->nb_points(), "dimensions mismatch");
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
                             -gmm::vect_sp(gmm::mat_col(baseN,l),
                                           gmm::mat_col(baseN,k))),
                 gmm::mat_col(baseN,k));
      }
      gmm::scale(gmm::mat_col(baseN,k),
                 1./gmm::vect_norm2(gmm::mat_col(baseN,k)));
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

  struct geometric_trans_naming_system
    : public dal::naming_system<geometric_trans> {
    geometric_trans_naming_system() :
      dal::naming_system<geometric_trans>("GT") {
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

  static std::string name_of_linear_qk_trans(size_type dim) {
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

  /* ********************************************************************* */
  /*       Precomputation on geometric transformations.                    */
  /* ********************************************************************* */

  DAL_DOUBLE_KEY(pre_geot_key_, pgeometric_trans, pstored_point_tab);

  geotrans_precomp_::geotrans_precomp_(pgeometric_trans pg,
                                       pstored_point_tab ps)
    : pgt(pg), pspt(ps) {}

  void geotrans_precomp_::init_val() const {
    c.clear();
    c.resize(pspt->size(), base_vector(pgt->nb_points()));
    for (size_type j = 0; j < pspt->size(); ++j)
      pgt->poly_vector_val((*pspt)[j], c[j]);
  }

  void geotrans_precomp_::init_grad() const {
    dim_type N = pgt->dim();
    pc.clear();
    pc.resize(pspt->size(), base_matrix(pgt->nb_points() , N));
    for (size_type j = 0; j < pspt->size(); ++j)
      pgt->poly_vector_grad((*pspt)[j], pc[j]);
  }

  void geotrans_precomp_::init_hess() const {
    base_poly P, Q;
    dim_type N = pgt->structure()->dim();
    hpc.clear();
    hpc.resize(pspt->size(), base_matrix(pgt->nb_points(), gmm::sqr(N)));
    for (size_type j = 0; j < pspt->size(); ++j)
      pgt->poly_vector_hess((*pspt)[j], hpc[j]);
  }

  base_node geotrans_precomp_::transform(size_type i,
                                         const base_matrix &G) const {
    if (c.empty()) init_val();
    size_type N = G.nrows(), k = pgt->nb_points();
    base_node P(N);
    base_matrix::const_iterator git = G.begin();
    for (size_type l = 0; l < k; ++l) {
      scalar_type a = c[i][l];
      base_node::iterator pit = P.begin(), pite = P.end();
      for (; pit != pite; ++git, ++pit) *pit += a * (*git);
    }
    return P;
  }

  pgeotrans_precomp geotrans_precomp(pgeometric_trans pg,
                                     pstored_point_tab pspt,
                                     dal::pstatic_stored_object dep) {
    dal::pstatic_stored_object o
      = dal::search_stored_object(pre_geot_key_(pg, pspt));
    if (o) return dal::stored_cast<geotrans_precomp_>(o);
    pgeotrans_precomp p = new geotrans_precomp_(pg, pspt);
    dal::add_stored_object(new pre_geot_key_(pg, pspt), p, pg, pspt,
                           dal::AUTODELETE_STATIC_OBJECT);
    if (dep) dal::add_dependency(p, dep);
    return p;
  }

  void delete_geotrans_precomp(pgeotrans_precomp pgp)
  { dal::del_stored_object(pgp); }

}  /* end of namespace bgeot.                                            */

