/*===========================================================================

 Copyright (C) 2002-2020 Yves Renard

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/


#include "getfem/bgeot_poly_composite.h"
#include "getfem/getfem_integration.h"
#include "getfem/getfem_mesh_fem.h"
#include "getfem/dal_naming_system.h"

namespace getfem {
 
  typedef const fem<bgeot::polynomial_composite> * ppolycompfem;

  struct polynomial_composite_fem : public fem<bgeot::polynomial_composite> {
    mesh m;
    mutable bgeot::mesh_precomposite mp;
    mesh_fem mf;
    mutable bgeot::pgeotrans_precomp pgp;
    mutable bgeot::pgeometric_trans pgt_stored;
    bgeot::pstored_point_tab pspt;

    
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
                           bgeot::pgeometric_trans pgt) const;

    
    polynomial_composite_fem(const mesh &m_, const mesh_fem &mf_)
      : m(m_), mp(m), mf(m), pgp(0), pgt_stored(0) {
      for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv)
        mf.set_finite_element(cv, mf_.fem_of_element(cv));
      mf.nb_dof();
      pspt = store_point_tab(m.points());
      // verification for the non-equivalent fems
      for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
        pfem pf1 = mf.fem_of_element(cv);
        if (!(pf1->is_equivalent())) {
          dal::bit_vector unshareable;
          for (const auto &i : mf.ind_basic_dof_of_element(cv))
            unshareable.add(i);
       
          for (dal::bv_visitor cv2(mf.convex_index()); !cv2.finished(); ++cv2) {
            if (cv2 != cv)
              for (const auto &i : mf.ind_basic_dof_of_element(cv2))
                GMM_ASSERT1(!(unshareable.is_in(i)), "Non equivalent elements "
                            "are allowed only if they do not share their dofs");
          }
        }
      }
    }
  };

  void polynomial_composite_fem::mat_trans(base_matrix &M, const base_matrix &G,
                                           bgeot::pgeometric_trans pgt) const {
    dim_type N = dim_type(G.nrows());
    gmm::copy(gmm::identity_matrix(), M);
    base_matrix G1, M1;

    if (pgt != pgt_stored) {
      pgt_stored = pgt;
      pgp = bgeot::geotrans_precomp(pgt, pspt, 0);
    }
    
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
      pfem pf1 = mf.fem_of_element(cv);
      if (!(pf1->is_equivalent())) {
        size_type npt=m.nb_points_of_convex(cv);
        size_type ndof=mf.nb_basic_dof_of_element(cv);
        GMM_ASSERT1(ndof == pf1->nb_base(0) && ndof == pf1->nb_dof(0),
                    "Sorry, limited implementation for the moment");
        // Compute the local G
        G1.resize(npt, N);
        M1.resize(ndof, ndof); 
        for (size_type i = 0; i < npt; ++i)
          gmm::copy(pgp->transform(m.ind_points_of_convex(cv)[i], G),
                    gmm::mat_col(G1, i));
        // Call for the local M
        pf1->mat_trans(M1, G1, m.trans_of_convex(cv));
        gmm::sub_index I(mf.ind_basic_dof_of_element(cv));
        // I is in fact an interval ...
        gmm::copy(M1, gmm::sub_matrix(M, I, I));
      }
    }
  }

  static pfem composite_fe_method(const getfem::mesh &m, 
				  const mesh_fem &mf, bgeot::pconvex_ref cr,
                                  bool ff = false) {
    
    GMM_ASSERT1(!mf.is_reduced(),
		"Sorry, does not work for reduced mesh_fems");
    auto p = std::make_shared<polynomial_composite_fem>(m, mf);
    
    p->mref_convex() = cr;
    p->dim() = cr->structure()->dim();
    p->is_polynomialcomp() = p->is_equivalent() = p->is_standard() = true;
    p->is_polynomial() = false;
    p->is_lagrange() = true;
    p->estimated_degree() = 0;
    p->init_cvs_node();

    std::vector<bgeot::polynomial_composite> base(mf.nb_basic_dof());
    std::fill(base.begin(), base.end(),
              bgeot::polynomial_composite(p->mp, true, ff));
    std::vector<pdof_description> dofd(mf.nb_basic_dof());
    
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
      pfem pf1 = mf.fem_of_element(cv);
      if (!pf1->is_lagrange()) p->is_lagrange() = false;
      if (!(pf1->is_equivalent())) p->is_equivalent() = false;
      
      GMM_ASSERT1(pf1->is_polynomial(), "Only for polynomial fems.");
      ppolyfem pf = ppolyfem(pf1.get());
      p->estimated_degree() = std::max(p->estimated_degree(),
				       pf->estimated_degree());
      for (size_type k = 0; k < pf->nb_dof(cv); ++k) {
	size_type igl = mf.ind_basic_dof_of_element(cv)[k];
	base_poly fu = pf->base()[k];
	base[igl].set_poly_of_subelt(cv, fu);
	dofd[igl] = pf->dof_types()[k];
      }
    }
    p->base().resize(mf.nb_basic_dof());
    for (size_type k = 0; k < mf.nb_basic_dof(); ++k) {  
      p->add_node(dofd[k], mf.point_of_basic_dof(k));
      p->base()[k] = base[k];
    }
    return p;
  }

  typedef dal::naming_system<virtual_fem>::param_list fem_param_list;

  pfem structured_composite_fem_method(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 2, "Bad number of parameters : "
		<< params.size() << " should be 2.");
    GMM_ASSERT1(params[0].type() == 1 && params[1].type() == 0,
		"Bad type of parameters");
    pfem pf = params[0].method();
    int k = int(::floor(params[1].num() + 0.01));
    GMM_ASSERT1(((pf->is_polynomial()) || !(pf->is_equivalent())) && k > 0
		&& k <= 150 && double(k) == params[1].num(), "Bad parameters");
    bgeot::pbasic_mesh pm;
    bgeot::pmesh_precomposite pmp;

    structured_mesh_for_convex(pf->ref_convex(0), short_type(k), pm, pmp);

    mesh m(*pm);
    mesh_fem mf(m);
    mf.set_finite_element(pm->convex_index(), pf);
    pfem p = composite_fe_method(m, mf, pf->ref_convex(0));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  pfem PK_composite_hierarch_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &) {
    GMM_ASSERT1(params.size() == 3, "Bad number of parameters : "
		<< params.size() << " should be 3.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0 &&
		params[2].type() == 0, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    int s = int(::floor(params[2].num() + 0.01)), t;
    GMM_ASSERT1(n > 0 && n < 100 && k > 0 && k <= 150 && s > 0 && s <= 150 &&
	(!(s & 1) || (s == 1)) && double(s) == params[2].num() &&
		double(n) == params[0].num() && double(k) == params[1].num(),
		"Bad parameters");
    std::stringstream name;
    if (s == 1) 
      name << "FEM_STRUCTURED_COMPOSITE(FEM_PK(" << n << "," << k << "),1)";
    else {
      for (t = 2; t <= s; ++t) if ((s % t) == 0) break;
      name << "FEM_GEN_HIERARCHICAL(FEM_PK_HIERARCHICAL_COMPOSITE(" << n
	   << "," << k << "," << s/t << "), FEM_STRUCTURED_COMPOSITE(FEM_PK("
	   << n << "," << k << ")," << s << "))";
    }
    return fem_descriptor(name.str());
  }

    pfem PK_composite_full_hierarch_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &) {
    GMM_ASSERT1(params.size() == 3, "Bad number of parameters : "
		<< params.size() << " should be 3.");
    GMM_ASSERT1(params[0].type() == 0 && params[1].type() == 0 &&
		params[2].type() == 0, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    int s = int(::floor(params[2].num() + 0.01)), t;
    GMM_ASSERT1(n > 0 && n < 100 && k > 0 && k <= 150 && s > 0 && s <= 150 &&
	(!(s & 1) || (s == 1)) && double(s) == params[2].num() &&
		double(n) == params[0].num() && double(k) == params[1].num(),
		"Bad parameters");
    std::stringstream name;
    if (s == 1) 
      name << "FEM_STRUCTURED_COMPOSITE(FEM_PK_HIERARCHICAL(" << n << ","
	   << k << "),1)";
    else {
      for (t = 2; t <= s; ++t) if ((s % t) == 0) break;
      name << "FEM_GEN_HIERARCHICAL(FEM_PK_FULL_HIERARCHICAL_COMPOSITE(" << n
	   << "," << k << "," << s/t
	   << "), FEM_STRUCTURED_COMPOSITE(FEM_PK_HIERARCHICAL("
	   << n << "," << k << ")," << s << "))";
    }
    return fem_descriptor(name.str());
  }
  

  /* ******************************************************************** */
  /*    P1 with piecewise linear bubble function on a triangle.           */
  /* ******************************************************************** */

  struct P1bubbletriangle__ : public fem<bgeot::polynomial_composite> {
    mesh m;
    bgeot::mesh_precomposite mp;
    P1bubbletriangle__(void);
  };

  P1bubbletriangle__::P1bubbletriangle__(void) {

    m.clear();
    size_type i0 = m.add_point(base_node(1.0/3.0, 1.0/3.0));
    size_type i1 = m.add_point(base_node(0.0, 0.0));
    size_type i2 = m.add_point(base_node(1.0, 0.0));
    size_type i3 = m.add_point(base_node(0.0, 1.0));
    m.add_triangle(i0, i2, i3);
    m.add_triangle(i0, i3, i1);
    m.add_triangle(i0, i1, i2);
    mp.initialise(m);

    std::stringstream s("1-x-y;1-x-y;1-x-y;x;x;x;y;y;y;3-3*x-3*y;3*x;3*y;");

    bgeot::pconvex_ref cr = bgeot::simplex_of_reference(2);
    mref_convex() = cr;
    dim() = cr->structure()->dim();
    is_polynomialcomp() = true;
    is_equivalent() = true;
    is_polynomial() = false;
    is_lagrange() = false;
    is_standard() = true;
    estimated_degree() = 3;
    init_cvs_node();

    base()=std::vector<bgeot::polynomial_composite>
      (4, bgeot::polynomial_composite(mp, false));
    for (size_type k = 0; k < 4; ++k)
      for (size_type ic = 0; ic < 3; ++ic)
        base()[k].set_poly_of_subelt(ic, bgeot::read_base_poly(2, s));

    for (size_type i = 0; i < 3; ++i) {
      base_node pt(0.0, 0.0);
      if (i) pt[i-1] = 1.0;
      add_node(lagrange_dof(2), pt);
    }

    add_node(bubble1_dof(2), base_node(1.0/3.0, 1.0/3.0));
  }


  pfem P1bubbletriangle_fem
  (fem_param_list &params,
   std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 0, "Bad number of parameters : "
		<< params.size() << " should be 0.");
    pfem p = std::make_shared<P1bubbletriangle__>();
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*    Hsieh-Clough-Tocher C^1 element (composite P3)                    */
  /* ******************************************************************** */

  struct HCT_triangle__ : public fem<bgeot::polynomial_composite> {
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
			   bgeot::pgeometric_trans pgt) const;
    mesh m;
    mutable bgeot::base_small_vector true_normals[3];
    mutable bgeot::mesh_precomposite mp;
    mutable bgeot::pgeotrans_precomp pgp;
    mutable pfem_precomp pfp;
    mutable bgeot::pgeometric_trans pgt_stored;
    mutable base_matrix K;

    HCT_triangle__(void);
  };

  void HCT_triangle__::mat_trans(base_matrix &M, const base_matrix &G,
				 bgeot::pgeometric_trans pgt) const {
    
    dim_type N = dim_type(G.nrows());
    
    GMM_ASSERT1(N == 2, "Sorry, this version of HCT "
		"element works only on dimension two.");
    if (pgt != pgt_stored) {
      pgt_stored = pgt;
      pgp = bgeot::geotrans_precomp(pgt, node_tab(0), 0);
      pfp = fem_precomp(std::make_shared<HCT_triangle__>(), node_tab(0), 0);
    }
    gmm::copy(gmm::identity_matrix(), M);
    
    gmm::mult(G, pgp->grad(0), K);
    for (size_type i = 0; i < 3; ++i) {
      if (i && !(pgt->is_linear())) gmm::mult(G, pgp->grad(3*i), K);  
      gmm::copy(K, gmm::sub_matrix(M, gmm::sub_interval(1+3*i, 2)));
    }

    // take the normal derivatives into account
    static base_matrix W(3, 12);
    base_small_vector norient(M_PI, M_PI * M_PI);
    if (pgt->is_linear()) gmm::lu_inverse(K); 
    for (unsigned i = 9; i < 12; ++i) {
      if (!(pgt->is_linear()))
	{ gmm::mult(G, pgp->grad(i), K); gmm::lu_inverse(K); }
      bgeot::base_small_vector n(2), v(2);
      gmm::mult(gmm::transposed(K), cvr->normals()[i-9], n);
      n /= gmm::vect_norm2(n);

      scalar_type ps = gmm::vect_sp(n, norient);
      if (ps < 0) n *= scalar_type(-1);
      true_normals[i-9] = n;

      if (gmm::abs(ps) < 1E-8)
	GMM_WARNING2("HCT_triangle : "
		     "The normal orientation may be not correct");
      gmm::mult(K, n, v);
      const bgeot::base_tensor &t = pfp->grad(i);
      // cout << "t = " << t << endl;
      for (unsigned j = 0; j < 12; ++j)
	W(i-9, j) = t(j, 0, 0) * v[0] + t(j, 0, 1) * v[1];
    }
    
    static base_matrix A(3, 3);
    static bgeot::base_vector w(3), coeff(3);
    static gmm::sub_interval SUBI(9, 3), SUBJ(0, 3);
    gmm::copy(gmm::sub_matrix(W, SUBJ, SUBI), A);
    gmm::lu_inverse(A);
    gmm::copy(gmm::transposed(A), gmm::sub_matrix(M, SUBI));

    for (unsigned j = 0; j < 9; ++j) {
      gmm::mult(W, gmm::mat_row(M, j), w);
      gmm::mult(A, gmm::scaled(w, -1.0), coeff);
      gmm::copy(coeff, gmm::sub_vector(gmm::mat_row(M, j), SUBI));
    }
  }

  HCT_triangle__::HCT_triangle__(void) : pgt_stored(0), K(2, 2) {

    m.clear();
    size_type i0 = m.add_point(base_node(1.0/3.0, 1.0/3.0));
    size_type i1 = m.add_point(base_node(0.0, 0.0));
    size_type i2 = m.add_point(base_node(1.0, 0.0));
    size_type i3 = m.add_point(base_node(0.0, 1.0));
    m.add_triangle(i0, i2, i3);
    m.add_triangle(i0, i3, i1);
    m.add_triangle(i0, i1, i2);
    mp.initialise(m);


 std::stringstream s
   ("-1 + 9*x + 9*y - 15*x^2 - 30*x*y - 15*y^2 + 7*x^3 + 21*x^2*y + 21*x*y^2 + 7*y^3;"
    "1 - 3*x^2 - 3*y^2 + 3*x^3 - 3*x^2*y + 2*y^3;"
    "1 - 3*x^2 - 3*y^2 + 2*x^3 - 3*x*y^2 + 3*y^3;"
    "-1/6 + 5/2*x - 9/2*x^2 - 4*x*y + 1/2*y^2 + 13/6*x^3 + 4*x^2*y + 3/2*x*y^2 - 1/3*y^3;"
    "x - 1/2*x^2 - 3*x*y - 7/6*x^3 + 2*x^2*y + 2*x*y^2;"
    "x - 2*x^2 - 3/2*y^2 + x^3 - 1/2*x*y^2 + 7/3*y^3;"
    "-1/6 + 5/2*y + 1/2*x^2 - 4*x*y - 9/2*y^2 - 1/3*x^3 + 3/2*x^2*y + 4*x*y^2 + 13/6*y^3;"
    "y - 3/2*x^2 - 2*y^2 + 7/3*x^3 - 1/2*x^2*y + y^3;"
    "y - 3*x*y - 1/2*y^2 + 2*x^2*y + 2*x*y^2 - 7/6*y^3;"
    "1 - 9/2*x - 9/2*y + 9*x^2 + 15*x*y + 6*y^2 - 9/2*x^3 - 21/2*x^2*y - 21/2*x*y^2 - 5/2*y^3;"
    "3*x^2 - 5/2*x^3 + 3/2*x^2*y;"
    "3*x^2 - 2*x^3 + 3/2*x*y^2 - 1/2*y^3;"
    "-1/6 + 3/4*x + 3/4*y - 2*x^2 - 5/2*x*y - y^2 + 17/12*x^3 + 7/4*x^2*y + 7/4*x*y^2 + 5/12*y^3;"
    "-x^2 + 13/12*x^3 - 1/4*x^2*y;"
    "-x^2 + x^3 - 1/4*x*y^2 + 1/12*y^3;"
    "2/3 - 13/4*x - 11/4*y + 9/2*x^2 + 19/2*x*y + 7/2*y^2 - 23/12*x^3 - 23/4*x^2*y - 25/4*x*y^2 - 17/12*y^3;"
    "-1/2*x^2 + 5/12*x^3 + 9/4*x^2*y;"
    "-x*y + 1/2*y^2 + 2*x^2*y + 7/4*x*y^2 - 13/12*y^3;"
    "1 - 9/2*x - 9/2*y + 6*x^2 + 15*x*y + 9*y^2 - 5/2*x^3 - 21/2*x^2*y - 21/2*x*y^2 - 9/2*y^3;"
    "3*y^2 - 1/2*x^3 + 3/2*x^2*y - 2*y^3;"
    "3*y^2 + 3/2*x*y^2 - 5/2*y^3;"
    "2/3 - 11/4*x - 13/4*y + 7/2*x^2 + 19/2*x*y + 9/2*y^2 - 17/12*x^3 - 25/4*x^2*y - 23/4*x*y^2 - 23/12*y^3;"
    "1/2*x^2 - x*y - 13/12*x^3 + 7/4*x^2*y + 2*x*y^2;"
    "-1/2*y^2 + 9/4*x*y^2 + 5/12*y^3;"
    "-1/6 + 3/4*x + 3/4*y - x^2 - 5/2*x*y - 2*y^2 + 5/12*x^3 + 7/4*x^2*y + 7/4*x*y^2 + 17/12*y^3;"
    "-y^2 + 1/12*x^3 - 1/4*x^2*y + y^3;"
    "-y^2 - 1/4*x*y^2 + 13/12*y^3;"
    "-sqrt(2)*2/3 + sqrt(2)*3*x + sqrt(2)*3*y - sqrt(2)*4*x^2 - sqrt(2)*10*x*y - sqrt(2)*4*y^2 + sqrt(2)*5/3*x^3 + sqrt(2)*7*x^2*y + sqrt(2)*7*x*y^2 + sqrt(2)*5/3*y^3;"
    "sqrt(2)*1/3*x^3 - sqrt(2)*x^2*y;"
    "-sqrt(2)*x*y^2 + sqrt(2)*1/3*y^3;"
    "2/3 - 2*x - 4*y + 2*x^2 + 8*x*y + 6*y^2 - 2/3*x^3 - 4*x^2*y - 6*x*y^2 - 8/3*y^3;"
    "2*x^2 - 4*x*y - 10/3*x^3 + 4*x^2*y + 4*x*y^2;"
    "-2*y^2 + 2*x*y^2 + 8/3*y^3;"
    "2/3 - 4*x - 2*y + 6*x^2 + 8*x*y + 2*y^2 - 8/3*x^3 - 6*x^2*y - 4*x*y^2 - 2/3*y^3;"
    "-2*x^2 + 8/3*x^3 + 2*x^2*y;"
    "-4*x*y + 2*y^2 + 4*x^2*y + 4*x*y^2 - 10/3*y^3;");

    bgeot::pconvex_ref cr = bgeot::simplex_of_reference(2);
    mref_convex() = cr;
    dim() = cr->structure()->dim();
    is_polynomialcomp() = true;
    is_equivalent() = false;
    is_polynomial() = false;
    is_lagrange() = false;
    is_standard() = false;
    estimated_degree() = 5;
    init_cvs_node();

    base()=std::vector<bgeot::polynomial_composite>
      (12, bgeot::polynomial_composite(mp, false));
    for (size_type k = 0; k < 12; ++k)
      for (size_type ic = 0; ic < 3; ++ic)
        base()[k].set_poly_of_subelt(ic, bgeot::read_base_poly(2, s));

    for (size_type i = 0; i < 3; ++i) {
      base_node pt(0.0, 0.0);
      if (i) pt[i-1] = 1.0;
      add_node(lagrange_dof(2), pt);
      add_node(derivative_dof(2, 0), pt);
      add_node(derivative_dof(2, 1), pt);
    }

    add_node(normal_derivative_dof(2), base_node(0.5, 0.5));
    add_node(normal_derivative_dof(2), base_node(0.0, 0.5));
    add_node(normal_derivative_dof(2), base_node(0.5, 0.0));
  }

  pfem HCT_triangle_fem
  (fem_param_list &params,
   std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 0, "Bad number of parameters : "
		<< params.size() << " should be 0.");
    pfem p = std::make_shared<HCT_triangle__>();
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }


  /* ******************************************************************** */
  /*    Reduced Hsieh-Clough-Tocher C^1 element (composite P3)            */
  /* ******************************************************************** */

  struct reduced_HCT_triangle__ : public fem<bgeot::polynomial_composite> {
    const HCT_triangle__ *HCT;
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
			   bgeot::pgeometric_trans pgt) const;
    virtual size_type nb_base(size_type) const { return 12; }
    mutable base_matrix P, Mhct;
    reduced_HCT_triangle__(void);
  };

  void reduced_HCT_triangle__::mat_trans(base_matrix &M, const base_matrix &G,
				 bgeot::pgeometric_trans pgt) const {
    HCT->mat_trans(Mhct, G, pgt);
    
    P(10, 1)=HCT->true_normals[1][0]*0.5; P(11, 1)=HCT->true_normals[2][0]*0.5;
    P(10, 2)=HCT->true_normals[1][1]*0.5; P(11, 2)=HCT->true_normals[2][1]*0.5;

    P( 9, 4)=HCT->true_normals[0][0]*0.5; P(11, 4)=HCT->true_normals[2][0]*0.5;
    P( 9, 5)=HCT->true_normals[0][1]*0.5; P(11, 5)=HCT->true_normals[2][1]*0.5;

    P( 9, 7)=HCT->true_normals[0][0]*0.5; P(10, 7)=HCT->true_normals[1][0]*0.5;
    P( 9, 8)=HCT->true_normals[0][1]*0.5; P(10, 8)=HCT->true_normals[1][1]*0.5;

    gmm::mult(gmm::transposed(P), Mhct, M);
  }

  reduced_HCT_triangle__::reduced_HCT_triangle__(void)
    : P(12, 9), Mhct(12, 12) {
    HCT = dynamic_cast<const HCT_triangle__ *>
      (&(*fem_descriptor("FEM_HCT_TRIANGLE")));

    bgeot::pconvex_ref cr = bgeot::simplex_of_reference(2);
    mref_convex() = cr;
    dim() = cr->structure()->dim();
    is_polynomialcomp() = true;
    is_equivalent() = false;
    is_polynomial() = false;
    is_lagrange() = false;
    is_standard() = false;
    estimated_degree() = 5;
    base() = HCT->base();

    gmm::copy(gmm::identity_matrix(), P);
    init_cvs_node();

    for (size_type i = 0; i < 3; ++i) {
      base_node pt(0.0, 0.0);
      if (i) pt[i-1] = 1.0;
      add_node(lagrange_dof(2), pt);
      add_node(derivative_dof(2, 0), pt);
      add_node(derivative_dof(2, 1), pt);
    }
  }


  pfem reduced_HCT_triangle_fem
  (fem_param_list &params,
   std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 0, "Bad number of parameters : "
		<< params.size() << " should be 0.");
    pfem p = std::make_shared<reduced_HCT_triangle__>();
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }


  /* ******************************************************************** */
  /*   C1 composite element on quadrilateral (piecewise P3, FVS element). */
  /* ******************************************************************** */

  struct quadc1p3__ : public fem<bgeot::polynomial_composite> {
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
			   bgeot::pgeometric_trans pgt) const;
    mesh m;
    mutable bgeot::mesh_precomposite mp;
    mutable bgeot::pgeotrans_precomp pgp;
    mutable pfem_precomp pfp;
    mutable bgeot::pgeometric_trans pgt_stored;
    mutable base_matrix K;
    mutable bgeot::base_small_vector true_normals[4];
    quadc1p3__(void);
  };

  void quadc1p3__::mat_trans(base_matrix &M, const base_matrix &G,
				 bgeot::pgeometric_trans pgt) const {
    
    dim_type N = dim_type(G.nrows());
    
    GMM_ASSERT1(N == 2, "Sorry, this version of reduced HCT "
		"element works only on dimension two.");
    if (pgt != pgt_stored) {
      pgt_stored = pgt;
      pgp = bgeot::geotrans_precomp(pgt, node_tab(0), 0);
      pfp = fem_precomp(std::make_shared<quadc1p3__>(), node_tab(0), 0);
    }
    gmm::copy(gmm::identity_matrix(), M);
    
    gmm::mult(G, pgp->grad(0), K);
    for (size_type i = 0; i < 4; ++i) {
      if (i && !(pgt->is_linear())) gmm::mult(G, pgp->grad(3*i), K);
      gmm::copy(K, gmm::sub_matrix(M, gmm::sub_interval(1+3*i, 2)));
    }

    // take the normal derivatives into account
    static base_matrix W(4, 16);
    base_small_vector norient(M_PI, M_PI * M_PI);
    if (pgt->is_linear()) gmm::lu_inverse(K); 
    for (unsigned i = 12; i < 16; ++i) {
      if (!(pgt->is_linear()))
	{ gmm::mult(G, pgp->grad(i), K); gmm::lu_inverse(K); }
      bgeot::base_small_vector n(2), v(2);
      gmm::mult(gmm::transposed(K), cvr->normals()[i-12], n);
      n /= gmm::vect_norm2(n);

      scalar_type ps = gmm::vect_sp(n, norient);
      if (ps < 0) n *= scalar_type(-1);
      true_normals[i-12] = n;
      if (gmm::abs(ps) < 1E-8)
	GMM_WARNING2("FVS_quadrilateral : "
		     "The normal orientation may be not correct");
      gmm::mult(K, n, v);
      const bgeot::base_tensor &t = pfp->grad(i);
      for (unsigned j = 0; j < 16; ++j)
	W(i-12, j) = t(j, 0, 0) * v[0] + t(j, 0, 1) * v[1];
    }
    
    static base_matrix A(4, 4);
    static bgeot::base_vector w(4), coeff(4);
    static gmm::sub_interval SUBI(12, 4), SUBJ(0, 4);
    gmm::copy(gmm::sub_matrix(W, SUBJ, SUBI), A);
    gmm::lu_inverse(A);
    gmm::copy(gmm::transposed(A), gmm::sub_matrix(M, SUBI));

    for (unsigned j = 0; j < 12; ++j) {
      gmm::mult(W, gmm::mat_row(M, j), w);
      gmm::mult(A, gmm::scaled(w, -1.0), coeff);
      gmm::copy(coeff, gmm::sub_vector(gmm::mat_row(M, j), SUBI));
    }
  }

  quadc1p3__::quadc1p3__(void) : pgt_stored(0), K(2, 2) {

    m.clear();
    size_type i0 = m.add_point(base_node(0.0, 0.0));
    size_type i1 = m.add_point(base_node(1.0, 0.0));
    size_type i2 = m.add_point(base_node(0.0, 1.0));
    size_type i3 = m.add_point(base_node(1.0, 1.0));
    size_type i4 = m.add_point(base_node(0.5, 0.5));
    m.add_triangle(i1, i3, i4);
    m.add_triangle(i2, i0, i4);
    m.add_triangle(i3, i2, i4);
    m.add_triangle(i0, i1, i4);
    mp.initialise(m);

 std::stringstream s
   ("2 - 3*x - 3*y + 6*x*y + x^3 - 3*x^2*y;"
    "1 - 3*x^2 - 3*y^2 + x^3 + 3*x^2*y + 2*y^3;"
    "2 - 3*x - 3*y + 6*x*y - 3*x*y^2 + y^3;"
    "1 - 3*x^2 - 3*y^2 + 2*x^3 + 3*x*y^2 + y^3;"
    "1/6 + 1/2*x - y - 3/2*x^2 + 2*x*y + 5/6*x^3 - x^2*y;"
    "x - 1/2*x^2 - 3*x*y + 1/6*x^3 + x^2*y + 2*x*y^2;"
    "1/6 + 1/2*x - y - x*y + 3/2*y^2 + 1/2*x*y^2 - 2/3*y^3;"
    "x - 2*x^2 - 3/2*y^2 + x^3 + 3/2*x*y^2 + 2/3*y^3;"
    "1/6 - x + 1/2*y + 3/2*x^2 - x*y - 2/3*x^3 + 1/2*x^2*y;"
    "y - 3/2*x^2 - 2*y^2 + 2/3*x^3 + 3/2*x^2*y + y^3;"
    "1/6 - x + 1/2*y + 2*x*y - 3/2*y^2 - x*y^2 + 5/6*y^3;"
    "y - 3*x*y - 1/2*y^2 + 2*x^2*y + x*y^2 + 1/6*y^3;"
    "-1 + 3*x + 3*y - 6*x*y - 3*y^2 - x^3 + 3*x^2*y + 2*y^3;"
    "3*x^2 - x^3 - 3*x^2*y;"
    "-1 + 3*x + 3*y - 6*x*y - 3*y^2 + 3*x*y^2 + y^3;"
    "3*x^2 - 2*x^3 - 3*x*y^2 + y^3;"
    "-2/3 + 1/2*x + 2*y - x*y - 2*y^2 + 1/6*x^3 - x^2*y + 2*x*y^2;"
    "-x^2 + 5/6*x^3 + x^2*y;"
    "-2/3 + 1/2*x + 2*y - x*y - 2*y^2 + 1/2*x*y^2 + 2/3*y^3;"
    "-x^2 + x^3 + 3/2*x*y^2 - 2/3*y^3;"
    "-5/6 + x + 5/2*y + 1/2*x^2 - 3*x*y - 2*y^2 - 2/3*x^3 + 3/2*x^2*y + y^3;"
    "-1/2*x^2 + 2/3*x^3 + 1/2*x^2*y;"
    "-5/6 + x + 5/2*y - 2*x*y - 5/2*y^2 + x*y^2 + 5/6*y^3;"
    "-x*y + 1/2*y^2 + 2*x^2*y - x*y^2 + 1/6*y^3;"
    "-1 + 3*x + 3*y - 3*x^2 - 6*x*y + x^3 + 3*x^2*y;"
    "3*y^2 + x^3 - 3*x^2*y - 2*y^3;"
    "-1 + 3*x + 3*y - 3*x^2 - 6*x*y + 2*x^3 + 3*x*y^2 - y^3;"
    "3*y^2 - 3*x*y^2 - y^3;"
    "-5/6 + 5/2*x + y - 5/2*x^2 - 2*x*y + 5/6*x^3 + x^2*y;"
    "1/2*x^2 - x*y + 1/6*x^3 - x^2*y + 2*x*y^2;"
    "-5/6 + 5/2*x + y - 2*x^2 - 3*x*y + 1/2*y^2 + x^3 + 3/2*x*y^2 - 2/3*y^3;"
    "-1/2*y^2 + 1/2*x*y^2 + 2/3*y^3;"
    "-2/3 + 2*x + 1/2*y - 2*x^2 - x*y + 2/3*x^3 + 1/2*x^2*y;"
    "-y^2 - 2/3*x^3 + 3/2*x^2*y + y^3;"
    "-2/3 + 2*x + 1/2*y - 2*x^2 - x*y + 2*x^2*y - x*y^2 + 1/6*y^3;"
    "-y^2 + x*y^2 + 5/6*y^3;"
    "1 - 3*x - 3*y + 3*x^2 + 6*x*y + 3*y^2 - x^3 - 3*x^2*y - 2*y^3;"
    "-x^3 + 3*x^2*y;"
    "1 - 3*x - 3*y + 3*x^2 + 6*x*y + 3*y^2 - 2*x^3 - 3*x*y^2 - y^3;"
    "3*x*y^2 - y^3;"
    "-2/3 + 3/2*x + 2*y - x^2 - 3*x*y - 2*y^2 + 1/6*x^3 + x^2*y + 2*x*y^2;"
    "5/6*x^3 - x^2*y;"
    "-2/3 + 3/2*x + 2*y - x^2 - 3*x*y - 2*y^2 + x^3 + 3/2*x*y^2 + 2/3*y^3;"
    "1/2*x*y^2 - 2/3*y^3;"
    "-2/3 + 2*x + 3/2*y - 2*x^2 - 3*x*y - y^2 + 2/3*x^3 + 3/2*x^2*y + y^3;"
    "-2/3*x^3 + 1/2*x^2*y;"
    "-2/3 + 2*x + 3/2*y - 2*x^2 - 3*x*y - y^2 + 2*x^2*y + x*y^2 + 1/6*y^3;"
    "-x*y^2 + 5/6*y^3;"
    "4/3 - 2*x - 4*y + 4*x*y + 4*y^2 + 2/3*x^3 - 4*x*y^2;"
    "-2/3*x^3;"
    "4/3 - 2*x - 4*y + 4*x*y + 4*y^2 - 2*x*y^2 - 4/3*y^3;"
    "-2*x*y^2 + 4/3*y^3;"
    "-2/3 + 2*x - 2*x^2 + 2/3*x^3;"
    "2*x^2 - 4*x*y - 2/3*x^3 + 4*x*y^2;"
    "-2/3 + 2*x - 4*x*y + 2*y^2 + 2*x*y^2 - 4/3*y^3;"
    "-2*y^2 + 2*x*y^2 + 4/3*y^3;"
    "4/3 - 4*x - 2*y + 4*x^2 + 4*x*y - 4/3*x^3 - 2*x^2*y;"
    "4/3*x^3 - 2*x^2*y;"
    "4/3 - 4*x - 2*y + 4*x^2 + 4*x*y - 4*x^2*y + 2/3*y^3;"
    "-2/3*y^3;"
    "-2/3 + 2*y + 2*x^2 - 4*x*y - 4/3*x^3 + 2*x^2*y;"
    "-2*x^2 + 4/3*x^3 + 2*x^2*y;"
    "-2/3 + 2*y - 2*y^2 + 2/3*y^3;"
    "-4*x*y + 2*y^2 + 4*x^2*y - 2/3*y^3;");

    bgeot::pconvex_ref cr = bgeot::parallelepiped_of_reference(2);
    mref_convex() = cr;
    dim() = cr->structure()->dim();
    is_polynomialcomp() = true;
    is_equivalent() = false;
    is_polynomial() = false;
    is_lagrange() = false;
    is_standard() = false;
    estimated_degree() = 5;
    init_cvs_node();

    base()=std::vector<bgeot::polynomial_composite>
      (16, bgeot::polynomial_composite(mp, false));
    for (size_type k = 0; k < 16; ++k)
      for (size_type ic = 0; ic < 4; ++ic)
        base()[k].set_poly_of_subelt(ic, bgeot::read_base_poly(2, s));

    for (size_type i = 0; i < 4; ++i) {
      base_node pt(0.0, 0.0);
      if (i & 1) pt[0] = 1.0;
      if (i & 2) pt[1] = 1.0;
      add_node(lagrange_dof(2), pt);
      add_node(derivative_dof(2, 0), pt);
      add_node(derivative_dof(2, 1), pt);
    }

    add_node(normal_derivative_dof(2), base_node(1.0, 0.5));
    add_node(normal_derivative_dof(2), base_node(0.0, 0.5));
    add_node(normal_derivative_dof(2), base_node(0.5, 1.0));
    add_node(normal_derivative_dof(2), base_node(0.5, 0.0));
  }


  pfem quadc1p3_fem
  (fem_param_list &params,
   std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 0, "Bad number of parameters : "
		<< params.size() << " should be 0.");
    pfem p = std::make_shared<quadc1p3__>();
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  /* ******************************************************************** */
  /*    Reduced C1 composite element on quadrilateral (piecewise P3).     */
  /* ******************************************************************** */

  struct reduced_quadc1p3__ : public fem<bgeot::polynomial_composite> {
    const quadc1p3__ *HCT;
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
			   bgeot::pgeometric_trans pgt) const;
    virtual size_type nb_base(size_type) const { return 16; }
    mutable base_matrix P, Mhct;
    reduced_quadc1p3__(void);
  };

  void reduced_quadc1p3__::mat_trans(base_matrix &M, const base_matrix &G,
				 bgeot::pgeometric_trans pgt) const {
    HCT->mat_trans(Mhct, G, pgt);
    
    P(13, 1)=HCT->true_normals[1][0]*0.5; P(15, 1)=HCT->true_normals[3][0]*0.5;
    P(13, 2)=HCT->true_normals[1][1]*0.5; P(15, 2)=HCT->true_normals[3][1]*0.5;

    P(12, 4)=HCT->true_normals[0][0]*0.5; P(15, 4)=HCT->true_normals[3][0]*0.5;
    P(12, 5)=HCT->true_normals[0][1]*0.5; P(15, 5)=HCT->true_normals[3][1]*0.5;

    P(13, 7)=HCT->true_normals[1][0]*0.5; P(14, 7)=HCT->true_normals[2][0]*0.5;
    P(13, 8)=HCT->true_normals[1][1]*0.5; P(14, 8)=HCT->true_normals[2][1]*0.5;

    P(12,10)=HCT->true_normals[0][0]*0.5; P(14,10)=HCT->true_normals[2][0]*0.5;
    P(12,11)=HCT->true_normals[0][1]*0.5; P(14,11)=HCT->true_normals[2][1]*0.5;

    gmm::mult(gmm::transposed(P), Mhct, M);
  }

  reduced_quadc1p3__::reduced_quadc1p3__(void)
    : P(16, 12), Mhct(16, 16) {
    HCT = dynamic_cast<const quadc1p3__ *>
      (&(*fem_descriptor("FEM_QUADC1_COMPOSITE")));

    bgeot::pconvex_ref cr = bgeot::parallelepiped_of_reference(2);
    mref_convex() = cr;
    dim() = cr->structure()->dim();
    is_polynomialcomp() = true;
    is_equivalent() = false;
    is_polynomial() = false;
    is_lagrange() = false;
    is_standard() = false;
    estimated_degree() = 5;
    base() = HCT->base();

    gmm::copy(gmm::identity_matrix(), P);
    init_cvs_node();

    for (size_type i = 0; i < 4; ++i) {
      base_node pt(0.0, 0.0);
      if (i & 1) pt[0] = 1.0;
      if (i & 2) pt[1] = 1.0;
      add_node(lagrange_dof(2), pt);
      add_node(derivative_dof(2, 0), pt);
      add_node(derivative_dof(2, 1), pt);
    }
  }


  pfem reduced_quadc1p3_fem
  (fem_param_list &params,
   std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() == 0, "Bad number of parameters : "
		<< params.size() << " should be 0.");
    pfem p = std::make_shared<reduced_quadc1p3__>();
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }


  /* ******************************************************************** */
  /*    HHO methods: First method for the interior of the elements and    */
  /*            a method for each face (or a single method for all faces) */
  /* ******************************************************************** */
  /* It is guaranted (and used) that the sub-element of index 0 is the    */
  /* element itself and the faces follows in their usual order.           */
  /* It has also to be guaranted that the internal degrees of freedom are */
  /* first. This is ensred by the dof enumeration of mesh_fem object      */
  /* since the interior element has the index 0.                          */
  
  pfem hho_method(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {
    GMM_ASSERT1(params.size() >= 2, "Bad number of parameters : "
		<< params.size() << " should be at least 2.");
    GMM_ASSERT1(params[0].type() == 1 && params[1].type() == 1,
		"Bad type of parameters");
    pfem pf = params[0].method();

    size_type nbf = pf->ref_convex(0)->structure()->nb_faces();
    GMM_ASSERT1(pf->is_polynomial(), "Only for polynomial elements");

    std::vector<pfem> pff(nbf);
    if (params.size() == 2)
      std::fill(pff.begin(), pff.end(), params[1].method());
    else {
      GMM_ASSERT1(params.size() == nbf+1, "Bad number of parameters : "
                  << params.size() << " a single method for all the faces or "
                  " a method for each face.");
      for (size_type i = 0; i < nbf; ++i) {
        GMM_ASSERT1(params[i+1].type() == 1, "Bad type of parameters");
        GMM_ASSERT1(params[i+1].method()->is_polynomial(),
                    "Only for polynomial elements");
        pff[i] = params[i+1].method();
      }
    }

    // Obtain the reference element
    bgeot::pbasic_mesh pm;
    bgeot::pmesh_precomposite pmp;
    structured_mesh_for_convex(pf->ref_convex(0), 1, pm, pmp);

    // Addition of faces
    bgeot::basic_mesh m0(*pm);
    for  (short_type i = 0; i < nbf; ++i) {
      bgeot::mesh_structure::ind_pt_face_ct
        ipts=m0.ind_points_of_face_of_convex(0,i);
      bgeot::pconvex_structure cvs
        = m0.structure_of_convex(0)->faces_structure()[i];
      m0.add_convex(bgeot::default_trans_of_cvs(cvs), ipts.begin());
    }

    // Build the mesh_fem
    mesh m1(m0);
    bgeot::mesh_precomposite mp; mp.initialise(m1);
    mesh_fem mf(m1);
    mf.set_finite_element(0, pf);
    for  (size_type i = 0; i < nbf; ++i)
      mf.set_finite_element(i+1, pff[i]);
    
    pfem p = composite_fe_method(m1, mf, pf->ref_convex(0), true);
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  pfem interior_fem_of_hho_method(pfem hho_method) {

    const polynomial_composite_fem *phho
      = dynamic_cast<const polynomial_composite_fem*>(hho_method.get());

    if (phho) {
      pfem pf0 = phho->mf.fem_of_element(0);
      pfem pf1 = phho->mf.fem_of_element(1);
      if (pf1 && (pf1->dim()+1 == pf0->dim()))
        return phho->mf.fem_of_element(0);
    }

    // GMM_WARNING2("probably not a HHO method");
    return hho_method;
  }








  


}  /* end of namespace getfem.                                            */
