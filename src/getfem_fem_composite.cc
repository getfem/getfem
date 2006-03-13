// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_fem_composite.cc : composite fem
//           
// Date    : August 26, 2002.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2002-2006 Yves Renard
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


#include <bgeot_poly_composite.h>
#include <getfem_integration.h>
#include <getfem_mesh_fem.h>
#include <dal_naming_system.h>

namespace getfem {
 
  typedef const fem<bgeot::polynomial_composite> * ppolycompfem;

  static ppolycompfem composite_fe_method(const bgeot::mesh_precomposite &mp, 
				   const mesh_fem &mf, bgeot::pconvex_ref cr) {
    
//     if (&(bgeot::basic_mesh(mf.linked_mesh()) != &(mp.linked_mesh()))
//       DAL_THROW(failure_error, "Meshes are different.");
    fem<bgeot::polynomial_composite> *p = new fem<bgeot::polynomial_composite>;

    p->mref_convex() = cr;
    p->dim() = cr->structure()->dim();
    p->is_polynomialcomp() = p->is_equivalent() = true;
    p->is_polynomial() = false;
    p->is_lagrange() = true;
    p->estimated_degree() = 0;
    p->init_cvs_node();

    std::vector<bgeot::polynomial_composite> base(mf.nb_dof());
    std::fill(base.begin(), base.end(), bgeot::polynomial_composite(mp));
    std::vector<pdof_description> dofd(mf.nb_dof());
    
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
      pfem pf1 = mf.fem_of_element(cv);
      if (!pf1->is_lagrange()) p->is_lagrange() = false;
      if (!(pf1->is_equivalent() && pf1->is_polynomial())) {
	delete p;
	DAL_THROW(failure_error, "Only for polynomial and equivalent fem.");
      }
      ppolyfem pf = ppolyfem(pf1.get());
      p->estimated_degree() = std::max(p->estimated_degree(),
				       pf->estimated_degree());
      for (size_type k = 0; k < pf->nb_dof(cv); ++k) {
	size_type igl = mf.ind_dof_of_element(cv)[k];
	base_poly fu = pf->base()[k];
	base[igl].poly_of_subelt(cv) = fu;
	dofd[igl] = pf->dof_types()[k];
      }
    }
    p->base().resize(mf.nb_dof());
    for (size_type k = 0; k < mf.nb_dof(); ++k) {  
      p->add_node(dofd[k], mf.point_of_dof(k));
      p->base()[k] = base[k];
    }
    return p;
  }

  typedef dal::naming_system<virtual_fem>::param_list fem_param_list;

  pfem structured_composite_fem_method(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &dependencies) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	  "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 1 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    pfem pf = params[0].method();
    int k = int(::floor(params[1].num() + 0.01));
    if (!(pf->is_polynomial()) || !(pf->is_equivalent()) || k <= 0
	|| k > 150 || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");

    bgeot::pbasic_mesh pm;
    bgeot::pmesh_precomposite pmp;

    structured_mesh_for_convex(pf->ref_convex(0), k, pm, pmp);

    mesh m(*pm);
    mesh_fem mf(m);
    mf.set_finite_element(pm->convex_index(), pf);
    const virtual_fem *p = composite_fe_method(*pmp, mf, pf->ref_convex(0));
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

  pfem PK_composite_hierarch_fem(fem_param_list &params,
	std::vector<dal::pstatic_stored_object> &) {
    if (params.size() != 3)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0 || params[2].type()!= 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    int s = int(::floor(params[2].num() + 0.01)), t;
    if (n <= 0 || n >= 100 || k <= 0 || k > 150 || s <= 0 || s > 150 ||
	((s & 1) && (s != 1)) || double(s) != params[2].num() ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");
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
    if (params.size() != 3)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0 || params[2].type()!= 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    int s = int(::floor(params[2].num() + 0.01)), t;
    if (n <= 0 || n >= 100 || k <= 0 || k > 150 || s <= 0 || s > 150 ||
	((s & 1) && (s != 1)) || double(s) != params[2].num() ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");
    std::stringstream name;
    if (s == 1) 
      name << "FEM_STRUCTURED_COMPOSITE(FEM_PK_HIERARCHICAL(" << n << "," << k << "),1)";
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
  /*    Hsieh-Clough-Tocher C^1 element (composite P3)                    */
  /* ******************************************************************** */

  struct HCT_triangle__ : public fem<bgeot::polynomial_composite> {
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
			   bgeot::pgeometric_trans pgt) const;
    mesh m;
    bgeot::mesh_precomposite mp;
    HCT_triangle__(void);
  };

  void HCT_triangle__::mat_trans(base_matrix &M, const base_matrix &G,
				 bgeot::pgeometric_trans pgt) const {
    
    static bgeot::pgeotrans_precomp pgp;
    static pfem_precomp pfp;
    static bgeot::pgeometric_trans pgt_stored = 0;
    static base_matrix K(2, 2);
    dim_type N = G.nrows();
    
    if (N != 2) DAL_THROW(failure_error, "Sorry, this version of HCT "
			  "element works only on dimension two.");
    if (pgt != pgt_stored) {
      pgt_stored = pgt;
      pgp = bgeot::geotrans_precomp(pgt, node_tab(0));
      pfp = fem_precomp(this, node_tab(0));
    }
    gmm::copy(gmm::identity_matrix(), M);
    
    gmm::mult(G, pgp->grad(0), K);
    for (size_type i = 0; i < 3; ++i) {
      if (i && !(pgt->is_linear())) gmm::mult(G, pgp->grad(i), K);      
      M(3+i, 3+i) = K(0,0); M(3+i, 6+i) = K(0,1);
      M(6+i, 3+i) = K(1,0); M(6+i, 6+i) = K(1,1);
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
      if (gmm::abs(ps) < 1E-8)
	DAL_WARNING2("HCT_triangle : "
		     "The normal orientation may be not correct");
      gmm::mult(K, n, v);
      const bgeot::base_tensor &t = pfp->grad(i);
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

  HCT_triangle__::HCT_triangle__(void) {

    m.clear();
    size_type i0 = m.add_point(base_node(1.0/3.0, 1.0/3.0));
    size_type i1 = m.add_point(base_node(0.0, 0.0));
    size_type i2 = m.add_point(base_node(1.0, 0.0));
    size_type i3 = m.add_point(base_node(0.0, 1.0));
    m.add_triangle(i0, i2, i3);
    m.add_triangle(i0, i3, i1);
    m.add_triangle(i0, i1, i2);
    mp = bgeot::mesh_precomposite(m);


 std::stringstream s
   ("-1 + 9*x + 9*y - 15*x^2 - 30*x*y - 15*y^2 + 7*x^3 + 21*x^2*y + 21*x*y^2 + 7*y^3;"
    "1 - 3*x^2 - 3*y^2 + 3*x^3 - 3*x^2*y + 2*y^3;"
    "1 - 3*x^2 - 3*y^2 + 2*x^3 - 3*x*y^2 + 3*y^3;"
    "1 - 9/2*x - 9/2*y + 9*x^2 + 15*x*y + 6*y^2 - 9/2*x^3 - 21/2*x^2*y - 21/2*x*y^2 - 5/2*y^3;"
    "3*x^2 - 5/2*x^3 + 3/2*x^2*y;"
    "3*x^2 - 2*x^3 + 3/2*x*y^2 - 1/2*y^3;"
    "1 - 9/2*x - 9/2*y + 6*x^2 + 15*x*y + 9*y^2 - 5/2*x^3 - 21/2*x^2*y - 21/2*x*y^2 - 9/2*y^3;"
    "3*y^2 - 1/2*x^3 + 3/2*x^2*y - 2*y^3;"
    "3*y^2 + 3/2*x*y^2 - 5/2*y^3;"
    "-1/6 + 5/2*x - 9/2*x^2 - 4*x*y + 1/2*y^2 + 13/6*x^3 + 4*x^2*y + 3/2*x*y^2 - 1/3*y^3;"
    "x - 1/2*x^2 - 3*x*y - 7/6*x^3 + 2*x^2*y + 2*x*y^2;"
    "x - 2*x^2 - 3/2*y^2 + x^3 - 1/2*x*y^2 + 7/3*y^3;"
    "-1/6 + 3/4*x + 3/4*y - 2*x^2 - 5/2*x*y - y^2 + 17/12*x^3 + 7/4*x^2*y + 7/4*x*y^2 + 5/12*y^3;"
    "-x^2 + 13/12*x^3 - 1/4*x^2*y;"
    "-x^2 + x^3 - 1/4*x*y^2 + 1/12*y^3;"
    "2/3 - 11/4*x - 13/4*y + 7/2*x^2 + 19/2*x*y + 9/2*y^2 - 17/12*x^3 - 25/4*x^2*y - 23/4*x*y^2 - 23/12*y^3;"
    "1/2*x^2 - x*y - 13/12*x^3 + 7/4*x^2*y + 2*x*y^2;"
    "-1/2*y^2 + 9/4*x*y^2 + 5/12*y^3;"
    "-1/6 + 5/2*y + 1/2*x^2 - 4*x*y - 9/2*y^2 - 1/3*x^3 + 3/2*x^2*y + 4*x*y^2 + 13/6*y^3;"
    "y - 3/2*x^2 - 2*y^2 + 7/3*x^3 - 1/2*x^2*y + y^3;"
    "y - 3*x*y - 1/2*y^2 + 2*x^2*y + 2*x*y^2 - 7/6*y^3;"
    "2/3 - 13/4*x - 11/4*y + 9/2*x^2 + 19/2*x*y + 7/2*y^2 - 23/12*x^3 - 23/4*x^2*y - 25/4*x*y^2 - 17/12*y^3;"
    "-1/2*x^2 + 5/12*x^3 + 9/4*x^2*y;"
    "-x*y + 1/2*y^2 + 2*x^2*y + 7/4*x*y^2 - 13/12*y^3;"
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
    estimated_degree() = 3;
    init_cvs_node();

    base()=std::vector<bgeot::polynomial_composite>
      (12, bgeot::polynomial_composite(mp, false));
    for (size_type k = 0; k < 12; ++k)
      for (size_type ic = 0; ic < 3; ++ic) {
	base()[k].poly_of_subelt(ic) = bgeot::read_base_poly(2, s);
	// cout << "poly read : " << base()[k].poly_of_subelt(ic) << endl;
      }

    pdof_description pdof = lagrange_dof(2);
    for (size_type i = 0; i < 3; ++i) {
      if (i == 1) pdof = derivative_dof(2, 0);
      if (i == 2) pdof = derivative_dof(2, 1);

      add_node(pdof, base_node(0.0, 0.0));
      add_node(pdof, base_node(1.0, 0.0));
      add_node(pdof, base_node(0.0, 1.0));
    }

    add_node(normal_derivative_dof(2), base_node(0.5, 0.5));
    add_node(normal_derivative_dof(2), base_node(0.0, 0.5));
    add_node(normal_derivative_dof(2), base_node(0.5, 0.0));
  }


  pfem HCT_triangle_fem
  (fem_param_list &params,
   std::vector<dal::pstatic_stored_object> &dependencies) {
    if (params.size() != 0)
      DAL_THROW(failure_error, "Bad number of parameters : " << params.size()
		<< " should be 0.");
    virtual_fem *p = new HCT_triangle__;
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }


  /* ******************************************************************** */
  /*    Reduced Hsieh-Clough-Tocher C^1 element (composite P3)            */
  /* ******************************************************************** */

  struct reduced_HCT_triangle__ : public fem<bgeot::polynomial_composite> {
    virtual void mat_trans(base_matrix &M, const base_matrix &G,
			   bgeot::pgeometric_trans pgt) const;
    mesh m;
    bgeot::mesh_precomposite mp;
    reduced_HCT_triangle__(void);
  };

  void reduced_HCT_triangle__::mat_trans(base_matrix &M, const base_matrix &G,
				 bgeot::pgeometric_trans pgt) const {
    
    static bgeot::pgeotrans_precomp pgp;
    static bgeot::pgeometric_trans pgt_stored = 0;
    static base_matrix K(2, 2);
    dim_type N = G.nrows();
    
    if (N != 2) 
      DAL_THROW(failure_error, "Sorry, this version of reduced HCT "
		"element works only on dimension two.");
    if (pgt != pgt_stored) {
      pgt_stored = pgt;
      pgp = bgeot::geotrans_precomp(pgt, node_tab(0));
    }
    gmm::copy(gmm::identity_matrix(), M);
    
    gmm::mult(G, pgp->grad(0), K);
    for (size_type i = 0; i < 3; ++i) {
      if (i && !(pgt->is_linear())) gmm::mult(G, pgp->grad(i), K);
      M(3+i, 3+i) = K(0,0); M(3+i, 6+i) = K(0,1);
      M(6+i, 3+i) = K(1,0); M(6+i, 6+i) = K(1,1);
    }
  }

  reduced_HCT_triangle__::reduced_HCT_triangle__(void) {

    m.clear();
    size_type i0 = m.add_point(base_node(1.0/3.0, 1.0/3.0));
    size_type i1 = m.add_point(base_node(0.0, 0.0));
    size_type i2 = m.add_point(base_node(1.0, 0.0));
    size_type i3 = m.add_point(base_node(0.0, 1.0));
    m.add_triangle(i0, i2, i3);
    m.add_triangle(i0, i3, i1);
    m.add_triangle(i0, i1, i2);
    mp = bgeot::mesh_precomposite(m);

 std::stringstream s
   ("-1 + 9*x + 9*y - 15*x^2 - 30*x*y - 15*y^2 + 7*x^3 + 21*x^2*y + 21*x*y^2 + 7*y^3;"
    "1 - 3*x^2 - 3*y^2 + 3*x^3 - 3*x^2*y + 2*y^3;"
    "1 - 3*x^2 - 3*y^2 + 2*x^3 - 3*x*y^2 + 3*y^3;"
    "1 - 9/2*x - 9/2*y + 9*x^2 + 15*x*y + 6*y^2 - 9/2*x^3 - 21/2*x^2*y - 21/2*x*y^2 - 5/2*y^3;"
    "3*x^2 - 5/2*x^3 + 3/2*x^2*y;"
    "3*x^2 - 2*x^3 + 3/2*x*y^2 - 1/2*y^3;"
    "1 - 9/2*x - 9/2*y + 6*x^2 + 15*x*y + 9*y^2 - 5/2*x^3 - 21/2*x^2*y - 21/2*x*y^2 - 9/2*y^3;"
    "3*y^2 - 1/2*x^3 + 3/2*x^2*y - 2*y^3;"
    "3*y^2 + 3/2*x*y^2 - 5/2*y^3;"
    "-1/2 + 7/2*x + 2*y - 11/2*x^2 - 8*x*y - 5/2*y^2 + 5/2*x^3 + 6*x^2*y + 9/2*x*y^2 + y^3;"
    "x - 3/2*x^2 - x*y + 1/2*x^3;"
    "x - 2*x^2 - 1/2*y^2 + x^3 - 3/2*x*y^2 + y^3;"
    "-1/2 + 9/4*x + 9/4*y - 4*x^2 - 15/2*x*y - 3*y^2 + 9/4*x^3 + 21/4*x^2*y + 21/4*x*y^2 + 5/4*y^3;"
    "-x^2 + 5/4*x^3 - 3/4*x^2*y;"
    "-x^2 + x^3 - 3/4*x*y^2 + 1/4*y^3;"
    "-1/4*x + 1/4*y + 1/2*x^2 + 1/2*x*y - 1/2*y^2 - 1/4*x^3 - 3/4*x^2*y + 3/4*x*y^2 + 1/4*y^3;"
    "-1/2*x^2 + x*y + 3/4*x^3 - 3/4*x^2*y;"
    "1/2*y^2 + 3/4*x*y^2 - 3/4*y^3;"
    "-1/2 + 2*x + 7/2*y - 5/2*x^2 - 8*x*y - 11/2*y^2 + x^3 + 9/2*x^2*y + 6*x*y^2 + 5/2*y^3;"
    "y - 1/2*x^2 - 2*y^2 + x^3 - 3/2*x^2*y + y^3;"
    "y - x*y - 3/2*y^2 + 1/2*y^3;"
    "1/4*x - 1/4*y - 1/2*x^2 + 1/2*x*y + 1/2*y^2 + 1/4*x^3 + 3/4*x^2*y - 3/4*x*y^2 - 1/4*y^3;"
    "1/2*x^2 - 3/4*x^3 + 3/4*x^2*y;"
    "x*y - 1/2*y^2 - 3/4*x*y^2 + 3/4*y^3;"
    "-1/2 + 9/4*x + 9/4*y - 3*x^2 - 15/2*x*y - 4*y^2 + 5/4*x^3 + 21/4*x^2*y + 21/4*x*y^2 + 9/4*y^3;"
    "-y^2 + 1/4*x^3 - 3/4*x^2*y + y^3;"
    "-y^2 - 3/4*x*y^2 + 5/4*y^3;");

    bgeot::pconvex_ref cr = bgeot::simplex_of_reference(2);
    mref_convex() = cr;
    dim() = cr->structure()->dim();
    is_polynomialcomp() = true;
    is_equivalent() = false;
    is_polynomial() = false;
    is_lagrange() = false;
    estimated_degree() = 3;
    init_cvs_node();

    base()=std::vector<bgeot::polynomial_composite>
      (9, bgeot::polynomial_composite(mp, false));
    for (size_type k = 0; k < 9; ++k)
      for (size_type ic = 0; ic < 3; ++ic)
	base()[k].poly_of_subelt(ic) = bgeot::read_base_poly(2, s);

    pdof_description pdof = lagrange_dof(2);
    for (size_type i = 0; i < 3; ++i) {
      if (i == 1) pdof = derivative_dof(2, 0);
      if (i == 2) pdof = derivative_dof(2, 1);

      add_node(pdof, base_node(0.0, 0.0));
      add_node(pdof, base_node(1.0, 0.0));
      add_node(pdof, base_node(0.0, 1.0));
    }
  }


  pfem reduced_HCT_triangle_fem
  (fem_param_list &params,
   std::vector<dal::pstatic_stored_object> &dependencies) {
    if (params.size() != 0)
      DAL_THROW(failure_error, "Bad number of parameters : " << params.size()
		<< " should be 0.");
    virtual_fem *p = new reduced_HCT_triangle__;
    dependencies.push_back(p->ref_convex(0));
    dependencies.push_back(p->node_tab(0));
    return p;
  }

}  /* end of namespace getfem.                                            */
