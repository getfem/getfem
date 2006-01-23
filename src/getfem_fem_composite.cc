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


#include <getfem_poly_composite.h>
#include <getfem_integration.h>
#include <getfem_mesh_fem.h>
#include <dal_naming_system.h>

namespace getfem {
 
  typedef const fem<polynomial_composite> * ppolycompfem;

  static ppolycompfem composite_fe_method(const mesh_precomposite &mp, 
				   const mesh_fem &mf, bgeot::pconvex_ref cr) {
    
    if (&(mf.linked_mesh()) != &(mp.linked_mesh()))
      DAL_THROW(failure_error, "Meshes are different.");
    fem<polynomial_composite> *p = new fem<polynomial_composite>;

    p->mref_convex() = cr;
    p->dim() = cr->structure()->dim();
    p->is_polynomialcomp() = p->is_equivalent() = true;
    p->is_polynomial() = false;
    p->is_lagrange() = true;
    p->estimated_degree() = 0;
    p->init_cvs_node();

    std::vector<polynomial_composite> base(mf.nb_dof());
    std::fill(base.begin(), base.end(), polynomial_composite(mp));
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

    pmesh pm;
    pmesh_precomposite pmp;

    structured_mesh_for_convex(pf->ref_convex(0), k, pm, pmp);

    mesh_fem mf(*pm);
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
  /*	Composite C1 P3 element on the triangle                           */
  /* ******************************************************************** */

//   struct composite_C1_triangle__ : public fem<polynomial_composite> {
//     virtual void mat_trans(base_matrix &M, const base_matrix &G,
// 			   bgeot::pgeometric_trans pgt) const;
//     composite_C1_triangle__(void);
//   };


//   void composite_C1_triangle__::mat_trans(base_matrix &M, const base_matrix &G,
// 					  bgeot::pgeometric_trans pgt) const {

//     static bgeot::pgeotrans_precomp pgp;
//     static bgeot::pgeometric_trans pgt_stored = 0;
//     static base_matrix K(2, 2);
//     dim_type N = G.nrows();
    
//     if (N != 2) DAL_THROW(failure_error, "Sorry, this version of hermite "
// 			  "element works only on dimension two.")
//       if (pgt != pgt_stored)
// 	{ pgt_stored = pgt; pgp = bgeot::geotrans_precomp(pgt, node_tab(0)); }
//     gmm::copy(gmm::identity_matrix(), M);
    
//     gmm::mult(G, pgp->grad(0), K);
//     for (size_type i = 0; i < 6; ++i) {
//       if (i && !(pgt->is_linear())) gmm::mult(G, pgp->grad(i), K);
//       M( 6+i, 6+i) = K(0,0); M( 6+i, 12+i) = K(0,1);
//       M(12+i, 6+i) = K(1,0); M(12+i, 12+i) = K(1,1);
//     }
//   }



//   composite_C1_triangle__::composite_C1_triangle__(void) {
    
//     std::stringstream s
//       ("1 - 3*x^2 - 6*x*y - 3*y^2 + 2*x^3 + 6*x^2*y + 6*x*y^2 + 2*y^3;"
//        "0;"
//        "0;"
//        "0;"
//        "3*x^2 + 13*x*y - 2*x^3 - 13*x^2*y - 13*x*y^2;"
//        "-7*x*y + 3*y^2 + 7*x^2*y + 7*x*y^2 - 2*y^3;"
//        "0;"
//        "1 - 3*x^2 - 6*x*y - 3*y^2 + 2*x^3 + 6*x^2*y + 6*x*y^2 + 2*y^3;"
//        "0;"
//        "0;"
//        "0;"
//        "3*x^2 - 2*x^3;"
//        "13*x*y + 3*y^2 - 13*x^2*y - 13*x*y^2 - 2*y^3;"
//        "3*x^2 - 7*x*y - 2*x^3 + 7*x^2*y + 7*x*y^2;"
//        "1 - 3*x^2 - 6*x*y - 3*y^2 + 2*x^3 + 6*x^2*y + 6*x*y^2 + 2*y^3;"
//        "0;"
//        "7*x*y - 7*x^2*y - 7*x*y^2;"
//        "1 - 3*x^2 - 13*x*y - 3*y^2 + 2*x^3 + 13*x^2*y + 13*x*y^2 + 2*y^3;"
//        "3*x^2 + 6*x*y - 2*x^3 - 6*x^2*y - 6*x*y^2;"
//        "6*x*y + 3*y^2 - 6*x^2*y - 6*x*y^2 - 2*y^3;"
//        "0;"
//        "0;"
//        "3*y^2 - 2*y^3;"
//        "0;"
//        "x - 2*x^2 - 2*x*y + x^3 + 2*x^2*y + x*y^2;"
//        "0;"
//        "0;"
//        "0;"
//        "-x^2 + 2*x*y + x^3 - 2*x^2*y - 2*x*y^2;"
//        "-2*x*y + 1*y^2 + 2*x^2*y + 2*x*y^2 - 1*y^3;"
//        "0;"
//        "1*x - 2*x^2 - 4*x*y + x^3 + 4*x^2*y + 3*x*y^2;"
//        "0;"
//        "0;"
//        "0;"
//        "-1*x^2 + 1*x^3;"
//        "-1*x*y + 1*x^2*y + 2*x*y^2;"
//        "1*x*y - 2*x^2*y - 1*x*y^2;"
//        "x - 2*x^2 + 1*x^3 - 1*x*y^2;"
//        "0;"
//        "-1*x*y + 1*x^2*y + 1*x*y^2;"
//        "-y + 3*x*y + 2*y^2 - 2*x^2*y - 3*x*y^2 - 1*y^3;"
//        "-1*x^2 + 1*x^3;"
//        "-2*x*y + 2*x^2*y + 3*x*y^2;"
//        "0;"
//        "0;"
//        "x*y^2;"
//        "0;"
//        "y - 2*x*y - 2*y^2 + x^2*y + 2*x*y^2 + y^3;"
//        "0;"
//        "0;"
//        "0;"
//        "-1*x*y + 2*x^2*y + 1*x*y^2;"
//        "1*x*y - 1*x^2*y - 2*x*y^2;"
//        "0;"
//        "y - 2*y^2 - x^2*y + 1*y^3;"
//        "0;"
//        "0;"
//        "0;"
//        "x^2*y;"
//        "2*x*y - y^2 - 2*x^2*y - 2*x*y^2 + y^3;"
//        "1*x^2 - 2*x*y - 1*x^3 + 2*x^2*y + 2*x*y^2;"
//        "1*y - 4*x*y - 2*y^2 + 3*x^2*y + 4*x*y^2 + y^3;"
//        "0;"
//        "-x*y + x^2*y + x*y^2;"
//        "-x + 2*x^2 + 3*x*y - 1*x^3 - 3*x^2*y - 2*x*y^2;"
//        "-2*x*y + 3*x^2*y + 2*x*y^2;"
//        "-y^2 + y^3;"
//        "0;"
//        "0;"
//        "-y^2 + y^3;"
//        "0;"
//        "-27*x*y + 27*x^2*y + 27*x*y^2;"
//        "27*x*y - 27*x^2*y - 27*x*y^2;"
//        "0;"
//        "0;");

//     bgeot::pconvex_ref cr = bgeot::simplex_of_reference(2);
//     pmesh pm;
//     pmesh_precomposite pmp;
//     structured_mesh_for_convex(cr, 2, pm, pmp);

//     pm->write_to_file(cout);

//     mref_convex() = cr;
//     dim() = cr->structure()->dim();
//     is_polynomialcomp() = true;
//     is_equivalent() = false;
//     is_polynomial() = false;
//     is_lagrange() = false;
//     estimated_degree() = 3;
//     init_cvs_node();

//     base()=std::vector<polynomial_composite>(19,polynomial_composite(*pmp));
//     for (size_type k = 0; k < 19; ++k)
//       for (size_type ic = 0; ic < 4; ++ic) {
// 	base()[k].poly_of_subelt(ic) = bgeot::read_base_poly(2, s);
// 	cout << "poly read : " << base()[k].poly_of_subelt(ic) << endl;
//       }
//     pdof_description pdof = lagrange_dof(2);
//     for (size_type i = 0; i < 3; ++i){
//       if (i == 1) pdof = derivative_dof(1, 0);
//       if (i == 2) pdof = derivative_dof(1, 1);
//       add_node(pdof, base_node(0.0, 0.0));
//       add_node(pdof, base_node(0.5, 0.0));
//       add_node(pdof, base_node(1.0, 0.0));
//       add_node(pdof, base_node(0.0, 0.5));
//       add_node(pdof, base_node(0.5, 0.5));
//       add_node(pdof, base_node(0.0, 1.0));
//     }
//     add_node(lagrange_dof(2), base_node(1.0/3.0, 1.0/3.0));
//   }


//   pfem composite_C1_triangle_fem
//   (fem_param_list &params,
//    std::vector<dal::pstatic_stored_object> &dependencies) {
//     if (params.size() != 0)
//       DAL_THROW(failure_error, "Bad number of parameters : " << params.size()
// 		<< " should be 0.");
//     virtual_fem *p = new composite_C1_triangle__;
//     dependencies.push_back(p->ref_convex(0));
//     dependencies.push_back(p->node_tab(0));
//     return p;
//   }








  
}  /* end of namespace getfem.                                            */
