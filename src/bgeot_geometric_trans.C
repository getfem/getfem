/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_geometric_trans.C : geometric transformations on convex*/
/*     									   */
/*                                                                         */
/* Date : December 20, 2000.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2000-2002  Yves Renard.                                   */
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



#include <bgeot_geometric_trans.h>
#include <dal_tree_sorted.h>
#include <ftool_naming.h>

namespace bgeot
{
  typedef ftool::naming_system<geometric_trans>::param_list gt_param_list;



  /* ******************************************************************** */
  /* transformation on simplex.                                           */
  /* ******************************************************************** */

  struct _simplex_trans : public geometric_trans
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

    _simplex_trans(dim_type nc, short_type k)
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
    return new _simplex_trans(n, k);
  }

  /* ******************************************************************** */
  /* direct product transformation                                        */
  /* ******************************************************************** */

  struct _cv_pr_t : public geometric_trans
  {
    _cv_pr_t(pgeometric_trans a, pgeometric_trans b)
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
    return new _cv_pr_t(a, b);
  }

  /* ******************************************************************** */
  /* linear direct product transformation.                                */
  /* ******************************************************************** */

  struct _cv_pr_tl : public geometric_trans
  {
    _cv_pr_tl(pgeometric_trans a, pgeometric_trans b)
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
    return new _cv_pr_tl(a, b);
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
  base_vector compute_normal(const base_matrix &G, size_type ir,
			     pgeometric_trans pgt, const base_node &pt) {
    dim_type P = pgt->structure()->dim(), N = G.nrows();
    short_type NP = pgt->nb_points();
    base_matrix K(N,P), CS(P,P), B(N,P), Grad(pgt->nb_points(),P), TMP1(P,P);
    base_vector un, up;
    base_poly Poly;
    
    if (G.ncols() != NP) DAL_THROW(dimension_error, "dimensions mismatch");
    
    un.resize(P); up.resize(N);
    un = pgt->normals()[ir];
    //cout << "un=" << un << endl;
    
    for (size_type i = 0; i < pgt->nb_points(); ++i) {
      for (dim_type n = 0; n < N; ++n) {
	Poly = pgt->poly_vector()[i];
	Poly.derivative(n);
	Grad(i,n) = Poly.eval(pt.begin());
      }
    }
    
    // on peut simplifier les calculs pour N = P
    // cout << "mat G : " << G << endl;
    // cout << "mat grad : " << Grad << endl;
    mat_product(G, Grad, K);
    mat_product_tn(K, K, CS);
    // cout << "CS = " << CS << endl;
    mat_inv_cholesky(CS, TMP1);
    mat_product(K, CS, B);
    mat_vect_product(B, un, up);
    
    return up;
  }

  /* ******************************************************************** */
  /*    Naming system                                                     */
  /* ******************************************************************** */

  static ftool::naming_system<geometric_trans> *_gt_naming_system = 0;
  
  static void init_gt_naming_system(void) {
    _gt_naming_system = new ftool::naming_system<geometric_trans>("GT");
    _gt_naming_system->add_suffix("PK", PK_gt);
    _gt_naming_system->add_suffix("QK", QK_gt);
    _gt_naming_system->add_suffix("PRISM", prism_gt);
    _gt_naming_system->add_suffix("PRODUCT", product_gt);
    _gt_naming_system->add_suffix("LINEAR_PRODUCT", linear_product_gt);
  }
  
  pgeometric_trans geometric_trans_descriptor(std::string name) {
    if (_gt_naming_system == 0) init_gt_naming_system();
    size_type i = 0;
    return _gt_naming_system->method(name, i);
  }

  std::string name_of_geometric_trans(pgeometric_trans p) {
    if (_gt_naming_system == 0) init_gt_naming_system();
    return _gt_naming_system->shorter_name_of_method(p);
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

  pgeometric_trans product_geotrans(pgeometric_trans pg1, pgeometric_trans pg2) {
    static pgeometric_trans pgt = 0;
    static pgeometric_trans _pg1 = 0;
    static pgeometric_trans _pg2 = 0;
    if (pg1 != _pg1 || pg2 != _pg2) {
      std::stringstream name;
      name << "GT_PRODUCT(" << name_of_geometric_trans(pg1) << "," 
	   << name_of_geometric_trans(pg2) << ")";
      pgt = geometric_trans_descriptor(name.str());
      _pg1 = pg1; _pg2 = pg2;
    }
    return pgt;    
  }

}  /* end of namespace bgeot.                                            */

