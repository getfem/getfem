/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_poly_integration.C : exact integration of polynomial   */
/*               on convexes of reference.                                 */
/*                                                                         */
/* Date : December 17, 2000.                                               */
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


#include <getfem_integration.h>
#include <ftool_naming.h>
#include <gmm.h>

#include <getfem_im_list.h>

namespace getfem
{
  typedef ftool::naming_system<integration_method>::param_list im_param_list;

  long_scalar_type poly_integration::int_poly(const base_poly &P) const
  {
    long_scalar_type res = 0.0;
    if (P.size() > int_monomials.size())
    {
      std::vector<long_scalar_type> *hum = &int_monomials;
      size_type i = P.size(), j = int_monomials.size();
      hum->resize(i);
      bgeot::power_index mi(P.dim()); mi[P.dim()-1] = P.degree();
      for (size_type k = i; k > j; --k, --mi)
	(*hum)[k-1] = int_monomial(mi);
    }
    base_poly::const_iterator it = P.begin(), ite = P.end();
    std::vector<long_scalar_type>::const_iterator itb = int_monomials.begin();
    for ( ; it != ite; ++it, ++itb) {
      res += long_scalar_type(*it) * long_scalar_type(*itb);
    }
    return res;
  }

  long_scalar_type
    poly_integration::int_poly_on_face(const base_poly &P, short_type f) const
  {
    long_scalar_type res = 0.0;
    std::vector<long_scalar_type> *hum = &(int_face_monomials[f]);
    if (P.size() > hum->size()) {
      size_type i = P.size(), j = hum->size();
      hum->resize(i);
      bgeot::power_index mi(P.dim()); mi[P.dim()-1] = P.degree();
      for (size_type k = i; k > j; --k, --mi)
	(*hum)[k-1] = int_monomial_on_face(mi, f);
    }
    base_poly::const_iterator it = P.begin(), ite = P.end();
    std::vector<long_scalar_type>::const_iterator itb = hum->begin();
    for ( ; it != ite; ++it, ++itb) res += long_scalar_type(*it) * long_scalar_type(*itb);
    return res;
  }

  /* ******************************************************************** */
  /* integration on simplex                                               */
  /* ******************************************************************** */

  struct _simplex_poly_integration : public poly_integration
  {
    long_scalar_type int_monomial(const bgeot::power_index &power) const;

    long_scalar_type int_monomial_on_face(const bgeot::power_index &power, 
				     short_type f) const;

    _simplex_poly_integration(bgeot::pconvex_structure c)
      { cvs = c;  int_face_monomials.resize(c->nb_faces()); }
  };


  long_scalar_type _simplex_poly_integration::int_monomial
  (const bgeot::power_index &power) const {
    long_scalar_type res = LONG_SCAL(1);
    short_type fa = 1;
    bgeot::power_index::const_iterator itm = power.begin(),
      itme = power.end();
    for ( ; itm != itme; ++itm)
      for (int k = 1; k <= *itm; ++k, ++fa)
	res *= long_scalar_type(k) / long_scalar_type(fa);
    
    for (int k = 0; k < cvs->dim(); k++) { res /= long_scalar_type(fa); fa++; }
    return res;
  }
  
  long_scalar_type _simplex_poly_integration::int_monomial_on_face
  (const bgeot::power_index &power, short_type f) const {
    long_scalar_type res = LONG_SCAL(0);
    
    if (f == 0 || power[f-1] == 0.0) {
      res = (f == 0) ? sqrt(long_scalar_type(cvs->dim())) : LONG_SCAL(1);
      short_type fa = 1;
      bgeot::power_index::const_iterator itm = power.begin(),
	itme = power.end();
      for ( ; itm != itme; ++itm)
	for (int k = 1; k <= *itm; ++k, ++fa)
	  res *= long_scalar_type(k) / long_scalar_type(fa);
      
      for (int k = 1; k < cvs->dim(); k++) { res /= long_scalar_type(fa); fa++; }
    }
    return res;
  }

  static pintegration_method exact_simplex(im_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n <= 0 || n >= 100 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");
    return new integration_method
      (new _simplex_poly_integration(bgeot::simplex_structure(n)));
  }

  /* ******************************************************************** */
  /* integration on direct product of convex structures                   */
  /* ******************************************************************** */

  struct _plyint_mul_structure : public poly_integration
  {
    ppoly_integration cv1, cv2;

    long_scalar_type int_monomial(const bgeot::power_index &power) const;

    long_scalar_type int_monomial_on_face(const bgeot::power_index &power, 
				     short_type f) const;

    _plyint_mul_structure(ppoly_integration a, ppoly_integration b);
  };

  long_scalar_type _plyint_mul_structure::int_monomial
  (const bgeot::power_index &power) const {
    bgeot::power_index mi1(cv1->dim()), mi2(cv2->dim());
    std::copy(power.begin(), power.begin() + cv1->dim(), mi1.begin());
    std::copy(power.begin() + cv1->dim(), power.end(), mi2.begin());
    return cv1->int_monomial(mi1) * cv2->int_monomial(mi2);
  }
  
  long_scalar_type _plyint_mul_structure::int_monomial_on_face
  (const bgeot::power_index &power, short_type f) const {
    bgeot::power_index mi1(cv1->dim()), mi2(cv2->dim());
    std::copy(power.begin(), power.begin() + cv1->dim(), mi1.begin());
    std::copy(power.begin() + cv1->dim(), power.end(), mi2.begin());
    short_type nfx = cv1->structure()->nb_faces();
    if (f < nfx)
      return cv1->int_monomial_on_face(mi1,f) * cv2->int_monomial(mi2);
    else
      return cv1->int_monomial(mi1) * cv2->int_monomial_on_face(mi2, f-nfx);
  }
  
  _plyint_mul_structure::_plyint_mul_structure(ppoly_integration a, 
					       ppoly_integration b) {
    cv1 = a; cv2 = b;
    cvs = bgeot::convex_product_structure(cv1->structure(),
					  cv2->structure());
    int_face_monomials.resize(cvs->nb_faces());
  }

  static pintegration_method product_exact(im_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	  "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 1 || params[1].type() != 1)
      DAL_THROW(failure_error, "Bad type of parameters");
    pintegration_method a = params[0].method();
    pintegration_method b = params[1].method();
    if (!(a->is_ppi && b->is_ppi))
      DAL_THROW(failure_error, "Bad parameters");
    return new integration_method(new _plyint_mul_structure(a->method.ppi,
							     b->method.ppi));
  }

  /* ******************************************************************** */
  /* integration on parallelepiped.                                       */
  /* ******************************************************************** */

  static pintegration_method exact_parallelepiped(im_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n <= 0 || n >= 100 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");

    std::stringstream name;
    if (n == 1)
      name << "IM_EXACT_SIMPLEX(1)";
    else 
      name << "IM_PRODUCT(IM_EXACT_PARALLELEPIPED(" << n-1
	   << "),IM_EXACT_SIMPLEX(1)))";
    return int_method_descriptor(name.str());
  }

  static pintegration_method exact_prism(im_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n <= 1 || n >= 100 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");

    std::stringstream name;
    name << "IM_PRODUCT(IM_EXACT_SIMPLEX(" << n-1
	 << "),IM_EXACT_SIMPLEX(1))";
    return int_method_descriptor(name.str());
  }

  /* ********************************************************************* */
  /* Approximated integration                                              */
  /* ********************************************************************* */

  void approx_integration::add_point(base_node pt,scalar_type w,short_type f) {
    if (valid) DAL_THROW(internal_error, 
			 "Impossible to modify a valid integration method.");
    if (dal::abs(w) > 1.0E-12) {
      ++f;
      if (f > cvr->structure()->nb_faces())
	DAL_THROW(internal_error, "Wrong argument.");
      size_type i = pt_to_store[f].search(pt);
      if (i == size_type(-1)) {
	i = pt_to_store[f].add(pt);
	int_coeffs.resize(int_coeffs.size() + 1); 
	for (size_type j = f; j <= cvr->structure()->nb_faces(); ++j)
	  repartition[j] += 1;
	for (size_type j = int_coeffs.size(); j > repartition[f]; --j)
	  int_coeffs[j-1] = int_coeffs[j-2];
	int_coeffs[repartition[f]-1] = 0.0;
      }
      int_coeffs[((f == 0) ? 0 : repartition[f-1]) + i] += w;
    }
  }

  void approx_integration::add_point_norepeat(base_node pt,scalar_type w,short_type f) {
    short_type f2 = f; f2++;
    if (pt_to_store[f2].search(pt) == size_type(-1)) add_point(pt,w,f);
  }

  void approx_integration::add_method_on_face(pintegration_method ppi,
					      short_type f) {
    papprox_integration pai = ppi->approx_method();
    if (valid) DAL_THROW(internal_error, 
			 "Impossible to modify a valid integration method.");
    if (pai->structure() != structure()->faces_structure()[f])
      DAL_THROW(internal_error, "structures missmatch");
    if (ppi->is_exact())
      DAL_THROW(internal_error, "Impossible with an exact method.");
    
    dim_type N = pai->structure()->dim();
    std::vector<base_node> pts(N);
    base_node pt(N+1);
    for (size_type i = 0; i < N; ++i)
      pts[i] = (cvr->dir_points_of_face(f))[i+1]
	- (cvr->dir_points_of_face(f))[0];

    base_matrix a(N+1, N), b(N, N), tmp(N, N);
    for (dim_type i = 0; i < N+1; ++i)
      for (dim_type j = 0; j < N; ++j)
	a(i, j) = pts[j][i];
    bgeot::mat_product_tn(a, a, b);
    scalar_type det = ::sqrt(dal::abs(bgeot::mat_gauss_det(b, tmp)));
    for (size_type i = 0; i < pai->nb_points_on_convex(); ++i) {
      pt = (cvr->dir_points_of_face(f))[0];
      for (dim_type j = 0; j < N; ++j)
	pt += pts[j] * (pai->integration_points()[i])[j];
      add_point(pt, pai->coeff(i) * det, f);
    }
  }
  

  void approx_integration::valid_method(void) {
    bgeot::stored_point_tab ptab(int_coeffs.size());
    size_type i = 0;
    for (short_type f = 0; f <= cvr->structure()->nb_faces(); ++f)
      for (size_type j = 0; j < pt_to_store[f].size(); ++j)
	{ ptab[i++] = pt_to_store[f][j]; }
    if (i != int_coeffs.size()) DAL_THROW(internal_error, "internal error.");
    pint_points = bgeot::store_point_tab(ptab);
    pt_to_store = std::vector<PT_TAB>();
    pt_to_store.clear();
    valid = true;
  }


  /* ********************************************************************* */
  /* methods stored in getfem_im_list.h                                    */
  /* ********************************************************************* */
  /// search a method in the getfem_im_list.

  pintegration_method im_list_integration(std::string name) {
    for (int i = 0; i < NB_IM; ++i)
      if (!strcmp(name.data(), im_desc_tab[i].method_name)) {
	bgeot::pgeometric_trans pgt
	  = bgeot::geometric_trans_descriptor(im_desc_tab[i].geotrans_name);
	dim_type N = pgt->structure()->dim();
	base_node pt(N);
	approx_integration *p = new approx_integration(pgt->convex_ref());
	size_type fr = im_desc_tab[i].firstreal;
	for (size_type j = 0; j < im_desc_tab[i].nb_points; ++j) {
	  for (dim_type k = 0; k < N; ++k)
	    pt[k] = im_desc_real[fr + j * (N+1) + k];
	  p->add_point(pt, im_desc_real[fr + j * (N+1) + N]);
	  // ajouter la gestion des symétries
	}

	for (short_type f = 0; f < pgt->structure()->nb_faces(); ++f) {
	  p->add_method_on_face
	    (int_method_descriptor
	     (im_desc_face_meth[im_desc_tab[i].firstface + f]), f);
	}

	p->valid_method();
	return new integration_method(p);
      }
    
    return 0;
  }


  /* ********************************************************************* */
  /* method de Gauss.                                                      */
  /* ********************************************************************* */

  static dal::dynamic_array<base_poly> *Legendre_polynomials;
  static dal::dynamic_array< std::vector<long_scalar_type> > *Legendres_roots;
  static void init_legendre(short_type de)
  {
    static int nb_lp = -1;

    if (nb_lp < 0)
    {
      Legendre_polynomials = new dal::dynamic_array<base_poly>();
      Legendres_roots = new dal::dynamic_array< std::vector<long_scalar_type> >();
      (*Legendre_polynomials)[0] = base_poly(1,0);
      (*Legendre_polynomials)[0].one();
      (*Legendre_polynomials)[1] = base_poly(1,1,0);
      (*Legendres_roots)[1].resize(1);
      (*Legendres_roots)[1][0] = 0.0;
      nb_lp = 1;
    }
    while (nb_lp < de)
    {
      ++nb_lp;
      (*Legendre_polynomials)[nb_lp] =
	(base_poly(1,1,0) * (*Legendre_polynomials)[nb_lp-1]
	 * ((2.0 * long_scalar_type(nb_lp) - 1.0) / long_scalar_type(nb_lp)))
	+ ((*Legendre_polynomials)[nb_lp-2]
	* ((1.0 - long_scalar_type(nb_lp)) / long_scalar_type(nb_lp)));
      (*Legendres_roots)[nb_lp].resize(nb_lp);
      (*Legendres_roots)[nb_lp][nb_lp/2] = 0.0;
      long_scalar_type a = -1.0, b, c, d, e, cv, ev, ecart, ecart2;
      for (int k = 0; k < nb_lp / 2; ++k) // + symetrie ...
      {
	b = (*Legendres_roots)[nb_lp-1][k];

	c = a, d = b;
	cv = (*Legendre_polynomials)[nb_lp].eval(&c);
	int imax = 10000;
	ecart = 2.0 * (d - c);
	while(c != d)
	{
	  --imax; if (imax == 0) break;
	  e = (c + d) / 2.0;
	  ecart2 = d - c; if (ecart2 >= ecart) break;
	  ecart = ecart2;
	  ev = (*Legendre_polynomials)[nb_lp].eval(&e);
	  if (ev * cv < 0.) { d = e; } else { c = e; cv = ev; }
	}

	(*Legendres_roots)[nb_lp][k] = c;
	(*Legendres_roots)[nb_lp][nb_lp-k-1] = -c;
	a = b;
      }
    } 
  }

  struct _gauss_approx_integration : public approx_integration
  {
    _gauss_approx_integration(short_type nbpt);
  };

  _gauss_approx_integration::_gauss_approx_integration(short_type nbpt) {
    if (nbpt > 32000) DAL_THROW(std::out_of_range, "too much points");
    
    cvr = bgeot::simplex_of_reference(1);
    bgeot::stored_point_tab int_points(nbpt+2);
    int_coeffs.resize(nbpt+2);
    repartition.resize(3);
    repartition[0] = nbpt; 
    repartition[1] = nbpt + 1;
    repartition[2] = nbpt + 2; 
    
    init_legendre(nbpt);
    
    for (short_type i = 0; i < nbpt; ++i) {
      int_points[i].resize(1);
      long_scalar_type lr = (*Legendres_roots)[nbpt][i];
      int_points[i][0] = 0.5 + 0.5 * lr;
      int_coeffs[i] = (1.0 - dal::sqr(lr))
	/ dal::sqr( long_scalar_type(nbpt) * ((*Legendre_polynomials)[nbpt-1].eval(&lr)));
    }
    
    int_points[nbpt].resize(1);
    int_points[nbpt][0] = 1.0; int_coeffs[nbpt] = 1.0;
    
    int_points[nbpt+1].resize(1);
    int_points[nbpt+1][0] = 0.0; int_coeffs[nbpt+1] = 1.0;
    pint_points = bgeot::store_point_tab(int_points);
    valid = true;
  }

  static pintegration_method gauss1d(im_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n < 0 || n >= 32000 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");
    if (n & 1) {
      std::stringstream name;
      name << "IM_GAUSS1D(" << n-1 << ")";
      return int_method_descriptor(name.str());
    }
    else
      return new integration_method(new _gauss_approx_integration(n/2 + 1));
  }

  /* ********************************************************************* */
  /* integration on simplexes                                              */
  /* ********************************************************************* */

  struct _Newton_Cotes_approx_integration : public approx_integration
  {
    // void calc_base_func(base_poly &p, short_type K, base_node &c) const;
    _Newton_Cotes_approx_integration(dim_type nc, short_type k);
  };

  _Newton_Cotes_approx_integration::_Newton_Cotes_approx_integration
  (dim_type nc, short_type k)
    : approx_integration(bgeot::simplex_of_reference(nc)) {
    size_type R = bgeot::alpha(nc,k);

    base_node c(nc); 
    if (nc == 0) {
      add_point(c, LONG_SCAL(1));
    }
    else {
      
      std::stringstream name;
      name << "IM_EXACT_SIMPLEX(" << int(nc) << ")";
      ppoly_integration ppi = int_method_descriptor(name.str())->method.ppi;
      
      size_type sum = 0, l;
      c.fill(scalar_type(0.0));
      if (k == 0) c.fill(1.0 / scalar_type(nc+1));
      
      bgeot::vsmatrix<long_scalar_type> M(R, R);
      bgeot::vsvector<long_scalar_type> F(R), U(R);
      std::vector<bgeot::power_index> base(R);
      std::vector<base_node> nodes(R);
      
      bgeot::power_index pi(nc);
      
      for (size_type r = 0; r < R; ++r, ++pi) {
	base[r] = pi; nodes[r] = c;
	if (k != 0 && nc > 0) {
	  l = 0; c[l] += 1.0 / scalar_type(k); sum++;
	  while (sum > k) {
	    sum -= int(floor(0.5+(c[l] * k)));
	    c[l] = 0.0; l++; if (l == nc) break;
	    c[l] += 1.0 / scalar_type(k); sum++;
	  }
	}
      }

//       if (nc == 1) {
// 	M = bgeot::vsmatrix<long_scalar_type>((R+1)/2, (R+1)/2);
// 	U = F = bgeot::vsvector<long_scalar_type>((R+1)/2);
// 	gmm::clear(M);
//       }

      for (size_type r = 0; r < R; ++r) {
// 	if (nc == 1) {
// 	  if (r < (R+1)/2) {
// 	    F[r] = ppi->int_monomial(base[R-1-r]);
// 	    cout << "F[" << r << "] = " << F[r] << endl; 
// 	  }
// 	}
// 	else {
	  F[r] = ppi->int_monomial(base[r]);
	  cout << "F[" << r << "] = " << F[r] << endl;
// 	}
	for (size_type q = 0; q < R; ++q) {
// 	  if (nc == 1) {
// 	    if (r < (R+1)/2) {
// 	      if (q < (R+1)/2) 
// 		M(r, q) += bgeot::eval_monomial(base[R-1-r], nodes[q].begin());
// 	      else
// 		M(r, R-1-q) += bgeot::eval_monomial(base[R-1-r], nodes[q].begin());
// 	    }
// 	  }
// 	  else
	    M(r, q) = bgeot::eval_monomial(base[r], nodes[q].begin());
	}
      }
      
      // gmm::iteration iter(1E-20, 1, 4000);
      // gmm::gmres(M, U, F, gmm::identity_matrix(), gmm::mat_nrows(M), iter);
      bgeot::mat_gauss_solve(M, F, U, LONG_SCALAR_EPS * 100);
      // bgeot::mat_gauss_solve(M, F, U, 1E-15);
//       if (nc == 1) {
// 	U.resize(R);
// 	for (size_type q = 0; q < R/2; ++q) U[R-q-1] = U[q];
//       }

      if (nc == 1)
	for (size_type r = 0; r < R; ++r)
	  cout << "node " << r << " : " << nodes[r] << " poids : " 
	       << U[r] << endl;
      
      for (size_type r = 0; r < R; ++r)
	add_point(nodes[r], U[r]);
      
      std::stringstream name2;
      name2 << "IM_NC(" << int(nc-1) << "," << int(k) << ")";
      for (short_type f = 0; f < structure()->nb_faces(); ++f)
	add_method_on_face(int_method_descriptor(name2.str()), f);
    }
    valid_method();
  }

  static pintegration_method Newton_Cotes(im_param_list &params) {
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
    return new integration_method(new _Newton_Cotes_approx_integration(n, k));
  }

  /* ********************************************************************* */
  /* integration on direct product of convex structures                    */
  /* ********************************************************************* */

  struct a_int_pro_integration : public approx_integration
  {
    a_int_pro_integration(papprox_integration a, papprox_integration b);
  };


  a_int_pro_integration::a_int_pro_integration(papprox_integration a,
					       papprox_integration b) {
    cvr = bgeot::convex_ref_product(a->ref_convex(), b->ref_convex());
    size_type n1 = a->nb_points_on_convex();
    size_type n2 = b->nb_points_on_convex();
    bgeot::stored_point_tab int_points;
    int_points.resize(n1 * n2);
    int_coeffs.resize(n1 * n2);
    repartition.resize(cvr->structure()->nb_faces()+1);
    repartition[0] = n1 * n2;
    for (size_type i1 = 0; i1 < n1; ++i1)
      for (size_type i2 = 0; i2 < n2; ++i2) {
	int_coeffs[i1 + i2 * n1] = a->coeff(i1) * b->coeff(i2);
	int_points[i1 + i2 * n1].resize(dim());
	std::copy(a->point(i1).begin(), a->point(i1).end(),
		  int_points[i1 + i2 * n1].begin());
	std::copy(b->point(i2).begin(), b->point(i2).end(),
		  int_points[i1 + i2 * n1].begin() + a->dim());
      }
    short_type f = 0;
    for (short_type f1 = 0; f1 < a->structure()->nb_faces(); ++f1, ++f) {
      n1 = a->nb_points_on_face(f1);
      n2 = b->nb_points_on_convex();
      size_type w = repartition[f];
      repartition[f+1] = w + n1 * n2;
      int_points.resize(repartition[f+1]);
      int_coeffs.resize(repartition[f+1]);
      for (size_type i1 = 0; i1 < n1; ++i1)
	for (size_type i2 = 0; i2 < n2; ++i2) {
	  int_coeffs[w + i1 + i2 * n1] = a->coeff_on_face(f1, i1)
	    * b->coeff(i2);
	  int_points[w + i1 + i2 * n1].resize(dim());
	  std::copy(a->point_on_face(f1, i1).begin(),
		    a->point_on_face(f1, i1).end(),
		    int_points[w + i1 + i2 * n1].begin());
	  std::copy(b->point(i2).begin(), b->point(i2).end(),
		    int_points[w + i1 + i2 * n1].begin() + a->dim());
	}
    }
    for (short_type f2 = 0; f2 < b->structure()->nb_faces(); ++f2, ++f) {
      n1 = a->nb_points_on_convex();
      n2 = b->nb_points_on_face(f2);
      size_type w = repartition[f];
      repartition[f+1] = w + n1 * n2;
      int_points.resize(repartition[f+1]);
      int_coeffs.resize(repartition[f+1]);
      for (size_type i1 = 0; i1 < n1; ++i1)
	for (size_type i2 = 0; i2 < n2; ++i2) {
	  int_coeffs[w + i1 + i2 * n1] = a->coeff(i1)
	    * b->coeff_on_face(f2, i2);
	  int_points[w + i1 + i2 * n1].resize(dim());
	  std::copy(a->point(i1).begin(), a->point(i1).end(),
		    int_points[w + i1 + i2 * n1].begin());
	  std::copy(b->point_on_face(f2, i2).begin(),
		    b->point_on_face(f2, i2).end(),
		    int_points[w + i1 + i2 * n1].begin() + a->dim());
	}
    }
    pint_points = bgeot::store_point_tab(int_points);
    valid = true;
  }

  static pintegration_method product_approx(im_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
       "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 1 || params[1].type() != 1)
      DAL_THROW(failure_error, "Bad type of parameters");
    pintegration_method a = params[0].method();
    pintegration_method b = params[1].method();
    if (a->is_ppi || b->is_ppi)
      DAL_THROW(failure_error, "Bad parameters");
    return new integration_method(new a_int_pro_integration(a->method.pai,
							    b->method.pai));
  }

  static pintegration_method product_which(im_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	  "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 1 || params[1].type() != 1)
      DAL_THROW(failure_error, "Bad type of parameters");
    pintegration_method a = params[0].method();
    pintegration_method b = params[1].method();
    if (a->is_ppi || b->is_ppi) return product_exact(params);
    else return product_approx(params);
  }


  /* ********************************************************************* */
  /* integration on parallelepiped with Newton Cotes formulae              */
  /* ********************************************************************* */

  static pintegration_method Newton_Cotes_para(im_param_list &params) {
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
      name << "IM_NC(1," << k << ")";
    else 
      name << "IM_PRODUCT(IM_NC_PARALLELEPIPED(" << n-1 << "," << k
	   << "),IM_NC(1," << k << "))";
    return int_method_descriptor(name.str());
  }

  static pintegration_method Newton_Cotes_prism(im_param_list &params) {
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
     name << "IM_PRODUCT(IM_NC(" << n-1 << "," << k << "),IM_NC(1,"
	 << k << "))";
    return int_method_descriptor(name.str());
  }

  /* ********************************************************************* */
  /* integration on parallelepiped with Gauss formulae                     */
  /* ********************************************************************* */

  static pintegration_method Gauss_paramul(im_param_list &params) {
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
      name << "IM_GAUSS1D(" << k << ")";
    else 
      name << "IM_PRODUCT(IM_GAUSS_PARALLELEPIPED(" << n-1 << "," << k
	   << "),IM_GAUSS1D(" << k << "))";
    return int_method_descriptor(name.str());
  }

  /* ********************************************************************* */
  /*                                                                       */
  /*   Particular integration methods on dimension 2.                      */
  /*                                                                       */
  /* ********************************************************************* */

  /* ********************************************************************* */
  /*    Integration on a triangle                                          */
  /* ********************************************************************* */

  static void triangle_add_point_full_symmetric(approx_integration *p, 
                                                long_scalar_type c1, 
                                                long_scalar_type c2,
                                                long_scalar_type w) {
    long_scalar_type c3 = 1.-c1-c2;
    p->add_point_norepeat(base_vector(c1,c2),w);
    p->add_point_norepeat(base_vector(c2,c1),w);
    p->add_point_norepeat(base_vector(c1,c3),w);
    p->add_point_norepeat(base_vector(c3,c1),w);
    p->add_point_norepeat(base_vector(c2,c3),w);
    p->add_point_norepeat(base_vector(c3,c2),w);
  }
  static void triangle_add_point_full_symmetric(approx_integration *p, 
                                                long_scalar_type c, 
                                                long_scalar_type w) {
    triangle_add_point_full_symmetric(p,c,c,w);
  }

  static pintegration_method approx_triangle(im_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n <= 0 || n >= 100 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");

    approx_integration *p
      = new approx_integration(bgeot::simplex_of_reference(2));

    switch (n) {
    case 1: 
      p->add_point(base_vector(1.0/3.0, 1.0/3.0), 0.5);
      break;
    case 2: 
      p->add_point(base_vector(1.0/6.0, 1.0/6.0), 1.0/6.0);
      p->add_point(base_vector(2.0/3.0, 1.0/6.0), 1.0/6.0);
      p->add_point(base_vector(1.0/6.0, 2.0/3.0), 1.0/6.0);
      break;
    case 3: 
      p->add_point(base_vector(1.0/3.0, 1.0/3.0), -27.0/96.0);
      p->add_point(base_vector(1.0/5.0, 1.0/5.0), 25.0/96.0);
      p->add_point(base_vector(3.0/5.0, 1.0/5.0), 25.0/96.0);
      p->add_point(base_vector(1.0/5.0, 3.0/5.0), 25.0/96.0);
      break;
    case 4: 
      {  
	double a = 0.445948490915965;
	double b = 0.091576213509771;
	double c = 0.111690794839005;
	double d = 0.054975871827661;
	p->add_point(base_vector(a, a),             c);
	p->add_point(base_vector(1.0 - 2.0 * a, a), c);
	p->add_point(base_vector(a, 1.0 - 2.0 * a), c);
	p->add_point(base_vector(b, b),             d);
	p->add_point(base_vector(1.0 - 2.0 * b, b), d);
	p->add_point(base_vector(b, 1.0 - 2.0 * b), d);
      }
      break;
    case 5:  
      {
	double a = 0.470142064105115;
	double b = 0.101286507323456;
	double c = 0.0661970763942530;
	double d = 0.0629695902724135;
	p->add_point(base_vector(a, a),             c);
	p->add_point(base_vector(1.0 - 2.0 * a, a), c);
	p->add_point(base_vector(a, 1.0 - 2.0 * a), c);
	p->add_point(base_vector(b, b),             d);
	p->add_point(base_vector(1.0 - 2.0 * b, b), d);
	p->add_point(base_vector(b, 1.0 - 2.0 * b), d);
	p->add_point(base_vector(1.0/3.0, 1.0/3.0), 9.0/80.0);
      }
      break;
    case 6: 
      {
	double a1 = 0.063089104491502;
	double a2 = 0.249286745170910;
	double aa = 0.310352451033785;
	double bb = 0.053145049844816;
	p->add_point(base_vector(a1, a1),             0.050844906370206 * 0.5);
	p->add_point(base_vector(1.0 - 2.0 * a1, a1), 0.050844906370206 * 0.5);
	p->add_point(base_vector(a1, 1.0 - 2.0 * a1), 0.050844906370206 * 0.5);
	p->add_point(base_vector(a2, a2),             0.116786275726378 * 0.5);
	p->add_point(base_vector(1.0 - 2.0 * a2, a2), 0.116786275726378 * 0.5);
	p->add_point(base_vector(a2, 1.0 - 2.0 * a2), 0.116786275726378 * 0.5);
	p->add_point(base_vector(aa, bb), 0.082851075618374 * 0.5);
	p->add_point(base_vector(aa, 1.0 - aa - bb), 0.082851075618374 * 0.5);
	p->add_point(base_vector(bb, aa), 0.082851075618374 * 0.5);
	p->add_point(base_vector(bb, 1.0 - aa - bb), 0.082851075618374 * 0.5);
	p->add_point(base_vector(1.0 - aa - bb, aa), 0.082851075618374 * 0.5);
	p->add_point(base_vector(1.0 - aa - bb, bb), 0.082851075618374 * 0.5);
      }
      break;
    case 7:
      {
	double r1  = 0.0651301029022;
	double r2  = 0.8697397941956;
	double r4  = 0.3128654960049;
	double r5  = 0.6384441885698;
	double r6  = 0.0486903154253;
	double r10 = 0.2603459660790;
	double r11 = 0.4793080678419;
	double r13 = 0.3333333333333;
	double w1  = 0.0533472356088;
	double w4  = 0.0771137608903;
	double w10 = 0.1756152574332;
	double w13 = -0.1495700444677;
	p->add_point(base_vector(r1, r1),   w1  * 0.5);
	p->add_point(base_vector(r2, r1),   w1  * 0.5);
	p->add_point(base_vector(r1, r2),   w1  * 0.5);
	p->add_point(base_vector(r4, r6),   w4  * 0.5);
	p->add_point(base_vector(r5, r4),   w4  * 0.5);
	p->add_point(base_vector(r6, r5),   w4  * 0.5);
	p->add_point(base_vector(r5, r6),   w4  * 0.5);
	p->add_point(base_vector(r4, r5),   w4  * 0.5);
	p->add_point(base_vector(r6, r4),   w4  * 0.5);
	p->add_point(base_vector(r10, r10), w10 * 0.5);
	p->add_point(base_vector(r11, r10), w10 * 0.5);
	p->add_point(base_vector(r10, r11), w10 * 0.5);
	p->add_point(base_vector(r13, r13), w13 * 0.5);
      }
      break;
    case 13:
      /*
	taken from the  Encyclopaedia of Cubature Formulas
	 http://www.cs.kuleuven.ac.be/~nines/research/ecf/ecf.html
      */
      {
        triangle_add_point_full_symmetric
	  (p, /* +3 */
	   LONG_SCAL(0.5),
	   LONG_SCAL(0.00267845189554543044455908674650066));
        triangle_add_point_full_symmetric
	  (p, /* +1 */
	   LONG_SCAL(0.333333333333333333333333333333333),
	   LONG_SCAL(0.0293480398063595158995969648597808));
        triangle_add_point_full_symmetric
	  (p, /* +3 */
	   LONG_SCAL(0.0246071886432302181878499494124643),
	   LONG_SCAL(0.00392538414805004016372590903990464));
        triangle_add_point_full_symmetric
	  (p, /* +3 */
	   LONG_SCAL(0.420308753101194683716920537182100),
	   LONG_SCAL(0.0253344765879434817105476355306468));
        triangle_add_point_full_symmetric
	  (p, /* +3 */
	   LONG_SCAL(0.227900255506160619646298948153592),
	   LONG_SCAL(0.0250401630452545330803738542916538));
        triangle_add_point_full_symmetric
	  (p, /* +3 */
	   LONG_SCAL(0.116213058883517905247155321839271),
	   LONG_SCAL(0.0158235572961491595176634480481793));
        triangle_add_point_full_symmetric
	  (p, /* +3 */
	   LONG_SCAL(0.476602980049079152951254215211496),
	   LONG_SCAL(0.0157462815379843978450278590138683));
        triangle_add_point_full_symmetric
	  (p, /* +6 */
	   LONG_SCAL(0.851775587145410469734660003794168),
	   LONG_SCAL(0.0227978945382486125477207592747430),
	   LONG_SCAL(0.00790126610763037567956187298486575));
        triangle_add_point_full_symmetric
	  (p, /* +6 */
	   LONG_SCAL(0.692797317566660854594116289398433),
	   LONG_SCAL(0.0162757709910885409437036075960413),
	   LONG_SCAL(0.00799081889046420266145965132482933));
        triangle_add_point_full_symmetric
	  (p, /* +6 */
	   LONG_SCAL(0.637955883864209538412552782122039),
	   LONG_SCAL(0.0897330604516053590796290561145196),
	   LONG_SCAL(0.0182757511120486476280967518782978));
        /* total : 37 points */
        assert(p->nb_points() == 37);
      }
      break;
    default : DAL_THROW(failure_error, "Method not implemented");
    }  
    
    std::stringstream name;
    name << "IM_GAUSS1D(" << n << ")";
    for (short_type f = 0; f < p->structure()->nb_faces(); ++f)
      p->add_method_on_face(int_method_descriptor(name.str()), f);
    p->valid_method();
    return new integration_method(p);
  }


  /* ********************************************************************* */
  /* Integration on quadrilaterals                                         */
  /* ********************************************************************* */

  static pintegration_method approx_quad(im_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n <= 0 || n >= 100 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");

    approx_integration *p
      = new approx_integration(bgeot::parallelepiped_of_reference(2));
 
    switch (n) {
    case 2: 
      p->add_point(base_vector(0.5+sqrt(1.0/6.0), 0.5),              1.0/3.0);
      p->add_point(base_vector(0.5-sqrt(1.0/24.0),0.5+sqrt(1.0/8.0)),1.0/3.0);
      p->add_point(base_vector(0.5-sqrt(1.0/24.0),0.5-sqrt(1.0/8.0)),1.0/3.0);
      break;
    case 3:
      p->add_point(base_vector(0.5 + ::sqrt(1.0/6.0), 0.5), 0.25);
      p->add_point(base_vector(0.5 - ::sqrt(1.0/6.0), 0.5), 0.25);
      p->add_point(base_vector(0.5, 0.5 + ::sqrt(1.0/6.0)), 0.25);
      p->add_point(base_vector(0.5, 0.5 - ::sqrt(1.0/6.0)), 0.25);
      break;
    case 5: 
      p->add_point(base_vector(0.5, 0.5), 2.0/7.0);
      p->add_point(base_vector(0.5, 0.5 + ::sqrt(14.0/60.0)), 5.0/63.0);
      p->add_point(base_vector(0.5, 0.5 - ::sqrt(14.0/60.0)), 5.0/63.0);
      p->add_point(base_vector(0.5 + ::sqrt(3.0/20.0),
			       0.5 + ::sqrt(3.0/20.0)), 5.0/36.0);
      p->add_point(base_vector(0.5 - ::sqrt(3.0/20.0),
			       0.5 + ::sqrt(3.0/20.0)), 5.0/36.0);
      p->add_point(base_vector(0.5 + ::sqrt(3.0/20.0),
			       0.5 - ::sqrt(3.0/20.0)), 5.0/36.0);
      p->add_point(base_vector(0.5 - ::sqrt(3.0/20.0),
			       0.5 - ::sqrt(3.0/20.0)), 5.0/36.0);
      break;
    default : DAL_THROW(failure_error, "Method not implemented");
    }
    std::stringstream name;
    name << "IM_GAUSS1D(" << n << ")";
    for (short_type f = 0; f < p->structure()->nb_faces(); ++f)
      p->add_method_on_face(int_method_descriptor(name.str()), f);
    p->valid_method();
    return new integration_method(p);
  }

  /* ********************************************************************* */
  /*                                                                       */
  /*   Particular integration methods on dimension 3.                      */
  /*                                                                       */
  /* ********************************************************************* */


  /* ********************************************************************* */
  /*  Integration on a tetrahedron                                         */
  /* ********************************************************************* */
  
  static pintegration_method approx_tetra(im_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n <= 0 || n >= 100 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");
    
    approx_integration *p = 
      new approx_integration(bgeot::simplex_of_reference(3));
    
    switch (n) {
    case 1: 
      p->add_point(base_vector(0.25, 0.25, 0.25), 1.0/6.0);
      break;
    case 2:
      {
	double c = 0.13819660112501052;       
	p->add_point(base_vector(c, c, c),           1.0/24.0);
	p->add_point(base_vector(1 - 3.0 * c, c, c), 1.0/24.0);
	p->add_point(base_vector(c, 1 - 3.0 * c, c), 1.0/24.0);
	p->add_point(base_vector(c, c, 1 - 3.0 * c), 1.0/24.0);
      }
      break;
    case 3: 
      p->add_point(base_vector(0.25, 0.25, 0.25), -4.0/30.0);
      p->add_point(base_vector(1.0/6.0, 1.0/6.0, 1.0/6.0), 9.0/120.0);
      p->add_point(base_vector(1.0/2.0, 1.0/6.0, 1.0/6.0), 9.0/120.0);
      p->add_point(base_vector(1.0/6.0, 1.0/2.0, 1.0/6.0), 9.0/120.0);
      p->add_point(base_vector(1.0/6.0, 1.0/6.0, 1.0/2.0), 9.0/120.0);
      break;
    case 5:  
      {
	double b1 = 0.31979362782962991;
	double b2 = 0.091971078052723033;
	double c1 = 0.040619116511110275;
	double c2 = 0.7240867658418309;
	double d  = 0.056350832689629156;
	double e  = 0.44364916731037084;
	double w0 = 0.019753086419753086;
	double w1 = 0.011511367871045398;
	double w2 = 0.01198951396316977;
	double w3 = 0.008818342151675485;
	p->add_point(base_vector(0.25, 0.25, 0.25), w0);
	p->add_point(base_vector(b1, b1, b1), w1);
	p->add_point(base_vector(c1, b1, b1), w1);
	p->add_point(base_vector(b1, c1, b1), w1);
	p->add_point(base_vector(b1, b1, c1), w1);
	p->add_point(base_vector(b2, b2, b2), w2);
	p->add_point(base_vector(c2, b2, b2), w2);
	p->add_point(base_vector(b2, c2, b2), w2);
	p->add_point(base_vector(b2, b2, c2), w2);
	p->add_point(base_vector(d, d, e), w3);
	p->add_point(base_vector(d, e, d), w3);
	p->add_point(base_vector(e, d, d), w3);
	p->add_point(base_vector(d, e, e), w3);
	p->add_point(base_vector(e, e, d), w3);    
	p->add_point(base_vector(e, d, e), w3);
      }
      break;
    default : DAL_THROW(failure_error, "Method not implemented");
    } 
    std::stringstream name;
    name << "IM_TRIANGLE(" << n << ")";
    for (short_type f = 0; f < p->structure()->nb_faces(); ++f)
      p->add_method_on_face(int_method_descriptor(name.str()), f);
    p->valid_method();
    return new integration_method(p);
  }

  /* ******************************************************************** */
  /*    Naming system                                                     */
  /* ******************************************************************** */
  
  pintegration_method structured_composite_int_method(im_param_list &);

  static ftool::naming_system<integration_method> *_im_naming_system = 0;
  
  static void init_im_naming_system(void) {
    _im_naming_system = new ftool::naming_system<integration_method>("IM");
    _im_naming_system->add_suffix("EXACT_SIMPLEX", exact_simplex);
    _im_naming_system->add_suffix("PRODUCT", product_which);
    _im_naming_system->add_suffix("EXACT_PARALLELEPIPED",exact_parallelepiped);
    _im_naming_system->add_suffix("EXACT_PRISM", exact_prism);
    _im_naming_system->add_suffix("GAUSS1D", gauss1d);
    _im_naming_system->add_suffix("NC", Newton_Cotes);
    _im_naming_system->add_suffix("NC_PARALLELEPIPED", Newton_Cotes_para);
    _im_naming_system->add_suffix("NC_PRISM", Newton_Cotes_prism);
    _im_naming_system->add_suffix("GAUSS_PARALLELEPIPED", Gauss_paramul);
    _im_naming_system->add_suffix("TRIANGLE", approx_triangle);
    _im_naming_system->add_suffix("QUAD", approx_quad);
    _im_naming_system->add_suffix("TETRAHEDRON", approx_tetra);
    _im_naming_system->add_suffix("STRUCTURED_COMPOSITE",
				  structured_composite_int_method);
  }
  
  pintegration_method int_method_descriptor(std::string name) {
    if (_im_naming_system == 0) init_im_naming_system();
    size_type i = 0;
    return _im_naming_system->method(name, i);
  }

  std::string name_of_int_method(pintegration_method p) {
    if (_im_naming_system == 0) init_im_naming_system();
    return _im_naming_system->shorter_name_of_method(p);
  }

  /* Fonctions pour la ref. directe.                                     */
  
  pintegration_method exact_simplex_im(size_type n) {
    static pintegration_method pim = 0;
    static size_type d = size_type(-2);
    if (d != n) {
      std::stringstream name;
      name << "IM_EXACT_SIMPLEX(" << n << ")";
      pim = int_method_descriptor(name.str());
      d = n;
    }
    return pim;
  }

  pintegration_method exact_parallelepiped_im(size_type n) {
    static pintegration_method pim = 0;
    static size_type d = size_type(-2);
    if (d != n) {
      std::stringstream name;
      name << "IM_EXACT_PARALLELEPIPED(" << n << ")";
      pim = int_method_descriptor(name.str());
      d = n;
    }
    return pim;
  }

  pintegration_method exact_prism_im(size_type n) {
    static pintegration_method pim = 0;
    static size_type d = size_type(-2);
    if (d != n) {
      std::stringstream name;
      name << "IM_EXACT_PRISM(" << n << ")";
      pim = int_method_descriptor(name.str());
      d = n;
    }
    return pim;
  }


  pintegration_method exact_classical_im(bgeot::pgeometric_trans pgt)
  {
    static bgeot::pgeometric_trans pgt_last = 0;
    static pintegration_method im_last = 0;
    bool found = false;

    if (pgt_last == pgt)
      return im_last;


    size_type n = pgt->structure()->dim();
    size_type nbp = pgt->basic_structure()->nb_points();
    std::stringstream name;

    /* Identifying P1-simplexes.                                          */

    if (nbp == n+1)
      if (pgt->basic_structure() == bgeot::simplex_structure(n))
    	{ name << "IM_EXACT_SIMPLEX("; found = true; }
    
    /* Identifying Q1-parallelepiped.                                     */

    if (!found && nbp == (size_type(1) << n))
      if (pgt->basic_structure() == bgeot::parallelepiped_structure(n))
    	{ name << "IM_EXACT_PARALLELEPIPED("; found = true; }

    /* Identifying Q1-prisms.                                             */
 
    if (!found && nbp == 2 * n)
      if (pgt->basic_structure() == bgeot::prism_structure(n))
     	{ name << "IM_EXACT_PRISM("; found = true; }
     
    // To be completed

    if (found) {
      name << int(n) << ')';
      im_last = int_method_descriptor(name.str());
      pgt_last = pgt;
      return im_last;
    }
 
    DAL_THROW(to_be_done_error,
	      "This element is not taken into account. Contact us");
  }

}  /* end of namespace getfem.                                           */

