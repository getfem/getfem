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
/* Copyright (C) 2000-2004  Yves Renard.                                   */
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
#include <bgeot_permutations.h>
#include <getfem_im_list.h>

namespace getfem
{

  typedef ftool::naming_system<integration_method>::param_list im_param_list;

  long_scalar_type poly_integration::int_poly(const base_poly &P) const {
    long_scalar_type res = 0.0;
    if (P.size() > int_monomials.size()) {
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
    poly_integration::int_poly_on_face(const base_poly &P,short_type f) const {
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
    for ( ; it != ite; ++it, ++itb) 
      res += long_scalar_type(*it) * long_scalar_type(*itb);
    return res;
  }

  /* ******************************************************************** */
  /* integration on simplex                                               */
  /* ******************************************************************** */

  struct simplex_poly_integration_ : public poly_integration {

    long_scalar_type int_monomial(const bgeot::power_index &power) const;

    long_scalar_type int_monomial_on_face(const bgeot::power_index &power, 
				     short_type f) const;

    simplex_poly_integration_(bgeot::pconvex_structure c)
      { cvs = c;  int_face_monomials.resize(c->nb_faces()); }
  };


  long_scalar_type
  simplex_poly_integration_::int_monomial(const bgeot::power_index &power) const{
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
  
  long_scalar_type simplex_poly_integration_::int_monomial_on_face
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
      
      for (int k = 1; k < cvs->dim(); k++) { res/=long_scalar_type(fa); fa++; }
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
      (new simplex_poly_integration_(bgeot::simplex_structure(n)));
  }

  /* ******************************************************************** */
  /* integration on direct product of convex structures                   */
  /* ******************************************************************** */

  struct plyint_mul_structure_ : public poly_integration
  {
    ppoly_integration cv1, cv2;

    long_scalar_type int_monomial(const bgeot::power_index &power) const;

    long_scalar_type int_monomial_on_face(const bgeot::power_index &power, 
				     short_type f) const;

    plyint_mul_structure_(ppoly_integration a, ppoly_integration b);
  };

  long_scalar_type plyint_mul_structure_::int_monomial
  (const bgeot::power_index &power) const {
    bgeot::power_index mi1(cv1->dim()), mi2(cv2->dim());
    std::copy(power.begin(), power.begin() + cv1->dim(), mi1.begin());
    std::copy(power.begin() + cv1->dim(), power.end(), mi2.begin());
    return cv1->int_monomial(mi1) * cv2->int_monomial(mi2);
  }
  
  long_scalar_type plyint_mul_structure_::int_monomial_on_face
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
  
  plyint_mul_structure_::plyint_mul_structure_(ppoly_integration a, 
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
    return new integration_method(new plyint_mul_structure_(a->method.ppi,
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

  void approx_integration::add_point(const base_node &pt,
				     scalar_type w,short_type f) {
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

  void approx_integration::add_point_norepeat(const base_node &pt,
					      scalar_type w,
					      short_type f) {
    short_type f2 = f; ++f2;
    if (pt_to_store[f2].search(pt) == size_type(-1)) add_point(pt,w,f);
  }

  void approx_integration::add_point_full_symmetric(base_node pt,
						    scalar_type w) {
    dim_type n = cvr->structure()->dim();
    dim_type k;
    base_node pt2(n);
    if (n+1 == cvr->structure()->nb_faces()) {
      // simplices
      base_node pt3(n+1);
      std::copy(pt.begin(), pt.end(), pt3.begin());
      pt3[n] = 1;  for (k = 0; k < n; ++k) pt3[n] -= pt[k];
      std::vector<int> ind(n); std::fill(ind.begin(), ind.end(), 0);
      std::vector<bool> ind2(n+1); 
      for (;;) {
	
	std::fill(ind2.begin(), ind2.end(), false);
	bool good = true;
	for (k = 0; k < n; ++k) 
	  if (ind2[ind[k]]) { good = false; break; } else ind2[ind[k]] = true;
	if (good) {
	  for (k = 0; k < n; ++k) pt2[k] = pt3[ind[k]];
	  add_point_norepeat(pt2, w);
	}
	ind[0]++; k = 0;
	while(ind[k] == n+1) { ind[k++] = 0; if (k == n) return; ind[k]++; }
      }
    }
    else if (cvr->structure()->nb_points() == (size_type(1) << n)) {
      // parallelepidedic
      for (size_type i = 0; i < (size_type(1) << n); ++i) {
	for (k = 0; k < n; ++k)
	  if (i & (size_type(1) << k)) pt2[k]=pt[k]; else pt2[k] = 1.0-pt[k];
	add_point_norepeat(pt2, w);
      }
    }
    else
      DAL_THROW(failure_error, "Fully symmetric option is only valid for"
		"simplices and parallelepipedic elements");

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
    scalar_type det = 1.0;
    base_node pt(N+1);
    std::vector<base_node> pts(N);  
      for (size_type i = 0; i < N; ++i)
	pts[i] = (cvr->dir_points_of_face(f))[i+1]
	  - (cvr->dir_points_of_face(f))[0];
    if (N) {
      base_matrix a(N+1, N), b(N, N), tmp(N, N);
      for (dim_type i = 0; i < N+1; ++i)
	for (dim_type j = 0; j < N; ++j)
	  a(i, j) = pts[j][i];
      
      gmm::mult(gmm::transposed(a), a, b);
      det = ::sqrt(dal::abs(gmm::lu_det(b)));
    }
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
  
  /// search a method in getfem_im_list.h
  static pintegration_method im_list_integration(std::string name) {
    // cerr << "searching " << name << endl;
    for (int i = 0; i < NB_IM; ++i)
      if (!name.compare(im_desc_tab[i].method_name)) {
	bgeot::pgeometric_trans pgt
	  = bgeot::geometric_trans_descriptor(im_desc_tab[i].geotrans_name);
	dim_type N = pgt->structure()->dim();
	base_node pt(N);
	approx_integration *pai = new approx_integration(pgt->convex_ref());
	size_type fr = im_desc_tab[i].firstreal;
	for (size_type j = 0; j < im_desc_tab[i].nb_points; ++j) {
	  for (dim_type k = 0; k < N; ++k)
	    pt[k] = atof(im_desc_real[fr+j*(N+1)+k]);
	  // pt[k] = LONG_SCALAR_ATOF(im_desc_real[fr+j*(N+1)+k]);
	  
	  switch (im_desc_node_type[im_desc_tab[i].firsttype + j]) {
	  case 2: {
            base_node pt2(pt.size());
            for (bgeot::permutation p(pt.size()); !p.finished(); ++p) {
              p.apply_to(pt,pt2);
              pai->add_point_full_symmetric(pt2, 
					    atof(im_desc_real[fr+j*(N+1)+N]));
	      // pai->add_point_full_symmetric(pt2,
	      //	        LONG_SCALAR_ATOF(im_desc_real[fr+j*(N+1)+N]));
            }
	  } break;
	  case 1: {
	    pai->add_point_full_symmetric(pt,atof(im_desc_real[fr+j*(N+1)+N]));
	    // pai->add_point_full_symmetric(pt,
	    //                   LONG_SCALAR_ATOF(im_desc_real[fr+j*(N+1)+N]));
	  } break;
	  case 0: {
	    pai->add_point(pt, atof(im_desc_real[fr+j*(N+1)+N]));
	    // pai->add_point(pt,LONG_SCALAR_ATOF(im_desc_real[fr+j*(N+1)+N]));
	  } break;
	  default: DAL_INTERNAL_ERROR("");
	  }
	}

	for (short_type f = 0; N > 0 && f < pgt->structure()->nb_faces(); ++f)
	  pai->add_method_on_face
	    (int_method_descriptor
	     (im_desc_face_meth[im_desc_tab[i].firstface + f]), f);

	pai->valid_method();
        // cerr << "finding " << name << endl;
	return new integration_method(pai);
      }
    return 0;
  }


  /* ********************************************************************* */
  /* Gauss method.                                                         */
  /* ********************************************************************* */

  static dal::dynamic_array<base_poly> *Legendre_polynomials;
  static dal::dynamic_array< std::vector<long_scalar_type> > *Legendres_roots;
  static void init_legendre(short_type de) {
    static int nb_lp = -1;

    if (nb_lp < 0) {
      Legendre_polynomials = new dal::dynamic_array<base_poly>();
      Legendres_roots =
	new dal::dynamic_array< std::vector<long_scalar_type> >();
      (*Legendre_polynomials)[0] = base_poly(1,0);
      (*Legendre_polynomials)[0].one();
      (*Legendre_polynomials)[1] = base_poly(1,1,0);
      (*Legendres_roots)[1].resize(1);
      (*Legendres_roots)[1][0] = 0.0;
      nb_lp = 1;
    }
    while (nb_lp < de) {
      ++nb_lp;
      (*Legendre_polynomials)[nb_lp] =
	(base_poly(1,1,0) * (*Legendre_polynomials)[nb_lp-1]
	 * ((2.0 * long_scalar_type(nb_lp) - 1.0) / long_scalar_type(nb_lp)))
	+ ((*Legendre_polynomials)[nb_lp-2]
	   * ((1.0 - long_scalar_type(nb_lp)) / long_scalar_type(nb_lp)));
      (*Legendres_roots)[nb_lp].resize(nb_lp);
      (*Legendres_roots)[nb_lp][nb_lp/2] = 0.0;
      long_scalar_type a = -1.0, b, c, d, e, cv, ev, ecart, ecart2;
      for (int k = 0; k < nb_lp / 2; ++k) { // + symetrie ...
	b = (*Legendres_roots)[nb_lp-1][k];

	c = a, d = b;
	cv = (*Legendre_polynomials)[nb_lp].eval(&c);
	int imax = 10000;
	ecart = 2.0 * (d - c);
	while(c != d) {
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

  struct gauss_approx_integration_ : public approx_integration {
    gauss_approx_integration_(short_type nbpt);
  };

  gauss_approx_integration_::gauss_approx_integration_(short_type nbpt) {
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
	/ dal::sqr( long_scalar_type(nbpt)
		    * ((*Legendre_polynomials)[nbpt-1].eval(&lr)));
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
      return new integration_method(new gauss_approx_integration_(n/2 + 1));
  }

  /* ********************************************************************* */
  /* integration on simplexes                                              */
  /* ********************************************************************* */

  struct Newton_Cotes_approx_integration_ : public approx_integration
  {
    // void calc_base_func(base_poly &p, short_type K, base_node &c) const;
    Newton_Cotes_approx_integration_(dim_type nc, short_type k);
  };

  Newton_Cotes_approx_integration_::Newton_Cotes_approx_integration_
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
      
      gmm::dense_matrix<long_scalar_type> M(R, R);
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
	  //cout << "F[" << r << "] = " << F[r] << endl;
// 	}
	for (size_type q = 0; q < R; ++q) {
// 	  if (nc == 1) {
// 	    if (r < (R+1)/2) {
// 	      if (q < (R+1)/2) 
// 		M(r, q) += bgeot::eval_monomial(base[R-1-r], nodes[q].begin());
// 	      else
// 		M(r, R-1-q) += bgeot::eval_monomial(base[R-1-r],
//	  nodes[q].begin());
// 	    }
// 	  }
// 	  else
	    M(r, q) = bgeot::eval_monomial(base[r], nodes[q].begin());
	}
      }
      
      gmm::lu_solve(M, U, F);
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
    return new integration_method(new Newton_Cotes_approx_integration_(n, k));
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

  /* ******************************************************************** */
  /*    Naming system                                                     */
  /* ******************************************************************** */
  
  pintegration_method structured_composite_int_method(im_param_list &);

  static ftool::naming_system<integration_method> *im_naming_system_ = 0;
  
  static void init_im_naming_system(void) {
    im_naming_system_ = new ftool::naming_system<integration_method>("IM");
    im_naming_system_->add_suffix("EXACT_SIMPLEX", exact_simplex);
    im_naming_system_->add_suffix("PRODUCT", product_which);
    im_naming_system_->add_suffix("EXACT_PARALLELEPIPED",exact_parallelepiped);
    im_naming_system_->add_suffix("EXACT_PRISM", exact_prism);
    im_naming_system_->add_suffix("GAUSS1D", gauss1d);
    im_naming_system_->add_suffix("NC", Newton_Cotes);
    im_naming_system_->add_suffix("NC_PARALLELEPIPED", Newton_Cotes_para);
    im_naming_system_->add_suffix("NC_PRISM", Newton_Cotes_prism);
    im_naming_system_->add_suffix("GAUSS_PARALLELEPIPED", Gauss_paramul);
    // im_naming_system_->add_suffix("TRIANGLE", approx_triangle);
    // im_naming_system_->add_suffix("QUAD", approx_quad);
    // im_naming_system_->add_suffix("TETRAHEDRON", approx_tetra);
    im_naming_system_->add_suffix("STRUCTURED_COMPOSITE",
				  structured_composite_int_method);
    im_naming_system_->add_generic_function(im_list_integration);
  }
  
  pintegration_method int_method_descriptor(std::string name) {
    if (im_naming_system_ == 0) init_im_naming_system();
    size_type i = 0;
    return im_naming_system_->method(name, i);
  }

  std::string name_of_int_method(pintegration_method p) {
    if (im_naming_system_ == 0) init_im_naming_system();
    return im_naming_system_->shorter_name_of_method(p);
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

