// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_geotrans_inv.cc : Allow to inverse geometric transf-
//           ormations and to localize a set of points.
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



#include <bgeot_geotrans_inv.h>
#include <gmm_solver_bfgs.h>
namespace bgeot
{ 
  /* inversion for linear geometric transformations */
  bool geotrans_inv_convex::invert_lin(const base_node& n, base_node& n_ref, scalar_type IN_EPS) {
    base_node y(n); for (size_type i=0; i < N; ++i) y[i] -= G(i,0);
    gmm::mult(gmm::transposed(B), y, n_ref);
    if (pgt->convex_ref()->is_in(n_ref) < IN_EPS) {
      if (P == N) return true;
      else {
	gmm::mult(K,gmm::scaled(n_ref,-1.0),y,y);
	//        y -= K * n_ref;
        if (vect_norm2(y) < IN_EPS) return true;
      }
    }
    return false;
  }
  
  void geotrans_inv_convex::update_B() {
    if (P != N) {
      gmm::mult(G,pc,K);
      gmm::mult(gmm::transposed(K), K, CS);
      gmm::lu_inverse(CS);
      gmm::mult(K, CS, B);
    }
    else {
      // L'inversion peut être optimisée par le non calcul global de B
      // et la resolution d'un système linéaire.
      gmm::mult(gmm::transposed(pc), gmm::transposed(G), K);
      gmm::copy(K,B);
      gmm::lu_inverse(K); B.swap(K); 
    }
  }

  class geotrans_inv_convex_bfgs {
    geotrans_inv_convex &gic;
    base_node xreal;
  public:
    geotrans_inv_convex_bfgs(geotrans_inv_convex &gic_, 
			     const base_node &xr) : gic(gic_), xreal(xr) {}
    scalar_type operator()(const base_node& x) const {
      base_node r = gic.pgt->transform(x, gic.cvpts) - xreal;
      return gmm::vect_norm2_sqr(r)/2.;
    }
    void operator()(const base_node& x, base_small_vector& gr) const {
      gic.pgt->gradient(x, gic.pc);
      gic.update_B();
      base_node r = gic.pgt->transform(x, gic.cvpts) - xreal;
      gr.resize(x.size());
      gmm::mult(gmm::transposed(gic.K), r, gr); 
    }
  };

  /* inversion for non-linear geometric transformations 
     (Newton on Grad(pgt)(y - pgt(x)) = 0 )
  */
  bool geotrans_inv_convex::invert_nonlin(const base_node& xreal, base_node& x, scalar_type IN_EPS) {
    base_node xn(P), y, z,x0;
    /* find an initial guess */
    x0 = pgt->geometric_nodes()[0]; y = cvpts[0];  
    scalar_type d = vect_dist2_sqr(y, xreal);
    for (size_type j = 1; j < pgt->nb_points(); ++j) { 
      scalar_type d2 = vect_dist2_sqr(cvpts[j], xreal);
      if (d2 < d)
        { d = d2; x0 = pgt->geometric_nodes()[j]; y = cvpts[j]; }
    }
    x = x0;
    base_node vres(N);
    base_node rn(xreal); rn -= y; 

    pgt->gradient(x, pc);
    update_B();
    gmm::mult(gmm::transposed(K), rn, vres);
    scalar_type res = vect_norm2(vres);

    //cerr << "DEBUT: res0=" << res << ", X=" << xreal << "\nB=" << B << ", K=" << K << "\n" << ", pc=" << pc << "\n";
    unsigned cnt = 50;
    while (res > EPS/10 && --cnt) {
      gmm::mult(gmm::transposed(B), rn, xn);
      scalar_type newres;
      for (unsigned i=1; i<=256; i*=2) {
	z = x + xn / scalar_type(i);
	y = pgt->transform(z, cvpts);
	/*y.fill(0.0);
	for (size_type k = 0; k < pgt->nb_points(); ++k) {
	  gmm::add(gmm::scaled(cvpts[k],
			       scalar_type(pgt->poly_vector()[k].eval(z.begin()))),y);
	}
	*/
	// cout << "Point : " << x << " : " << y << " ptab : " << ptab[i] << endl; getchar();
	
	rn = xreal - y; 
	
	pgt->gradient(z, pc);
	update_B();
	
	if (P != N) {
	  gmm::mult(gmm::transposed(K), rn, vres);
	  newres = vect_norm2(vres); 
	} else {
	  newres = vect_norm2(rn); // "better" residu
	}
	if (newres < 1.5*res) break;
      }
      x = z; res = newres;
      //cout << "cnt=" << cnt << ", x=" << x << ", res=" << res << "\n";
    }
    //cout << " invert_nonlin done\n";
    //cerr << "cnt=" << cnt << ", P=" << P << ", N=" << N << ", G=" << G << "\nX=" << xreal << " Xref=" << x << "\nresidu=" << res << "\nB=" << B << ", K=" << K << "\n" << ", pc=" << pc << "\n-------------------^^^^^^^^\n";
    if (cnt == 0) {
      //cout << "BFGS in geotrans_inv_convex!\n";
      geotrans_inv_convex_bfgs b(*this, xreal);
      gmm::iteration iter(EPS,0);
      x = x0;
      gmm::bfgs(b,b,x,10,iter);
      rn = pgt->transform(x,cvpts) - xreal; 
      
      if (pgt->convex_ref()->is_in(x) < IN_EPS &&
	  N==P && vect_norm2(rn) > IN_EPS)
	DAL_THROW(dal::failure_error, 
		  "inversion of non-linear geometric transformation "
		  "failed ! (too much iterations)");
    }
    // Test un peu sevère peut-être en ce qui concerne rn.
    if (pgt->convex_ref()->is_in(x) < IN_EPS
        && (P == N || vect_norm2(rn) < IN_EPS)) {
      //cout << "point " << x << "in IN (" << pgt->convex_ref()->is_in(x) << ")\n";
      return true;
    } //else cout << "point IS OUT\n";
    return false;
  }

}  /* end of namespace bgeot.                                             */
