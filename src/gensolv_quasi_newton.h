/* *********************************************************************** */
/*                                                                         */
/* Library :  GENeric SOLVers (gensolv) version 1.0                        */
/* File    :  gensolv_quasi_newton.h : Quasi Newton algorithms.            */
/*     									   */
/* Date : February 4, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; version 2 of the License.                 */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software Foundation, */
/* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.         */
/*                                                                         */
/* *********************************************************************** */


#ifndef __GENSOLV_QUASI_NEWTON_H
#define __GENSOLV_QUASI_NEWTON_H

#include <gensolv_line_search.h>

namespace gensolv
{

  /* ********************************************************************* */
  /*                                                                       */
  /* Stockage des approximation avec mise à jour de rang un.               */
  /*                                                                       */
  /* ********************************************************************* */

  /** Represents a matrix by a sequence off rank one actualizations.      **/
  template <class VECT> class approx_mat_rank_one
  {
    public :
      typedef typename VECT::size_type  size_type;
      typedef typename VECT::value_type value_type;  

    protected :

      dal::dynamic_array<VECT> s, z;
      dal::dynamic_array<value_type> a;
      size_type nb;
      size_type limited;

    public :

      void clear(void) { s.clear(); z.clear(); a.clear(); nb = 0; }

      approx_mat_rank_one(void) { nb = 0; limited = 0; }

      void mult(const VECT &v, VECT &w)
      {
	copy(v, w);
	for (size_type k = 0; k < nb; ++k)
	  add_mult(w, dot(w, s[k]) / a[k], z[k]);
      }

      void add(const VECT &as, const VECT &az, value_type aa)
      { 
	s[nb] = VECT(as.size()); z[nb] = VECT(az.size());
	copy(as, s[nb]); copy(az, z[nb]); a[nb++] = aa;
      }

  };


  /* ********************************************************************* */
  /*                                                                       */
  /* Actualisations methods.                                               */
  /*                                                                       */
  /* ********************************************************************* */

  template <class VECT> class broyden_inverse_rank_one
  {
    protected :
      approx_mat_rank_one<VECT> H;

    public :
      typedef typename VECT::size_type  size_type;
      typedef typename VECT::value_type value_type;  


      void actualization(const VECT &s, const VECT &y)
      { /* s = x2 - x1; y = f(x2) - f(x1)                                  */
	VECT z(y.size()); H.mult(y, z);
	value_type a = dot(s, z);
        gensolv::mult(z, value_type(-1)); add(z, s);
	H.add(s, z, a);
      }

      template <class VECT2> inline void mult(const VECT2 &v, VECT2 &w)
      { H.mult(v, w); }

      template<class FUNC> void init(FUNC &f, const VECT &x)
      { 
	H.clear();
      }

  };



  // faire le stockage des inverses de gradient ... stocker les actualisation
  // pour pouvoir faire rapidement des produits matrices vecteurs.


  /* ********************************************************************* */
  /*                                                                       */
  /* Algorithme de quasi-Newton générique.                                 */
  /*                                                                       */
  /* ********************************************************************* */


  /** Generic Quasi-Newton algorithm with parametrizable line search
   *  and actualisation.
   */
  template <class FUNC, class VECT, class LINE_SEARCH, class ACTU>
    void quasi_newton(const FUNC &F, VECT &x1, LINE_SEARCH linesearch, 
		      ACTU actu, typename VECT::value_type RES_LS,
		      typename VECT::value_type RESIDU,
		      typename VECT::value_type EPS, int noisy = 0)
  {
    typedef typename VECT::value_type value_type;
    typedef typename VECT::size_type  size_type;
    // - Estimation de H0, pour le moment H = Id,
    // il faudrait estimer H0 au moins sur la diagonale ?
    actu.init(F, x1);
    VECT y1(x1.size()), y2(x1.size()), x2(x1.size()), y3(x1.size());
    VECT x3(x1.size()), p(x1.size());
    F.value(x1, y1);
    value_type residu = ::sqrt(dot(y1, y1) / y1.size()), a, aa; 
    if (noisy > 0) cout << "residu initial " << residu << endl;
    size_type iter, reinit = 100, adding = 250, stay = 0;

    for (iter = 0; residu > RESIDU; ++iter)
    {
      // - calcul de la direction de descente. 
      actu.mult(y1, p);  // p = -H F(x1) = B^{-1} F(x1)
      // copy(y1, p);
      mult(p, value_type(-1));
      // - Appel du line search avec estimation du résidu et test d'arrêt.
      norm2_dir_function<FUNC> norm2_F(F, x1, p);
      a = linesearch.search(norm2_F, 0.0, 1.0, RES_LS, noisy-1);
      aa = (dal::abs(a) < EPS) ? 1.0 : a;
      copy(x1, x2); add_mult(x2, aa, p); // x2 = x1 + a * p
      F.value(x2, y2); // y2 = F(x2) calcul en trop,on peut recup de norm2_F
      // Actualisation 
      copy(x2, x3); add_mult(x3, value_type(-1), x1);
      copy(y2, y3); add_mult(y3, value_type(-1), y1);
      actu.actualization(x3, y3);            // Actualization of H 
      if (a == aa) { copy(x2, x1); copy(y2, y1); if (aa > 0) stay = 0; };
      residu = ::sqrt(dot(y1, y1) / y1.size());
      if (noisy > 0) cout << "iter " << iter << " residu " << residu
			  << " coefficient a = " << a << endl;
      if (iter == reinit || stay == 10)
	{ actu.init(F, x1); reinit=iter + adding; adding += 50; stay=0; cout << "reinit \n"; }
      
    }
  }
}  /* end of namespace gensolv.                                           */


#endif /* __GENSOLV_FUNCTION_H                                            */
