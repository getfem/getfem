/* *********************************************************************** */
/*                                                                         */
/* Library :  GENeric SOLVers (gensolv) version 1.0                        */
/* File    :  gensolv_gauss_seidel.h : Gauss Seidel algorithm for non      */
/*            linear systems.                                              */
/*     									   */
/* Date : December 18, 2001.                                               */
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


#ifndef __GENSOLV_GAUSS_SEIDEL_H
#define __GENSOLV_GAUSS_SEIDEL_H

#include <gensolv_secant.h>

namespace gensolv
{

  /** Basic Gauss-Seidel Algorithm for non-linear systems.
   *  (with a one variable line-search)
   *  correct only for the gradient of a function
   */
  template <class FUNC>
    void gauss_seidel(const FUNC &f, typename FUNC::base_vector &x, 
		      double EPS, double PREC, int noisy = 0)
  {
    typedef typename FUNC::size_type size_type;
    size_type n = f.size(), i;
    norm_type norm;
    size_type count = 0;

    do
    {
      for (i = 0, norm = 0; i < n; ++i)
      {
	norm += dal::sqr(f(x, i));
	secant(component_function<FUNC>(f, i, i, x), x[i], x[i] + 1.0,
	       EPS, PREC, noisy - 1);

      }
      norm /= n; norm = ::sqrt(norm);
      if (noisy > 0)
      { cout << "Iteration " << count++ << " residu " << norm << endl; }
    } while (norm > PREC);
  }
  
 

}  /* end of namespace gensolv.                                           */


#endif /* __GENSOLV_GAUSS_SEIDEL_H                                        */
