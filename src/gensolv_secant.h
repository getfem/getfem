/* *********************************************************************** */
/*                                                                         */
/* Library :  GENeric SOLVers (gensolv) version 1.0                        */
/* File    :  gensolv_secant.h : secant algorithms.                        */
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


#ifndef __GENSOLV_SECANT_H
#define __GENSOLV_SECANT_H

#include <gensolv_function.h>

namespace gensolv
{

  /// secant algorithm with a control on the bounds.
  template <class FUNC>
    typename FUNC::base_scalar secant(const FUNC &f,
				      typename FUNC::base_scalar a,
				      typename FUNC::base_scalar b,
				      double EPS, double PREC,
				      int noisy = 0)
  {
    typedef typename FUNC::base_scalar base_scalar;
    base_scalar c, fa, fb, fc;
    base_scalar neg, pos;
    bool bneg = false, bpos = false, both = false;
    int count = 0;

    if (dal::abs(a - b) < PREC) a = b + 1.0;

    fa = f(a); fb = f(b);
    if (fa < 0) { neg = a; bneg = true; } else { pos = a; bpos = true; }
    if (fb < 0) { neg = b; bneg = true; } else { pos = b; bpos = true; }
    both = (bneg && bpos);

    // test d'arret à évaluer ...
    while ((dal::abs(a - b) > PREC) && (dal::abs(fb) > EPS))
    {
      /* compute the next position with a security.                       */
      if ( dal::abs(fb - fa) < EPS ) 
      { 
	if (both) c = 0.5 * (neg + pos);
	else if (dal::abs(a-b) < PREC) c = a + 1.0;
	else c = 0.5 * (a + b); 
      }
      else c = (a * fb - b * fa) / (fb - fa);

      /* control with the bounds.                                         */
      if (both && ((c > neg && c > pos) || (c < pos && c < neg)))
	c = 0.5 * (neg + pos);

      fc = f(c);
      /* bounds updating.                                                 */
      if (both)
      { if (fc < 0) neg = c; else pos = c; }
      else
      {
	if (fc < 0 && !bneg) { neg = c; bneg = true; }
	if (fa > 0 && !bpos) { pos = c; bpos = true; }
	both = (bneg && bpos);
      }
      a = b; fa = fb; b = c; fb = fc;

      if (noisy > 0)
      {
	cout << "iteration " << count << " a = " << a << " b = " << b
	     << " c = " << c << " fa = " << fa << " fb = " << fb << " fc = "
	     << fc << " neg = " << neg << " pos = " << pos << endl;
      }
      count++; assert(count < 1000);
    }
    return b;
  }
 

}  /* end of namespace gensolv.                                           */


#endif /* __GENSOLV_SECANT_H                                              */
