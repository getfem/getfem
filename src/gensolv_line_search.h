/* *********************************************************************** */
/*                                                                         */
/* Library :  GENeric SOLVers (gensolv) version 1.0                        */
/* File    :  gensolv_line_search.h : 1D optimization algorithms.          */
/*     									   */
/* Date : December 21, 2001.                                               */
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


#ifndef __GENSOLV_LINE_SEARCH_H
#define __GENSOLV_LINE_SEARCH_H

#include <gensolv_function.h>

namespace gensolv
{

  /* ********************************************************************* */
  /*    Algorithms without explicit gradient computation                   */
  /* ********************************************************************* */

  struct line_search_dicho
  { // a tester ...
  
    template <class FUNC> typename FUNC::base_scalar search(FUNC &f,
					     typename FUNC::base_scalar a,
					     typename FUNC::base_scalar b,
					     typename FUNC::base_scalar PREC,
					     int noisy = 0)
    {
      typename FUNC::base_scalar t[5], u[5];
      t[0] = a;    t[2] = 0.5 * (a+b); t[4] = b;
      u[0] = f(a); u[2] = f(t[2]);     u[4] = f(b);
      
      while (dal::abs(t[0] - t[4]) > PREC)
      {
	t[1] = 0.5 * (t[0] + t[2]); u[1] = f(t[1]);   
	t[3] = 0.5 * (t[4] + t[2]); u[3] = f(t[3]);
	int i, j;
	for (j = 0, i = 0; j < 5; ++j) if (u[j] < u[i]) i = j;
	switch (i)
	{
          case 0 : t[4] = t[1]; u[4] = u[1]; 
	           t[2] = 0.5 * (t[4]+t[0]); u[2] = f(t[2]); break;
          case 1 : t[4] = t[2]; u[4] = u[2]; t[2] = t[1]; u[2] = u[1]; break;
          case 2 : t[0] = t[1]; u[0] = u[1]; t[4] = t[3]; u[4] = u[3]; break;
          case 3 : t[0] = t[2]; u[0] = u[2]; t[2] = t[3]; u[2] = u[3]; break;
          case 4 : t[0] = t[3]; u[0] = u[3]; 
	           t[2] = 0.5 * (t[4]+t[0]); u[2] = f(t[2]); break;
	}
      }   
      return t[2];
    }
  };
  
  /** Line-search of dichotomie type without computation of the derivative and
   *  optimized with a gold number rule. The algorihm can extend the initial
   *  interval if needed.
   *  See Michel Minoux, Programmation mathématique, Dunod 1983
   */
  struct line_search_gold_number_with_extension
  { // a tester ...
  
    template <class FUNC> typename FUNC::base_scalar search(FUNC &f,
					     typename FUNC::base_scalar a,
					     typename FUNC::base_scalar b,
					     typename FUNC::base_scalar PREC,
					     int noisy = 0)
    {
      typename FUNC::base_scalar t[4], u[4];
      const typename FUNC::base_scalar ratio = 0.381966; // # 1-2/(1+sqrt(5))
      t[0] = a; t[1] = a + (b-a)*ratio; t[2] = b - (b-a)*ratio; t[3] = b;
      u[0] = f(t[0]); u[1] = f(t[1]); u[2] = f(t[2]); u[3] = f(t[3]);
      int i, j;
      while (dal::abs(t[0] - t[3]) > PREC)
      {
	if (noisy > 0) std::cout << t[0] << " : " << t[1] << " : " << t[2]
			    << " : " << t[3] << endl << "  " << u[0]
			    << " : " << u[1] << " : " << u[2] << " : "
			    << u[3] << endl;
	for (j = 0, i = 0; j < 4; ++j) if (u[j] < u[i]) i = j;
	if (noisy > 0) std::cout << "i = " << i << endl;
	switch (i)
	{
          case 0 : t[3] = t[1]; u[3] = u[1]; t[1] = t[0]; u[1] = u[0];
                   t[0] = (t[1] - ratio * t[3])/(1.0 - ratio); u[0] = f(t[0]);
		   t[2] = t[0] + t[3] - t[1]; u[2] = f(t[2]); break;
          case 1 : t[3] = t[2]; u[3] = u[2]; t[2] = t[1]; u[2] = u[1];
	           t[1] = t[0] + (t[2] - t[0]) * ratio / (1.0 - ratio);
		   u[1] = f(t[1]); break;
          case 2 : t[0] = t[1]; u[0] = u[1]; t[1] = t[2]; u[1] = u[2];
                   t[2] = t[3] - (t[3] - t[1]) * ratio / (1.0 - ratio);
		   u[2] = f(t[2]); break;
          case 3 : t[0] = t[2]; u[0] = u[2]; t[2] = t[3]; u[2] = u[3];
	           t[3] = (t[2] - ratio * t[0])/(1.0 - ratio); u[3] = f(t[3]);
		   t[1] = t[3] + t[0] - t[2]; u[1] = f(t[1]); break;
	}
      }
      for (j = 0, i = 0; j < 4; ++j) if (u[j] < u[i]) i = j;
      return t[i];
    }
  };



}  /* end of namespace gensolv.                                            */


#endif /* ___GENSOLV_LINE_SEARCH_H                                         */
