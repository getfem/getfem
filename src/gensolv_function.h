/* *********************************************************************** */
/*                                                                         */
/* Library :  GENeric SOLVers (gensolv) version 1.0                        */
/* File    :  gensolv_function.h : definitions of functions.               */
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


#ifndef __GENSOLV_FUNCTION_H
#define __GENSOLV_FUNCTION_H

#include <dal_std.h>
#include <functional>

namespace gensolv
{

  /* ********************************************************************* */
  /*                                                                       */
  /* function definitions                                                  */
  /*                                                                       */
  /* ********************************************************************* */



  /// represents a vector function
  template <class BVECT>
    struct vect_function : public std::unary_function<BVECT, BVECT>
  {
    public :
      typedef typename BVECT::size_type size_type;
      typedef typename BVECT::value_type base_scalar;
      typedef BVECT base_vector;
  };

  /// represents a scalar function
  template <class T>
    struct scal_function : public std::unary_function<T, T>
  {
    public :
      typedef T base_scalar;
  };


  /// represents a component of a vector function as a scalar function
  template <class FUNC>
    struct component_function // non ... 
      : public scal_function<typename FUNC::base_scalar>
  {
    typedef typename FUNC::size_type size_type;
    typedef typename FUNC::base_vector base_vector;
    typedef typename FUNC::base_scalar base_scalar;

    protected :
      const FUNC *f;
      base_vector *v;
      size_type i, j;

    public :
      base_scalar operator()(base_scalar x) const
      { (*v)[j] = x; return (*f)(*v, i); }
      component_function(const FUNC &func, size_type ii,
			 size_type jj, base_vector &vv)
      { f = &func; i = ii; j = jj; v = &vv; }
  };


  /** represents the norm^2 d'une fonction vectorielle dans une direction
   * comme une fonction d'une variable
   */
  template <class FUNC>
    struct norm2_dir_function
      : public scal_function<typename FUNC::base_scalar>
  {
    typedef typename FUNC::size_type size_type;
    typedef typename FUNC::base_vector base_vector;
    typedef typename FUNC::base_scalar base_scalar;

    protected :
      const FUNC *f;
      const base_vector *v, *x;
      base_vector z, y;

    public :
      base_scalar operator()(base_scalar a)
      { 
	copy(*x, z); add_mult(z, a, *v); (*f).value(z, y); 
	// print_vector(*x); print_vector(*v); print_vector(z); print_vector(y);
	// cout << "a = " << a << endl;
	return dot(y, y);
      }
      norm2_dir_function(const FUNC &func, const base_vector &px,
			                   const base_vector &pv)
      { f = &func; v = &pv; x = &px; z = base_vector(x->size()); y = base_vector(x->size()); }
  };


  /// interface for a matrix
  template <class MAT, class BVECT>
    class mat_function : public vect_function<BVECT>
  {
    protected :
      
      MAT *m;

    public :

      typedef typename vect_function<BVECT>::size_type size_type;
      typedef typename vect_function<BVECT>::base_vector base_vector;
      typedef typename vect_function<BVECT>::base_scalar base_scalar;
      

      // operation à revoir
      base_vector operator()(const base_vector &x) const
        { return (*m) * x; }
      base_scalar operator()(const base_vector &x, size_type i) const
	{ return dot(m->line(i), x); }
      void value(const base_vector &x, base_vector &y) const
	{ y = (*m) * x; }
      size_type size(void) const { return m.nrows(); }
      mat_function(const MAT &mat)
	{ m = &mat; assert(mat.nrows() == mat.ncols()); }
  };
 

  /* ********************************************************************* */
  /*                                                                       */
  /* basic type definitions                                                */
  /*                                                                       */
  /* ********************************************************************* */

  typedef double norm_type;


  /* ********************************************************************* */
  /*                                                                       */
  /*    Convention pour les types de matrices et de vecteurs.              */
  /*    On doit ensuite interfacer les types de MTL et Blitz ...           */
  /*                                                                       */
  /* ********************************************************************* */
  /* vectors :                                                             */
  /*   copy(w, v)           v <--- w                                       */
  /*   add(w, v)            v <--- v + w                                   */
  /*   add_mult(v, a, w)    v <--- v + a*w                                 */
  /*   dot(v, w)            v.w                                            */
  /*   mult(v, a)           v <--- a*v                                     */
  /*                                                                       */
  /* matrices :                                                            */
  /*                                                                       */
  /* ********************************************************************* */





}  /* end of namespace gensolv.                                            */


#endif /* __GENSOLV_FUNCTION_H                                             */
