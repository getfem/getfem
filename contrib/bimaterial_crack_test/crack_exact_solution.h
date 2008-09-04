// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2007-2008 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// As a special exception, you  may use  this file  as it is a part of a free
// software  library  without  restriction.  Specifically,  if   other  files
// instantiate  templates  or  use macros or inline functions from this file,
// or  you compile this  file  and  link  it  with other files  to produce an
// executable, this file  does  not  by itself cause the resulting executable
// to be covered  by the GNU Lesser General Public License.  This   exception
// does not  however  invalidate  any  other  reasons why the executable file
// might be covered by the GNU Lesser General Public License.
//
//===========================================================================

#include "getfem/getfem_mesh_fem_global_function.h"


/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

#define VALIDATE_XFEM

#ifdef VALIDATE_XFEM

struct crack_exact_solution_function : public getfem::abstract_xy_function {
  unsigned function_num;
  unsigned component_num; /* 0 -> x component, 1 -> y component */
  scalar_type lambda, mu;
  virtual scalar_type val(scalar_type x, scalar_type y) const;
  virtual base_small_vector grad(scalar_type x, scalar_type y) const;
  virtual base_matrix hess(scalar_type, scalar_type) const
  { GMM_ASSERT1(false, "Sorry, to be done ..."); }
  crack_exact_solution_function(unsigned fnum, 
				unsigned cnum,
				scalar_type l, scalar_type m) {
    function_num = fnum; 
    component_num = cnum;
    lambda = l; mu = m;
  }
  base_small_vector eval(const base_node &x, base_matrix *pgrad) const;
};


struct crack_exact_solution {
  getfem::mesh_fem_global_function mf;
  getfem::base_vector U;

  crack_exact_solution(getfem::mesh &me) : mf(me) {}
  
  void init(int function_num, scalar_type lambda, scalar_type mu,
	    getfem::level_set &ls);
};

inline base_small_vector sol_f(const base_node &x) {
  int N = x.size();
  base_small_vector res(N);
  return res;
}

#else

inline base_small_vector sol_f(const base_node &x) {
  int N = x.size();
  base_small_vector res(N); res[N-1] = x[N-1];
  return res;
}

#endif
