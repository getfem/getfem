// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2006-2006 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

#include <getfemint.h>
#include <getfemint_integ.h>

using namespace getfemint;

/*MLABCOM
  @TEXT INTEG:INIT('INTEG_init')
MLABCOM*/

/*@TEXT INTEG:INIT('INTEG_init')
    INTEG:INIT(@string method)
    Return a FEM Integration Method from a string description.<Par>

    Possible descriptions are (non exhaustive list):<par>
    * IM_EXACT_SIMPLEX(n)<par>
      Exact integration on simplices (works only with linear geometric
      transformations and PK FEMs). Note that 'exact integration'
      should be avoided in general, since they only apply to linear
      geometric transformations, are quite slow, and subject to
      numerical stability problems for high degree FEMs.

      <par>
    * IM_PRODUCT(a, b)
      Product of two integration methods.<Par>

    * IM_EXACT_PARALLELEPIPED(n)<par>
      Exact integration on parallelepipeds.<Par>

    * IM_EXACT_PRISM(n)      <par>
      Exact integration on prisms.<Par>

    * IM_GAUSS1D(K)       <par>
      Gauss method on the segment, order K (K=1,3,...99).<Par>

    * IM_NC(N,K)          <par>
      Newton-Cotes approximative integration on simplexes, order K.<Par>

    * IM_NC_PARALLELEPIPED(N,K) <par>
      Product of Newton-Cotes integration on parallelepipeds.<Par>

    * IM_NC_PRISM(N,K) <par>
      Product of Newton-Cotes integration on prisms.<Par>

    * IM_GAUSS_PARALLELEPIPED(N,K) <par>
      Product of Gauss1D integration on parallelepipeds.<Par>

    * IM_TRIANGLE(K)<par>
      Gauss methods on triangles (K=1,3,5,6,7,8,9,10,13,17,19).<Par>

    * IM_QUAD(K)<par>
      Gauss methods on quadrilaterons (K=2, 3, 5, .. 17). Note that IM_GAUSS_PARALLELEPIPED should be prefered for QK FEMs.<Par>

    * IM_TETRAHEDRON(K)<par>
      Gauss methods on tetrahedrons (K=1, 2, 3, 5, 6 or 8).<Par>

    * IM_SIMPLEX4D(3)<par>
      Gauss method on a 4-dimensional simplex.<Par>

    * IM_STRUCTURED_COMPOSITE(IM1, K)<par> 
      Composite method on a grid with K divisions.<Par>

    * IM_HCT_COMPOSITE(IM1)<par>
      Composite integration suited to the HCT composite finite element.<Par>

    example:<par>
      INTEG:INIT('IM_PRODUCT(IM_GAUSS1D(5),IM_GAUSS1D(5))')<par>
      is the same as:<par>
      INTEG:INIT('IM_GAUSS_PARALLELEPIPED(2,5)')<par>
      @*/

void gf_integ(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  std::string cmd = in.pop().to_string();
  out.pop().from_object_id(getfemint::ind_integ(getfem::int_method_descriptor(cmd)),
			   INTEG_CLASS_ID);
}

