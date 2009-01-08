// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2008 Yves Renard, Julien Pommier.
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
//===========================================================================

#include <getfemint.h>
#include <getfemint_integ.h>

using namespace getfemint;

/*MLABCOM
  @TEXT INTEG:INIT('INTEG_init')
MLABCOM*/

/*@TEXT INTEG:INIT('INTEG_init')
@tinteg = INTEG:INIT(@str method)<Par>

Here is a list of some integration methods defined in getfem++<par>
(see the description of finite element and integration methods<par>
for a complete reference):<Par>

* IM_EXACT_SIMPLEX(n)<par>
   Exact integration on simplices (works only with linear geometric<par>
   transformations and PK @tfem's).<par>
* IM_PRODUCT(A,B)<par>
   Product of two integration methods.<par>
* IM_EXACT_PARALLELEPIPED(n)<par>
   Exact integration on parallelepipeds.<par>
* IM_EXACT_PRISM(n)<par>
   Exact integration on prisms.<par>
* IM_GAUSS1D(k)<par>
   Gauss method on the segment, order `k=1,3,...99`.<par>
* IM_NC(n,k)<par>
   Newton-Cotes approximative integration on simplexes, order `k`.<par>
* IM_NC_PARALLELEPIPED(n,k)<par>
   Product of Newton-Cotes integration on parallelepipeds.<par>
* IM_NC_PRISM(n,k)<par>
   Product of Newton-Cotes integration on prisms.<par>
* IM_GAUSS_PARALLELEPIPED(n,k) <par>
   Product of Gauss1D integration on parallelepipeds.<par>
* IM_TRIANGLE(k)<par>
   Gauss methods on triangles `k=1,3,5,6,7,8,9,10,13,17,19`.<par>
* IM_QUAD(k)<par>
   Gauss methods on quadrilaterons `k=2, 3, 5, .. 17`. Note that
   IM_GAUSS_PARALLELEPIPED should be prefered for QK @tfem's.<par>
* IM_TETRAHEDRON(k)<par>
   Gauss methods on tetrahedrons `k=1, 2, 3, 5, 6 or 8`.<par>
* IM_SIMPLEX4D(3)<par>
   Gauss method on a 4-dimensional simplex.<par>
* IM_STRUCTURED_COMPOSITE(im,k)<par>
   Composite method on a grid with `k` divisions.<par>
* IM_HCT_COMPOSITE(im)<par>
   Composite integration suited to the HCT composite finite element.<Par>

Example:<par>
- INTEG:INIT('IM_PRODUCT(IM_GAUSS1D(5),IM_GAUSS1D(5))')<par>
is the same as:<par>
- INTEG:INIT('IM_GAUSS_PARALLELEPIPED(2,5)')<Par>

Note that 'exact integration' should be avoided in general, since they<par>
only apply to linear geometric transformations, are quite slow, and<par>
subject to numerical stability problems for high degree @tfem's.<par>
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
