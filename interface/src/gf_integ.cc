/*===========================================================================

 Copyright (C) 2006-2020 Yves Renard, Julien Pommier.

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/

#include <getfemint.h>
#include <getfem/getfem_integration.h>

using namespace getfemint;

/*@GFDOC
  General object for obtaining handles to various integrations methods on
  convexes (used when the elementary matrices are built).
@*/



void gf_integ(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  /*@INIT I = ('.init', @str method)
    Here is a list of some integration methods defined in GetFEM (see the
    description of finite element and integration methods for a complete
    reference):

     - IM_EXACT_SIMPLEX(n) :
       Exact integration on simplices (works only with linear geometric
       transformations and PK @tfem's).
     - IM_PRODUCT(A,B) :
       Product of two integration methods.
     - IM_EXACT_PARALLELEPIPED(n) :
       Exact integration on parallelepipeds.
     - IM_EXACT_PRISM(n) :
       Exact integration on prisms.
     - IM_GAUSS1D(k) :
       Gauss method on the segment, order `k=1,3,...,99`.
     - IM_NC(n,k) :
       Newton-Cotes approximative integration on simplexes, order `k`.
     - IM_NC_PARALLELEPIPED(n,k) :
       Product of Newton-Cotes integration on parallelepipeds.
     - IM_NC_PRISM(n,k) :
       Product of Newton-Cotes integration on prisms.
     - IM_GAUSS_PARALLELEPIPED(n,k) :
       Product of Gauss1D integration on parallelepipeds.
     - IM_TRIANGLE(k) :
       Gauss methods on triangles `k=1,3,5,6,7,8,9,10,13,17,19`.
     - IM_QUAD(k) :
       Gauss methods on quadrilaterons `k=2,3,5, ...,17`. Note that
       IM_GAUSS_PARALLELEPIPED should be prefered for QK @tfem's.
     - IM_TETRAHEDRON(k) :
       Gauss methods on tetrahedrons `k=1,2,3,5,6 or 8`.
     - IM_SIMPLEX4D(3) :
       Gauss method on a 4-dimensional simplex.
     - IM_STRUCTURED_COMPOSITE(im,k) :
       Composite method on a grid with `k` divisions.
     - IM_HCT_COMPOSITE(im) :
       Composite integration suited to the HCT composite finite element.

    Example:

     - I = INTEG:INIT('IM_PRODUCT(IM_GAUSS1D(5),IM_GAUSS1D(5))')

    is the same as:

     - I = INTEG:INIT('IM_GAUSS_PARALLELEPIPED(2,5)')

    Note that 'exact integration' should be avoided in general, since they
    only apply to linear geometric transformations, are quite slow, and
    subject to numerical stability problems for high degree @tfem's. @*/
  std::string cmd = in.pop().to_string();
  out.pop().from_object_id
    (store_integ_object(getfem::int_method_descriptor(cmd)), INTEG_CLASS_ID);
}
