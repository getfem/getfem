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
#include <getfem/bgeot_geometric_trans.h>

using namespace getfemint;

/*@GFDOC
   The geometric transformation must be used when you are building a custom
   mesh convex by convex (see the add_convex() function of @tmesh): it also
   defines the kind of convex (triangle, hexahedron, prism, etc..)
  @*/

/*@INIT GT = ('.name', @str name)

The name argument contains the specification of the geometric transformation
as a string, which may be:

  - GT_PK(n,k) :
    Transformation on simplexes, dim `n`, degree `k`.
  - GT_QK(n,k) :
    Transformation on parallelepipeds, dim `n`, degree `k`.
  - GT_PRISM(n,k) :
    Transformation on prisms, dim `n`, degree `k`.
  - GT_PRODUCT(A,B) :
    Tensorial product of two transformations.
  - GT_LINEAR_PRODUCT(@tgt gt1,@tgt gt2) :
    Linear tensorial product of two transformations
@*/

void gf_geotrans(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  std::string cmd = in.pop().to_string();
  id_type id = store_geotrans_object(bgeot::geometric_trans_descriptor(cmd));
  out.pop().from_object_id(id, GEOTRANS_CLASS_ID);
}
