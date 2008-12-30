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
#include <getfemint_pgt.h>

using namespace getfemint;

/*MLABCOM
  FUNCTION I = gf_geotrans(name)
    @TEXT GEOTRANS:INIT('GEOTRANS_init')
MLABCOM*/

/*@TEXT GEOTRANS:INIT('GEOTRANS_init')
GEOTRANS:INIT(@str descriptors)<Par>
Possible descriptions are (non exhaustive list):<Par>
* GT_PK(N,K)<par>
   Transformation on simplexes, dim N, degree K.<par>
* GT_QK(N,K)<par>
   Transformation on parallelepipeds, dim N, degree K.<par>
* GT_PRISM(N,K)<par>
   Transformation on prisms, dim N, degree K.<par>
* GT_PRODUCT(a,b)<par>
   Tensorial product of two transformations.<par>
* GT_LINEAR_PRODUCT(a,b)<par>
   Linear tensorial product of two transformations
@*/

void gf_geotrans(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  std::string cmd = in.pop().to_string();
  out.pop().from_object_id(getfemint::ind_pgt(bgeot::geometric_trans_descriptor(cmd)),
			   GEOTRANS_CLASS_ID);
}
