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

#include <getfemint_misc.h>
#include <getfemint_mesh_im.h>

using namespace getfemint;


static void gf_mesh_im_set_integ_(getfem::mesh_im *mim, getfemint::mexargs_in& in)
{
  getfem::pintegration_method pim = 0;
  pim = in.pop().to_integration_method();

  /* check or build the convex list */
  dal::bit_vector bv;
  if (in.remaining() == 1) {
    bv = in.pop().to_bit_vector(&mim->linked_mesh().convex_index(), -1);
  } else {
    bv = mim->linked_mesh().convex_index();
  }
  
  /* check for the validity of the operation */
  for (dal::bv_visitor cv(bv); !cv.finished(); ++cv) {
    if (!mim->linked_mesh().convex_index().is_in(cv))
      THROW_ERROR("Convex " << cv+config::base_index() << " was not found in mesh");
    if (pim->structure() != 
	mim->linked_mesh().structure_of_convex(cv)->basic_structure())
      infomsg() << "Warning: structure of the Integration Method seems to be incompatible with the structure of the convex\n";
  }
    
  /* all the work done here */
  mim->set_integration_method(bv, pim);
}

/* set the classical integ of order IM_DEGREE on the mesh_im, with a classical integration
   method */
static void gf_mesh_im_set_classical_integ(getfem::mesh_im *mim, getfemint::mexargs_in& in) {
  dim_type IM_DEGREE = dim_type(-1);
  if (in.remaining()) IM_DEGREE = in.pop().to_integer(-1,255);
  dal::bit_vector bv;
  if (in.remaining() == 1) {
    bv = in.pop().to_bit_vector(&mim->linked_mesh().convex_index(), -1);
  } else {
    bv = mim->linked_mesh().convex_index();
  }
  mim->set_integration_method(bv,IM_DEGREE);
}

/* WARNING: gf_mesh_im.cc also uses this function! do not change its interface! */
void gf_mesh_im_set_integ(getfem::mesh_im *mim, getfemint::mexargs_in& in) {
  if (in.front().is_object_id())
    gf_mesh_im_set_integ_(mim, in);
  else 
    gf_mesh_im_set_classical_integ(mim, in);
}

/*MLABCOM
  FUNCTION [x] = gf_mesh_im_set(meshim MIM, operation [, args])

  General function for modifying mesh_im objects
  
  @SET MESHIM:SET('integ')

  $Id$
MLABCOM*/

void gf_mesh_im_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }

  getfem::mesh_im *mim = in.pop().to_mesh_im();
  std::string cmd        = in.pop().to_string();
  if (check_cmd(cmd, "integ", in, out, 1, 2, 0, 0)) {
    /*@SET MESHIM:SET('integ', {@integ IM|@int IM_DEGREE} [, @ivec CVIDX])
      Set the integration method.
      
      Assign an integration method to all convexes whose #ids are
      listed in CVIDX. If CVIDX is not given, the integration is
      assigned to all convexes. It is possible to assign a specific
      integration method with an integration method handle IM obtained
      via INTEG:INIT('IM_SOMETHING'), or to let getfem choose a
      suitable integration method with IM_DEGREE (choosen such that
      polynomials of degree <= IM_DEGREE are exactly integrated. If
      IM_DEGREE=-1, then the dummy integration method IM_NONE will be
      used.) @*/
    gf_mesh_im_set_integ(mim, in);
  } else bad_cmd(cmd);
}
