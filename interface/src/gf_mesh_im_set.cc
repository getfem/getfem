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

#include <getfemint_misc.h>
#include <getfem/getfem_mesh_im.h>
#include <getfem/getfem_mesh_im_level_set.h>

using namespace getfemint;


static void gf_mesh_im_set_integ_(getfem::mesh_im *mim,
				  getfemint::mexargs_in& in) {
  getfem::pintegration_method pim = to_integ_object(in.pop());

  bool all_cv = false;
  /* check or build the convex list */
  dal::bit_vector bv;
  if (in.remaining() == 1)
    bv = in.pop().to_bit_vector(&mim->linked_mesh().convex_index(),
				-config::base_index());
  else
    all_cv = true;

  /* check for the validity of the operation */
  for (dal::bv_visitor cv(bv); !cv.finished(); ++cv) {
    //if (!mim->linked_mesh().convex_index().is_in(cv))
    //  THROW_ERROR("Convex " << cv+config::base_index()
    //                        << " was not found in mesh");
    if (pim->structure() !=
	bgeot::basic_structure(mim->linked_mesh().structure_of_convex(cv)))
      infomsg() << "Warning: structure of the Integration Method seems "
	"to be incompatible with the structure of the convex\n";
  }

  /* all the work done here */
  if (!all_cv) {
    mim->set_integration_method(bv, pim);
  }  else {
    mim->set_integration_method(pim);
  }

}

/* set the classical integ of order IM_DEGREE on the mesh_im, with a classical integration
   method */
static void gf_mesh_im_set_classical_integ(getfem::mesh_im *mim,
					   getfemint::mexargs_in& in) {
  dim_type IM_DEGREE = dim_type(-1);
  if (in.remaining()) IM_DEGREE = dim_type(in.pop().to_integer(-1,255));
  bool all_cv = false;
  dal::bit_vector bv;
  if (in.remaining() == 1) {
    bv = in.pop().to_bit_vector(&mim->linked_mesh().convex_index(),
				-config::base_index());
  } else
    all_cv = true;
  if (!all_cv) {
    mim->set_integration_method(bv,IM_DEGREE);
  } else {
    mim->set_integration_method(IM_DEGREE);
  }
  
}

/* WARNING: gf_mesh_im.cc also uses this function! do not change its interface! */
void gf_mesh_im_set_integ(getfem::mesh_im *mim, getfemint::mexargs_in& in) {
  if (in.front().is_object_id())
    gf_mesh_im_set_integ_(mim, in);
  else
    gf_mesh_im_set_classical_integ(mim, in);
}

/*@GFDOC
  General function for modifying mesh_im objects
  @*/

void gf_mesh_im_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }

  getfem::mesh_im *mim = to_meshim_object(in.pop());
  std::string cmd = in.pop().to_string();
  if (check_cmd(cmd, "integ", in, out, 1, 2, 0, 0)) {
    /*@SET ('integ',{@tinteg im|@int im_degree}[, @ivec CVids])
    Set the integration method.

    Assign an integration method to all convexes whose #ids are
    listed in `CVids`. If `CVids` is not given, the integration is
    assigned to all convexes. It is possible to assign a specific
    integration method with an integration method handle `im` obtained
    via INTEG:INIT('IM_SOMETHING'), or to let getfem choose a suitable
    integration method with `im_degree` (choosen such that polynomials
    of :math:`\text{degree} \leq \text{im\_degree}` are exactly integrated.
    If `im_degree=-1`, then the dummy integration method IM_NONE will 
    be used.)@*/
    gf_mesh_im_set_integ(mim, in);
  } else if (check_cmd(cmd, "adapt", in, out, 0, 0, 0, 0)) {
    /*@SET ('adapt')
    For a @tmim levelset object only. Adapt the integration methods to a
    change of the levelset function.@*/
    getfem::mesh_im_level_set *mimls
      = dynamic_cast<getfem::mesh_im_level_set *>(mim);
    if (!mimls) THROW_BADARG("The command 'adapt' can only be "
			     "applied to a mesh_im_level_set object");
    mimls->adapt();
  } else bad_cmd(cmd);
}
