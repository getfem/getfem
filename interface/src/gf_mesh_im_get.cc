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
#include <getfemint_integ.h>
/*
  $Id$

  ChangeLog:
  $Log: gf_mesh_im_get.cc,v $
  Revision 1.4  2006/03/28 10:06:35  pommier
  *** empty log message ***

  Revision 1.3  2006/02/14 17:57:17  pommier
  *** empty log message ***

  Revision 1.2  2006/01/18 11:21:52  pommier
  *** empty log message ***

  Revision 1.1  2005/03/08 16:50:12  pommier
  added meshim, many doc updates

 */


using namespace getfemint;

static void
get_integ_of_convexes(const getfem::mesh_im& mim, mexargs_in& in, mexargs_out& out)
{  
  dal::bit_vector cvlst;
  if (in.remaining()) cvlst = in.pop().to_bit_vector(&mim.linked_mesh().convex_index());
  else { cvlst = mim.linked_mesh().convex_index(); }
  std::vector<id_type> ids; ids.reserve(cvlst.card());
  for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv) {
    if (mim.convex_index().is_in(cv))
      ids.push_back(ind_integ(mim.int_method_of_element(cv)));
    else ids.push_back(id_type(-1));
  }
  out.return_packed_obj_ids(ids, INTEG_CLASS_ID);
}

/*MLABCOM
  FUNCTION [x] = gf_mesh_im_get(meshim MIM, operation [, args])

  General function extracting information from mesh_im objects.

  @GET    MESHIM:GET('integ')
  Example:
     cvid=gf_mesh_get(mim,'cvid');
     [f,c2f]=gf_mesh_im_get(mim, 'integ');
     for i=1:size(f), sf{i}=gf_integ_get('char',f(i)); end;
     for i=1:size(c2f),
       disp(sprintf('the integration of convex %d is %s',...
            cvid(i),sf{i}));
     end;
  @GET    MESHIM:GET('eltm')
  @GET    MESHIM:GET('save')
  @GET    MESHIM:GET('char')
  @GET    MESHIM:GET('linked mesh')
  @GET    MESHIM:GET('memsize')

  $Id$
MLABCOM*/

void gf_mesh_im_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfemint_mesh_im *mi_mim = in.pop().to_getfemint_mesh_im();
  getfem::mesh_im *mim   = &mi_mim->mesh_im();
  std::string cmd        = in.pop().to_string();
  if (check_cmd(cmd, "integ", in, out, 0, 1, 0, 2)) {
    /*@GET [INTEG, CV2I] = MESHIM:GET('integ' [, CVLIST])
      Return a list of integration methods used by the @tmim.

      INTEG is an array of all @tinteg objects found in the convexes
      given in CVLST. If CV2F was supplied as an output argument, it
      contains, for each convex listed in CVLST, the index of its
      correspounding integration method in INTEG.

      Convexes which are not part of the mesh, or convexes which do
      not have any integration method have their correspounding entry
      in CV2I set to -1.
      @*/
    get_integ_of_convexes(*mim, in, out);
  } else if (check_cmd(cmd, "convex_index", in, out, 0, 0, 0, 1)) {
    /*@GET CVLST = MESHFEM:GET('convex_index')
      Return the list of convexes who have a integration
      method. Convexes who have the dummy IM_NONE method are not
      listed.
      @*/
    dal::bit_vector bv = mim->convex_index();
    for (dal::bv_visitor ic(mim->convex_index()); !ic.finished(); ++ic) {
      if (mim->int_method_of_element(ic)->type() == getfem::IM_NONE)
	bv.sup(ic);
    }
    out.pop().from_bit_vector(bv);
  } else if (check_cmd(cmd, "eltm", in, out, 2, 3, 0, 1)) {
    /*@GET M = MESHIM:GET('eltm', @eltm MET, @int CV [@int F])
      Return the elementary matrix (or tensor) integrated on the convex CV.

      !!WARNING!! Be sure that the fem used for the construction of
      MET is compatible with the fem assigned to element CV ! This is
      not checked by the function ! If the argument F is given, then
      the elementary tensor is integrated on the face F of CV instead
      of the whole convex.
      @*/
    getfem::pmat_elem_type pmet = in.pop().to_mat_elem_type();
    size_type cv = in.pop().to_convex_number(mim->linked_mesh());
    /* one should check that the fem given to the MET is 
       compatible with the fem of the element (not easy ..) */
    getfem::base_tensor t;
    /* if the convex has a IM, then it has been added to the convex index 
       of the mesh_fem 
    */
    check_cv_im(*mim, cv);
    getfem::pmat_elem_computation pmec = 
      getfem::mat_elem(pmet, 
		       mim->int_method_of_element(cv) , 
		       mim->linked_mesh().trans_of_convex(cv));
    if (!in.remaining()) {
      pmec->gen_compute(t, mim->linked_mesh().points_of_convex(cv), cv);
    } else {
      unsigned nbf = 
	mim->linked_mesh().structure_of_convex(cv)->nb_faces();
      size_type f = in.pop().to_face_number(nbf);
      pmec->gen_compute_on_face(t, mim->linked_mesh().points_of_convex(cv), f, cv);
    }
    out.pop().from_tensor(t);
  } else if (check_cmd(cmd, "im_nodes", in, out, 0, 1, 0, 1)) {
    /*@GET MESHIM:GET('im_nodes' [, CVLST]) 
      Return the coordinates of the integration points, with their
      weights.

      CVLST may be a list of convexes, or a list of convex faces, such
      as returned by MESH:GET('region') @*/
    getfem::base_vector tmp;
    unsigned N = mim->linked_mesh().dim();
    getfem::mesh_region mr;
    if (in.remaining()) mr = in.pop().to_mesh_region();
    else mr = getfem::mesh_region(mim->linked_mesh().convex_index());
    for (getfem::mr_visitor ir(mr); !ir.finished(); ++ir) {
      size_type cv = ir.cv();
      if (!mim->convex_index().is_in(cv)) continue;
      getfem::pintegration_method pim = mim->int_method_of_element(cv);
      if (pim->type() != getfem::IM_APPROX) continue;
      getfem::papprox_integration pai = pim->approx_method();
      bgeot::pgeometric_trans pgt = mim->linked_mesh().trans_of_convex(cv);
      size_type nbpt = 
	(ir.is_face() ? pai->nb_points_on_face(ir.f()) : pai->nb_points_on_convex());
      for (unsigned ii=0; ii < nbpt; ++ii) {
	getfem::base_node Pref;
	scalar_type w;
	if (ir.is_face()) {
	  Pref = pai->point_on_face(ir.f(), ii);
	  w = pai->coeff_on_face(ir.f(), ii);
	} else {
	  Pref = pai->point(ii);
	  w = pai->coeff(ii);
	}
	getfem::base_node P = pgt->transform(Pref, mim->linked_mesh().points_of_convex(cv));
	for (unsigned j=0; j < N; ++j) tmp.push_back(P[j]);
	tmp.push_back(w);
      }
    }
    darray ww = out.pop().create_darray(N+1, tmp.size() / (N+1));
    std::copy(tmp.begin(), tmp.end(), &ww[0]);
  } else if (check_cmd(cmd, "save", in, out, 1, 2, 0, 0)) {
    /*@GET MESHIM:GET('save', filename [,'with mesh'])
      Saves a @tmim in a text file (and optionaly its linked mesh object).
      @*/
    std::string s = in.pop().to_string();
    bool with_mesh = false;
    if (in.remaining()) {
      if (cmd_strmatch(in.pop().to_string(), "with mesh")) {
	with_mesh = true;
      } else THROW_BADARG("expecting string 'with mesh'");
    }
    std::ofstream o(s.c_str());
    if (!o) THROW_ERROR("impossible to write in file '" << s << "'");
    o << "% GETFEM MESH+FEM FILE " << endl;
    o << "% GETFEM VERSION " << GETFEM_VERSION << endl;
    if (with_mesh) mim->linked_mesh().write_to_file(o);
    mim->write_to_file(o);
    o.close();
  } else if (check_cmd(cmd, "char", in, out, 0, 0, 0, 1)) {
    /*@GET MESHIM:GET('char' [,'with mesh'])
      Output a string description of the mesh_im. 

      By default, it does not include the description of the linked
      mesh object. @*/
    std::stringstream s;
    if (in.remaining() && cmd_strmatch(in.pop().to_string(),"with mesh"))
      mim->linked_mesh().write_to_file(s);
    mim->write_to_file(s);
    out.pop().from_string(s.str().c_str());
  } else if (check_cmd(cmd, "linked mesh", in, out, 0, 0, 0, 1)) {
    /*@GET M=MESHIM:GET('linked mesh')
      Returns a reference to the mesh object linked to MIM.
      @*/
    out.pop().from_object_id(mi_mim->linked_mesh_id(), MESH_CLASS_ID);
  } else if (check_cmd(cmd, "memsize", in, out, 0, 0, 0, 1)) {
    /*@GET MESHIM:GET('memsize')
      Return the amount of memory (in bytes) used by the mesh_im object.

      The result does not take into account the linked mesh object.
      @*/
    out.pop().from_integer(mim->memsize());
  } else bad_cmd(cmd);
}
