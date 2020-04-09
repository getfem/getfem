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

#include <getfemint_workspace.h>
#include <getfem/getfem_mesh_im.h>
#include <getfem/getfem_integration.h>
#include <getfem/getfem_mat_elem.h>
/*
  $Id$
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
      ids.push_back(store_integ_object(mim.int_method_of_element(cv)));
    else ids.push_back(id_type(-1));
  }
  out.return_packed_obj_ids(ids, INTEG_CLASS_ID);
}

/*@GFDOC
  General function extracting information from mesh_im objects.
  @*/




// Object for the declaration of a new sub-command.

struct sub_gf_mim_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
                   getfemint::mexargs_out& out,
                   getfem::mesh_im *mim) = 0;
};

typedef std::shared_ptr<sub_gf_mim_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_mim_get {				\
      virtual void run(getfemint::mexargs_in& in,			\
                       getfemint::mexargs_out& out,			\
                       getfem::mesh_im *mim)				\
      { dummy_func(in); dummy_func(out);  dummy_func(mim); code }	\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }





void gf_mesh_im_get(getfemint::mexargs_in& m_in,
                    getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {

    /*@GET @CELL{I, CV2I} = ('integ'[, @mat CVids])
    Return a list of integration methods used by the @tmim.

    `I` is an array of all @tinteg objects found in the convexes
    given in `CVids`. If `CV2I` was supplied as an output argument, it
    contains, for each convex listed in `CVids`, the index of its
    correspounding integration method in `I`.

    Convexes which are not part of the mesh, or convexes which do
    not have any integration method have their correspounding entry
    in `CV2I` set to -1.
    @MATLAB{
    Example::

      cvid=gf_mesh_get(mim,'cvid');
      [f,c2f]=gf_mesh_im_get(mim, 'integ');
      for i=1:size(f), sf{i}=gf_integ_get('char',f(i)); end;
      for i=1:size(c2f),
        disp(sprintf('the integration of convex %d is %s',...
        cvid(i),sf{i}));
      end;
    }
@*/
    sub_command
      ("integ", 0, 1, 0, 2,
       get_integ_of_convexes(*mim, in, out);
       );


    /*@GET CVids = ('convex_index')
    Return the list of convexes who have a integration method.

    Convexes who have the dummy IM_NONE method are not listed.@*/
    sub_command
      ("convex_index", 0, 0, 0, 1,
       dal::bit_vector bv = mim->convex_index();
       for (dal::bv_visitor ic(mim->convex_index()); !ic.finished(); ++ic) {
         if (mim->int_method_of_element(ic)->type() == getfem::IM_NONE)
           bv.sup(ic);
       }
       out.pop().from_bit_vector(bv);
       );


    /*@GET M = ('eltm', @teltm em, @int cv [, @int f])
    Return the elementary matrix (or tensor) integrated on the convex `cv`.

    **WARNING**

    Be sure that the fem used for the construction of `em` is compatible
    with the fem assigned to element `cv` ! This is not checked by the
    function ! If the argument `f` is given, then the elementary tensor
    is integrated on the face `f` of `cv` instead of the whole convex.@*/
    sub_command
      ("eltm", 2, 3, 0, 1,
       getfem::pmat_elem_type pmet = to_eltm_object(in.pop());
       size_type cv = in.pop().to_convex_number(mim->linked_mesh());
       /* one should check that the fem given to the MET is
          compatible with the fem of the element (not easy ..) */
       getfem::base_tensor t;
       /* if the convex has a IM, then it has been added to the convex index
          of the mesh_fem
       */
       if (!mim->convex_index()[cv])
	 THROW_ERROR("convex " << cv+config::base_index()
		     << " has no integration method!");
       getfem::pmat_elem_computation pmec =
       getfem::mat_elem(pmet,
                        mim->int_method_of_element(cv) ,
                        mim->linked_mesh().trans_of_convex(cv));
       if (!in.remaining()) {
         pmec->gen_compute(t, mim->linked_mesh().points_of_convex(cv), cv);
       } else {
         short_type nbf =
           mim->linked_mesh().structure_of_convex(cv)->nb_faces();
         short_type f = in.pop().to_face_number(nbf);
         pmec->gen_compute_on_face(t, mim->linked_mesh().points_of_convex(cv), f, cv);
       }
       out.pop().from_tensor(t);
       );


    /*@GET Ip = ('im_nodes'[, @mat CVids])
    Return the coordinates of the integration points, with their weights.

    `CVids` may be a list of convexes, or a list of convex faces, such
    as returned by MESH:GET('region')

    **WARNING**

    Convexes which are not part of the mesh, or convexes which
    do not have an approximate integration method do not have
    their corresponding entry (this has no meaning for exact
    integration methods!).@*/
    sub_command
      ("im_nodes", 0, 1, 0, 1,
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
       darray ww = out.pop().create_darray(N+1, unsigned(tmp.size() / (N+1)));
       if (tmp.size()) std::copy(tmp.begin(), tmp.end(), &ww[0]);
       );


    /*@GET ('save',@str filename[, 'with mesh'])
      Saves a @tmim in a text file (and optionally its linked mesh object).@*/
    sub_command
      ("save", 1, 2, 0, 0,
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
       );


    /*@GET ('char'[,'with mesh'])
      Output a string description of the @tmim.
      
      By default, it does not include the description of the linked
      @tmesh object.@*/
    sub_command
      ("char", 0, 0, 0, 1,
       std::stringstream s;
       if (in.remaining() && cmd_strmatch(in.pop().to_string(),"with mesh"))
         mim->linked_mesh().write_to_file(s);
       mim->write_to_file(s);
       out.pop().from_string(s.str().c_str());
       );


    /*@GET ('display')
      displays a short summary for a @tmim object.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       infomsg() << "gfMeshIm object in dimension "
       << int(mim->linked_mesh().dim())
       << " with " << mim->linked_mesh().nb_points() << " points and "
       << mim->linked_mesh().convex_index().card() << " elements\n";
       );


    /*@GET m = ('linked mesh')
    Returns a reference to the @tmesh object linked to `mim`.@*/
    sub_command
      ("linked mesh", 0, 0, 0, 1,
       id_type id = workspace().object((const void *)(&mim->linked_mesh()));
       if (id == id_type(-1)) {
	 auto pst = workspace().hidden_object(workspace().object(mim),
					      &mim->linked_mesh());
	 if (!pst.get()) THROW_INTERNAL_ERROR;
	 std::shared_ptr<getfem::mesh> pm = 
	   std::const_pointer_cast<getfem::mesh>
	   (std::dynamic_pointer_cast<const getfem::mesh>(pst));
	 id = store_mesh_object(pm);
       }
       out.pop().from_object_id(id, MESH_CLASS_ID);
       );


    /*@GET z = ('memsize')
    Return the amount of memory (in bytes) used by the @tmim object.

    The result does not take into account the linked @tmesh object.@*/
    sub_command
      ("memsize", 0, 0, 0, 1,
       out.pop().from_integer(int(mim->memsize()));
       );

  }


  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  getfem::mesh_im *mim   = to_meshim_object(m_in.pop());
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
              it->second->arg_in_max, it->second->arg_out_min,
              it->second->arg_out_max);
    it->second->run(m_in, m_out, mim);
  }
  else bad_cmd(init_cmd);

}
