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
#include <getfem/bgeot_convex_structure.h>

using namespace getfemint;

/*@GFDOC
  General function for querying information about convex_structure objects.

  The convex structures are internal structures of GetFEM. They do not
  contain points positions. These structures are recursive, since the faces
  of a convex structures are convex structures.
@*/



// Object for the declaration of a new sub-command.

struct sub_gf_cvstruct_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
                   getfemint::mexargs_out& out,
                   const bgeot::pconvex_structure &cs) = 0;
};

typedef std::shared_ptr<sub_gf_cvstruct_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_cvstruct_get {				\
      virtual void run(getfemint::mexargs_in& in,			\
                       getfemint::mexargs_out& out,			\
                       const bgeot::pconvex_structure &cs)		\
      { dummy_func(in); dummy_func(out); dummy_func(cs); code }		\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }





void gf_cvstruct_get(getfemint::mexargs_in& m_in,
                     getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {

    /*@RDATTR n = ('nbpts')
      Get the number of points of the convex structure.@*/
    sub_command
      ("nbpts", 0, 0, 0, 1,
       out.pop().from_scalar(cs->nb_points());
       );


    /*@RDATTR d = ('dim')
      Get the dimension of the convex structure.@*/
    sub_command
      ("dim", 0, 0, 0, 1,
       out.pop().from_scalar(cs->dim());
       );


    /*@RDATTR cs = ('basic structure')
    Get the simplest convex structure.

    For example, the 'basic structure' of the 6-node triangle, is the
    canonical 3-noded triangle.@*/
    sub_command
      ("basic_structure", 0, 0, 0, 1,
       out.pop().from_object_id
       (store_cvstruct_object(bgeot::basic_structure(cs)),
        CVSTRUCT_CLASS_ID);
       );


    /*@RDATTR cs = ('face', @int F)
      Return the convex structure of the face `F`.@*/
    sub_command
      ("face", 1, 1, 0, 1,
       short_type f = in.pop().to_face_number(cs->nb_faces());
       out.pop().from_object_id
       (store_cvstruct_object(cs->faces_structure()[f]),
        CVSTRUCT_CLASS_ID);
       );


    /*@GET I = ('facepts', @int F)
      Return the list of point indices for the face `F`.@*/
    sub_command
      ("facepts", 1, 1, 0, 1,
       short_type f = short_type(in.pop().to_face_number(cs->nb_faces()));
       iarray w = out.pop().create_iarray_h(cs->nb_points_of_face(f));
       for (size_type i=0; i < w.size(); ++i)
         w[i] = cs->ind_points_of_face(f)[i]+config::base_index();
       );

    /*@GET s = ('char')
      Output a string description of the @tcvstruct.@*/
    sub_command
      ("char", 0, 0, 0, 1,
       GMM_ASSERT1(false,
                   "No output format for a convex structure. To be done");
       );

    /*@GET ('display')
      displays a short summary for a @tcvstruct object.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       infomsg() << "gfCvStruct (convex structure) in dimension "
       << int(cs->dim()) << " with " << cs->nb_points() << "points. \n";
       );


  }

  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  bgeot::pconvex_structure cs = to_cvstruct_object(m_in.pop());
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);


  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
              it->second->arg_in_max, it->second->arg_out_min,
              it->second->arg_out_max);
    it->second->run(m_in, m_out, cs);
  }
  else bad_cmd(init_cmd);
}
