/*===========================================================================

 Copyright (C) 2014-2020 Konstantinos Poulios.

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
#include <getfemint_workspace.h>
#include <getfem/getfem_im_data.h>

/*
  $Id: gf_mesh_im_data_get.cc 4705 2014-07-09 07:31:24Z logari81 $
 */


using namespace getfemint;


/*@GFDOC
  General function extracting information from mesh_im_data objects.
  @*/


// Object for the declaration of a new sub-command.

struct sub_gf_mimd_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
                   getfemint::mexargs_out& out,
                   getfem::im_data *mimd) = 0;
};

typedef std::shared_ptr<sub_gf_mimd_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_mimd_get {				\
      virtual void run(getfemint::mexargs_in& in,			\
                       getfemint::mexargs_out& out,			\
                       getfem::im_data *mimd)				\
      { dummy_func(in); dummy_func(out);  dummy_func(mimd); code }	\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }





void gf_mesh_im_data_get(getfemint::mexargs_in& m_in,
                         getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {

    /*@GET ('region')
      Output the region that the @tmimd is restricted to.
    @*/
    sub_command
      ("region", 0, 0, 0, 1,
       out.pop().from_integer(int(mimd->filtered_region()));
       );

    /*@GET ('nbpts')
      Output the number of integration points (filtered in the considered region).
    @*/
    sub_command
      ("nbpts", 0, 0, 0, 1,
       out.pop().from_integer(int(mimd->nb_filtered_index()));
       );

    /*@GET ('nb tensor elements')
      Output the size of the stored data (per integration point).
    @*/
    sub_command
      ("nb tensor elements", 0, 0, 0, 1,
       if (mimd->tensor_size().size()) {
         out.pop().from_integer(int(mimd->nb_tensor_elem()));
       }
       );

    /*@GET ('tensor size')
      Output the dimensions of the stored data (per integration point).
    @*/
    sub_command
      ("tensor size", 0, 0, 0, 1,
       if (mimd->tensor_size().size()) {
         iarray oidx = out.pop().create_iarray_h(unsigned(mimd->tensor_size().size()));
         std::copy(mimd->tensor_size().begin(),
                   mimd->tensor_size().end(), &oidx[0]);
       }
       );

    /*@GET ('display')
      displays a short summary for a @tmimd object.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       infomsg() << "gfMeshImData object containing data of size "
       << mimd->tensor_size()
       << " on a mesh in dimension "
       << int(mimd->linked_mesh_im().linked_mesh().dim())
       << " with " << mimd->linked_mesh_im().linked_mesh().nb_points()
       << " points and "
       << mimd->linked_mesh_im().linked_mesh().convex_index().card()
       << " elements\n";
       );


    /*@GET m = ('linked mesh')
    Returns a reference to the @tmesh object linked to `mim`.@*/
    sub_command
      ("linked mesh", 0, 0, 0, 1,
       id_type id = workspace().object(&mimd->linked_mesh_im().linked_mesh());
       if (id == id_type(-1)) {
	 auto pst = workspace().hidden_object
	   (workspace().object(&mimd->linked_mesh_im()),
	    &mimd->linked_mesh_im().linked_mesh());
	 if (!pst.get()) THROW_INTERNAL_ERROR;
	 std::shared_ptr<getfem::mesh> pm = 
	   std::const_pointer_cast<getfem::mesh>
	   (std::dynamic_pointer_cast<const getfem::mesh>(pst));
	 id = store_mesh_object(pm);
       }
       out.pop().from_object_id(id, MESH_CLASS_ID);
       );

  }


  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  getfem::im_data *mimd = to_meshimdata_object(m_in.pop());
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
              it->second->arg_in_max, it->second->arg_out_min,
              it->second->arg_out_max);
    it->second->run(m_in, m_out, mimd);
  }
  else bad_cmd(init_cmd);

}
