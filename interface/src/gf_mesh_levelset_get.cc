/*===========================================================================

 Copyright (C) 2006-2020 Julien Pommier.

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
// $Id$
#include <getfem/getfem_mesh_level_set.h>
#include <getfemint.h>
#include <getfemint_workspace.h>
#include <getfemint_levelset.h>

using namespace getfemint;

/*@GFDOC
  General function for querying information about @tmesh_levelset objects.
@*/



// Object for the declaration of a new sub-command.

struct sub_gf_lset_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   getfem::mesh_level_set &mls) = 0;
};

typedef std::shared_ptr<sub_gf_lset_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_lset_get {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       getfem::mesh_level_set &mls)			\
      { dummy_func(in); dummy_func(out); dummy_func(mls); code }	\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }                           



void gf_mesh_levelset_get(getfemint::mexargs_in& m_in,
                          getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {


    /*@GET M = ('cut_mesh')
      Return a @tmesh cut by the linked @tls's.@*/
    sub_command
      ("cut_mesh", 0, 0, 0, 1,
       auto mm = std::make_shared<getfem::mesh>();
       mls.global_cut_mesh(*mm);
       out.pop().from_object_id(store_mesh_object(mm), MESH_CLASS_ID);
       );


    /*@GET LM = ('linked_mesh')
      Return a reference to the linked @tmesh.@*/
    sub_command
      ("linked_mesh", 0, 0, 0, 1,
       id_type id = workspace().object(&mls.linked_mesh());
       if (id == id_type(-1)) THROW_INTERNAL_ERROR;
       out.pop().from_object_id(id, MESH_CLASS_ID);
       );


    /*@GET nbls = ('nb_ls')
      Return the number of linked @tls's.@*/
    sub_command
      ("nb_ls", 0, 0, 0, 1,
       out.pop().from_integer(int(mls.nb_level_sets()));
       );


    /*@GET LS = ('levelsets')
      Return a list of references to the linked @tls's.@*/
    sub_command
      ("levelsets", 0, 0, 0, 1,
       std::vector<id_type> ids;
       for (unsigned i=0; i < mls.nb_level_sets(); ++i) {
	 id_type id = workspace().object((const void *)(mls.get_level_set(i)));
	 GMM_ASSERT1(id != id_type(-1), "Unknown levelset !");
	 ids.push_back(id);
       }
       out.pop().from_object_id(ids, LEVELSET_CLASS_ID);
       );


    /*@GET CVIDs = ('crack_tip_convexes')
      Return the list of convex #id's of the linked @tmesh on
      which have a tip of any linked @tls's.@*/
    sub_command
      ("crack_tip_convexes", 0, 0, 0, 1,
       out.pop().from_bit_vector(mls.crack_tip_convexes());
       );


    /*@GET SIZE = ('memsize')
      Return the amount of memory (in bytes) used by the @tmls.@*/
    sub_command
      ("memsize", 0, 0, 0, 1,
       out.pop().from_integer(int(mls.memsize()));
       );


    /*@GET s = ('char')
      Output a (unique) string representation of the @tmlsn.

      This can be used to perform comparisons between two
      different @tmls objects.
      This function is to be completed.
      @*/
    sub_command
      ("char", 0, 0, 0, 1,
       GMM_ASSERT1(false, "Sorry, function to be done");
       // std::string s = ...;
       // out.pop().from_string(s.c_str());
       );


    /*@GET ('display')
      displays a short summary for a @tmls object.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       infomsg() << "gfMeshLevelSet object in dimension "
       << int(mls.linked_mesh().dim())
       << " with " << mls.linked_mesh().nb_points() << " points and "
       << mls.linked_mesh().convex_index().card() << " elements\n";
       );


  }



  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  getfem::mesh_level_set &mls = *(to_mesh_levelset_object(m_in.pop()));

  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out, mls);
  }
  else bad_cmd(init_cmd);

}
