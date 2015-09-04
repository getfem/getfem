/*===========================================================================

 Copyright (C) 2006-2015 Yves Renard, Julien Pommier.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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
#include <getfemint_workspace.h>
#include <getfemint_mesh_fem.h>
#include <getfemint_mesh.h>
#include <getfem/getfem_mesh_fem_sum.h>
#include <getfem/getfem_mesh_fem_product.h>
#include <getfemint_levelset.h>
#include <getfemint_mesh_levelset.h>
#include <getfem/getfem_mesh_fem_level_set.h>
#include <getfem/getfem_partial_mesh_fem.h>
#include <getfem/getfem_mesh_fem_global_function.h>

using namespace getfemint;

/*@GFDOC
  This object represents a finite element method defined on a whole mesh.
@*/


// Object for the declaration of a new sub-command.

typedef getfemint_mesh_fem *pgetfemint_mesh_fem;

struct sub_gf_mf : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   getfemint_mesh *mm,
		   pgetfemint_mesh_fem &mmf,
		   unsigned q_dim) = 0;
};

typedef boost::intrusive_ptr<sub_gf_mf> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_mf {					\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       getfemint_mesh *mm,				\
		       pgetfemint_mesh_fem &mmf,			\
		       unsigned q_dim)					\
      { dummy_func(in); dummy_func(out); dummy_func(mm);		\
	dummy_func(q_dim); code }					\
    };									\
    psub_command psubc = new subc;					\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }

/*@INIT MF = ('.mesh', @tmesh m[, @int Qdim1=1[, @int Qdim2=1, ...]])
  Build a new @tmf object.
  
  The `Qdim` parameters specifies the dimension of the field represented
  by the finite element method. Qdim1 = 1 for a scalar field,
  Qdim1 = n for a vector field off size n, Qdim1=m, Qdim2=n for
  a matrix field of size mxn ...
  Returns the handle of the created object. @*/

void gf_mesh_fem(getfemint::mexargs_in& m_in,
		 getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {

    /*@INIT MF = ('load', @str fname[, @tmesh m])
      Load a @tmf from a file.

      If the mesh `m` is not supplied (this kind of file does not store the
      mesh), then it is read from the file `fname` and its descriptor is
      returned as the second output argument.@*/
    sub_command
      ("load", 1, 2, 0, 1,
       std::string fname = in.pop().to_string();
       if (in.remaining()) mm = in.pop().to_getfemint_mesh();
       else {
	 getfem::mesh *m = new getfem::mesh();
	 m->read_from_file(fname);
	 mm = getfemint_mesh::get_from(m);
       }
       mmf = getfemint_mesh_fem::new_from(mm,q_dim);
       mmf->mesh_fem().read_from_file(fname);
       );

    /*@INIT MF = ('from string', @str s[, @tmesh m])
      Create a @tmf object from its string description.

      See also ``MESH_FEM:GET('char')``@*/
    sub_command
      ("from string", 1, 2, 0, 1,
       std::stringstream ss(in.pop().to_string());
       if (in.remaining()) mm = in.pop().to_getfemint_mesh();
       else {
	 getfem::mesh *m = new getfem::mesh();
	 m->read_from_file(ss);
	 mm = getfemint_mesh::get_from(m);
       }
       mmf = getfemint_mesh_fem::new_from(mm,q_dim);
       mmf->mesh_fem().read_from_file(ss);
       );


    /*@INIT MF = ('clone', @tmf mf)
      Create a copy of a @tmf.@*/
    sub_command
      ("clone", 1, 1, 0, 1,
       getfemint_mesh_fem *mmf2 = in.pop().to_getfemint_mesh_fem();
       mm = object_to_mesh(workspace().object(mmf2->linked_mesh_id()));
       mmf = getfemint_mesh_fem::new_from(mm,q_dim);
       std::stringstream ss;
       mmf2->mesh_fem().write_to_file(ss);
       mmf->mesh_fem().read_from_file(ss);
       );

    
    /*@INIT MF = ('sum', @tmf mf1, @tmf mf2[, @tmf mf3[, ...]])
      Create a @tmf that spans two (or more) @tmf's.

      All @tmf must share the same mesh.

      After that, you should not modify the FEM of `mf1`, `mf2` etc.@*/
    sub_command
      ("sum", 1, -1, 0, 1,
       std::vector<const getfem::mesh_fem*> mftab;
       getfem::mesh_fem_sum *msum = 0;
       while (in.remaining()) {
	 getfemint_mesh_fem *gfimf = in.pop().to_getfemint_mesh_fem();
	 if (msum == 0) {
	   mm = object_to_mesh(workspace().object(gfimf->linked_mesh_id()));
	   msum = new getfem::mesh_fem_sum(gfimf->linked_mesh());
	
	   mmf = getfemint_mesh_fem::get_from(msum);
	 }
	 workspace().set_dependance(mmf, gfimf);
	 mftab.push_back(&gfimf->mesh_fem());
       }
       msum->set_mesh_fems(mftab);
       msum->adapt();
       );

    /*@INIT MF = ('product', @tmf mf1, @tmf mf2)
      Create a @tmf that spans all the product of a selection of shape
      functions of `mf1` by all shape functions of `mf2`.
      Designed for Xfem enrichment.

      `mf1` and `mf2` must share the same mesh.

      After that, you should not modify the FEM of `mf1`, `mf2`.@*/
    sub_command
      ("product", 2, 2, 0, 1,
       getfemint_mesh_fem *gfimf1 = in.pop().to_getfemint_mesh_fem();
       getfemint_mesh_fem *gfimf2 = in.pop().to_getfemint_mesh_fem();
       getfem::mesh_fem_product *mprod = 
         new getfem::mesh_fem_product(gfimf1->mesh_fem(), gfimf2->mesh_fem());
       mmf = getfemint_mesh_fem::get_from(mprod);
       workspace().set_dependance(mmf, gfimf1);
       workspace().set_dependance(mmf, gfimf2);
       );


    /*@INIT MF = ('levelset', @tmls mls, @tmf mf)
      Create a @tmf that is conformal to implicit surfaces defined in
      @tmls.@*/
    sub_command
      ("levelset", 2, 2, 0, 1,
       getfemint_mesh_levelset *gmls = in.pop().to_getfemint_mesh_levelset();
       getfemint_mesh_fem *gmf = in.pop().to_getfemint_mesh_fem();
       getfem::mesh_fem_level_set *mfls =
       new getfem::mesh_fem_level_set(gmls->mesh_levelset(),
				      gmf->mesh_fem());
       mmf = getfemint_mesh_fem::get_from(mfls);
       workspace().set_dependance(mmf, gmf);
       workspace().set_dependance(mmf, gmls);
       mfls->adapt();
       );


    /*@INIT MF = ('global function', @tmesh m, @tls ls, @CELL{@tgf GF1,...}[, @int Qdim_m])
      Create a @tmf whose base functions are global function given by the
      user in the system of coordinate defined by the iso-values of the two
      level-set function of `ls`. @*/
    sub_command
      ("global function", 3, 4, 0, 1,
       mm = in.pop().to_getfemint_mesh();
       getfemint_levelset *gls = in.pop().to_getfemint_levelset();
       mexargs_in *in_gf = new mexargs_in(1, &in.pop().arg, true);
       if (in.remaining() && in.front().is_integer())
	 q_dim = in.pop().to_integer(1,256);

       std::vector<getfem::pglobal_function> vfunc(size_type(in_gf->narg()));
       for (size_type i = 0; i < vfunc.size(); ++i) {
	 getfem::abstract_xy_function *s = in_gf->pop().to_global_function();
	 vfunc[i] = getfem::global_function_on_level_set(gls->levelset(), *s);
       }

       getfem::mesh_fem_global_function *mfgf = new getfem::mesh_fem_global_function(mm->mesh());
       mfgf->set_qdim(dim_type(q_dim));
       mfgf->set_functions(vfunc);

       mmf = getfemint_mesh_fem::get_from(mfgf);
       );


    /*@INIT MF = ('partial', @tmf mf, @ivec DOFs[, @ivec RCVs])
      Build a restricted @tmf by keeping only a subset of the degrees of
      freedom of `mf`.

      If `RCVs` is given, no FEM will be put on the convexes listed in
      `RCVs`.@*/
    sub_command
      ("partial", 2, 3, 0, 1,
       getfemint_mesh_fem *gmf = in.pop().to_getfemint_mesh_fem();
       dal::bit_vector doflst = in.pop().to_bit_vector();
       dal::bit_vector rcvlst;
       if (in.remaining()) rcvlst = in.pop().to_bit_vector();

       getfem::partial_mesh_fem *ppmf =
       new getfem::partial_mesh_fem(gmf->mesh_fem());
       ppmf->adapt(doflst, rcvlst);

       mmf = getfemint_mesh_fem::get_from(ppmf);
       workspace().set_dependance(mmf, gmf);
       );
  }


  if (m_in.narg() < 1) THROW_BADARG("Wrong number of input arguments");
  getfemint_mesh *mm = NULL;
  getfemint_mesh_fem *mmf = NULL;
  unsigned q_dim = 1;

  if (m_in.front().is_string()) {

    std::string init_cmd   = m_in.pop().to_string();
    std::string cmd        = cmd_normalize(init_cmd);


    SUBC_TAB::iterator it = subc_tab.find(cmd);
    if (it != subc_tab.end()) {
      check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
		it->second->arg_in_max, it->second->arg_out_min,
		it->second->arg_out_max);
      it->second->run(m_in, m_out, mm, mmf, q_dim);
    }
    else bad_cmd(init_cmd);

  } else if (check_cmd("MeshFem", "MeshFem", m_in, m_out, 1, 7, 0, 1)) {
    /* Documentation iof the commande moved in first position to appear first
       in the documentation. */
    mm = m_in.pop().to_getfemint_mesh();
    bgeot::multi_index mi;
    dim_type qdim = 1;
    while (m_in.remaining()) {
      dim_type q = dim_type(m_in.pop().to_integer(1,65536));
      mi.push_back(q);
      qdim = dim_type(qdim*q);
    }
    if (mi.size() == 0) mi.push_back(qdim);
    mmf = getfemint_mesh_fem::new_from(mm,qdim);
    mmf->mesh_fem().set_qdim(mi);

    workspace().set_dependance(mmf, mm);
  }
  m_out.pop().from_object_id(mmf->get_id(), MESHFEM_CLASS_ID);
}
