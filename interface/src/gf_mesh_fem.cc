/*===========================================================================

 Copyright (C) 2006-2020 Yves Renard, Julien Pommier.

 This file is a part of GetFEM

 GetFEM is free software;  you can  redistribute it  and/or modify it under
 the  terms  of the  GNU  Lesser General Public License as published by the
 Free Software Foundation;  either version 3  of  the License,  or (at your
 option) any  later  version  along with  the GCC Runtime Library Exception
 either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/

#include <getfem/getfem_mesh_fem_sum.h>
#include <getfem/getfem_mesh_fem_product.h>
#include <getfem/getfem_mesh_level_set.h>
#include <getfem/getfem_mesh_fem_level_set.h>
#include <getfem/getfem_partial_mesh_fem.h>
#include <getfem/getfem_mesh_fem_global_function.h>
#include <getfemint_misc.h>
#include <getfemint_workspace.h>
#include <getfemint_levelset.h>

using namespace getfemint;

/*@GFDOC
  This object represents a finite element method defined on a whole mesh.
@*/


// Object for the declaration of a new sub-command.

struct sub_gf_mf : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
                   getfemint::mexargs_out& out,
                   const getfem::mesh *mm,
                   std::shared_ptr<getfem::mesh_fem> &mmf,
                   unsigned q_dim) = 0;
};

typedef std::shared_ptr<sub_gf_mf> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_mf {                                        \
      virtual void run(getfemint::mexargs_in& in,                           \
                       getfemint::mexargs_out& out,                         \
                       const getfem::mesh *mm,                              \
                       std::shared_ptr<getfem::mesh_fem> &mmf,              \
                       unsigned q_dim)                                      \
      { dummy_func(in); dummy_func(out); dummy_func(mm);                    \
        dummy_func(q_dim); code }                                           \
    };                                                                      \
    psub_command psubc = std::make_shared<subc>();                          \
    psubc->arg_in_min = arginmin;   psubc->arg_in_max = arginmax;           \
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;         \
    subc_tab[cmd_normalize(name)] = psubc;                                  \
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
  static std::map<std::string, psub_command > subc_tab;

  if (subc_tab.empty()) {

    /*@INIT MF = ('load', @str fname[, @tmesh m])
      Load a @tmf from a file.

      If the mesh `m` is not supplied (this kind of file does not store the
      mesh), then it is read from the file `fname` and its descriptor is
      returned as the second output argument.@*/
    sub_command
      ("load", 1, 2, 0, 1,
       std::string fname = in.pop().to_string();
       if (in.remaining())  {
         mm = extract_mesh_object(in.pop());
         mmf = std::make_shared<getfem::mesh_fem>(*mm, q_dim);
       } else {
         auto m = std::make_shared<getfem::mesh>();
         m->read_from_file(fname);
         store_mesh_object(m);
         mm = m.get();
         mmf = std::make_shared<getfem::mesh_fem>(*mm, q_dim);
         workspace().add_hidden_object(store_meshfem_object(mmf), m);
       }
       mmf->read_from_file(fname);
       );

    /*@INIT MF = ('from string', @str s[, @tmesh m])
      Create a @tmf object from its string description.

      See also ``MESH_FEM:GET('char')``@*/
    sub_command
      ("from string", 1, 2, 0, 1,
       std::stringstream ss(in.pop().to_string());
       if (in.remaining()) {
         mm = extract_mesh_object(in.pop());
         mmf = std::make_shared<getfem::mesh_fem>(*mm, q_dim);

       } else {
         auto m = std::make_shared<getfem::mesh>();
         m->read_from_file(ss);
         store_mesh_object(m);
         mm = m.get();
         mmf = std::make_shared<getfem::mesh_fem>(*mm, q_dim);
         workspace().add_hidden_object(store_meshfem_object(mmf), m);
       }
       mmf->read_from_file(ss);
       );


    /*@INIT MF = ('clone', @tmf mf)
      Create a copy of a @tmf.@*/
    sub_command
      ("clone", 1, 1, 0, 1,
       getfem::mesh_fem *mmf_in = to_meshfem_object(in.pop());
       mm = &mmf_in->linked_mesh();
       getfem::mesh_fem_sum *mfsum
         = dynamic_cast<getfem::mesh_fem_sum *>(mmf_in);
       getfem::mesh_fem_product *mfprod
         = dynamic_cast<getfem::mesh_fem_product *>(mmf_in);
       getfem::mesh_fem_level_set *mfls
         = dynamic_cast<getfem::mesh_fem_level_set *>(mmf_in);
       getfem::partial_mesh_fem *mfpart
         = dynamic_cast<getfem::partial_mesh_fem *>(mmf_in);
       getfem::mesh_fem_global_function *mfglob
         = dynamic_cast<getfem::mesh_fem_global_function *>(mmf_in);
       if (mfsum)
         mmf = std::make_shared<getfem::mesh_fem_sum>(*mfsum);
       else if (mfprod)
         mmf = std::make_shared<getfem::mesh_fem_product>(*mfprod);
       else if (mfls) {
         std::shared_ptr<getfem::mesh_fem_level_set> mmfls
           = std::make_shared<getfem::mesh_fem_level_set>(mfls->linked_mesh_level_set(),
                                                          mfls->linked_mesh_fem());
         mmfls->adapt();
         mmf = mmfls;
       } else if (mfpart) {
         GMM_WARNING1("Cloning a partial_mesh_fem simply clones the underlying"
                      " adapted mesh_fem");
         mmf = std::make_shared<getfem::mesh_fem>(mfpart->linked_mesh_fem());
       } else if (mfglob)
         mmf = std::make_shared<getfem::mesh_fem_global_function>(*mfglob);
       else
         mmf = std::make_shared<getfem::mesh_fem>(*mmf_in);
       );


    /*@INIT MF = ('sum', @tmf mf1, @tmf mf2[, @tmf mf3[, ...]])
      Create a @tmf that spans two (or more) @tmf's.

      All @tmf must share the same mesh.

      After that, you should not modify the FEM of `mf1`, `mf2` etc.@*/
    sub_command
      ("sum", 1, -1, 0, 1,
       std::vector<const getfem::mesh_fem*> mftab;
       std::shared_ptr<getfem::mesh_fem_sum> mfsum;
       while (in.remaining()) {
         getfem::mesh_fem *gfimf = to_meshfem_object(in.pop());
         if (mmf.get() == 0) {
           mm = &gfimf->linked_mesh();
           mfsum = std::make_shared<getfem::mesh_fem_sum>(*mm);
           mmf = mfsum;
           store_meshfem_object(mmf);
         }
         workspace().set_dependence(mmf.get(), gfimf);
         mftab.push_back(gfimf);
       }
       mfsum->set_mesh_fems(mftab);
       mfsum->adapt();
       mmf = mfsum;
       );

    /*@INIT MF = ('product', @tmf mf1, @tmf mf2)
      Create a @tmf that spans all the product of a selection of shape
      functions of `mf1` by all shape functions of `mf2`.
      Designed for Xfem enrichment.

      `mf1` and `mf2` must share the same mesh.

      After that, you should not modify the FEM of `mf1`, `mf2`.@*/
    sub_command
      ("product", 2, 2, 0, 1,
       getfem::mesh_fem *gfimf1 = to_meshfem_object(in.pop());
       getfem::mesh_fem *gfimf2 = to_meshfem_object(in.pop());
       mmf = std::make_shared<getfem::mesh_fem_product>(*gfimf1, *gfimf2);
       store_meshfem_object(mmf);
       workspace().set_dependence(mmf.get(), gfimf1);
       workspace().set_dependence(mmf.get(), gfimf2);
       );


    /*@INIT MF = ('levelset', @tmls mls, @tmf mf)
      Create a @tmf that is conformal to implicit surfaces defined in
      @tmls.@*/
    sub_command
      ("levelset", 2, 2, 0, 1,
       getfem::mesh_level_set &mls = *(to_mesh_levelset_object(in.pop()));
       getfem::mesh_fem *gmf = to_meshfem_object(in.pop());
       auto mfls = std::make_shared<getfem::mesh_fem_level_set>(mls, *gmf);
       mfls->adapt();
       mmf = mfls;
       store_meshfem_object(mmf);
       workspace().set_dependence(mmf.get(), gmf);
       workspace().set_dependence(mmf.get(), &mls);
       );


    /*@INIT MF = ('global function', @tmesh m, @tls ls, @CELL{@tgf GF1,...}[, @int Qdim_m])
      Create a @tmf whose base functions are global function given by the
      user in the system of coordinate defined by the iso-values of the two
      level-set function of `ls`. @*/
    sub_command
      ("global function", 3, 4, 0, 1,
       mm = extract_mesh_object(in.pop());
       auto pls = to_levelset_object(in.pop());
       mexargs_in in_gf(1, &in.pop().arg, true);
       if (in.remaining() && in.front().is_integer())
         q_dim = in.pop().to_integer(1,256);

       std::vector<getfem::pglobal_function> vfuncs(size_type(in_gf.narg()));
       for (auto &vfunc : vfuncs) {
         getfem::pxy_function s = to_global_function_object(in_gf.pop());
         vfunc = getfem::global_function_on_level_set(*pls, s);
       }

       auto mfgf = std::make_shared<getfem::mesh_fem_global_function>(*mm);
       mfgf->set_qdim(dim_type(q_dim));
       mfgf->set_functions(vfuncs);
       mmf = mfgf;
       );


    /*@INIT MF = ('bspline_uniform', @tmesh m, @int NX[, @int NY[, @int NZ]], @int order[, @str bcX_low[, @str bcY_low[, @str bcZ_low]][, @str bcX_high[, @str bcY_high[, @str bcZ_high]]]])
      Create a @tmf on mesh `m`, whose base functions are global functions
      corresponding to bspline basis of order `order`, in an NX x NY x NZ
      grid (just NX in 1D or NX x NY in 2D) that spans the entire bounding
      box of `m`.
      Optionally boundary conditions at the edges of the domain can be
      defined with `bcX_low`, `bcY_low`, `bcZ_low`, `bcX_high`, `bcY_high`,
      and `bcZ_high` set to 'free' (default) or 'periodic' or 'symmetry'. @*/
    sub_command
      ("bspline_uniform", 3, 11, 0, 1,
       mm = extract_mesh_object(in.pop());
       dim_type dim = mm->dim();
       if (dim > 3)
         THROW_ERROR("Uniform bspline only supported for dim = 1,2,3");
       size_type NX = in.pop().to_integer(1,1000);
       size_type NY = (dim >= 2) ? in.pop().to_integer(1,1000) : 0;
       size_type NZ = (dim == 3) ? in.pop().to_integer(1,1000) : 0;
       if (!in.remaining() || !in.front().is_integer())
         THROW_ERROR("One integer was expected for bspline order");
       size_type order = in.pop().to_integer(3,5);
       std::string bcx_low("free");
       std::string bcy_low("free");
       std::string bcz_low("free");
       std::string bcx_high("");
       std::string bcy_high("");
       std::string bcz_high("");
       if (in.remaining())             bcx_low = in.pop().to_string();
       if (dim >= 2 && in.remaining()) bcy_low = in.pop().to_string();
       if (dim == 3 && in.remaining()) bcz_low = in.pop().to_string();
       if (in.remaining())             bcx_high = in.pop().to_string();
       if (dim >= 2 && in.remaining()) bcy_high = in.pop().to_string();
       if (dim == 3 && in.remaining()) bcz_high = in.pop().to_string();
       if (in.remaining())
         THROW_ERROR("Too many arguments for bspline mesh_fem");
       getfem::bspline_boundary bcX_low(getfem::bspline_boundary::FREE);
       getfem::bspline_boundary bcY_low(getfem::bspline_boundary::FREE);
       getfem::bspline_boundary bcZ_low(getfem::bspline_boundary::FREE);
       getfem::bspline_boundary bcX_high(getfem::bspline_boundary::FREE);
       getfem::bspline_boundary bcY_high(getfem::bspline_boundary::FREE);
       getfem::bspline_boundary bcZ_high(getfem::bspline_boundary::FREE);
       if (bcx_low == "periodic")
         bcX_high = bcX_low = getfem::bspline_boundary::PERIODIC;
       else if (bcx_low == "symmetry")
         bcX_high = bcX_low = getfem::bspline_boundary::SYMMETRY;
       else if (bcx_low != "free")
         THROW_ERROR("Unknown boundary condition " << bcx_low);

       if (bcy_low == "periodic")
         bcY_high = bcY_low = getfem::bspline_boundary::PERIODIC;
       else if (bcy_low == "symmetry")
         bcY_high = bcY_low = getfem::bspline_boundary::SYMMETRY;
       else if (bcy_low != "free")
         THROW_ERROR("Unknown boundary condition " << bcy_low);

       if (bcz_low == "periodic")
         bcZ_high = bcZ_low = getfem::bspline_boundary::PERIODIC;
       else if (bcz_low == "symmetry")
         bcZ_high = bcZ_low = getfem::bspline_boundary::SYMMETRY;
       else if (bcz_low != "free")
         THROW_ERROR("Unknown boundary condition " << bcz_low);

       if (!bcx_high.empty()) {
         if (bcx_high == "periodic")
           bcX_high = getfem::bspline_boundary::PERIODIC;
         else if (bcx_high == "symmetry")
           bcX_high = getfem::bspline_boundary::SYMMETRY;
         else if (bcx_high == "free")
           bcX_high = getfem::bspline_boundary::FREE;
         else
           THROW_ERROR("Unknown boundary condition " << bcx_high);
       }

       if (!bcy_high.empty()) {
         if (bcy_high == "periodic")
           bcY_high = getfem::bspline_boundary::PERIODIC;
         else if (bcy_high == "symmetry")
           bcY_high = getfem::bspline_boundary::SYMMETRY;
         else if (bcy_high == "free")
           bcY_high = getfem::bspline_boundary::FREE;
         else
           THROW_ERROR("Unknown boundary condition " << bcy_high);
       }

       if (!bcz_high.empty()) {
         if (bcz_high == "periodic")
           bcZ_high = getfem::bspline_boundary::PERIODIC;
         else if (bcz_high == "symmetry")
           bcZ_high = getfem::bspline_boundary::SYMMETRY;
         else if (bcz_high == "free")
           bcZ_high = getfem::bspline_boundary::FREE;
         else
           THROW_ERROR("Unknown boundary condition " << bcz_high);
       }

       auto mfgf = std::make_shared<getfem::mesh_fem_global_function>(*mm);
       mfgf->set_qdim(1.);
       if (dim == 1)
         define_uniform_bspline_basis_functions_for_mesh_fem
           (*mfgf, NX, order, bcX_low, bcX_high);
       else if (dim == 2)
         define_uniform_bspline_basis_functions_for_mesh_fem
           (*mfgf, NX, NY, order, bcX_low, bcY_low, bcX_high, bcY_high);
       else
         define_uniform_bspline_basis_functions_for_mesh_fem
           (*mfgf, NX, NY, NZ, order, bcX_low, bcY_low, bcZ_low,
                                      bcX_high, bcY_high, bcZ_high);
       mmf = mfgf;
       );


    /*@INIT MF = ('partial', @tmf mf, @ivec DOFs[, @ivec RCVs])
      Build a restricted @tmf by keeping only a subset of the degrees of
      freedom of `mf`.

      If `RCVs` is given, no FEM will be put on the convexes listed in
      `RCVs`.@*/
    sub_command
      ("partial", 2, 3, 0, 1,
       getfem::mesh_fem *gmf = to_meshfem_object(in.pop());
       dal::bit_vector doflst = in.pop().to_bit_vector();
       dal::bit_vector rcvlst;
       if (in.remaining()) rcvlst = in.pop().to_bit_vector();

       auto ppmf = std::make_shared<getfem::partial_mesh_fem>(*gmf);
       ppmf->adapt(doflst, rcvlst);
       mmf = ppmf;
       store_meshfem_object(mmf);
       workspace().set_dependence(mmf.get(), gmf);
       );
  }


  if (m_in.narg() < 1) THROW_BADARG("Wrong number of input arguments");
  const getfem::mesh *mm = NULL;
  std::shared_ptr<getfem::mesh_fem> mmf;
  unsigned q_dim = 1;

  if (m_in.front().is_string()) {

    std::string init_cmd   = m_in.pop().to_string();
    std::string cmd        = cmd_normalize(init_cmd);

    auto it = subc_tab.find(cmd);
    if (it != subc_tab.end()) {
      check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
                it->second->arg_in_max, it->second->arg_out_min,
                it->second->arg_out_max);
      it->second->run(m_in, m_out, mm, mmf, q_dim);
    }
    else bad_cmd(init_cmd);

  } else if (check_cmd("MeshFem", "MeshFem", m_in, m_out, 1, 7, 0, 1)) {
    /* Documentation of the command moved to the beginning so that it appears
       first in the documentation. */
    mm = extract_mesh_object(m_in.pop());
    bgeot::multi_index mi;
    dim_type qdim = 1;
    while (m_in.remaining()) {
      dim_type q = dim_type(m_in.pop().to_integer(1,65536));
      mi.push_back(q);
      qdim = dim_type(qdim*q);
    }
    if (mi.size() == 0) mi.push_back(qdim);
    mmf = std::make_shared<getfem::mesh_fem>(*mm, qdim);
    mmf->set_qdim(mi);
     store_meshfem_object(mmf);
     workspace().set_dependence(mmf.get(), mm);
  }

  id_type id = store_meshfem_object(mmf);
  m_out.pop().from_object_id(id, MESHFEM_CLASS_ID);
}
