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

#include <getfemint_misc.h>
#include <getfemint_gsparse.h>
#include <getfem/getfem_partial_mesh_fem.h>
#include <gmm/gmm_range_basis.h>
#include <getfem/getfem_mesh_fem_level_set.h>
#include <getfem/getfem_mesh_fem_product.h>
#include <getfem/getfem_fem.h>

using namespace getfemint;


static void set_fem(getfem::mesh_fem *mf, getfemint::mexargs_in& in)
{
  getfem::pfem fem = to_fem_object(in.pop());

  /* check or build the convex list */
  dal::bit_vector bv;
  bool all_cv = false;
  if (in.remaining() == 1)
    bv = in.pop().to_bit_vector(&mf->linked_mesh().convex_index(), -config::base_index());
  else
    all_cv = true;

  /* check for the validity of the operation */
  for (dal::bv_visitor cv(bv); !cv.finished(); ++cv) {
    if (!mf->linked_mesh().convex_index().is_in(cv))
      THROW_ERROR("Convex " << cv+config::base_index()
                  << " was not found in mesh");
    if (fem->basic_structure(cv) != bgeot::basic_structure(mf->linked_mesh().structure_of_convex(cv)))
      infomsg() << "Warning: structure of the FEM seems to be incompatible "
        "with the structure of the convex (if you are using high degree "
        "geom. transf. ignore this)\n";
  }

  /* all the work done here */
  if (!all_cv)
    mf->set_finite_element(bv, fem);
  else
    mf->set_finite_element(fem);
}

/* set the classical fem of order on the mesh_fem, with a classical integration
   method */
static void set_classical_fem(getfem::mesh_fem *mf, getfemint::mexargs_in& in,
                              bool discontinuous) {
  dim_type K = dim_type(in.pop().to_integer(0,255));

  bool complete(false);
  if (in.remaining() && in.front().is_string()) {
    std::string s = in.pop().to_string();
    if (cmd_strmatch(s, "complete"))
      complete = true;
    else
      { THROW_BADARG("Invalid option" << s); }
  }

  scalar_type alpha = 0.0;
  if (discontinuous && in.remaining()) alpha = in.pop().to_scalar();

  dal::bit_vector bv;
  if (in.remaining()) {
    bv = in.pop().to_bit_vector(&mf->linked_mesh().convex_index(),
                                -config::base_index());
    if (!discontinuous) {
      mf->set_classical_finite_element(bv, K, complete);
    } else {
      mf->set_classical_discontinuous_finite_element(bv, K, alpha, complete);
    }
  } else {
    if (!discontinuous) {
      mf->set_classical_finite_element(K, complete);
    } else {
      mf->set_classical_discontinuous_finite_element(K, alpha, complete);
    }
  }
}

/*@GFDOC
  General function for modifying mesh_fem objects.
  @*/





// Object for the declaration of a new sub-command.

struct sub_gf_mf_set : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
                   getfemint::mexargs_out& out,
                   getfem::mesh_fem *mf) = 0;
};

typedef std::shared_ptr<sub_gf_mf_set> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_mf_set {                                    \
      virtual void run(getfemint::mexargs_in& in,                           \
                       getfemint::mexargs_out& out,                         \
                       getfem::mesh_fem *mf)                                \
      { dummy_func(in); dummy_func(out); code }                             \
    };                                                                      \
    psub_command psubc = std::make_shared<subc>();                          \
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;             \
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;         \
    subc_tab[cmd_normalize(name)] = psubc;                                  \
  }



void gf_mesh_fem_set(getfemint::mexargs_in& m_in,
                     getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;
  
  if (subc_tab.size() == 0) {

    
    /*@SET ('fem', @tfem f[, @ivec CVids])
      Set the Finite Element Method.
      
      Assign an FEM `f` to all convexes whose #ids are listed in `CVids`.
      If `CVids` is not given, the integration is assigned to all convexes.
      
      See the help of FEM:INIT to obtain a list of available FEM methods.@*/
    sub_command
      ("fem", 1, 2, 0, 0,
       set_fem(mf, in);
       );


    /*@SET ('classical fem', @int k[[, 'complete'], @ivec CVids])
    Assign a classical (Lagrange polynomial) fem of order `k` to the @tmf.
    The option 'complete' requests complete Langrange polynomial elements,
    even if the element geometric transformation is an incomplete one
    (e.g. 8-node quadrilateral or 20-node hexahedral).

    Uses FEM_PK for simplexes, FEM_QK for parallelepipeds etc.@*/
    sub_command
      ("classical fem", 1, 3, 0, 0,
       set_classical_fem(mf, in, false);
       );


    /*@SET ('classical discontinuous fem', @int k[[, 'complete'], @tscalar alpha[, @ivec CVIDX]])
    Assigns a classical (Lagrange polynomial) discontinuous fem of order k.

    Similar to MESH_FEM:SET('set classical fem') except that
    FEM_PK_DISCONTINUOUS is used. Param `alpha` the node inset,
    :math:`0 \leq alpha < 1`, where 0 implies usual dof nodes, greater values
    move the nodes toward the center of gravity, and 1 means that all
    degrees of freedom collapse on the center of gravity.
    The option 'complete' requests complete Langrange polynomial elements,
    even if the element geometric transformation is an incomplete one
    (e.g. 8-node quadrilateral or 20-node hexahedral).@*/
    sub_command
      ("classical discontinuous fem", 1, 4, 0, 0,
       set_classical_fem(mf, in, true);
       );


    /*@SET ('qdim', @int Q)
      Change the `Q` dimension of the field that is interpolated by the @tmf.
      
      `Q = 1` means that the @tmf describes a scalar field, `Q = N` means
      that the @tmf describes a vector field of dimension N.@*/
    sub_command
      ("qdim", 1, 1, 0, 0,
       size_type q_dim = in.pop().to_integer(1,255);
       mf->set_qdim(dim_type(q_dim));
       );


    /*@SET ('reduction matrices', @mat R, @mat E)
      Set the reduction and extension matrices and valid their use.@*/
    sub_command
      ("reduction matrices", 2, 2, 0, 0,
       std::shared_ptr<gsparse> R = in.pop().to_sparse();
       std::shared_ptr<gsparse> E = in.pop().to_sparse();
       if (R->is_complex() || E->is_complex())
         THROW_BADARG("Reduction and extension matrices should be real matrices");
       if (R->storage()==gsparse::CSCMAT && E->storage()==gsparse::CSCMAT)
         mf->set_reduction_matrices(R->real_csc(), E->real_csc());
       else if (R->storage()==gsparse::CSCMAT && E->storage()==gsparse::WSCMAT)
         mf->set_reduction_matrices(R->real_csc(), E->real_wsc());
       else if (R->storage()==gsparse::WSCMAT && E->storage()==gsparse::CSCMAT)
         mf->set_reduction_matrices(R->real_wsc(), E->real_csc());
       else if (R->storage()==gsparse::WSCMAT && E->storage()==gsparse::WSCMAT)
         mf->set_reduction_matrices(R->real_wsc(), E->real_wsc());
       else
         THROW_BADARG("Reduction and extension matrices should be "
                      "sparse matrices");
       );


    /*@SET ('reduction', @int s)
      Set or unset the use of the reduction/extension matrices.@*/
    sub_command
      ("reduction", 1, 1, 0, 0,
       size_type s = in.pop().to_integer(0,255);
       mf->set_reduction(s != size_type(0));
       );

    /*@SET ('reduce meshfem', @mat RM)
      Set reduction mesh fem
      This function selects the degrees of freedom of the finite element
      method by selecting a set of independent vectors of the matrix RM.
      The numer of columns of RM should corresponds to the number of degrees
      of freedom of the finite element method.  @*/
    sub_command
      ("reduce meshfem", 1, 1, 0, 0,
       std::shared_ptr<gsparse>  RM = in.pop().to_sparse();
       std::set<size_type> cols;
       cols.clear();
       gmm::range_basis(RM->real_csc(), cols, 1e-12);
       mf->reduce_to_basic_dof(cols);
       );

    /*@SET ('dof partition', @ivec DOFP)
      Change the 'dof_partition' array.
      
      `DOFP` is a vector holding a integer value for each convex of the @tmf.
      See MESH_FEM:GET('dof partition') for a description of "dof partition".@*/
    sub_command
      ("dof partition", 1, 1, 0, 0,
       iarray v =
       in.pop().to_iarray(int(mf->linked_mesh().convex_index().last_true()+1));
       for (unsigned i=0; i < v.size(); ++i)
         mf->set_dof_partition(i, v[i]);
       );


    /*@SET ('set partial', @ivec DOFs[, @ivec RCVs])
      Can only be applied to a partial @tmf. Change the subset of the
      degrees of freedom of `mf`.

      If `RCVs` is given, no FEM will be put on the convexes listed
      in `RCVs`.@*/
    sub_command
      ("set partial", 1, 2, 0, 0,
       dal::bit_vector doflst = in.pop().to_bit_vector();
       dal::bit_vector rcvlst;
       if (in.remaining()) rcvlst = in.pop().to_bit_vector();
       
       getfem::partial_mesh_fem *ppmf
       = dynamic_cast<getfem::partial_mesh_fem *>(mf);
       if (!ppmf) THROW_BADARG("The command 'set partial' can only be "
                               "applied to a partial mesh_fem object");
       ppmf->adapt(doflst, rcvlst);
       );

    /*@SET ('adapt')
      For a @tmf levelset object only. Adapt the mesh_fem object to a
      change of the levelset function. @*/
    sub_command
      ("adapt", 0, 0, 0, 0,
       getfem::mesh_fem_level_set *mfls
       = dynamic_cast<getfem::mesh_fem_level_set *>(mf);
       if (!mfls) THROW_BADARG("The command 'adapt' can only be "
                               "applied to a mesh_fem_level_set object");
       mfls->adapt();
       );

    /*@SET ('set enriched dofs', @ivec DOFs)
      For a @tmf product object only. Set te enriched dofs and adapt the @tmf product.
      @*/
    sub_command
      ("set enriched dofs", 1, 1, 0, 0,
       getfem::mesh_fem_product *mfprod
       = dynamic_cast<getfem::mesh_fem_product *>(mf);
       if (!mfprod) THROW_BADARG("The command 'set enriched dofs' can only be "
                                 "applied to a mesh_fem_product object");
       dal::bit_vector doflst = in.pop().to_bit_vector();
       mfprod->set_enrichment(doflst);
       );

  }


  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  getfem::mesh_fem *mf   = to_meshfem_object(m_in.pop());
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
              it->second->arg_in_max, it->second->arg_out_min,
              it->second->arg_out_max);
    it->second->run(m_in, m_out, mf);
  }
  else bad_cmd(init_cmd);

}
