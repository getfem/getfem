// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2009 Julien Pommier.
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
#include <getfemint_mesh_fem.h>
#include <getfemint_gsparse.h>

using namespace getfemint;


static void set_fem(getfem::mesh_fem *mf, getfemint::mexargs_in& in)
{
  getfem::pfem               fem = in.pop().to_fem();

  /* check or build the convex list */
  dal::bit_vector bv;
  bool all_cv = false;
  if (in.remaining() == 1)
    bv = in.pop().to_bit_vector(&mf->linked_mesh().convex_index(), -1);
  else
    all_cv = true;

  /* check for the validity of the operation */
  for (dal::bv_visitor cv(bv); !cv.finished(); ++cv) {
    if (!mf->linked_mesh().convex_index().is_in(cv))
      THROW_ERROR("Convex " << cv+config::base_index()
		  << " was not found in mesh");
    if (fem->basic_structure(cv) != mf->linked_mesh().structure_of_convex(cv)->basic_structure())
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
static void set_classical_fem(getfem::mesh_fem *mf, getfemint::mexargs_in& in, bool discontinuous) {
  dim_type K = dim_type(in.pop().to_integer(0,255)); //, IM_DEGREE = dim_type(-1);
  dal::bit_vector bv;
  if (in.remaining() == 1) {
    bv = in.pop().to_bit_vector(&mf->linked_mesh().convex_index(), -1);
  } else {
    bv = mf->linked_mesh().convex_index();
  }
  if (!discontinuous) {
    mf->set_classical_finite_element(bv,K);
  } else {
    mf->set_classical_discontinuous_finite_element(bv,K);
  }
}

/*MLABCOM
  FUNCTION [x] = gf_mesh_fem_set(meshfem MF, operation [, args])

  General function for modifying mesh_fem objects.

  @SET MESHFEM:SET('fem')
  @SET MESHFEM:SET('classical fem')
  @SET MESHFEM:SET('classical discontinuous fem')
  @SET MESHFEM:SET('qdim')
  @SET MESHFEM:SET('reduction')
  @SET MESHFEM:SET('reduction matrices')
  @SET MESHFEM:SET('dof partition')

  $Id$
MLABCOM*/

void gf_mesh_fem_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }

  getfem::mesh_fem *mf = in.pop().to_mesh_fem();
  std::string cmd        = in.pop().to_string();
  if (check_cmd(cmd, "fem", in, out, 1, 2, 0, 0)) {
    /*@SET MESHFEM:SET('fem',@tfem f[, @ivec CVids])
    Set the Finite Element Method.

    Assign a FEM `f` to all convexes whose #ids are listed in `CVids`.
    If `CVids` is not given, the integration is assigned to all convexes.

    See the help of FEM:INIT to obtain a list of available FEM methods.@*/
    set_fem(mf, in);
  } else if (check_cmd(cmd, "classical fem", in, out, 1, 2, 0, 0)) {
    /*@SET MESHFEM:SET('classical fem',@int k[, @ivec CVids])
    Assign a classical (Lagrange polynomial) fem of order `k` to the @tmf.

    Uses FEM_PK for simplexes, FEM_QK for parallelepipeds etc.@*/
    set_classical_fem(mf, in, false);
  } else if (check_cmd(cmd, "classical discontinuous fem", in, out,
		       1, 2, 0, 0)) {
    /*@SET MESHFEM:SET('classical discontinuous fem',@int K, [@int IM_DEGREE [,@ivec CVIDX]])
    Assigns a classical (Lagrange polynomial) discontinuous fem or order K.

    Similar to MESHFEM:SET('classical fem') except that
    FEM_PK_DISCONTINUOUS is used.@*/
    set_classical_fem(mf, in, true);
  } else if (check_cmd(cmd, "qdim", in, out, 1, 1, 0, 0)) {
    /*@SET MESHFEM:SET('qdim',@int Q)
    Change the `Q` dimension of the field that is interpolated by the @tmf.

    `Q = 1` means that the @tmf describes a scalar field, `Q = N` means
    that the @tmf describes a vector field of dimension N.@*/
    size_type q_dim = in.pop().to_integer(1,255);
    mf->set_qdim(dim_type(q_dim));
  } else if (check_cmd(cmd, "reduction matrices", in, out, 2, 2, 0, 0)) {
    /*@SET MESHFEM:SET('reduction matrices',@mat R,@mat E)
    Set the reduction and extension matrices and valid their use.@*/
    dal::shared_ptr<gsparse> R = in.pop().to_sparse();
    dal::shared_ptr<gsparse> E = in.pop().to_sparse();
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
  } else if (check_cmd(cmd, "reduction", in, out, 1, 1, 0, 0)) {
    /*@SET MESHFEM:SET('reduction',@int s)
    Set or unset the use of the reduction/extension matrices.@*/
    size_type s = in.pop().to_integer(0,255);
    mf->set_reduction(s != size_type(0));
  } else if (check_cmd(cmd, "dof partition", in, out, 1, 1, 0, 0)) {
    /*@SET MESHFEM:SET('dof partition',@ivec DOFP)
    Change the 'dof_partition' array.

    `DOFP` is a vector holding a integer value for each convex of the @tmf.
    See MESHFEM:GET('dof partition') for a description of "dof partition".@*/
    iarray v =
      in.pop().to_iarray(int(mf->linked_mesh().convex_index().last_true()+1));
    for (unsigned i=0; i < v.size(); ++i)
      mf->set_dof_partition(i, v[i]);
  } else bad_cmd(cmd);
}
