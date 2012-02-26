/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2001-2012 Y. Renard, J. Pommier.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/

/**\file gf_mdbrick_set.cc
   \brief getfemint_mdbrick setter.
*/

#include <getfemint.h>
#include <getfemint_mdbrick.h>
#include <getfemint_mesh_fem.h>
#include <getfemint_workspace.h>
#include <getfemint_gsparse.h>

using namespace getfemint;


bool is_constraints_brick(getfemint_mdbrick *b) {
  if (!b->is_complex())
    return b->cast0<getfem::mdbrick_constraint<real_model_state> >() != 0;
  else return b->cast0<getfem::mdbrick_constraint<cplx_model_state> >() != 0;
}

getfem::mdbrick_constraint<real_model_state> *
to_constraints_brick(getfemint_mdbrick *b, scalar_type) {
  return b->cast<getfem::mdbrick_constraint<real_model_state> >();
}

getfem::mdbrick_constraint<cplx_model_state> *
to_constraints_brick(getfemint_mdbrick *b, complex_type) {
  return b->cast<getfem::mdbrick_constraint<cplx_model_state> >();
}


/*@GFDOC
  Modify a model brick object.
@*/
void gf_mdbrick_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfemint_mdbrick *b   = in.pop().to_getfemint_mdbrick(true);
  std::string cmd        = in.pop().to_string();
  if (check_cmd(cmd, "param", in, out, 2, 3, 0, 0)) {
    /*@SET ('param', @str name, {@tmf mf,V | V})
    Change the value of a brick parameter.

    `name` is the name of the parameter. `V` should contain the
    new parameter value (vector or float). If a @tmf is given,
    `V` should hold the field values over that @tmf (i.e. its
    last dimension should be MESH_FEM:GET('nbdof') or 1 for
    constant field).@*/
    std::string pname = in.pop().to_string();
    for (unsigned i=0; i < pname.size(); ++i)
      if (pname[i] == ' ') pname[i] = '_';

    getfem::mdbrick_abstract_parameter *p = b->param(pname);
    if (!p) THROW_BADARG("wrong parameter name for this brick: " << pname);

    real_mdbrick_parameter *rp = dynamic_cast<real_mdbrick_parameter*>(p);
    cplx_mdbrick_parameter *cp = dynamic_cast<cplx_mdbrick_parameter*>(p);

    const getfem::mesh_fem *mf = &p->mf();
    array_dimensions d, d_expected1, d_expected2;

    darray rw;
    carray cw;

    if (in.front().is_mesh_fem()) {
      getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
      mf = &gfi_mf->mesh_fem();
      if (mf->get_qdim() != 1)
	THROW_BADARG("cannot use a mesh_fem whose qdim is not 1 "
		     "for a brick parameter");
      workspace().set_dependance(b, gfi_mf);
    }

    d_expected1.assign(p->fsizes());
    d_expected2.assign(p->fsizes());
    d_expected2.push_back(unsigned(mf->nb_dof()));

    if (cp) {
      cw = in.pop().to_carray(); d = cw;
    } else {
      rw = in.pop().to_darray(); d = rw;
    }

    bool sz_ok = (d.ndim() >= d_expected1.ndim());
    for (unsigned i=0; sz_ok && i < d.ndim(); ++i) {
      if (d.dim(i) != d_expected1.dim(i)) sz_ok = false;
    }
    if (!sz_ok) {
      sz_ok = (d.ndim() >= d_expected2.ndim());
      for (unsigned i=0; sz_ok && i < d.ndim(); ++i) {
	if (d.dim(i) != d_expected2.dim(i)) sz_ok = false;
      }
    }
    if (!sz_ok && !config::has_1D_arrays() &&
	d_expected2.ndim() == 1 &&
	d.dim(0)==1 && d.dim(1) == d_expected2.dim(0))
      sz_ok = true;

    /* for python, a 1d of size n is ok if we were expecting a 1 x n matrix */
    if (!sz_ok && config::has_1D_arrays() &&
        d_expected2.ndim() == 2 && d_expected2.dim(0) == 1 &&
        d.ndim() == 1 && d.dim(0) == d_expected2.dim(1))
      sz_ok = true;
    /* for python, a float is ok if we were expecting a constant field*/
    if (!sz_ok && config::has_1D_arrays() &&
        d.ndim() == 0)
      sz_ok = true;

    if (!sz_ok) THROW_BADARG("wrong size for the parameter " << pname
			     << ", expected an array of size "
			     << d_expected2 << " ( or " << d_expected1 <<
			     " for a constant field), got an array of size "
			     << d);

    if (cp) cp->set(*mf, cw);
    else    rp->set(*mf, rw);
  } else if (check_cmd(cmd, "penalization_epsilon", in, out, 1, 1, 0, 0)) {
    /*@SET ('penalization_epsilon', @scalar eps)
    Change the penalization coefficient of a constraint brick.

    This is only applicable to the bricks which inherit from the
    constraint brick, such as the Dirichlet ones. And of course it
    is not effective when the constraint is enforced via direct
    elimination or via Lagrange multipliers. The default value of
    `eps` is 1e-9.@*/
    if (!is_constraints_brick(b))
      THROW_BADARG("this is only applicable to constraint bricks, and it sub-classes");
    scalar_type eps = in.pop().to_scalar();
    if (eps <= 0)
      THROW_BADARG("Penalization parameter  should be small, but positive!");
    if (!b->is_complex())
      to_constraints_brick(b, scalar_type())->set_penalization_parameter(eps);
    else
      to_constraints_brick(b, complex_type())->set_penalization_parameter(eps);
  } else if (check_cmd(cmd, "constraints", in, out, 2, 2, 0, 0)) {
    /*@SET ('constraints', @mat H, @vec R)
    Set the constraints imposed by a constraint brick.

    This is only applicable to the bricks which inherit from the
    constraint brick, such as the Dirichlet ones. Imposes `H.U=R`.@*/
    if (!is_constraints_brick(b))
      THROW_BADARG("this is only applicable to constraint bricks, and it sub-classes");
    dal::shared_ptr<gsparse> H;
    if (in.front().is_sparse())
      H = in.pop().to_sparse();
    else {
      darray dH = in.pop().to_darray();
      if (dH.ndim() != 2)
	THROW_BADARG("Constraint matrix should be a matrix!");
      H.reset(new gsparse(dH.getm(), dH.getn(),
			  gsparse::WSCMAT, gsparse::REAL));
      for (unsigned i=0; i < dH.getm(); ++i)
	for (unsigned j=0; j < dH.getn(); ++j)
	  if (dH(i,j))
	    H->real_wsc()(i,j) = dH(i,j);
    }
    if (b->is_complex() != H->is_complex()) {
      THROW_BADARG("incompatible constraint matrix (mixing complexes with reals");
    }
    H->to_csc();
    if (!b->is_complex()) {
      gf_real_sparse_csc_const_ref C(H->real_csc()); //; in.pop().to_sparse(C);
      darray rhs = in.pop().to_darray(int(gmm::mat_nrows(C)));
      to_constraints_brick(b, scalar_type())->set_constraints(C, rhs);
    } else {
      gf_cplx_sparse_csc_const_ref C(H->cplx_csc()); //; in.pop().to_sparse(C);
      carray rhs = in.pop().to_carray(int(gmm::mat_nrows(C)));
      to_constraints_brick(b, complex_type())->set_constraints(C, rhs);
    }
  } else if (check_cmd(cmd, "constraints_rhs", in, out, 1, 1, 0, 0)) {
    /*@SET ('constraints_rhs', @mat H, @vec R)
    Set the right hand side of the constraints imposed by a constraint brick.

    This is only applicable to the bricks which inherit from the
    constraint brick, such as the Dirichlet ones.@*/
    if (!is_constraints_brick(b))
      THROW_BADARG("this is only applicable to constraint bricks, and its derivatives");
    if (!b->is_complex())
      to_constraints_brick(b, scalar_type())->set_constraints_rhs(in.pop().to_darray());
    else
      to_constraints_brick(b, complex_type())->set_constraints_rhs(in.pop().to_carray());
  } else bad_cmd(cmd);
}
