// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2005-2008 Julien Pommier.
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

/**\file gf_mdbrick.cc
   \brief mdbrick construction wrapper.
*/
#include <getfemint.h>
#include <getfemint_workspace.h>
#include <getfemint_misc.h>
#include <getfemint_mdbrick.h>
#include <getfemint_mesh_im.h>
#include <getfemint_mesh_fem.h>
#include <getfem/getfem_nonlinear_elasticity.h>
#include <getfem/getfem_plasticity.h>
#include <getfem/getfem_fourth_order.h>
#include <getfem/getfem_Navier_Stokes.h>
#include <getfem/getfem_linearized_plates.h>

using namespace getfemint;


getfem::mesh_im &
pop_mesh_im(mexargs_in &in, getfemint_mdbrick *b) {
  getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
  workspace().set_dependance(b, gfi_mim);
  return gfi_mim->mesh_im();
}


getfem::mesh_fem &
pop_mesh_fem(mexargs_in &in, getfemint_mdbrick *b) {
  getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
  workspace().set_dependance(b, gfi_mf);
  return gfi_mf->mesh_fem();
}

getfemint_mdbrick &
pop_mdbrick(mexargs_in &in, getfemint_mdbrick *b) {
  getfemint_mdbrick *gfi_mdb = in.pop().to_getfemint_mdbrick();
  workspace().set_dependance(b, gfi_mdb);
  return *gfi_mdb;
}

bool get_complexity(mexargs_in &in, bool default_v = false) {
  if (in.remaining() && in.front().is_string()) {
    std::string s = in.front().to_string();
    if (cmd_strmatch(s, "complex")) { in.pop(); return true; }
    else if (cmd_strmatch(s, "real")) { in.pop(); return false; }
  }
  return default_v;
}

size_type get_num_fem(mexargs_in &in, getfemint_mdbrick &b) {
  size_type num_fem = 0; if (in.remaining()) num_fem = in.pop().to_integer();
  if (num_fem >= b.mdbrick().nb_mesh_fems())
    THROW_BADARG("wrong mesh_fem number :" << num_fem);
  return num_fem;
}

getfem::constraints_type get_constraints_type(mexargs_in &in) {
  if (!in.remaining())
    THROW_BADARG("missing argument: expected a constraints policy: "
		 "'augmented', 'penalized' or 'eliminated'");
  std::string dtype = in.pop().to_string();
  if (cmd_strmatch(dtype, "augmented")) {
    return getfem::AUGMENTED_CONSTRAINTS;
  } else if (cmd_strmatch(dtype, "penalized")) {
    return getfem::PENALIZED_CONSTRAINTS;
  } else if (cmd_strmatch(dtype, "eliminated")) {
    return getfem::ELIMINATED_CONSTRAINTS;
  } else 
    THROW_BADARG("expected a constraints policy: 'augmented', 'penalized' or 'eliminated'");
}

#define SET_BRICK_R(class_, name_, args_)				\
  b->set_brick(new getfem::class_<real_model_state> args_, name_);	\
  
#define SET_BRICK_C(class_, name_, args_)				\
  b->set_brick(new getfem::class_<cplx_model_state> args_, name_);	\
  
#define SET_BRICK2(class_, name_, args1_, args2_)			\
  if (!is_complex)							\
    { SET_BRICK_R(class_, name_, args1_); }				\
  else									\
    { SET_BRICK_C(class_, name_, args2_); }				\
  
#define SET_BRICK(class_, name_, args_)	SET_BRICK2(class_,name_,args_,args_)

/*MLABCOM

  FUNCTION M = gf_mdbrick(brick_name, [, args])
  General constructor for mdbrick object. Returns a getfem handle to the newly
  created object.

  Many of the bricks take a "numfem" optional parameter, which
  is the meshfem number in the stack of parent bricks (by default
  numfem=0, i.e. it refers to the first meshfem in the stack of
  bricks).
  
  @INIT MDBRICK:INIT ('constraint')
  @INIT MDBRICK:INIT ('dirichlet')
  @INIT MDBRICK:INIT ('dirichlet on normal component')
  @INIT MDBRICK:INIT ('dirichlet on normal derivative')
  @INIT MDBRICK:INIT ('generalized dirichlet')
  @INIT MDBRICK:INIT ('source term')
  @INIT MDBRICK:INIT ('normal source term')
  @INIT MDBRICK:INIT ('normal derivative source term')
  @INIT MDBRICK:INIT ('neumann KirchhoffLove source term')
  @INIT MDBRICK:INIT ('qu term')
  @INIT MDBRICK:INIT ('mass matrix')
  @INIT MDBRICK:INIT ('generic elliptic')
  @INIT MDBRICK:INIT ('helmholtz')
  @INIT MDBRICK:INIT ('isotropic linearized elasticity')
  @INIT MDBRICK:INIT ('linear incompressibility term')
  @INIT MDBRICK:INIT ('nonlinear elasticity')
  @INIT MDBRICK:INIT ('nonlinear elasticity incompressibility term')
  @INIT MDBRICK:INIT ('small deformations plasticity')
  @INIT MDBRICK:INIT ('dynamic')
  @INIT MDBRICK:INIT ('navier stokes')
  @INIT MDBRICK:INIT ('bilaplacian')
  @INIT MDBRICK:INIT ('isotropic_linearized_plate')
  @INIT MDBRICK:INIT ('mixed_isotropic_linearized_plate')
  @INIT MDBRICK:INIT ('plate_source_term')
  @INIT MDBRICK:INIT ('plate_simple_support')
  @INIT MDBRICK:INIT ('plate_clamped_support')
  @INIT MDBRICK:INIT ('plate_closing')
  $Id$
MLABCOM*/
void gf_mdbrick(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  std::string cmd    = in.pop().to_string();
  getfemint_mdbrick *b = new getfemint_mdbrick();
  out.pop().from_object_id(workspace().push_object(b), MDBRICK_CLASS_ID);
  if (check_cmd(cmd, "constraint", in, out, 2, -1, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('constraint', @tbrick parent, @str CTYPE [, @int numfem])
      Build a generic constraint brick. 

      It may be useful in some situations, such as the Stokes problem
      where the pressure in defined modulo a constant. In such a
      situation, this brick can be used to add an additional
      constraint on the pressure value.      
      CTYPE has to be chosen among 'augmented', 'penalized', and
      'eliminated'.  The constraint can be specified with
      MDBRICK:SET('constraints'). Note that Dirichlet bricks (except the
      'generalized Dirichlet' one) are also specializations of the
      'constraint' brick. @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    bool is_complex = parent.is_complex();
    getfem::constraints_type ctype = get_constraints_type(in);
    size_type num_fem = get_num_fem(in, parent);
    SET_BRICK2(mdbrick_constraint, "Constraint", 
	       (parent.real_mdbrick(), num_fem),
	       (parent.cplx_mdbrick(), num_fem));
    b->set_constraints_type(ctype);
  } else if (check_cmd(cmd, "dirichlet", in, out, 4, -1, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('dirichlet', @tbrick parent, @int BNUM, @tmf MFMULT, @str CTYPE [, @int numfem])
      Build a Dirichlet condition brick which impose the value of a field along a mesh boundary.
      
      The BNUM parameter selects on which mesh region the Dirichlet
      condition is imposed. CTYPE has to be chosen among 'augmented',
      'penalized', and 'eliminated'. The MFMULT may generally be taken
      as the @tmf of the unknown, but for 'augmented' Dirichlet
      conditions, you may have to respect the Inf-Sup condition and
      choose an adequate @tmf.
      @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    bool is_complex = parent.is_complex();
    size_type bound = in.pop().to_integer();
    const getfem::mesh_fem &mf_mult = pop_mesh_fem(in, b);
    getfem::constraints_type ctype = get_constraints_type(in);
    size_type num_fem = get_num_fem(in, parent);
    SET_BRICK2(mdbrick_Dirichlet, "Dirichlet", 
	       (parent.real_mdbrick(), bound, mf_mult, num_fem),
	       (parent.cplx_mdbrick(), bound, mf_mult, num_fem));
    b->set_constraints_type(ctype);
  } else if (check_cmd(cmd, "dirichlet on normal component", in, out, 4, -1, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('dirichlet on normal component', @tbrick parent, @int BNUM, @tmf MFMULT, @str CTYPE [, @int numfem])
      Build a Dirichlet condition brick which imposes the value of the normal component of a vector field.
      @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    bool is_complex = parent.is_complex();
    size_type bound = in.pop().to_integer();
    const getfem::mesh_fem &mf_mult = pop_mesh_fem(in, b);
    getfem::constraints_type ctype = get_constraints_type(in);
    size_type num_fem = get_num_fem(in, parent);
    SET_BRICK2(mdbrick_normal_component_Dirichlet, 
	       "DirichletOnNormalComponent", 
	       (parent.real_mdbrick(), bound, mf_mult, num_fem),
	       (parent.cplx_mdbrick(), bound, mf_mult, num_fem));
    b->set_constraints_type(ctype);
  } else if (check_cmd(cmd, "dirichlet on normal derivative", 
		       in, out, 4, -1, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('dirichlet on normal derivative', @tbrick parent, @int BNUM, @tmf MFMULT, @str CTYPE [, @int numfem])
      Build a Dirichlet condition brick which imposes the value of the normal derivative of the unknown.
      @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    //bool is_complex = parent.is_complex();
    const getfem::mesh_fem &mf_mult = pop_mesh_fem(in, b);
    size_type bound = in.pop().to_integer();
    getfem::constraints_type ctype = get_constraints_type(in);
    size_type num_fem = get_num_fem(in, parent);
    SET_BRICK_R(mdbrick_normal_derivative_Dirichlet, 
		"DirichletOnNormalDerivative", 
		(parent.real_mdbrick(), bound, mf_mult, num_fem));
    b->set_constraints_type(ctype);
  } else if (check_cmd(cmd, "generalized dirichlet", in, out, 2, -1, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('generalized dirichlet', @tbrick parent, @int BNUM [, @int numfem])
      This is the "old" Dirichlet brick of getfem.
      
      This brick can be used to impose general Dirichlet conditions
      'h(x)u(x) = r(x)' , however it may have some issues with elaborated FEM (such as Argyris, etc). It should be avoided when possible.
      @*/

    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    bool is_complex = parent.is_complex();
    size_type bound = in.pop().to_integer();
    size_type num_fem = get_num_fem(in, parent);
    SET_BRICK2(mdbrick_generalized_Dirichlet, "GeneralizedDirichlet", 
	       (parent.real_mdbrick(), bound, num_fem),
	       (parent.cplx_mdbrick(), bound, num_fem));

  } else if (check_cmd(cmd, "source term", in, out, 1, 3, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('source term', @tbrick parent, [, @int BNUM=-1[, @int numfem]])
      Add a boundary or volumic source term ( \int B.v ).

      If BNUM is omitted (or set to -1) , the brick adds a volumic
      source term on the whole mesh. For BNUM >= 0, the source term is
      imposed on the mesh region BNUM. Use MDBRICK:SET('param','source
      term',mf,B) to set the source term field. The source term is
      expected as a vector field of size Q (with Q = qdim). @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    bool is_complex = parent.is_complex();
    
    size_type bound = size_type(-1); 
    if (in.remaining()) bound = in.pop().to_integer();
    size_type num_fem = get_num_fem(in, parent);
    const getfem::mesh_fem &mf_u = 
      parent.mdbrick().get_mesh_fem(num_fem);
    const getfem::mesh_fem &mf_d = getfem::classical_mesh_fem(mf_u.linked_mesh(), 0);
    size_type n = mf_u.get_qdim() * mf_d.nb_dof();
    SET_BRICK2(mdbrick_source_term, "SourceTerm",
	       (parent.real_mdbrick(), mf_d, 
		real_model_state::vector_type(n), bound, num_fem),
	       (parent.cplx_mdbrick(), mf_d, 
		cplx_model_state::vector_type(n), bound, num_fem));
  } else if (check_cmd(cmd, "normal source term", in, out, 2, 3, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('normal source term', @tbrick parent, @int BNUM [, @int numfem])
      Add a boundary source term ( \int (Bn).v ).

      The source term is imposed on the mesh region BNUM (which of
      course is not allowed to be a volumic region, only boundary
      regions are allowed). Use MDBRICK:SET('param','source term',mf,B)
      to set the source term field. The source term B is expected as
      tensor field of size QxN (with Q = qdim, N = mesh dim). For
      example, if you consider an elasticity problem, this brick may
      be used to impose a force on the boundary with B as the stress
      tensor.
    @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    bool is_complex = parent.is_complex();
    
    size_type bound = size_type(-1); 
    if (in.remaining()) bound = in.pop().to_integer();
    size_type num_fem = get_num_fem(in, parent);
    const getfem::mesh_fem &mf_u = 
      parent.mdbrick().get_mesh_fem(num_fem);
    const getfem::mesh_fem &mf_d = getfem::classical_mesh_fem(mf_u.linked_mesh(), 0);
    size_type n = mf_u.get_qdim() * mf_u.linked_mesh().dim() * mf_d.nb_dof();
    SET_BRICK2(mdbrick_normal_source_term, "NormalSourceTerm",
	       (parent.real_mdbrick(), mf_d, 
		real_model_state::vector_type(n), bound, num_fem),
	       (parent.cplx_mdbrick(), mf_d, 
		cplx_model_state::vector_type(n), bound, num_fem));
  } else if (check_cmd(cmd, "normal derivative source term", in, out, 2, 3, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('normal derivative source term', @tbrick parent, @int BNUM [, @int numfem])
      Add a boundary source term ( \int (\partial_n B).v ).

      The source term is imposed on the mesh region BNUM.
      Use MDBRICK:SET('param','source term',mf,B)
      to set the source term field, which is expected as a vector field 
      of size Q (with Q = qdim).
    @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    bool is_complex = parent.is_complex();
    size_type bound = in.pop().to_integer();
    size_type num_fem = get_num_fem(in, parent);
    const getfem::mesh_fem &mf_u = parent.mdbrick().get_mesh_fem(num_fem);
    const getfem::mesh_fem &mf_d = getfem::classical_mesh_fem(mf_u.linked_mesh(), 0);
    size_type n = mf_u.get_qdim() * mf_d.nb_dof();
    SET_BRICK2(mdbrick_normal_derivative_source_term, "NormalDerivativeSourceTerm",
	       (parent.real_mdbrick(), mf_d, 
		real_model_state::vector_type(n), bound, num_fem),
	       (parent.cplx_mdbrick(), mf_d, 
		cplx_model_state::vector_type(n), bound, num_fem));
  } else if (check_cmd(cmd, "neumann Kirchhoff-Love source term", in, out, 2, 3, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('neumann KirchhoffLove source term', @tbrick parent, @int BNUM [, @int numfem])

      Add a boundary source term for neumann Kirchhoff-Love plate
      problems (should be used with the Kirchhoff-Love flavour of the
      bilaplacian brick).
    @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    //bool is_complex = parent.is_complex();
    size_type bound = in.pop().to_integer();
    size_type num_fem = get_num_fem(in, parent);
    const getfem::mesh_fem &mf_u = parent.mdbrick().get_mesh_fem(num_fem);
    const getfem::mesh_fem &mf_d = getfem::classical_mesh_fem(mf_u.linked_mesh(), 0);
    size_type mdim = mf_u.linked_mesh().dim();
    SET_BRICK_R(mdbrick_neumann_KL_term, "NeumannKirchhoffLoveSourceTerm",
	       (parent.real_mdbrick(), mf_d, 
		real_model_state::vector_type(mdim*mdim*mf_d.nb_dof()), 
		real_model_state::vector_type(mdim*mf_d.nb_dof()), 
		bound, num_fem));
  } else if (check_cmd(cmd, "qu term", in, out, 1, 4, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('qu term', @tbrick parent, [, @int BNUM [, @int numfem]])
      Update the tangent matrix with a \int (Qu).v term.

      The Q(x) parameter is a matrix field of size qdim x qdim. An
      example of use is for the "iku" part of Robin boundary
      conditions \partial_n u + iku = ... 
      @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    bool is_complex = parent.is_complex();
    scalar_type rQ = 0;
    complex_type cQ = 0;
    /*if (in.remaining()) {
      if (!is_complex) rQ = in.pop().to_scalar();
      else             cQ = in.pop().to_scalar(complex_type());
      }*/
    size_type bound = size_type(-1); 
    if (in.remaining()) bound = in.pop().to_integer();
    size_type num_fem = get_num_fem(in, parent);

    SET_BRICK2(mdbrick_QU_term, "QUTerm",
	       (parent.real_mdbrick(), rQ, bound, num_fem),
	       (parent.cplx_mdbrick(), cQ, bound, num_fem));

  } else if (check_cmd(cmd, "mass matrix", in, out, 2, -1, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('mass matrix', @tmim mim, @tmf mf_u [,'real'|'complex'])
      Build a mass-matrix brick.
      @*/
    getfem::mesh_im &mim = pop_mesh_im(in, b);
    getfem::mesh_fem &mf_u = pop_mesh_fem(in, b);
    bool is_complex = get_complexity(in);
    
    //darray rho = in.pop().to_darray(mf_d.nb_dof());
    SET_BRICK(mdbrick_mass_matrix, "MassMatrix", (mim, mf_u));

  } else if (check_cmd(cmd, "generic elliptic", in, out, 2, -1, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('generic elliptic', @tmim MIM, @tmf mfu [,'scalar'|'matrix'|'tensor'][,'real'|'complex'])
      Setup a generic elliptic problem ( (A*grad(U)).grad(V) ).

      The brick parameter 'A' may be a scalar field, a matrix field, or a tensor field (default is scalar).
      @*/
    getfem::mesh_im &mim = pop_mesh_im(in, b);
    getfem::mesh_fem &mf_u = pop_mesh_fem(in, b);
    std::string s = "scalar";
    if (in.remaining()) s = in.pop().to_string();
    bool is_complex = get_complexity(in);    
    SET_BRICK(mdbrick_generic_elliptic, "ScalarElliptic", (mim, mf_u));
    if (!is_complex) {
      getfem::mdbrick_generic_elliptic<real_model_state> *bb = 
	b->cast<getfem::mdbrick_generic_elliptic<real_model_state> >();
      if (cmd_strmatch(s, "scalar")) bb->set_coeff_dimension(0);
      if (cmd_strmatch(s, "matrix")) bb->set_coeff_dimension(2);
      if (cmd_strmatch(s, "tensor")) bb->set_coeff_dimension(4);
    } else {
      getfem::mdbrick_generic_elliptic<cplx_model_state> *bb = 
	b->cast<getfem::mdbrick_generic_elliptic<cplx_model_state> >();
      if (cmd_strmatch(s, "scalar")) bb->set_coeff_dimension(0);
      if (cmd_strmatch(s, "matrix")) bb->set_coeff_dimension(2);
      if (cmd_strmatch(s, "tensor")) bb->set_coeff_dimension(4);
    }
  } else if (check_cmd(cmd, "helmholtz", in, out, 2, -1, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('helmholtz', @tmim MIM, @tmf mfu [,'real'|'complex'])
      Setup a Helmholtz problem.

      The brick has one parameter, 'wave_number'.
      @*/
    getfem::mesh_im &mim = pop_mesh_im(in, b);
    getfem::mesh_fem &mf_u = pop_mesh_fem(in, b);
    bool is_complex = get_complexity(in, true);
    SET_BRICK(mdbrick_Helmholtz, "Helmholtz", (mim, mf_u));
  } else if (check_cmd(cmd, "isotropic linearized elasticity", in, out, 2, -1, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('isotropic linearized elasticity', @tmim MIM, @tmf mfu)
      Setup a linear elasticity problem.

      The brick has two scalar parameter, 'lambda' and 'mu' (the Lame coefficients).
    @*/
    getfem::mesh_im &mim = pop_mesh_im(in, b);
    getfem::mesh_fem &mf_u = pop_mesh_fem(in, b);
    SET_BRICK_R(mdbrick_isotropic_linearized_elasticity, 
		"IsotropicLinearizedElasticity", (mim, mf_u));

  } else if (check_cmd(cmd, "linear incompressibility term", in, out, 2, -1, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('linear incompressibility term', @tbrick parent, @tmf mf_p [, @int numfem])
      Add an incompressibily constraint (div u = 0).
      @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    bool is_complex = parent.is_complex();
    getfem::mesh_fem &mf_p = pop_mesh_fem(in, b);
    size_type num_fem = get_num_fem(in, parent);
    SET_BRICK2(mdbrick_linear_incomp,
	      "LinearIncompressibilityTerm", 
	       (parent.real_mdbrick(), mf_p, num_fem),
	       (parent.cplx_mdbrick(), mf_p, num_fem));
  } else if (check_cmd(cmd, "nonlinear elasticity", in, out, 3, -1, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('nonlinear elasticity', @tmim MIM, @tmf mfu, @str lawname)
      Setup a nonlinear elasticity (large deformations) problem.

      The material law can be chosen among
      - 'SaintVenant Kirchhoff' (linearized material law)
      - 'Mooney Rivlin' (to be used with the nonlinear incompressibily term)
      - 'Ciarlet Geymonat'
      @*/
    getfem::mesh_im &mim = pop_mesh_im(in, b);
    getfem::mesh_fem &mf_u = pop_mesh_fem(in, b);
    std::string lawname = in.pop().to_string();
    std::auto_ptr<getfem::abstract_hyperelastic_law> l = 
      abstract_hyperelastic_law_from_name(lawname);
    
    real_model_state::vector_type P(l->nb_params());
    std::fill(P.begin(), P.end(), 1.);
    if (dynamic_cast<getfem::Ciarlet_Geymonat_hyperelastic_law*>(l.get())) 
      P.back() = -1; // for ciarlet geymonat, a good parameter set is [1, 1, -1] 
    SET_BRICK_R(mdbrick_nonlinear_elasticity, 
		"NonlinearElasticity", (*l,mim,mf_u,P));
    b->hyperelastic_law = l; /* won't leak .. */
  } else if (check_cmd(cmd, "nonlinear elasticity incompressibility term", in, out, 2, 3, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('nonlinear elasticity incompressibility term', @tbrick parent, @tmf mf_p [, @int numfem])
      Add an incompressibily constraint to a large strain elasticity problem.
      @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    getfem::mesh_fem &mf_p = pop_mesh_fem(in, b);
    size_type num_fem = get_num_fem(in, parent);
    
    
    SET_BRICK_R(mdbrick_nonlinear_incomp,
		"NonlinearIncompressibilityTerm", 
		(parent.real_mdbrick(), mf_p, num_fem));
  } else if (check_cmd(cmd, "small deformations plasticity", in, out, 3, 3, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('small deformations plasticity', @tmim MIM, @tmf mfu, @scalar THRESHOLD)
      Setup a plasticity problem (with small deformations).

      The THRESHOLD parameter is the maximum value of the Von Mises
      stress before 'plastification' of the material.
      @*/
    getfem::mesh_im &mim = pop_mesh_im(in, b);
    getfem::mesh_fem &mf_u = pop_mesh_fem(in, b);
    scalar_type stress_threshold = in.pop().to_scalar();
    b->plasticity_stress_projection.reset(new getfem::VM_projection());
    SET_BRICK_R(mdbrick_plasticity, 
		"SmallDeformationsPlasticity", 
		(mim, mf_u, 100, 40, stress_threshold, 
		 *b->plasticity_stress_projection));
  } else if (check_cmd(cmd, "dynamic", in, out, 2, 3, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('dynamic', @tbrick parent, @scalar rho [, @int numfem])
      Dynamic brick. This brick is not ready.
      @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    bool is_complex = parent.is_complex();
    scalar_type rho = in.pop().to_scalar();
    size_type num_fem = get_num_fem(in, parent);
    SET_BRICK2(mdbrick_dynamic, "Dynamic", 
	       (parent.real_mdbrick(),rho,num_fem),
	       (parent.cplx_mdbrick(),rho,num_fem));
  } else if (check_cmd(cmd, "bilaplacian", in, out, 2, 3, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('bilaplacian', @tmim MIM, @tmf mfu, ['Kirchhoff-Love'])
      Setup a bilaplacian problem.

      If the Kirchhoff-Love option is specified, the Kirchhoff-Love
      plate model is used.
      @*/
    getfem::mesh_im &mim = pop_mesh_im(in, b);
    getfem::mesh_fem &mf_u = pop_mesh_fem(in, b);
    bool use_KL = false;
    if (in.remaining() && in.front().is_string()) {
      std::string opt = in.pop().to_string();
      if (cmd_strmatch(opt, "KL") ||
	  cmd_strmatch(opt, "Kirchhoff-Love"))
	use_KL = true;
      else THROW_BADARG("wrong option: " << opt);
    }
    //bool is_complex = get_complexity(in);    

    SET_BRICK_R(mdbrick_bilaplacian, 
		"Bilaplacian",
		(mim, mf_u, use_KL));
  } else if (check_cmd(cmd, "navier stokes", in, out, 4, 4, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('navier stokes', @tmim MIM, @tmf mfu, @tmf mfp)
      Setup a Navier-Stokes problem (this brick is not ready, do not use it).
      @*/
    getfem::mesh_im &mim = pop_mesh_im(in, b);
    getfem::mesh_fem &mf_u = pop_mesh_fem(in, b);
    getfem::mesh_fem &mf_p = pop_mesh_fem(in, b);
    scalar_type nu = in.pop().to_scalar();
    SET_BRICK_R(mdbrick_navier_stokes, "IncompressibleNavierStokes",
		(mim, mf_u, mf_p, nu));
  } else if (check_cmd(cmd, "isotropic_linearized_plate", in, out, 6, 6, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('isotropic_linearized_plate', @tmim MIM, @tmim MIMSUB, @tmf MF_UT, @tmf MF_U3, @tmf MF_THETA, @scalar EPSILON)
      
    Setup a linear plate model brick (for moderately thick plates,
    using the Reissner-Mindlin model). EPSILON is the plate thinkness,
    the @mf MF_UT and MF_U3 are used respectively for the membrane
    displacement and the transverse displacement of the plate. The @mf
    MF_THETA is the rotation of the normal ("section rotations").

    The second integration method MIMSUB can be chosen equal to MIM,
    or different if you want to perform sub-integration on the
    transverse shear term (mitc4 projection).

    This brick has two parameters "lambda" and "mu" (the Lamé coefficients)
    @*/
    getfem::mesh_im &mim = pop_mesh_im(in, b);
    getfem::mesh_im &mim_subint = pop_mesh_im(in, b);
    getfem::mesh_fem &mf_ut = pop_mesh_fem(in, b);
    getfem::mesh_fem &mf_u3 = pop_mesh_fem(in, b);
    getfem::mesh_fem &mf_theta = pop_mesh_fem(in, b);
    scalar_type epsilon = in.pop().to_scalar();
    SET_BRICK_R(mdbrick_isotropic_linearized_plate, "IsotropicLinearizedPlate",
		(mim, mim_subint, mf_ut, mf_u3, mf_theta, 100., 40., epsilon));
  } else if (check_cmd(cmd, "mixed_isotropic_linearized_plate", in, out, 5, 5, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('mixed_isotropic_linearized_plate', @tmim MIM, @tmf MF_UT, @tmf MF_U3, @tmf MF_THETA, @scalar EPSILON)
    
    Setup a mixed linear plate model brick (for thin plates, using
    Kirchhoff-Love model).

    For a non-mixed version, use the bilaplacian brick. 
    @*/
    getfem::mesh_im &mim = pop_mesh_im(in, b);
    getfem::mesh_fem &mf_ut = pop_mesh_fem(in, b);
    getfem::mesh_fem &mf_u3 = pop_mesh_fem(in, b);
    getfem::mesh_fem &mf_theta = pop_mesh_fem(in, b);
    scalar_type epsilon = in.pop().to_scalar();
    SET_BRICK_R(mdbrick_mixed_isotropic_linearized_plate, "MixedIsotropicLinearizedPlate",
		(mim, mf_ut, mf_u3, mf_theta, 100., 40., epsilon));

  } else if (check_cmd(cmd, "plate_source_term", in, out, 1, 3, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('plate_source_term', @tbrick parent, [, @int BNUM=-1[, @int numfem]])

      Add a boundary or a volumic source term to a plate problem. This
      brick has two parameters: "B" is the displacement (ut and u3)
      source term, "M" is the moment source term (i.e. the source term
      on the rotation of the normal).
      @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    size_type bound = size_type(-1);
    if (in.remaining()) bound = in.pop().to_integer();
    size_type num_fem = get_num_fem(in, parent);
    const getfem::mesh_fem &mf_u = parent.mdbrick().get_mesh_fem(num_fem);
    const getfem::mesh_fem &mf_d = getfem::classical_mesh_fem(mf_u.linked_mesh(), 0);
    size_type n = mf_d.nb_dof();
    SET_BRICK_R(mdbrick_plate_source_term, "PlateSourceTerm",
		(parent.real_mdbrick(), mf_d, 
		 real_model_state::vector_type(3*n), 
		 real_model_state::vector_type(2*n), 
		 bound, num_fem));

  } else if (check_cmd(cmd, "plate_simple_support", in, out, 3, 4, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('plate_simple_support', @tbrick parent, @int BNUM, @str CTYPE [, @int numfem])

    Add a "simple support" boundary condition to a plate problem
    (homogeneous Dirichlet condition on the displacement, free
    rotation). CTYPE specifies how the constraint is enforced
    ('penalized', 'augmented' or 'eliminated'). @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    size_type bound = in.pop().to_integer();
    getfem::constraints_type ctype = get_constraints_type(in);
    size_type num_fem = get_num_fem(in, parent);
    SET_BRICK_R(mdbrick_plate_simple_support, "PlateSimpleSupport",
		(parent.real_mdbrick(), bound, num_fem, ctype));

  } else if (check_cmd(cmd, "plate_clamped_support", in, out, 3, 4, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('plate_clamped_support', @tbrick parent, @int BNUM, @str CTYPE [, @int numfem])

    Add a "clamped support" boundary condition to a plate problem
    (homogeneous Dirichlet condition on the displacement and on the
    rotation). CTYPE specifies how the constraint is enforced
    ('penalized', 'augmented' or 'eliminated'). @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    size_type bound = in.pop().to_integer();
    getfem::constraints_type ctype = get_constraints_type(in);
    size_type num_fem = get_num_fem(in, parent);
    SET_BRICK_R(mdbrick_plate_clamped_support, "PlateClampedSupport",
		(parent.real_mdbrick(), bound, num_fem, ctype));
  } else if (check_cmd(cmd, "plate_closing", in, out, 1, 2, 0, 1)) {
    /*@INIT B=MDBRICK:INIT('plate_closing', @tbrick parent [, @int numfem])

    Add a free edges condition for the mixed plate model brick.

    This brick is required when the mixed linearized plate brick is
    used. It must be inserted after all other boundary conditions
    (the reason is that the brick has to inspect all other boundary
    conditions to determine the number of disconnected boundary parts
    which are free edges). @*/
    getfemint_mdbrick &parent = pop_mdbrick(in, b);
    size_type num_fem = get_num_fem(in, parent);
    SET_BRICK_R(mdbrick_plate_closing, "PlateClosing",
		(parent.real_mdbrick(), num_fem));
  } else bad_cmd(cmd);

  if (in.remaining()) THROW_BADARG("too many arguments");
}
