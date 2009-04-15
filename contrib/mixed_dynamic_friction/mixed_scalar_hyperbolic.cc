// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2008 Yves Renard.
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

/**
 * Dynamic friction in linear elasticity.
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
 */


#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_import.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_Coulomb_friction.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_derivatives.h"
#include "gmm/gmm.h"
#include <fstream>

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::short_type;   /* = unsigned long */
using bgeot::dim_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. 
 * default types of getfem_model_solvers.h
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/*
 * structure for the friction problem
 */
struct hyperbolic_problem {

  enum {
    DIRICHLET_BOUNDARY, CONTACT_BOUNDARY
  };
  getfem::mesh mesh;  /* the mesh */
  getfem::mesh_im  mim;      /* integration methods.                         */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the friction solution     */
  getfem::mesh_fem mf_v;     /* main mesh_fem, for the friction solution     */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_mult;  /* mesh_fem for the Dirichlet condition.        */

  scalar_type residual;      /* max residual for the iterative solvers       */
  size_type N, noisy, scheme;
  scalar_type T, dt, r;
  scalar_type dirichlet_val, dtexport;
  bool  dxexport, compare;

  std::string datafilename, refrootname, INTEGRATION;
  bgeot::md_param PARAM;

  void solve(void);
  void init(void);
  hyperbolic_problem(void) : mim(mesh), mf_u(mesh), mf_v(mesh),
			     mf_rhs(mesh), mf_mult(mesh) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void hyperbolic_problem::init(void) {
  std::string FEM_TYPE_U  = PARAM.string_value("FEM_TYPE_U","FEM name");
  std::string FEM_TYPE_V  = PARAM.string_value("FEM_TYPE_V","FEM name");
  INTEGRATION = PARAM.string_value("INTEGRATION","Name of integration method");
  cout << "FEM_TYPE_U = "  << FEM_TYPE_U << "\n";
  cout << "FEM_TYPE_V = "  << FEM_TYPE_V << "\n";
  cout << "INTEGRATION = " << INTEGRATION << "\n";

  std::string meshname
    (PARAM.string_value("MESHNAME", "Nom du fichier de maillage"));

  getfem::import_mesh(meshname, mesh);
  N = mesh.dim();
  mesh.optimize_structure();

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  refrootname = PARAM.string_value("REFROOTNAME",
				   "Base name of reference files.");
  compare = (PARAM.int_value("COMPARE") != 0);
  residual = PARAM.real_value("RESIDUAL");
  if (residual == 0.) residual = 1e-10;

  dirichlet_val = PARAM.real_value("DIRICHLET_VAL","Dirichlet condition");
  T = PARAM.real_value("T", "from [0,T] the time interval");
  dt = PARAM.real_value("DT", "time step");
  dxexport = (PARAM.int_value("DX_EXPORT", "Exporting on OpenDX format")
	      != 0);
  dtexport = PARAM.real_value("DT_EXPORT", "time step for the export");
  dtexport = dt * double(int(dtexport / dt + 0.5));
 
  r = PARAM.real_value("R", "augmentation parameter");
  noisy = PARAM.int_value("NOISY", "verbosity of iterative methods");
  scheme = PARAM.int_value("SCHEME", "scheme");

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE_U);
  getfem::pfem pf_v = getfem::fem_descriptor(FEM_TYPE_V);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_v.set_finite_element(mesh.convex_index(), pf_v);
  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    if (!pf_u->is_lagrange()) {
      GMM_ASSERT1(false, "You are using a non-lagrange FEM. "
		  << "In that case you need to set "
		  << "DATA_FEM_TYPE in the .param file");
    }
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }

  std::string mult_fem_name = PARAM.string_value("MULT_FEM_TYPE");
  if (mult_fem_name.size() == 0) {
    if (!pf_u->is_lagrange()) {
      GMM_ASSERT1(false, "You are using a non-lagrange FEM. "
		  << "In that case you need to set "
		  << "DATA_FEM_TYPE in the .param file");
    }
    mf_mult.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_mult.set_finite_element(mesh.convex_index(), 
			       getfem::fem_descriptor(mult_fem_name));
  }
  
  /* set boundary conditions */
  base_node center(0.,0.,20.);
  std::cout << "Reperage des bord de contact et Dirichlet\n";  
  for (dal::bv_visitor cv(mesh.convex_index()); !cv.finished(); ++cv) {
    short_type nf = mesh.structure_of_convex(cv)->nb_faces();
    for (short_type f = 0; f < nf; f++) {
      if (!mesh.is_convex_having_neighbour(cv, f)) {
	base_small_vector un = mesh.normal_of_face_of_convex(cv, f);
	un /= gmm::vect_norm2(un);	
	base_node pt = mesh.points_of_face_of_convex(cv,f)[0];
	mesh.region(DIRICHLET_BOUNDARY).add(cv, f);

	// if (un[N-1] < -0.000001 && (N != 3 || (gmm::vect_dist2(pt, center)
	//		   > .99*sqrt(25. + 15*15) && pt[N-1] < 20.1)))
	//  mesh.region(CONTACT_BOUNDARY).add(cv, f);
	// if (un[N-1] > 0.1)
	//  mesh.region(DIRICHLET_BOUNDARY).add(cv, f);
      }
    }
  }
}


/**************************************************************************/
/*  Brique dynamique mixte.                                               */
/**************************************************************************/

namespace getfem {

# define MDBRICK_MIXED_DYNAMIC 738053

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_mixed_dynamic : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    const mesh_fem *mf_u;
     const mesh_fem &mf_v;
    mdbrick_parameter<VECTOR> RHO_;
    VECTOR DFU, DFV;
    T_MATRIX B_, C_;
    size_type num_fem;
    value_type Bcoef, Kcoef, Ccoef;
    std::set<size_type> boundary_sup;
    bool M_uptodate;

    virtual void proper_update(void) {
      mf_u = this->mesh_fems[num_fem];
      M_uptodate = false;
    }

    void proper_update_M(void) {
      GMM_TRACE2("Assembling mass matrices for mdbrick_mixed_dynamic");
      gmm::clear(B_);
      gmm::resize(B_, mf_v.nb_dof(), mf_u->nb_dof());
      asm_mass_matrix_param(B_, *(this->mesh_ims[0]), mf_v, *mf_u,
			    RHO_.mf(), RHO_.get());
      gmm::clear(C_);
      gmm::resize(C_, mf_v.nb_dof(), mf_v.nb_dof());
      asm_mass_matrix_param(C_, *(this->mesh_ims[0]), mf_v,
			    RHO_.mf(), RHO_.get());
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), mf_v.nb_dof());
      gmm::sub_interval SUBJ(i0+this->mesh_fem_positions[num_fem],
			     mf_u->nb_dof());

      if (Kcoef != value_type(1)) gmm::scale(MS.tangent_matrix(), Kcoef);
      gmm::copy(gmm::scaled(get_B(), Bcoef),
		gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
      gmm::copy(gmm::transposed(gmm::scaled(get_B(), Bcoef)),
		gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
      gmm::copy(gmm::scaled(get_C(), -Ccoef),
		gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
    }
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), mf_v.nb_dof());
      gmm::sub_interval SUBJ(i0+this->mesh_fem_positions[num_fem],
			     mf_u->nb_dof());

      if (Kcoef != value_type(1))  gmm::scale(MS.residual(), Kcoef);
      gmm::mult(get_B(), gmm::scaled(gmm::sub_vector(MS.state(), SUBJ), Bcoef),
		gmm::sub_vector(MS.residual(), SUBI));
      gmm::mult_add(gmm::transposed(get_B()),
		    gmm::scaled(gmm::sub_vector(MS.state(), SUBI), Bcoef),
		    gmm::sub_vector(MS.residual(), SUBJ));
      gmm::mult_add(get_C(),
		    gmm::scaled(gmm::sub_vector(MS.state(), SUBI), -Ccoef),
		    gmm::sub_vector(MS.residual(), SUBI));
      if (gmm::vect_size(DFU) > 0)
	gmm::add(gmm::scaled(DFU, -value_type(1)),
		 gmm::sub_vector(MS.residual(), SUBJ));
      if (gmm::vect_size(DFV) > 0)
	gmm::add(gmm::scaled(DFV, -value_type(1)),
		 gmm::sub_vector(MS.residual(), SUBI));
    }

    void set_dynamic_coeff(value_type a, value_type b, value_type c)
    { Bcoef=a; Kcoef=b; Ccoef = c; }
    template <class VEC> void set_DFU(const VEC &DF_)
    { gmm::resize(DFU, gmm::vect_size(DF_)); gmm::copy(DF_, DFU); }
    template <class VEC> void set_DFV(const VEC &DF_)
    { gmm::resize(DFV, gmm::vect_size(DF_)); gmm::copy(DF_, DFV); }

    const T_MATRIX &get_B(void) {
      this->context_check();
      if (!M_uptodate || this->parameters_is_any_modified()) {
	proper_update_M();
	M_uptodate = true;
	this->parameters_set_uptodate();
      }
      return B_; 
    }

    const T_MATRIX &get_C(void) {
      this->context_check();
      if (!M_uptodate || this->parameters_is_any_modified()) {
	proper_update_M();
	M_uptodate = true;
	this->parameters_set_uptodate();
      }
      return C_; 
    }

    /** provide access to the value of the solution corresponding to
	the local mesh_fem mf_v. */
    SUBVECTOR get_V(MODEL_STATE &MS) {
      gmm::sub_interval SUBU = gmm::sub_interval(this->first_index()
				     + sub_problem.nb_dof(), mf_v.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    /**
       @param num_fem_ the mesh_fem number on which this brick is is applied.
    */
    mdbrick_mixed_dynamic(mdbrick_abstract<MODEL_STATE> &problem,
			  const mesh_fem &mf_v_, value_type RHO__, 
			  size_type num_fem_=0)
      : sub_problem(problem), mf_v(mf_v_), RHO_("rho", this),
	num_fem(num_fem_) {
      
      Bcoef = Kcoef = Ccoef = value_type(1);
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = false;
      this->add_proper_mesh_fem(mf_v, MDBRICK_MIXED_DYNAMIC);
      this->force_update();
      
      RHO_.set(classical_mesh_fem(mf_u->linked_mesh(), 0), RHO__);
    }
  };



}



/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

void hyperbolic_problem::solve(void) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  N = mesh.dim();
  cout << "Number of dof for u: " << mf_u.nb_dof() << endl;
  cout << "Number of dof for v: " << mf_v.nb_dof() << endl;

  // Choosing a reference dof.
  GMM_ASSERT1(!mf_u.is_reduced(), "To be adapted");
  size_type ref_dof = 0, ref_mult = size_type(-1);
  base_node P(N); gmm::fill(P, 0.5);
  for (size_type i = 1; i < mf_u.nb_dof(); ++i)
    if (gmm::vect_dist2(mf_u.point_of_basic_dof(i), P)
	< gmm::vect_dist2(mf_u.point_of_basic_dof(ref_dof), P))
      ref_dof = i;
  cout << "ref_dof = " << ref_dof << " point = "
       << mf_u.point_of_basic_dof(ref_dof) << endl;

  // Linearized elasticity brick.
  getfem::mdbrick_generic_elliptic<> ELAS(mim, mf_u);

  scalar_type h = mesh.minimal_convex_radius_estimate();
  cout << "minimal convex radius estimate : " << h << endl;
  cout << "CFL estimate = " << h << endl;

  // Defining the volumic source term.
  scalar_type vol_term = PARAM.real_value("VOLUMIC_TERM",
					  "Volumic source term");
  plain_vector F(nb_dof_rhs, vol_term);
  // Volumic source term brick.
  getfem::mdbrick_source_term<> VOL_F(ELAS, mf_rhs, F);

  // Dirichlet condition brick.
  gmm::fill(F, dirichlet_val);
//   for (size_type i = 0;i < mf_rhs.nb_dof(); ++i)
//     F[i] = dirichlet_val * (1.0 - gmm::sqr(0.5 - mf_rhs.point_of_dof(i)[0]))
//       * (1.0 - 0.2*mf_rhs.point_of_dof(i)[1]);
  
  getfem::mdbrick_Dirichlet<> DIRICHLET(VOL_F, DIRICHLET_BOUNDARY, mf_mult);
  DIRICHLET.rhs().set(mf_rhs, F);
  
  // Contact condition for Lagrange elements
  dal::bit_vector cn, dn = mf_u.basic_dof_on_region(DIRICHLET_BOUNDARY);
  for (size_type i = 0; i < mf_u.nb_dof(); ++i)
    if (!dn.is_in(i)) {

      if (mesh.points().search_node(mf_u.point_of_basic_dof(i)) != size_type(-1))
	cn.add(i);
    }

 
  //  cout << "cn = " << cn << endl;
  cout << "Number of contacting nodes : " << cn.card() << endl;
  sparse_matrix BN(cn.card(), mf_u.nb_dof());
  plain_vector gap(cn.card());
  size_type jj = 0;
  for (dal::bv_visitor i(cn); !i.finished(); ++i) {
    BN(jj, i) = -1.;
    if (i == ref_dof) ref_mult = jj;
    ++jj;
  }
  GMM_ASSERT1(ref_mult != size_type(-1), "No ref_mult !");

  getfem::mdbrick_Coulomb_friction<> FRICTION(DIRICHLET, BN, gap);
  FRICTION.set_r(r);
  
  // Dynamic brick.
  getfem::mdbrick_mixed_dynamic<> DYNAMIC(FRICTION, mf_v, 1.0);
  
  cout << "Total number of variables: " << DYNAMIC.nb_dof() << endl;
  getfem::standard_model_state MS(DYNAMIC);
  
  // cout << "C = " << DYNAMIC.get_C() << endl;

  plain_vector WN(mf_u.nb_dof());
  plain_vector DFU(mf_u.nb_dof()), DFV(mf_v.nb_dof());
  plain_vector U0(mf_u.nb_dof()), V0(mf_v.nb_dof());
  plain_vector U1(mf_u.nb_dof()), V1(mf_v.nb_dof());
  plain_vector Udemi(mf_u.nb_dof()), Vdemi(mf_v.nb_dof());
  plain_vector LN1(gmm::mat_nrows(BN));
  scalar_type t(0), t_export(dtexport), beta_(0);
  scalar_type alpha_(0);

  // Initial conditions (U0, V0)
  gmm::clear(V0);
  {  // projection of F on U0.
    sparse_matrix BB(mf_u.nb_dof(), mf_rhs.nb_dof());
    sparse_matrix MM(mf_u.nb_dof(), mf_u.nb_dof());
    plain_vector BBF(mf_u.nb_dof());
    getfem::asm_mass_matrix(MM, mim, mf_u);
    getfem::asm_mass_matrix(BB, mim, mf_u, mf_rhs);
    gmm::mult(BB, F, BBF);
    double rcond;
    gmm::SuperLU_solve(MM, U0, BBF, rcond);
  }
  
  gmm::iteration iter(residual, 0, 40000);
  iter.set_noisy(int(noisy));

  std::auto_ptr<getfem::dx_export> exp;
  getfem::stored_mesh_slice sl;
  if (dxexport) {
    exp.reset(new getfem::dx_export(datafilename + ".dx", false));
    if (N <= 2)
      sl.build(mesh, getfem::slicer_none(),4);
    else
      sl.build(mesh, getfem::slicer_boundary(mesh),4);
    exp->exporting(sl,true);
    exp->exporting_mesh_edges();
    exp->write_point_data(mf_u, U0, "stepinit"); 
    exp->serie_add_object("deformationsteps");
  }

  std::ofstream fileout((datafilename+".data").c_str(), std::ios::out);
 
  cout << "t=0, initial energy: " << 0.5*gmm::vect_sp(ELAS.get_K(), U0, U0)
    + 0.5 * gmm::vect_sp(DYNAMIC.get_C(), V0, V0)
    - gmm::vect_sp(VOL_F.get_F(), U0) << endl;
  

  size_type nitexp = 0;
  while (t <= T) {

    switch(scheme) {
    case 0 : beta_ = 1./dt; alpha_ = 1.; break;
    case 1 : beta_ = 2./dt; alpha_ = 1.; break;
    }
    
    gmm::mult(DYNAMIC.get_B(), gmm::scaled(U0, beta_), DFV);    
    gmm::mult(gmm::transposed(DYNAMIC.get_B()), gmm::scaled(V0, beta_), DFU);
    gmm::clear(WN);

    FRICTION.set_WN(WN); FRICTION.set_r(r); 
    FRICTION.set_alpha(1.); 
    DYNAMIC.set_dynamic_coeff(beta_, alpha_, 1.);
    DYNAMIC.set_DFU(DFU);
    DYNAMIC.set_DFV(DFV);
    
    iter.init();
    gmm::default_newton_line_search ls(size_type(-1), 4.0/3.0,
				       1.0/20.0, 9.0/10.0, 1.1);
    getfem::standard_solve(MS, DYNAMIC, iter,
			   getfem::default_linear_solver(DYNAMIC), ls);
    
    
    gmm::copy(ELAS.get_solution(MS), Udemi);
    gmm::copy(DYNAMIC.get_V(MS), Vdemi);
    gmm::copy(FRICTION.get_LN(MS), LN1);

    switch(scheme) {
    case 0 :
      gmm::copy(Vdemi, V1);
      gmm::copy(Udemi, U1);
      break;
    case 1 :
      gmm::add(gmm::scaled(V0, -1.), gmm::scaled(Vdemi, 2.), V1);
      gmm::add(gmm::scaled(U0, -1.), gmm::scaled(Udemi, 2.), U1);
      break;
    }

    scalar_type J1 =  0.5*gmm::vect_sp(ELAS.get_K(), U1, U1)
      + 0.5 * gmm::vect_sp(DYNAMIC.get_C(), V1, V1)
      - gmm::vect_sp(VOL_F.get_F(), U1);

    scalar_type Jdemi = 0.5*gmm::vect_sp(ELAS.get_K(), Udemi, Udemi)
      + 0.5 * gmm::vect_sp(DYNAMIC.get_C(), Vdemi, Vdemi)
      - gmm::vect_sp(VOL_F.get_F(), Udemi);

    t += dt; 
    scalar_type umin = dirichlet_val;
    size_type nbco = 0;
    for (dal::bv_visitor i(cn); !i.finished(); ++i)
      { umin = std::min(Udemi[i], umin); if (Udemi[i] < 1E-4) nbco++; }
    cout << "t = " << t << " J1 = " << J1 << " Jdemi = " << Jdemi
	 << " nbcontact : " << nbco << " umin = " << umin << endl;
    
    gmm::copy(U1, U0); gmm::copy(V1, V0);

    fileout << t << "   \t" << J1 << "   \t" << Udemi[ref_dof] << "   \t"
	    << LN1[ref_mult] << "\n";

    
    if (dxexport && t >= t_export-dt/20.0) {
      
      t_export += dtexport; ++nitexp;

      exp->write_point_data(mf_u, Udemi);
      exp->serie_add_object("deformationsteps");

      if (PARAM.int_value("VTK_EXPORT")) {
	char numit[100]; sprintf(numit, "%d", int(nitexp));
	cout << "export to " << datafilename + numit + ".vtk" << "..\n";
	getfem::vtk_export vexp(datafilename + numit + ".vtk",
			   PARAM.int_value("VTK_EXPORT")==1);
	vexp.exporting(mf_u); 
	vexp.write_point_data(mf_u, Udemi, "solution");
	cout << "export done, you can view the data file with (for example)\n"
	  "mayavi -d " << datafilename << numit << ".vtk -f "
	  "WarpScalar -m BandedSurfaceMap -m Outline\n";
	
      }
    }

  }

  // save the last time step for convergence test.
  mf_u.write_to_file(datafilename + ".mfu", true);
  gmm::vecsave(datafilename + ".U", U0);
  
  if (compare) {
    getfem::mesh m_ref;
    m_ref.read_from_file(refrootname + ".mfu");
    getfem::mesh_fem mf_ref(m_ref);
    mf_ref.read_from_file(refrootname + ".mfu");
    plain_vector Uref(mf_ref.nb_dof());
    gmm::vecload(refrootname + ".U", Uref);
    
    plain_vector U(mf_ref.nb_dof());
    
    getfem::interpolation(mf_u, mf_ref, U0, U, true);
    getfem::mesh_im mim_ref(m_ref);
    getfem::pintegration_method ppi = 
      getfem::int_method_descriptor(INTEGRATION);
    mim_ref.set_integration_method(m_ref.convex_index(), ppi);
    
    cout << "To ref L2 ERROR:"
	 << getfem::asm_L2_dist(mim_ref, mf_ref, U, mf_ref, Uref) << endl;
    
    cout << "To ref H1 ERROR:"
	 << getfem::asm_H1_dist(mim_ref, mf_ref, U, mf_ref, Uref) << endl;
    
  }
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  hyperbolic_problem p;
  p.PARAM.read_command_line(argc, argv);
  p.init();
  p.solve();
  
  cout << "To see the simulation, you have to set DX_EXPORT to 1 in "
    "dynamic_friction.param and DT_EXPORT to a suitable value (for "
    "instance equal to DT). Then you can use Open_DX (type just \"dx\" "
    "if it is installed on your system) with the Visual Program "
    "mixed_scalar_hyperbolic.net (use for instance \"Edit Visual Programs ...\" "
    "with dynamic_friction.net, then \"execute once\" in Execute menu and "
    "use the sequencer to start the animation).\n";
  

  return 0; 
}
