/*===========================================================================
 
 Copyright (C) 2002-2012 Yves Renard, Julien Pommier.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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

/**
 * Linear Elastostatic problem with a growing delaminated crack.
 *
 * Research program.
 */

#if 0  // Broken program which uses old getfem brick. Should be adapted.
       // Kept because it is the only example of 3D delamination.

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_mesh_fem_product.h"
#include "getfem/getfem_mesh_fem_global_function.h"
#include "getfem/getfem_spider_fem.h"
#include "getfem/getfem_mesh_fem_sum.h"
#include "getfem/getfem_import.h"
#include "gmm/gmm.h"



using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;


/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_vector;
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::short_type;  /* = short */
using bgeot::dim_type;
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

// #define VALIDATE_XFEM

#ifdef VALIDATE_XFEM

struct exact_solution_3D {
  getfem::mesh_fem_global_function mf;
  getfem::base_vector U;

  exact_solution_3D(getfem::mesh &me) : mf(me) {}
  
  void init(int mode, scalar_type lambda, scalar_type mu,
	    getfem::level_set &ls) {
    std::vector<getfem::pglobal_function> cfun(4);
    for (unsigned j=0; j < 4; ++j) {
      getfem::crack_singular_xy_function *s = 
	new getfem::crack_singular_xy_function(j);
      cfun[j] = getfem::global_function_on_level_set(ls, *s);
    }
    mf.set_functions(cfun);
    
    mf.set_qdim(1);
    assert(mf.linked_mesh().dim() == 3);
    U.resize(4*mf.linked_mesh().dim()); assert(mf.nb_dof() == 4);
    getfem::base_vector::iterator it = U.begin();
    scalar_type coeff=0.;
    switch(mode) {
      case 1: {
	scalar_type A=2+2*mu/(lambda+2*mu), B=-2*(lambda+mu)/(lambda+2*mu);
	/* "colonne" 1: ux, colonne 2: uy */
	*it++ = 0;       *it++=0; *it++ = A-B; /* sin(theta/2) */
	*it++ = -(A+B);  *it++=0; *it++ = 0;   /* cos(theta/2) */
	*it++ = B;       *it++=0; *it++ = 0;   /* sin(theta/2)*sin(theta) */ 
	*it++ = 0;       *it++=0; *it++ = B;   /* cos(theta/2)*cos(theta) */
	coeff = 1/sqrt(2*M_PI);
      } break;
      case 2: {
	scalar_type C1 = (lambda+3*mu)/(lambda+mu); 
	*it++ = 1-C1-2; *it++=0; *it++ = 0;
	*it++ = 0;      *it++=0; *it++ = -(C1-2+1);
	*it++ = 0;      *it++=0; *it++ = 1;
	*it++ = -1;     *it++=0; *it++ = 0;
	coeff = 2*(mu+lambda)/(lambda+2*mu)/sqrt(2*M_PI);
      } break;
      default:
	assert(0);
	break;
    }
    U *= coeff;
  }
};

base_small_vector sol_f(const base_node &x) {
  return base_small_vector(x.size());
}

#else

base_small_vector sol_f(const base_node &x) {
  return base_small_vector(x.size());
}

#endif

scalar_type cutoff_C2(scalar_type r, scalar_type r1, scalar_type r2) {
  if (r < r1) return 1.0;
  if (r < r2) return pow(r2 - r, 3.) *
    (r2*r2 - 5*r1*r2 + 3*r*r2 + 10*r1*r1 + 6*r*r - 15*r*r1) / pow(r2 - r1, 5.);
  return 0.0;
}

/**************************************************************************/
/*  Structure for the crack problem.                                      */
/**************************************************************************/

struct crack_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1};
  getfem::mesh mesh;  /* the mesh */
  getfem::mesh_level_set mls;       /* the level set aware mesh.             */
  getfem::mesh_im standard_mim; 
  getfem::mesh_im_level_set mim;    /* the integration methods.              */
  getfem::mesh_im_level_set mim_crack;    /* the integration on the crack.   */
  getfem::mesh_fem mf_pre_u;
  getfem::mesh_fem mf_mult;
  getfem::mesh_fem_level_set mfls_u; 
  getfem::mesh_fem_global_function mf_sing_u;
  getfem::mesh_fem mf_partition_of_unity;
  getfem::mesh_fem_product mf_product;
  getfem::mesh_fem_sum mf_u_sum;
  getfem::mesh_fem mf_SD;

  getfem::mesh_fem& mf_u() { return mf_u_sum; }
  

  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
#ifdef VALIDATE_XFEM
  exact_solution_3D exact_sol;
#endif
  scalar_type lambda, mu;    /* Lamé coefficients.                           */
  scalar_type Gc;
  scalar_type neumann_force;
  bool is_prescribed_disp;
  scalar_type prescribed_disp;

  getfem::level_set ls;      /* The two level sets defining the crack.       */
  
  scalar_type residual;       /* max residual for the iterative solvers      */
  unsigned dir_with_mult;
  scalar_type enr_area_radius, min_h;
  int enrichment_option;
  
  std::string datafilename;
  bgeot::md_param PARAM;

  scalar_type fracture_energy(void);
  void shape_derivative(const plain_vector &U, plain_vector &SD);
  void update_level_set(const plain_vector &SD);
  bool solve(plain_vector &U);
  void save_U(plain_vector &U);
  void init(void);
  crack_problem(void)
    : mls(mesh), standard_mim(mesh), mim(mls), 
      mim_crack(mls, getfem::mesh_im_level_set::INTEGRATE_BOUNDARY),
      mf_pre_u(mesh), mf_mult(mesh), mfls_u(mls, mf_pre_u),
      mf_sing_u(mesh), mf_partition_of_unity(mesh),
      mf_product(mf_partition_of_unity, mf_sing_u), 
      mf_u_sum(mesh), mf_SD(mesh), mf_rhs(mesh),
#ifdef VALIDATE_XFEM	
			exact_sol(mesh),
#endif
			ls(mesh, 1, true) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */

static base_node Pmin, Pmax;

void crack_problem::init(void) {
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
					 "Name of simplex integration method");
  std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION");


  enrichment_option = int(PARAM.int_value("ENRICHMENT_OPTION",
					  "Enrichment option"));
  prescribed_disp = PARAM.real_value("PRESCRIBED_DISPLACEMENT",
				     "prescribed displacement");
  is_prescribed_disp = (PARAM.int_value("IS_PRESCRIBED_DISPLACEMENT",
					"prescribed displacement") != 0);
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";
  
  /* First step : build the mesh */
  std::string MESH_FILE = PARAM.string_value("MESH_FILE", "Mesh file");
  getfem::import_mesh(MESH_FILE, mesh);
  unsigned N = mesh.dim();
  base_node Ptrans(N);
  Pmin.resize(N); Pmax.resize(N); 
  mesh.bounding_box(Pmin, Pmax);
  Ptrans[N-1] = -(Pmin[N-1] + Pmax[N-1])/2.0;
  mesh.translation(Ptrans);
  min_h = mesh.minimal_convex_radius_estimate();

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;
  enr_area_radius = PARAM.real_value("RADIUS_ENR_AREA",
				     "radius of the enrichment area");

  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  Gc = PARAM.real_value("GC", "Fracture energy density");
  neumann_force = PARAM.real_value("NEUMANN_FORCE", "Neumann force");

  mf_u().set_qdim(dim_type(N));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);
  getfem::pintegration_method simp_ppi = 
    getfem::int_method_descriptor(SIMPLEX_INTEGRATION);
  getfem::pintegration_method sing_ppi = 
    (SINGULAR_INTEGRATION.size() ? 
     getfem::int_method_descriptor(SINGULAR_INTEGRATION) : 0);
  
  standard_mim.set_integration_method(mesh.convex_index(), ppi);
  mim.set_integration_method(mesh.convex_index(), ppi);
  mim_crack.set_integration_method(mesh.convex_index(), ppi);
  mls.add_level_set(ls);
  mim.set_simplex_im(simp_ppi, sing_ppi);
  mim_crack.set_simplex_im(simp_ppi, sing_ppi);
  mf_pre_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_mult.set_finite_element(mesh.convex_index(), pf_u);
  mf_mult.set_qdim(dim_type(N));
  mf_partition_of_unity.set_classical_finite_element(1);
  mf_SD.set_qdim(dim_type(N));
  mf_SD.set_classical_finite_element(ls.degree());

  dir_with_mult = unsigned(PARAM.int_value("DIRICHLET_WITH_MULTIPLIERS", "Well"));

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }

  /* set boundary conditions
   * (Neuman on the upper face, Dirichlet elsewhere) */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
#ifdef VALIDATE_XFEM
    base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);
    if (gmm::abs(un[1]) <= 1e-7) {
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
    }
#else
    base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);
    base_node barycenter
      = gmm::mean_value(mesh.points_of_face_of_convex (i.cv(), i.f()).begin(),
			mesh.points_of_face_of_convex(i.cv(), i.f()).end());
    if (gmm::abs(un[2] - 1.0) <= 1.0E-7) // new Neumann face
      if (is_prescribed_disp && barycenter[0] > Pmax[0]*0.9)
	mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
      else
	mesh.region(NEUMANN_BOUNDARY_NUM).add(i.cv(), i.f());
    else if (is_prescribed_disp && gmm::abs(un[2] + 1.0) <= 1.0E-7 && barycenter[0] > Pmax[0]*0.9)
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
    else if (gmm::abs(un[0] + 1.0) <= 1.0E-7) // new Dirichlet face
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
#endif
  }

#ifdef VALIDATE_XFEM
  exact_sol.init(1, lambda, mu, ls);
#endif
}



base_small_vector ls_function(const base_node P) {
  scalar_type x = P[0] /*, y = P[1]*/, z = P[2];
  base_small_vector res(2);
  res[0] = z;
  res[1] = Pmax[0]*0.6 - x;
  return res;
}

/**************************************************************************/
/*  Computation of the shape derivative.                                  */
/**************************************************************************/

template<typename VECT1> class shape_der_nonlinear_term 
  : public getfem::nonlinear_elem_term {
  
  const getfem::mesh_fem &mf;
  std::vector<scalar_type> U;
  size_type N;
  base_vector coeff;
  base_matrix gradU, E, Sigma;
  bgeot::multi_index sizes_;
  scalar_type lambda, mu;
  
public:
  shape_der_nonlinear_term(const getfem::mesh_fem &mf_, const VECT1 &U_,
			  scalar_type lambda_, scalar_type mu_) 
    : mf(mf_), U(mf_.nb_basic_dof()), N(mf_.get_qdim()),
      gradU(N, N), E(N, N), Sigma(N,N), sizes_(N,N),
      lambda(lambda_), mu(mu_) {
    assert(N == mf_.linked_mesh().dim());
    mf.extend_vector(U_, U);
  }
  
  const bgeot::multi_index &sizes(size_type) const { return sizes_; }
  
  virtual void compute(getfem::fem_interpolation_context& ctx,
		       bgeot::base_tensor &t) {
    size_type cv = ctx.convex_num();
    coeff.resize(mf.nb_basic_dof_of_element(cv));
    gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf.ind_basic_dof_of_element(cv))),
	      coeff);
    ctx.pf()->interpolation_grad(ctx, coeff, gradU, mf.get_qdim());
    gmm::copy(gradU, E); gmm::add(gmm::transposed(gradU), E);
    gmm::scale(E, 0.5);
    gmm::copy(gmm::identity_matrix(), Sigma);
    gmm::scale(Sigma, lambda * gmm::mat_trace(E));
    gmm::add(gmm::scaled(E, 2.*mu), Sigma);
    
    for (size_type i = 0; i < N; ++i) 
      for (size_type j = 0; j < N; ++j) {
	t(i,j) = 0.0;
	for (size_type k = 0; k < N; ++k)
	  t(i,j) += Sigma(i,k) * gradU(k,j); 
      }
  }
};


scalar_type crack_problem::fracture_energy(void) {
  plain_vector V(1);
  getfem::generic_assembly assem1("V()+=comp()");
  assem1.push_mi(mim_crack);
  assem1.push_vec(V);
  assem1.assembly();
  return V[0]*Gc;
}

void crack_problem::shape_derivative(const plain_vector &U, plain_vector &SD) {
  gmm::clear(SD);
  
  GMM_TRACE2("Computing shape derivative, part 1");

  shape_der_nonlinear_term<plain_vector> nl(mf_u(), U, lambda, mu);

  // derivative of fracture energy
  getfem::generic_assembly assem1("V(#1)+=comp(vGrad(#1))(:,i,i)");
  assem1.push_mi(mim_crack);
  assem1.push_mf(mf_SD);
  assem1.push_vec(SD);
  assem1.assembly();

  gmm::scale(SD, Gc);


  GMM_TRACE2("Computing shape derivative, part 2");
#if 1
  // derivative of elastic energy
  getfem::generic_assembly assem2
    ("t=comp(NonLin(#1).vGrad(#2));"
     //"t2=comp(vGrad(#1).vGrad(#1).vGrad(#2))(k,:,:,l,:,:,:,:,:).u(k).u(l);"
     //"print u; print(t(:,:,1,1,1)); print(t2(l,l,:,:,1,1,1));"
     "V(#2) += 0.5*t(i,i,:,j,j) - t(i,j,:,j,i);");
  assem2.push_mi(mim);
  assem2.push_mf(mf_u());
  assem2.push_nonlinear_term(&nl);
  assem2.push_mf(mf_SD);
  assem2.push_data(U);
  assem2.push_vec(SD);

  double t0 = gmm::uclock_sec();
  assem2.assembly();
  cerr << " done : " << gmm::uclock_sec() - t0 << "\n";

#else 

  // derivative of elastic energy
  char s[500];
  sprintf(s, "u=data(#1);"
	  // "t=comp(vGrad(#1).vGrad(#1).vGrad(#2))(k,:,: ,l,:,:,:,:,:)
	  // ".u(k).u(l);"
	  "t=comp(vGrad(#1)(k,:,:).vGrad(#1)(l,:,:).vGrad(#2).u(k).u(l));"
	  "w=t(:,:,:,:,:,m,m)*0.5; z=t(:,:,:,l,:,l,:);"
	  "V(#2)+= %g*(w(i,j,i,j,:)+w(j,i,i,j,:)-z(j,i,i,:,j)"
	  "           -z(i,j,i,:,j)) + %g*(w(i,i,j,j,:) - z(k,k,i,:,i))",
	  mu, lambda);
  
  getfem::generic_assembly assem2(s);
  assem2.push_mi(mim);
  assem2.push_mf(mf_u());
  assem2.push_mf(mf_SD);
  assem2.push_data(U);
  assem2.push_vec(SD);

  double t0 = gmm::uclock_sec();
  assem2.assembly();
  cerr << " done : " << gmm::uclock_sec() - t0 << "\n";

#endif


  size_type N = mesh.dim();
  // The third component is set to zero.
  gmm::clear(gmm::sub_vector(SD, gmm::sub_slice(N-1,gmm::vect_size(SD)/N, N)));

  for (size_type i = 0; i < mf_SD.nb_dof(); ++i)
    SD[i] *= cutoff_C2(gmm::abs(ls.values(1)[i/N]), 0.1, 0.2);

  {
    GMM_TRACE2("Exporting shape derivative");
    getfem::stored_mesh_slice sl; 
    sl.build(mesh, getfem::slicer_half_space(base_node(0, 0, 0), base_node(0, 0, 1), 0), 3);
    getfem::vtk_export exp(datafilename + "_grad.vtk",
			   PARAM.int_value("VTK_EXPORT")==1);
    exp.exporting(sl); 
    exp.write_point_data(mf_SD, SD, "Shape Gradient");
    cout << "export done, you can view the data file with (for example)\n"
      "mayavi -d " << datafilename << "_grad.vtk -f "
      "ExtractVectorNorm -m BandedSurfaceMap -m Outline -m Glyph\n";
  }
}

/**************************************************************************/
/*  Level-set update.                                                     */
/**************************************************************************/

namespace getfem {

  template<typename MAT, typename VECT>
  void asm_transport_dc1
  (const MAT &M, const mesh_im &mim, const mesh_fem &mf,
   const mesh_fem &mf_data,
   const VECT &A, const mesh_region &rg = mesh_region::all_convexes()) {

    generic_assembly 
      assem("a=data$1(#2); M$1(#1,#1)+="
	    "comp(Base(#1).vBase(#2).Grad(#1))(:,i,j,:,j).a(i)");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mf_data);
    assem.push_data(A);
    assem.push_mat(const_cast<MAT&>(M));
    assem.assembly(rg);
  }

# define MDBRICK_TRANSPORT_DC1 932845
  
  /* "Transport" brick ( @f$ K = \int \theta.\nabla u v @f$ ). */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_transport_dc1
    : public mdbrick_abstract_linear_pde<MODEL_STATE> {
    
    TYPEDEF_MODEL_STATE_TYPES;
    
    mdbrick_parameter<VECTOR> theta_; /* vector field parameter */

    void proper_update_K(void) {
      asm_transport_dc1
	(this->K, this->mim, this->mf_u, theta_.mf(), theta_.get(),
	 this->mf_u.linked_mesh().get_mpi_region());
    }

  public :
    mdbrick_parameter<VECTOR> &theta(void) { return theta_; }
    const mdbrick_parameter<VECTOR> &theta(void) const { return theta_; }

    mdbrick_transport_dc1(const mesh_im &mim_, const mesh_fem &mf_u_)
      : mdbrick_abstract_linear_pde<MODEL_STATE>(mim_, mf_u_,
						 MDBRICK_TRANSPORT_DC1),
	theta_("theta", mf_u_.linked_mesh(), this) { }
  };


  /* compute 
      version == 0 : tanh(phi0) * grad(phin)/norm(grad(phin))
      version == 1 : tanh(phi0)
      version == 2 : (3tanh(phi0) - tanh(phin))/2
      version == 3 : abs(tanh(phi0))
  */
  template<typename VECT1> class transportdc_nonlinear_term 
    : public getfem::nonlinear_elem_term {
    
    const mesh_fem &mf;
    std::vector<scalar_type> phi0, phin;
    size_type N;
    base_vector coeff, Phi;
    base_matrix gradPhin;
    bgeot::multi_index sizes_;
    int version; 
    
  public:
    transportdc_nonlinear_term(const mesh_fem &mf_, const VECT1 &phi0_,
			       const VECT1 &phin_, int version_) 
      : mf(mf_), phi0(mf_.nb_basic_dof()), phin(mf_.nb_basic_dof()),
	N(mf_.linked_mesh().dim()), Phi(1),
	gradPhin(1,N), sizes_(1), version(version_) {
      mf.extend_vector(phi0_, phi0);
      mf.extend_vector(phin_, phin);
      sizes_[0] = short_type((version == 1) ? 1 : N);
    }
    
    const bgeot::multi_index &sizes(size_type) const { return sizes_; }
    
    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num();
      coeff.resize(mf.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector(phi0,
				gmm::sub_index(mf.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation(ctx, coeff, Phi, 1);
      gmm::copy(gmm::sub_vector(phin,
				gmm::sub_index(mf.ind_basic_dof_of_element(cv))), coeff);
      
      scalar_type sgnphi = gmm::sgn(Phi[0]);
      scalar_type no = ::tanh(Phi[0]);
     
      switch (version) {
      case 0 :
	ctx.pf()->interpolation_grad(ctx, coeff, gradPhin, 1);
	// gradPhin(0,0) = -0.7; gradPhin(0,1) = -0.7;  gradPhin(0,2) = 0.0; // for test to be suppressed
	no /= std::max(1e-6, gmm::mat_euclidean_norm(gradPhin));
	for (size_type i = 0; i < N; ++i) 
	  t[i] = no * gradPhin(0, i);
	break;
      case 1 : t[0] = no; break;
      case 2 :
	ctx.pf()->interpolation(ctx, coeff, Phi, 1);
	t[0] = no + gmm::abs(no)*0.5*(sgnphi - gmm::sgn(Phi[0]));
	break;
      case 3 : t[0] = gmm::abs(no); break;
      default : GMM_ASSERT1(false, "Oups...");
      }
    }
  };

  template<typename MAT, typename VECT>
  void asm_transport_dc2
  (const MAT &M, const mesh_im &mim, const mesh_fem &mf,
   const VECT &phi0, const VECT &phin,
   const mesh_region &rg = mesh_region::all_convexes()) {

    transportdc_nonlinear_term<VECT> nterm1(mf, phi0, phin, 0); 

    generic_assembly 
      assem("M$1(#1,#1)+=comp(Base(#1).NonLin$1(#1).Grad(#1))(:,i,:,i)");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_nonlinear_term(&nterm1);
    assem.push_mat(const_cast<MAT&>(M));
    assem.assembly(rg);
  }

  template<typename VECT1, typename VECT>
  void asm_transport_dc2_st
  (const VECT1 &V, const mesh_im &mim, const mesh_fem &mf,
   const VECT &phi0, const VECT &phin,
   const mesh_region &rg = mesh_region::all_convexes()) {

    transportdc_nonlinear_term<VECT> nterm(mf, phi0, phin, 2); 

    generic_assembly 
      assem("V$1(#1)+=comp(Base(#1).NonLin(#1))(:,1)");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(const_cast<VECT1&>(V));
    assem.assembly(rg);
  }


# define MDBRICK_TRANSPORT_DC2 3423456

  /* "Level-set regularization" brick. */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_transport_dc2
    : public mdbrick_abstract_linear_pde<MODEL_STATE> {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_parameter<VECTOR> phi0_, phin_; /* vector field parameter */

    T_MATRIX K0;
    bool K0init;
    scalar_type h;

    void proper_update_K(void) {
      if (!K0init) { // viscous term for regularization
	gmm::resize(K0, this->mf_u.nb_dof(), this->mf_u.nb_dof());
	gmm::clear(K0);
	generic_assembly 
	  assem("M$1(#1,#1)+=comp(Grad(#1).Grad(#1))(:,i,:,i)"); 
	assem.push_mi(this->mim);
	assem.push_mf(this->mf_u);
	assem.push_mat(K0);
	assem.assembly(this->mf_u.linked_mesh().get_mpi_region());
	gmm::scale(K0, h*0.1);
      }
      gmm::copy(K0, this->K);
      asm_transport_dc2
	(this->K, this->mim, this->mf_u, phi0_.get(), phin_.get(),
	 this->mf_u.linked_mesh().get_mpi_region());
    }

  public :
    mdbrick_parameter<VECTOR> &phi0(void) { return phi0_; }
    const mdbrick_parameter<VECTOR> &phi0(void) const { return phi0_; }
    mdbrick_parameter<VECTOR> &phin(void) { return phin_; }
    const mdbrick_parameter<VECTOR> &phin(void) const { return phin_; }

    mdbrick_transport_dc2(const mesh_im &mim_, const mesh_fem &mf_u_,
			  scalar_type h_)
      : mdbrick_abstract_linear_pde<MODEL_STATE>(mim_, mf_u_,
						 MDBRICK_TRANSPORT_DC2),
        phi0_("phi0", mf_u_.linked_mesh(), this),
	phin_("phin", mf_u_.linked_mesh(), this), 
	K0init(false), h(h_) { }
  };

}

void crack_problem::update_level_set(const plain_vector &SD) {


  for (size_type ntt = 0; ntt < 3; ++ntt) {
    
    // Level-set advance
    
    scalar_type h = min_h;
    cout << "min_h = " << h << endl;
    cout << "Norminf(SD) = " << gmm::vect_norminf(SD) << endl;
    scalar_type dt = 0.8 * h / gmm::vect_norminf(SD);
    
    plain_vector phi0 = ls.values(1), DF(gmm::vect_size(phi0));
    plain_vector DF0(gmm::vect_size(phi0));
    
    gmm::copy(ls.values(1), phi0);
    
    getfem::mdbrick_transport_dc1<> transport(standard_mim, ls.get_mesh_fem());
    transport.theta().set(mf_SD, gmm::scaled(SD, -1.0));
    
    getfem::mdbrick_dynamic<> dynamic(transport, 1.0);
    dynamic.set_dynamic_coeff(1.0/dt, 1.0);
    
    gmm::mult(dynamic.get_M(), gmm::scaled(phi0, 1.0/dt), DF);
    dynamic.set_DF(DF);
    
    getfem::standard_model_state MS(dynamic);
    gmm::iteration iter(residual, 1, 40000);
    getfem::standard_solve(MS, dynamic, iter);
    
    for (size_type i = 0; i < ls.get_mesh_fem().nb_dof(); ++i)
      ls.values(1)[i] = std::min(ls.values(1)[i], MS.state()[i]);
    
    // Level-set reinitialization -> get back a signed distance
    
    dt = 0.8 * h;
    
    gmm::copy(ls.values(1), phi0);
    plain_vector phin = phi0;
    
    getfem::mdbrick_transport_dc2<> transport2(standard_mim, ls.get_mesh_fem(),h);
    getfem::mdbrick_dynamic<> dynamic2(transport2, 1.0);
    transport2.phi0().set(ls.get_mesh_fem(), phi0);
    
    for (size_type nt = 0; nt < 20; ++nt) {
      
      gmm::clear(DF0);
      asm_transport_dc2_st
	(DF0, standard_mim, ls.get_mesh_fem(), phi0, phin,
	 ls.get_mesh_fem().linked_mesh().get_mpi_region());
      transport2.phin().set(ls.get_mesh_fem(), phin);
      dynamic2.set_dynamic_coeff(1.0/dt, 1.0);
      gmm::mult(dynamic2.get_M(), gmm::scaled(phin, 1.0/dt), DF0, DF);
      dynamic2.set_DF(DF);
      
      if (0) {
	// explicit scheme
	gmm::iteration iter2(residual, 0, 40000);
	plain_vector DF1(gmm::vect_size(phi0));
	gmm::mult(transport2.get_K(), gmm::scaled(phin, -1.0), DF1);
	gmm::add(DF1, DF);
	gmm::cg(dynamic2.get_M(), phin, gmm::scaled(DF, dt),
		gmm::identity_matrix(), iter2);
      }
      else {
	// implicit scheme
	gmm::iteration iter2(residual, 1, 40000);
	getfem::standard_model_state MS2(dynamic2);
	getfem::standard_solve(MS2, dynamic2, iter2);
	gmm::copy(MS2.state(), phin);
      }
    }
    gmm::copy(phin, ls.values(1));
  
  }


  {
    GMM_TRACE2("Exporting level set");
    getfem::stored_mesh_slice sl; 
    sl.build(mesh, getfem::slicer_half_space(base_node(0, 0, 0),
					     base_node(0, 0, 1), 0), 3);
    getfem::vtk_export exp(datafilename + "_ls.vtk",
			   PARAM.int_value("VTK_EXPORT")==1);
    exp.exporting(sl); 
    exp.write_point_data(ls.get_mesh_fem(), ls.values(1), "Level_set");
    cout << "export done, you can view the data file with (for example)\n"
      "mayavi -d " << datafilename << "_ls.vtk "
      "-m BandedSurfaceMap -m Outline -f WarpScalar -m IsoSurface\n";
  }

  ls.touch();
  mls.adapt();
  mim.adapt();
  mim_crack.adapt();
 
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool crack_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();

//   cout << "testing mims..\n";
//   for (dal::bv_visitor cv(mim.linked_mesh().convex_index()); !cv.finished(); ++cv) {
//     base_node G = gmm::mean_value(mim.linked_mesh().points_of_convex(cv));
//     cerr << "on cv\t" << cv << " nb_integ_pts = " << mim.int_method_of_element(cv)->approx_method()->nb_points_on_convex() << "\t , at " << G << "\t";
//     scalar_type err = getfem::test_integration_error(mim.int_method_of_element(cv)->approx_method(), 2);
//     cout << " max integ error = " << err << "\n";
//     assert(err < 1e-10);
//   }


//   getfem::mesh_fem mf(mim.linked_mesh()); mf.set_classical_finite_element(0);
//   scalar_type vol1 = gmm::sqr(getfem::asm_L2_norm(mim, mf, std::vector<double>(mf.nb_dof(), 1.0)));
//   scalar_type surf1 = gmm::sqr(getfem::asm_L2_norm(mim, mf, std::vector<double>(mf.nb_dof(), 1.0), NEUMANN_BOUNDARY_NUM));
//   scalar_type surfcrack =  gmm::sqr(getfem::asm_L2_norm(mim_crack, mf, std::vector<double>(mf.nb_dof(), 1.0)));


//   cout << "surface of the crack : " << surfcrack << endl;
//   cout << "volume of mesh is " << vol1 << " surf = " << surf1 << "\n";
//   getfem::mesh_im mim2(mim.linked_mesh()); mim2.set_integration_method(mim.linked_mesh().convex_index(), 6);
//   scalar_type vol2 = gmm::sqr(getfem::asm_L2_norm(mim2, mf, std::vector<double>(mf.nb_dof(), 1.0)));
//   scalar_type surf2 = gmm::sqr(getfem::asm_L2_norm(mim2, mf, std::vector<double>(mf.nb_dof(), 1.0), NEUMANN_BOUNDARY_NUM));
//   cout << "volume of mesh is " << vol2 << " surf = " << surf2 << "\n";
//   assert(gmm::abs(vol1-vol2) < 1e-5);
//   assert(gmm::abs(surf1-surf2) < 1e-5);

  
  mfls_u.adapt();
  std::vector<getfem::pglobal_function> vfunc(4);
  for (unsigned j=0; j < 4; ++j) {
    getfem::crack_singular_xy_function *s = 
      new getfem::crack_singular_xy_function(j);
    vfunc[j] = getfem::global_function_on_level_set(ls, *s);
  }
  
  mf_sing_u.set_functions(vfunc);


//   for (dal::bv_visitor i(mim_crack.convex_index()); !i.finished(); ++i)
//     if (mls.is_convex_cut(i)) {
//       cerr << i << " primary zone: " << mls.primary_zone_of_convex(i) << " sz.size=" << mls.zoneset_of_convex(i).size() << " " << gmm::mean_value(mls.linked_mesh().points_of_convex(i)) << "\n";
//       if (mls.zoneset_of_convex(i).size() == 1) {
//       }
//     }
  
  dal::bit_vector crack_tip_convexes = mls.crack_tip_convexes();
  cerr << "crack_tip_convexes = " << crack_tip_convexes << "\n";

  switch (enrichment_option) {
  case 1 :
    {
      dal::bit_vector enriched_dofs;
      plain_vector X(mf_partition_of_unity.nb_dof());
      plain_vector Y(mf_partition_of_unity.nb_dof());
      getfem::interpolation(ls.get_mesh_fem(), mf_partition_of_unity,
			    ls.values(1), X);
      getfem::interpolation(ls.get_mesh_fem(), mf_partition_of_unity,
			    ls.values(0), Y);
      for (size_type j = 0; j < mf_partition_of_unity.nb_dof(); ++j) {
	if (gmm::sqr(X[j]) + gmm::sqr(Y[j]) <= gmm::sqr(enr_area_radius))
	  enriched_dofs.add(j);
      }


      for (dal::bv_visitor cv(crack_tip_convexes); !cv.finished(); ++cv) {
	//cerr << "inspecting convex centered at " << gmm::mean_value(mf_partition_of_unity.linked_mesh().points_of_convex(cv)) << "\n";
	for (unsigned j=0; j < mf_partition_of_unity.nb_basic_dof_of_element(cv); ++j) {
	  size_type d = mf_partition_of_unity.ind_basic_dof_of_element(cv)[j];
	  if (!enriched_dofs.is_in(d)) {
	    enriched_dofs.add(d);
	    /*cerr << "added supplementary cracktip dof " << d << ", distance is:"
	      << sqrt(gmm::sqr(X[d]) + gmm::sqr(Y[d])) << "\n";*/
	  }
	}
      }

      /*for (dal::bv_visitor ii(enriched_dofs); !ii.finished(); ++ii) {
	cerr << "enriched_dof " << ii << "  @ " << mf_partition_of_unity.point_of_dof(ii) << "\n";
	}*/
      cerr << "total enriched cracktip dof: " << enriched_dofs.card() << "\n\n";


      if (enriched_dofs.card() < 3)
	GMM_WARNING0("There is " << enriched_dofs.card() <<
		     " enriched dofs for the crack tip");
      mf_product.set_enrichment(enriched_dofs);
      mf_u_sum.set_mesh_fems(mf_product, mfls_u);
    }
    break;
    default : mf_u_sum.set_mesh_fems(mfls_u); break;
  }

  U.resize(mf_u().nb_dof());

  cout << "Number of dof for u: " << mf_u().nb_dof() << endl;

  // Linearized elasticity brick.
  getfem::mdbrick_isotropic_linearized_elasticity<> ELAS(mim, mf_u(),
							 lambda, mu);

  // Defining the volumic source term.
  plain_vector F(nb_dof_rhs * N);
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    gmm::copy(sol_f(mf_rhs.point_of_basic_dof(i)),
	      gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  // Volumic source term brick.
  getfem::mdbrick_source_term<> VOL_F(ELAS, mf_rhs, F);

  // Defining the Neumann condition right hand side.
  gmm::clear(F);

  // Neumann condition brick.
  if (!is_prescribed_disp) {
    base_small_vector f(N); f[N-1] = neumann_force;
    for (size_type i = 0; i < nb_dof_rhs; ++i)
      if (mf_rhs.point_of_basic_dof(i)[0] > 0.8 * Pmax[0])
	gmm::copy(f, gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  }
  getfem::mdbrick_source_term<> NEUMANN(VOL_F, mf_rhs, F, NEUMANN_BOUNDARY_NUM);
  
  gmm::clear(F);
  if (is_prescribed_disp) {
    for (size_type i = 0; i < nb_dof_rhs; ++i)
      if (mf_rhs.point_of_basic_dof(i)[0] > 0.4 * Pmax[0])
	if (mf_rhs.point_of_basic_dof(i)[2] >= Pmax[2]/2.0)
	  F[i*N+2] = prescribed_disp;
  }
  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<> final_model(NEUMANN,
					  DIRICHLET_BOUNDARY_NUM, mf_mult);
  final_model.rhs().set(mf_rhs, F);  
  final_model.set_constraints_type(getfem::constraints_type(dir_with_mult));

  // Generic solve.
  size_type nnb = final_model.nb_dof();
  cout << "Total number of variables : " << nnb << endl;
  getfem::standard_model_state MS(final_model);
  gmm::iteration iter(residual, 1, 40000);
  getfem::standard_solve(MS, final_model, iter);

  // Solution extraction
  gmm::copy(ELAS.get_solution(MS), U);

  // Energy computation
  double rupt_e = fracture_energy();
  double elas_e = gmm::vect_sp(ELAS.get_K(), U, U) * 0.5;
  double rhs_e = - gmm::vect_sp(VOL_F.get_F(), U)
    - gmm::vect_sp(NEUMANN.get_F(), U);
  cout << "Potential energy = " <<  elas_e + rhs_e
       << " Fracture energy = " << rupt_e
       << " Total energy = " << elas_e + rhs_e + rupt_e << endl;

  return (iter.converged());
}
  
/**************************************************************************/
/*  Export U.                                                             */
/**************************************************************************/

void crack_problem::save_U(plain_vector &U) {
	
  cout << "Post processing...\n";
  getfem::mesh mcut;
  mls.global_cut_mesh(mcut);
  //       getfem::mesh_fem mf(mcut, mf_u().get_qdim());
  //       mf.set_classical_discontinuous_finite_element(2, 0.001);
  
  //       plain_vector V(mf.nb_dof());
  //       getfem::interpolation(mf_u(), mf, U, V);
  
  getfem::stored_mesh_slice sl;
  getfem::mesh mcut_refined;
  //sl.build(mcut, /*getfem::slicer_boundary(mcut), */getfem::slicer_build_mesh(mcut_refined), 2);
  
  
  getfem::mesh_slicer slicer(mcut);
  getfem::slicer_build_stored_mesh_slice sbuild(sl);
  //getfem::slicer_build_mesh sbmesh(mcut_refined);
  slicer.push_back_action(sbuild);
  //slicer.push_back_action(sbmesh);
  
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mcut, border_faces);
  // merge outer_faces with faces of the crack
  for (getfem::mr_visitor i(mcut.region(0)); !i.finished(); ++i)
    border_faces.add(i.cv(), i.f());
  
  
  getfem::mesh_fem mfcut(mcut, mf_u().get_qdim());
  getfem::mesh_fem mfcut_vm(mcut, 1);
  mfcut.set_classical_discontinuous_finite_element(3, 0.001);
  mfcut_vm.set_classical_discontinuous_finite_element(2, 0.001);
  
  dal::bit_vector reflst;
  
  /* choose an adequate slice refinement based on the distance to the
     crack tip */
  unsigned nn = 1;
  std::vector<short_type> nrefine(mcut.convex_index().last_true()+1);
  for (dal::bv_visitor cv(mcut.convex_index()); !cv.finished(); ++cv) {
    scalar_type dmin=0, d;
    base_node PPmin,P;
    for (unsigned i=0; i < mcut.nb_points_of_convex(cv); ++i) {
      P = mcut.points_of_convex(cv)[i];
      d = gmm::vect_norm2(ls_function(P));
      if (d < dmin || i == 0) { dmin = d; PPmin = P; }
    }
    
    if (dmin < 1e-5) { nrefine[cv] = short_type(nn*4); reflst.add(cv); }
    else if (dmin < .1) nrefine[cv] = short_type(nn*2);
    else nrefine[cv] = short_type(nn);
  }
  
  mfcut.set_classical_discontinuous_finite_element(reflst, 6, 0.001);
  mfcut_vm.set_classical_discontinuous_finite_element(reflst, 5, 0.001);
  
  slicer.exec(nrefine, border_faces);
  
  plain_vector V(mfcut.nb_dof()), VonMises(mfcut_vm.nb_dof());
  getfem::interpolation(mf_u(), mfcut, U, V);
  getfem::interpolation_von_mises(mfcut, mfcut_vm, V, VonMises, mu);
  
  //getfem::interpolation(mfcut, sl, U, W);
  
  if (PARAM.int_value("VTK_EXPORT")) {
    cout << "export to " << datafilename + ".vtk" << "..\n";
    getfem::vtk_export exp(datafilename + ".vtk",
			   PARAM.int_value("VTK_EXPORT")==1);
    exp.exporting(sl); 
    exp.write_point_data(mfcut_vm, VonMises, "Von Mises Stress");
    exp.write_point_data(mfcut, V, "elastostatic_displacement");
    cout << "export done, you can view the data file with (for example)\n"
      "mayavi -d " << datafilename << ".vtk -f "
      "WarpVector -m BandedSurfaceMap -m Outline\n";
  }
  
  /*getfem::mesh_fem mf_refined(mcut_refined, mf_u().get_qdim());
    mf_refined.set_classical_discontinuous_finite_element(1, 0.001);
    plain_vector W(mf_refined.nb_dof());
    getfem::interpolation(mf_u(), mf_refined, U, W);
    
    if (PARAM.int_value("VTK_EXPORT")) {
    cout << "export to " << datafilename + ".vtk" << "..\n";
    getfem::vtk_export exp(datafilename + ".vtk",
    PARAM.int_value("VTK_EXPORT")==1);
    exp.exporting(mf_refined); 
    exp.write_point_data(mf_refined, W, "elastostatic_displacement");
    cout << "export done, you can view the data file with (for example)\n"
    "mayavi -d " << datafilename << ".vtk -f ExtractVectorNorm -f "
    "WarpVector -m BandedSurfaceMap -m Outline\n";
    }
  */
  
}

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

#endif


int main(
#if 0     
         int argc, char *argv[]
#endif
         ) {
#if 0

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  //getfem::getfem_mesh_level_set_noisy();

  try {
    crack_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh");
    plain_vector U(p.mf_u().nb_dof());

    // Level-set initialization
    p.ls.reinit();
    for (size_type d = 0; d < p.ls.get_mesh_fem().nb_dof(); ++d) {
      const base_node P = p.ls.get_mesh_fem().point_of_basic_dof(d);
      p.ls.values(0)[d] = ls_function(P)[0];
      p.ls.values(1)[d] = ls_function(P)[1];
      //+ 0.3 - gmm::sqr((ls.get_mesh_fem().point_of_dof(d))[1])*2.0;
    }
    p.ls.touch();
    p.mls.adapt();
    p.mim.adapt();
    p.mim_crack.adapt();

    p.solve(U);
    p.save_U(U);
    cout << " U saved" << endl;

    for (size_type it = 0; it < 100; ++it) {
      
      plain_vector SD(p.mf_SD.nb_dof());
      p.shape_derivative(U, SD);
      p.update_level_set(SD);
      p.solve(U);
      p.save_U(U);

      getchar();
    }
  }
  GMM_STANDARD_CATCH_ERROR;

#endif

  return 0; 
}
