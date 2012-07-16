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
 * Goal : scalar Dirichlet problem with Xfem.
 *
 * Research program.
 */

#include "getfem/getfem_assembling.h" /* import assembly methods      */
#include "getfem/getfem_export.h"     /* export functions             */
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_partial_mesh_fem.h"
#include "getfem/getfem_Coulomb_friction.h"
#include "getfem/getfem_import.h"
#include "getfem/getfem_inter_element.h"
#include "gmm/gmm.h"

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;


/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_vector;
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::short_type;  /* = short */
using bgeot::size_type;   /* = unsigned long */
using bgeot::dim_type;
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

typedef gmm::row_matrix<sparse_vector> sparse_row_matrix;

int Nxy;

int DIRICHLET_BOUNDARY_NUM = 0;
int NEUMANN_BOUNDARY_NUM = 1;

scalar_type nu;
double Radius;
int u_version;


/*
 * Exact solution for the velocity u
 */

  base_small_vector u_exact(const base_node &p){  
	  return base_small_vector	(							 
			p[0]*p[0]*(p[0]-1.0)*(p[0]-1.0)*(2.0*p[1]-6.0*p[1]*p[1]+4.0*p[1]*p[1]*p[1]),
		   -p[1]*p[1]*(p[1]-1.0)*(p[1]-1.0)*(2.0*p[0]-6.0*p[0]*p[0]+4.0*p[0]*p[0]*p[0]) ) ;  
  }
  

/* 
 * Exact solution for the pressure p
 */

  double p_exact(const base_node &p) {
      return p[0]*p[0]*p[1]*p[1];
  }

/*
 * Value of the normal derivative on the level-set. The fluid is outside
 */

  base_small_vector g_exact_ls(const base_node &p){
      double R = Radius;
      return base_small_vector( (-1.0/(2.0*R))*2.0*p[0]*(2.0*p[0]-1.0)*(2.0*p[0]-1.0)*(p[0]-1.0)*2.0*p[1]*(2.0*p[1]-1.0)*(p[1]-1.0)-(1.0/(2.0*R))*(2.0*p[1]-1.0)*p[0]*p[0]*(1.0-p[0])*(1.0-p[0])*(2.0-12.0*p[1]+12.0*p[1]*p[1]),
          (1.0/(2.0*R))*(2*p[0]-1.0)*p[1]*p[1]*(p[1]-1.0)*(p[1]-1.0)*(2.0-12.0*p[0]+12.0*p[0]*p[0])+(1.0/(2.0*R))*2.0*p[1]*(2.0*p[1]-1.0)*(2.0*p[1]-1.0)*(p[1]-1.0)*2.0*p[0]*(2.0*p[0]-1.0)*(p[0]-1.0) ) ;
  }

/*
 * Value of the normal derivative on the box boundary.
 */
  base_small_vector g_exact(const base_node &p){
      
      if (p[0] < 1e-7) {
        return base_small_vector( 0.0, -2.0*p[1]*p[1]*(p[1]-1.0)*(p[1]-1.0)
               ) ;
      }
      if (p[0] > 1.0 - 1e-7) {
        return base_small_vector( 0.0, 2.0*p[1]*p[1]*(p[1]-1.0)*(p[1]-1.0)
               ) ;
      }
      if (p[1] > 1.0 - 1e-7) {
        return base_small_vector( -2.0*p[0]*p[0]*(p[0]-1.0)*(p[0]-1.0), 0.0
               ) ;
      }
      if (p[1] < 1e-7) {
        return base_small_vector( 2.0*p[0]*p[0]*(p[0]-1.0)*(p[0]-1.0), 0.0
               ) ;
      }
      else {
        return base_small_vector( 0.0, 0.0
               ) ;
      }   
  }

/*
 * Right-hand-side corresponding to the exact solution
 */

  base_small_vector rhs(const base_node &p){
    return  base_small_vector( -nu*( (4*(2*p[1]-1))*(6*p[1]*p[1]*p[0]*p[0]+p[1]*p[1]-6*p[1]*p[1]*p[0]-6*p[0]*p[0]*p[1]-p[1]+6*p[0]*p[1]-6*p[0]*p[0]*p[0]+3*p[0]*p[0]+3*p[0]*p[0]*p[0]*p[0])   ) + 2.0*p[0]*p[1]*p[1],
							  -nu*( -(4*(-1+2*p[0]))*(6*p[1]*p[1]*p[0]*p[0]+p[0]*p[0]-6*p[0]*p[0]*p[1]-6*p[1]*p[1]*p[0]+6*p[0]*p[1]-p[0]+3*p[1]*p[1]*p[1]*p[1]-6*p[1]*p[1]*p[1]+3*p[1]*p[1])  ) + 2.0*p[1]*p[0]*p[0] );
  }

/*
 * The level-set : a circle of radius R = Radius
 */

  double ls_value(const base_node &p) {
      double R = Radius;
      return  (p[0]-0.5)*(p[0]-0.5) + (p[1]-0.5)*(p[1]-0.5)- R*R;
  }
  

/*
 * Test procedure
 */

void test_mim(getfem::mesh_im_level_set &mim, getfem::mesh_fem &mf_rhs,
	      bool bound) {
  if (!u_version) {
    unsigned N = unsigned(mim.linked_mesh().dim());
    size_type nbdof = mf_rhs.nb_dof();
    plain_vector V(nbdof), W(1);
    std::fill(V.begin(), V.end(), 1.0);
    
    getfem::generic_assembly assem("u=data(#1); V()+=comp(Base(#1))(i).u(i);");
    assem.push_mi(mim); assem.push_mf(mf_rhs); assem.push_data(V);
    assem.push_vec(W);
    assem.assembly(getfem::mesh_region::all_convexes());
    double exact(0), R2 = Radius*Radius, R3 = R2*Radius;
    switch (N) {
      case 1: exact = bound ? 1.0 : 2.0*Radius; break;
      case 2: exact = bound ? Radius*M_PI : R2*M_PI; break;
      case 3: exact = bound ? 2.0*M_PI*R2 : 4.0*M_PI*R3/3.0; break;
      default: assert(N <= 3);
    }
    if (bound) cout << "Boundary length: "; else cout << "Area: ";
    cout << W[0] << " should be " << exact << endl;
    assert(gmm::abs(exact-W[0])/exact < 0.01); 
  }
}

/*
 * Assembly of stabilization terms
 */


template<typename VECT1> class level_set_unit_normal 
  : public getfem::nonlinear_elem_term {
  const getfem::mesh_fem &mf;
  const VECT1 &U;
  size_type N;
  base_matrix gradU;
  bgeot::base_vector coeff;
  bgeot::multi_index sizes_;
public:
  level_set_unit_normal(const getfem::mesh_fem &mf_, const VECT1 &U_) 
    : mf(mf_), U(U_), N(mf_.linked_mesh().dim()), gradU(1, N)
  { sizes_.resize(1); sizes_[0] = short_type(N); }
  const bgeot::multi_index &sizes() const {  return sizes_; }
  virtual void compute(getfem::fem_interpolation_context& ctx,
		       bgeot::base_tensor &t) {
    size_type cv = ctx.convex_num();
    coeff.resize(mf.nb_dof_of_element(cv));
    gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf.ind_dof_of_element(cv))),
	      coeff);
    ctx.pf()->interpolation_grad(ctx, coeff, gradU, 1);
    scalar_type norm = gmm::vect_norm2(gmm::mat_row(gradU, 0));
    for (size_type i = 0; i < N; ++i) t[i] = - gradU(0, i) / norm;
    // cout << "point " << ctx.xreal() << " norm = " << t << endl;
  }
};


template<class MAT>
void asm_stabilization_mixed_term
(const MAT &RM_, const getfem::mesh_im &mim, const getfem::mesh_fem &mf, 
 const getfem::mesh_fem &mf_mult, getfem::level_set &ls,
 const getfem::mesh_region &rg = getfem::mesh_region::all_convexes()) {
  MAT &RM = const_cast<MAT &>(RM_);

  level_set_unit_normal<std::vector<scalar_type> >
    nterm(ls.get_mesh_fem(), ls.values());

  getfem::generic_assembly assem("t=comp(vBase(#2).vGrad(#1).NonLin(#3));"
				 "M(#2, #1)+= t(:,j,:,j,i,i)");
  assem.push_mi(mim);
  assem.push_mf(mf);
  assem.push_mf(mf_mult);
  assem.push_mf(ls.get_mesh_fem());
  assem.push_mat(RM);
  assem.push_nonlinear_term(&nterm);
  assem.assembly(rg);
}


template<class MAT>
void asm_stabilization_symm_term
(const MAT &RM_, const getfem::mesh_im &mim, const getfem::mesh_fem &mf, 
 getfem::level_set &ls,
 const getfem::mesh_region &rg = getfem::mesh_region::all_convexes()) {
  MAT &RM = const_cast<MAT &>(RM_);

  level_set_unit_normal<std::vector<scalar_type> >
    nterm(ls.get_mesh_fem(), ls.values());

  getfem::generic_assembly
    assem("t=comp(vGrad(#1).NonLin(#2).vGrad(#1).NonLin(#2));"
	  "M(#1, #1)+= sym(t(:,k,i,i,:,k,j,j))");
  assem.push_mi(mim);
  assem.push_mf(mf);
  assem.push_mf(ls.get_mesh_fem());
  assem.push_mat(RM);
  assem.push_nonlinear_term(&nterm);
  assem.assembly(rg);
}





template<class MAT>
void asm_mass_matrix_pl(const MAT &RM_, const getfem::mesh_im &mim,  const getfem::mesh_fem &mf_p, const getfem::mesh_fem &mf, 
						 getfem::level_set &ls,
						 const getfem::mesh_region &rg = getfem::mesh_region::all_convexes()) {
	MAT &RM = const_cast<MAT &>(RM_);
	
	level_set_unit_normal<std::vector<scalar_type> >
    nterm(ls.get_mesh_fem(), ls.values());
	
	getfem::generic_assembly
    assem("t=comp(Base(#1).vBase(#2).NonLin(#3));"
		  "M(#1, #2)+= t(:,:,i,i)");
	assem.push_mi(mim);
	assem.push_mf(mf_p);
	assem.push_mf(mf);
	assem.push_mf(ls.get_mesh_fem());
	assem.push_mat(RM);
	assem.push_nonlinear_term(&nterm);
	assem.assembly(rg);
}
					


template<class MAT>
void asm_mass_matrix_up(const MAT &RM_, const getfem::mesh_im &mim, const getfem::mesh_fem &mf_p, const getfem::mesh_fem &mf, 
					getfem::level_set &ls,
					const getfem::mesh_region &rg = getfem::mesh_region::all_convexes()) {
	MAT &RM = const_cast<MAT &>(RM_);
	
	level_set_unit_normal<std::vector<scalar_type> >
    nterm(ls.get_mesh_fem(), ls.values());
	
	getfem::generic_assembly
    assem("t=comp(Base(#1).NonLin(#3).NonLin(#3).vGrad(#2));"
		  "M(#1, #2)+= t(:,i,j,:,i,j)");
	assem.push_mi(mim);
	assem.push_mf(mf_p);
	assem.push_mf(mf);
	assem.push_mf(ls.get_mesh_fem());
	assem.push_mat(RM);
	assem.push_nonlinear_term(&nterm);
	assem.assembly(rg);
}
				   
				   
























/*
 * Elementary extrapolation matrices
 */

void compute_mass_matrix_extra_element(base_matrix &M, const getfem::mesh_im &mim, const getfem::mesh_fem &mf,
 size_type cv1, size_type cv2) {

  getfem::pfem pf1_old = 0;
  static getfem::pfem_precomp pfp1 = 0;
  static getfem::papprox_integration pai1_old = 0;
  bgeot::geotrans_inv_convex gic;
  bgeot::base_tensor t1, t2;
  getfem::base_matrix G1, G2;
  
  const getfem::mesh &m(mf.linked_mesh());
  
  GMM_ASSERT1(mf.convex_index().is_in(cv1) && mim.convex_index().is_in(cv1) &&
	      mf.convex_index().is_in(cv2) && mim.convex_index().is_in(cv2),
	      "Bad element");
    
  bgeot::pgeometric_trans pgt1 = m.trans_of_convex(cv1);
  getfem::pintegration_method pim1 = mim.int_method_of_element(cv1);
  getfem::papprox_integration pai1 =
    getfem::get_approx_im_or_fail(pim1);
  getfem::pfem pf1 = mf.fem_of_element(cv1);
  size_type nbd1 = pf1->nb_dof(cv1);
  
  if (pf1 != pf1_old || pai1 != pai1_old) {
    pfp1 = fem_precomp(pf1, &pai1->integration_points(), pim1);
    pf1_old = pf1; pai1_old = pai1;
  }
  
  bgeot::vectors_to_base_matrix(G1, m.points_of_convex(cv1));
  getfem::fem_interpolation_context ctx1(pgt1, pfp1, 0, G1, cv1,size_type(-1));
  
  getfem::pfem pf2 = mf.fem_of_element(cv2);
  size_type nbd2 = pf1->nb_dof(cv2);
  base_node xref2(pf2->dim());
  bgeot::pgeometric_trans pgt2 = m.trans_of_convex(cv2);
  gic.init(m.points_of_convex(cv2), pgt2);

  gmm::resize(M, nbd1, nbd2); gmm::clear(M);

  bgeot::vectors_to_base_matrix(G2, m.points_of_convex(cv2));
  
  getfem::fem_interpolation_context ctx2(pgt2, pf2, base_node(pgt2->dim()),
					 G2, cv2, size_type(-1));

  for (unsigned ii=0; ii < pai1->nb_points_on_convex(); ++ii) {
    ctx1.set_ii(ii);
    scalar_type coeff = pai1->integration_coefficients()[ii] * ctx1.J();
    bool converged;
    gic.invert(ctx1.xreal(), xref2, converged);
    GMM_ASSERT1(converged, "geometric transformation not well inverted ... !");
    // cout << "xref2 = " << xref2 << endl;
    ctx2.set_xref(xref2);

    pf1->real_base_value(ctx1, t1);
    pf2->real_base_value(ctx2, t2);
    
    for (size_type i = 0; i < nbd1; ++i)
      for (size_type j = 0; j < nbd2; ++j)
	M(i,j) += t1[i] * t2[j] * coeff;
  }
  
}



/* 
 * Main program 
 */

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.
    
  // Read parameters.
  bgeot::md_param PARAM;
  PARAM.read_command_line(argc, argv);
  u_version = int(PARAM.int_value("EXACT_SOL", "Which exact solution"));
  nu = double(PARAM.real_value("NU", "Viscosité"));
  
  
  // Load the mesh
  getfem::mesh mesh;
 
   std::string MESH_FILE = PARAM.string_value("MESH_FILE", "Mesh file");
   getfem::import_mesh(MESH_FILE, mesh);

  /* //to make tables for results
  Nxy  = PARAM.int_value("Nxy","Number of subdivision for the mesh in x and y direction");
  std::ostringstream oss;
  oss << Nxy;
  std::string NxyStr = oss.str();
  std::string meshname1 = "structured:GT=\"GT_PK(2,1)\";SIZES=[1,1];NOISED=0;";
    //(PARAM.string_value("MESHNAME1", "Nom du fichier de maillage"));
  std::string meshname2
    (PARAM.string_value("MESHNAME2", "Nom du fichier de maillage"));
  std::string meshname = meshname1 +  "NSUBDIV=[" + NxyStr + "," + NxyStr + "];" +  meshname2;
  getfem::import_mesh(meshname, mesh);
  */
  unsigned N = mesh.dim();
  
  // center the mesh in (0, 0).
  base_node Pmin(N), Pmax(N);
  mesh.bounding_box(Pmin, Pmax);
  Pmin += Pmax; Pmin /= -2.0;
  // Pmin[N-1] = -Pmax[N-1];
  
  // mesh.translation(Pmin);
  
  scalar_type h = mesh.minimal_convex_radius_estimate();
  cout << "h = " << h << endl;
  
/*
 * Level set definition
 */

  unsigned lsdeg = unsigned(PARAM.int_value("LEVEL_SET_DEGREE", "level set degree"));
  bool simplify_level_set = 
    (PARAM.int_value("SIMPLIFY_LEVEL_SET",
		     "simplification or not of the level sets") != 0);
  Radius = PARAM.real_value("RADIUS", "Domain radius");
  getfem::level_set ls(mesh, dim_type(lsdeg));


  //getfem::level_set lsup(mesh, dim_type(lsdeg), true), lsdown(mesh, dim_type(lsdeg), true);
  const getfem::mesh_fem &lsmf = ls.get_mesh_fem();
  for (unsigned i = 0; i < lsmf.nb_dof(); ++i) {
    ls.values()[i] = ls_value(lsmf.point_of_dof(i));
	  }
 
	
  if (simplify_level_set) {
    scalar_type simplify_rate = std::min(0.03, 0.05 * sqrt(h));
    cout << "Simplification of level sets with rate: " <<
      simplify_rate << endl;
    ls.simplify(simplify_rate);
 //   lsup.simplify(simplify_rate);
 //   lsdown.simplify(simplify_rate); 
  }
  
	getfem::mesh_level_set mls(mesh); //, mlsup(mesh), mlsdown(mesh);
  mls.add_level_set(ls);
  mls.adapt();
 // mlsup.add_level_set(lsup);
 // mlsup.adapt();
 // mlsdown.add_level_set(lsdown);
 // mlsdown.adapt();

  
  getfem::mesh mcut;
  mls.global_cut_mesh(mcut);
  mcut.write_to_file("cut.mesh");
	
/*	getfem::mesh mcutdown;
	mlsdown.global_cut_mesh(mcutdown);
	mcutdown.write_to_file("cutdown.mesh");
	getfem::mesh mcutup;
	mlsup.global_cut_mesh(mcutup);
	mcutup.write_to_file("cutup.mesh");
	
*/	
	
	
	
  
  // Integration method on the domain
  std::string IM = PARAM.string_value("IM", "Mesh file");
  std::string IMS = PARAM.string_value("IM_SIMPLEX", "Mesh file");
  int intins = getfem::mesh_im_level_set::INTEGRATE_OUTSIDE;
  getfem::mesh_im uncutmim(mesh);
  uncutmim.set_integration_method(mesh.convex_index(),
				  getfem::int_method_descriptor(IM));
  getfem::mesh_im_level_set mim(mls, intins,
				getfem::int_method_descriptor(IMS));
  mim.set_integration_method(mesh.convex_index(),
			     getfem::int_method_descriptor(IM));
  mim.adapt();
  
  
  // Integration methods on the boudary
  int intbound = getfem::mesh_im_level_set::INTEGRATE_BOUNDARY;
  getfem::mesh_im_level_set mimbound(mls, intbound,
					 getfem::int_method_descriptor(IMS));
  mimbound.set_integration_method(mesh.convex_index(),
				      getfem::int_method_descriptor(IM));
  mimbound.adapt();
 // getfem::mesh_im_level_set mimboundup(mlsup, intbound,
//				       getfem::int_method_descriptor(IMS));
 // mimboundup.set_integration_method(mesh.convex_index(),
//				    getfem::int_method_descriptor(IM));
 // mimboundup.adapt();

  
  // Finite element method for the unknown velocity u
  getfem::mesh_fem pre_mf(mesh);
  std::string FEM = PARAM.string_value("FEM", "finite element method");
  pre_mf.set_qdim(N);
  pre_mf.set_finite_element(mesh.convex_index(),
			    getfem::fem_descriptor(FEM));
	
  getfem::partial_mesh_fem mf(pre_mf);
	
  dal::bit_vector kept_dof = select_dofs_from_im(pre_mf, mim);
	
  dal::bit_vector rejected_elt;
  for (dal::bv_visitor cv(mim.convex_index()); !cv.finished(); ++cv)
    if (mim.int_method_of_element(cv) == getfem::im_none())
      rejected_elt.add(cv);
  mf.adapt(kept_dof, rejected_elt);
  size_type nb_dof = mf.nb_dof();
  cout << "nb_dof = " << nb_dof << endl;
  
  // Finite element method for the pressure p
  getfem::mesh_fem pre_mf_p(mesh);
  std::string FEM_p = PARAM.string_value("FEM_p", "fem for pressure");
  pre_mf_p.set_finite_element(mesh.convex_index(),
				 getfem::fem_descriptor(FEM_p));
  getfem::partial_mesh_fem mf_p(pre_mf_p);
  dal::bit_vector kept_dof_p = select_dofs_from_im(pre_mf_p, mim);
  dal::bit_vector rejected_elt_p;
  for (dal::bv_visitor cv(mim.convex_index()); !cv.finished(); ++cv)
    if (mim.int_method_of_element(cv) == getfem::im_none())
      rejected_elt_p.add(cv);
  mf_p.adapt(kept_dof_p, rejected_elt_p);
  size_type nb_dof_p = mf_p.nb_dof();
  cout << "nb_dof_p = " << nb_dof_p << endl;
  
  // Finite element method for the rhs
  getfem::mesh_fem mf_rhs(mesh);
  std::string FEMR = PARAM.string_value("FEM_RHS", "finite element method");
  mf_rhs.set_qdim(N);
  mf_rhs.set_finite_element(mesh.convex_index(),
			    getfem::fem_descriptor(FEMR));
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  cout << "nb_dof_rhs = " << nb_dof_rhs << endl;
  
  // A P0 finite element
  const getfem::mesh_fem &mf_P0 = getfem::classical_mesh_fem(mesh, 0);
  
  // Finite element method for the multipliers
  getfem::mesh_fem pre_mf_mult(mesh);
  std::string FEMM = PARAM.string_value("FEM_MULT", "fem for multipliers");
  pre_mf_mult.set_qdim(N);
  pre_mf_mult.set_finite_element(mesh.convex_index(),
				 getfem::fem_descriptor(FEMM));
  getfem::partial_mesh_fem mf_mult(pre_mf_mult);
	
	
	
	
	
  dal::bit_vector kept_dof_mult
    = select_dofs_from_im(pre_mf_mult, mimbound, N-1);
  mf_mult.adapt(kept_dof_mult, rejected_elt);
	
	
  size_type nb_dof_mult = mf_mult.nb_dof();
  cout << "nb_dof_mult = " << nb_dof_mult << endl;

 /*
	// Alternative : RANGE_BASIS ............. gmm 
	sparse_matrix BBB(pre_mf.nb_dof(), pre_mf.nb_dof());
	getfem::asm_mass_matrix(BBB, mim, pre_mf, pre_mf);
	std::set<size_type> columns;
	gmm::range_basis(BBB,columns,1.e-12);
	dal::bit_vector kept_dof_mult;
	for (std::set<size_type>:: iterator it=columns.begin();it!=columns.end();++it){
		kept_dof_mult.add(*it);	
	};
	*/
	
	
		
	
	
	
	
	
  // Mass matrix on the boundary
  sparse_matrix B2(mf_rhs.nb_dof(), pre_mf.nb_dof());
  getfem::asm_mass_matrix(B2, mimbound, mf_rhs, pre_mf);
	
  sparse_matrix B(nb_dof_mult, nb_dof);
  getfem::asm_mass_matrix(B, mimbound, mf_mult, mf);
	
  int stabilized_dirichlet =
    int(PARAM.int_value("STABILIZED_DIRICHLET", "Stabilized version of "
			"Dirichlet condition or not"));
  scalar_type dir_gamma0(0);
  sparse_matrix MA(nb_dof_mult, nb_dof_mult), KA(nb_dof, nb_dof);
  sparse_matrix BA(nb_dof_mult, nb_dof);
	
	
	
	// Stabilization of the pressure
	sparse_matrix Mup(nb_dof_p, nb_dof);
	sparse_matrix Mpp(nb_dof_p, nb_dof_p);
	sparse_matrix Mpl(nb_dof_p, nb_dof_mult);
	
	
	getfem::asm_mass_matrix(Mpp,mimbound,mf_p,mf_p);
	asm_mass_matrix_pl(Mpl,mimbound,mf_p,mf_mult,ls);
	asm_mass_matrix_up(Mup,mimbound,mf_p,mf,ls);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
  if (stabilized_dirichlet > 0) {

    sparse_row_matrix E1(nb_dof, nb_dof);

    if (stabilized_dirichlet == 2) {
      // Computation of the extrapolation operator
      scalar_type min_ratio =
	PARAM.real_value("MINIMAL_ELT_RATIO",
			 "Threshold ratio for the fully stab Dirichlet");

      cout << "Computation of the extrapolation operator" << endl;
      dal::bit_vector elt_black_list, dof_black_list;
      size_type nbe = mf_P0.nb_dof();
      plain_vector ratios(nbe);
      sparse_matrix MC1(nbe, nbe), MC2(nbe, nbe);
      getfem::asm_mass_matrix(MC1, mim, mf_P0);
      getfem::asm_mass_matrix(MC2, uncutmim, mf_P0);
      for (size_type i = 0; i < nbe; ++i) {
	size_type cv = mf_P0.first_convex_of_dof(i);
	ratios[cv] = gmm::abs(MC1(i,i)) / gmm::abs(MC2(i,i));
	if (ratios[cv] > 0 && ratios[cv] < min_ratio) elt_black_list.add(cv);
      }
      
	
      sparse_matrix EO(nb_dof, nb_dof);
      sparse_row_matrix T1(nb_dof, nb_dof), EX(nb_dof, nb_dof);
      asm_mass_matrix(EO, uncutmim, mf);

      for (size_type i = 0; i < nb_dof; ++i) {
	bool found = false;
	getfem::mesh::ind_cv_ct ct = mf.convex_to_dof(i);
	getfem::mesh::ind_cv_ct::const_iterator it;
	for (it = ct.begin(); it != ct.end(); ++it)
	  if (!elt_black_list.is_in(*it)) found = true;
	if (found)
	  { gmm::clear(gmm::mat_col(EO, i)); EO(i,i) = scalar_type(1); }
	else
	  dof_black_list.add(i);
      }

      bgeot::mesh_structure::ind_set is;
      base_matrix Mloc;
      for (dal::bv_visitor i(elt_black_list); !i.finished(); ++i) {
	mesh.neighbours_of_convex(i, is);
	size_type cv2 = size_type(-1);
	scalar_type ratio = scalar_type(0);
	for (size_type j = 0; j < is.size(); ++j) {
	  scalar_type r = ratios[is[j]];
	  if (r > ratio) { ratio = r; cv2 = is[j]; }
	}
	GMM_ASSERT1(cv2 != size_type(-1), "internal error");
	compute_mass_matrix_extra_element(Mloc, uncutmim, mf, i, cv2);
	for (size_type ii = 0; ii < gmm::mat_nrows(Mloc); ++ii) 
	  for (size_type jj = 0; jj < gmm::mat_ncols(Mloc); ++jj)
	    EX(mf.ind_dof_of_element(i)[ii], mf.ind_dof_of_element(cv2)[jj])
	      += Mloc(ii, jj);
      }

      gmm::copy(gmm::identity_matrix(), E1);
      gmm::copy(gmm::identity_matrix(), T1);
      for (dal::bv_visitor i(dof_black_list); !i.finished(); ++i)
	gmm::copy(gmm::mat_row(EX, i), gmm::mat_row(T1, i));

      plain_vector BE(nb_dof), BS(nb_dof);
      for (dal::bv_visitor i(dof_black_list); !i.finished(); ++i) {
	BE[i] = scalar_type(1);
	// TODO: store LU decomp.
	double rcond; 
	gmm::SuperLU_solve(EO, BS, BE, rcond);
	gmm::mult(gmm::transposed(T1), BS, gmm::mat_row(E1, i));
	BE[i] = scalar_type(0);
      }
      gmm::clean(E1, 1e-13);

//       gmm::clean(EO, 1e-13); cout << "E0 = " << gmm::transposed(EO) << endl; getchar();
//       cout << "T1 = " << T1 << endl; getchar();
      

      cout << "Extrapolation operator computed" << endl;

//       sparse_row_matrix A1(nb_dof, nb_dof);
//       gmm::mult(E1, gmm::transposed(EO), A1);
//       gmm::clean(A1, 1e-13);
//       cout << "A1 = " << A1 << endl;

    }

    dir_gamma0 = PARAM.real_value("DIRICHLET_GAMMA0",
				  "Stabilization parameter for "
				  "Dirichlet condition");
    getfem::asm_mass_matrix(MA, mimbound, mf_mult);
    asm_stabilization_mixed_term(BA, mimbound, mf, mf_mult, ls);
    asm_stabilization_symm_term(KA, mimbound, mf, ls);
    gmm::scale(MA, -dir_gamma0 * h);
    gmm::scale(BA, -dir_gamma0 * h);
    gmm::scale(KA, -dir_gamma0 * h);

	gmm::scale(Mup,  dir_gamma0 * h);
	gmm::scale(Mpp, -dir_gamma0 * h);  
	gmm::scale(Mpl, dir_gamma0 * h);  
	  
	  
    if (stabilized_dirichlet == 2) {
      sparse_matrix A1(nb_dof_mult, nb_dof);
      gmm::copy(BA, A1);
      gmm::mult(gmm::transposed(E1), gmm::transposed(A1), gmm::transposed(BA));
      sparse_matrix A2(nb_dof, nb_dof);
      gmm::mult(gmm::transposed(E1), KA, A2);
      gmm::mult(gmm::transposed(E1), gmm::transposed(A2), gmm::transposed(KA));
    }
  
     gmm::add(BA, B);

  }

  // Tests
  test_mim(mim, mf_rhs, false);
  test_mim(mimbound, mf_rhs, true);

  /* 
   * Choose boundaries
   */
  getfem::mesh_region r;
  getfem::outer_faces_of_mesh(mesh,r);

  base_node BBmin, BBmax;
  mesh.bounding_box(BBmin, BBmax);
  cout << "mesh bounding box: " << BBmin << " ... " << BBmax << endl;

  for (getfem::mr_visitor i(r); !i.finished(); ++i) {
//  base_node G = gmm::mean_value(mesh.points_of_face_of_convex(i.cv(),i.f()));

/*    if ((gmm::abs(G[0]-BBmin[0]) < 1e-7)||((gmm::abs(G[0]-BBmax[0]) < 1e-7))||((gmm::abs(G[1]-BBmin[1]) < 1e-7 && gmm::abs(G[0]-BBmin[0]) > h && gmm::abs(G[0]-BBmax[0]) > h))||((gmm::abs(G[1]-BBmax[1]) < 1e-7 && gmm::abs(G[0]-BBmax[0]) > h && gmm::abs(G[0]-BBmin[0]) > h))) {
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());
 compteur = compteur +1;
cout << compteur << " " << G[0] << " " << G[1] << endl;
    }
    else {
      mesh.region(NEUMANN_BOUNDARY_NUM).add(i.cv(),i.f());
    }
*/
     mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(),i.f());
  }

	
  // Model system

  getfem::model MS;


 // Linearized elasticity brick
 // MS.add_initialized_fixed_size_data("lambada", plain_vector(1, 0.0));    // lambada vu comme un scalaire egal à 0
  MS.add_initialized_fixed_size_data("nu", plain_vector(1, nu));
  MS.add_fem_variable("u", mf);
  getfem::add_generic_elliptic_brick(MS, mim, "u", "nu");
	
	
//	getfem::add_isotropic_linearized_elasticity_brick(MS, mim, "u", "lambada", "nu");	

  // Linearized incompressibility condition brick
  MS.add_fem_variable("p", mf_p);  // Adding the pressure as variable
  add_linear_incompressibility(MS, mim, "u", "p");	
	
 /*
  // Neumann condition on the box boundary
  plain_vector G(nb_dof_rhs);
  getfem::interpolation_function(mf_rhs, G, g_exact);

  MS.add_initialized_fem_data("NeumannData", mf_rhs, G);
  getfem::add_source_term_brick(MS, mim, "u", "NeumannData", NEUMANN_BOUNDARY_NUM);
*/
	
/*  // Neumann condition on the level-set lsup
  MS.add_multiplier("LAMBDAN", mf_mult, "u");
  plain_vector G(nb_dof_rhs);
  getfem::interpolation_function(mf_rhs, G, g_exact);
  plain_vector G2(nb_dof);
  gmm::mult(gmm::transposed(B2), G, G2);
  MS.add_initialized_fem_data("Neumanndata", mf, G2);
  //MS.add_initialized_fixed_size_data("B2data", gmm::transposed(B2));
  //getfem::add_constraint_with_multipliers(MS, "u", "LAMBDAN", gmm::transposed(B2), G2);
  getfem::add_source_term_brick(MS, mim, "u", "Neumanndata");
*/
	
	
	/*plain_vector G(nb_dof_rhs);
	getfem::interpolation_function(mf_rhs, G, g_exact_ls);
	plain_vector G2(mf_mult.nb_dof());
	getfem::interpolation(mf_rhs, mf_mult,G,G2);

	
	MS.add_initialized_fem_data("Neumanndata", mf_mult, G2);
	add_source_term_brick(MS,mimboundup,"u","Neumanndata", NEUMANN_BOUNDARY_NUM);
	*/

  // Dirichlet on the box boundary
 
	
//	cout << B << endl;
	

  // Dirichlet condition on the box boundary
  MS.add_initialized_fem_data("Dirichletdata", mf, plain_vector(nb_dof));
  getfem::add_Dirichlet_condition_with_multipliers(MS, mim, "u", mf, DIRICHLET_BOUNDARY_NUM, "Dirichletdata");


  // Dirichlet condition on the level-set lsdown
  MS.add_fem_variable("LAMBDA", mf_mult);
  plain_vector H(nb_dof_rhs);
  getfem::interpolation_function(mf_rhs, H, u_exact);
 	
	plain_vector H3(nb_dof_mult);
	getfem :: asm_source_term(H3,mimbound,mf_mult,mf_rhs,H);
	
	cout << mf_mult.nb_dof() << endl;
  getfem::add_constraint_with_multipliers(MS, "u", "LAMBDA", B, H3);

//	getfem::add_constraint_with_penalization(MS, "u", 1.e9, B, H3);
	
  if (stabilized_dirichlet > 0){
       getfem::add_explicit_matrix(MS, "LAMBDA", "LAMBDA", MA);
       getfem::add_explicit_matrix(MS, "u", "u", KA);
	  
	  getfem::add_explicit_matrix(MS, "u", "p", gmm::transposed(Mup));
	  getfem::add_explicit_matrix(MS, "p", "p", Mpp);
	  getfem::add_explicit_matrix(MS, "p", "LAMBDA", Mpl);
	  getfem::add_explicit_matrix(MS, "LAMBDA", "p", gmm::transposed(Mpl));	  
  }


  // Volumic source term
  plain_vector F(nb_dof_rhs); 
  getfem::interpolation_function(mf_rhs, F, rhs);
  MS.add_initialized_fem_data("VolumicData", mf_rhs, F);
  getfem::add_source_term_brick(MS, mim, "u", "VolumicData");
//cout << F2 << endl;  

	

 /* MS.assembly(getfem::model::BUILD_ALL);
	MS.interval_of_variables("u");
	MS.real_tangent_matrix();
	*/
	
	
	
	
  // Solving the problem
  // cout << "Total number of unknown: " << final_brick->nb_dof() << endl;
  gmm::iteration iter(1e-9, 1, 40000);
  getfem::standard_solve(MS, iter);
  plain_vector U(nb_dof);
  gmm::copy(MS.real_variable("u"), U);
	
//  plain_vector LAMBDA(nb_dof_rhs);
//  gmm::copy(MS.real_variable("LAMBDA"), LAMBDA);

  // interpolation of the solution on mf_rhs
  plain_vector Uint(nb_dof_rhs), Vint(nb_dof_rhs), Eint(nb_dof_rhs), Vex(nb_dof_rhs), Uex(nb_dof);
  
  getfem::interpolation_function(mf_rhs, Vex, u_exact);
  getfem::interpolation(mf_rhs, mf, Vex, Uex);

  getfem::interpolation_function(mf_rhs, Vint, u_exact);
  getfem::interpolation(mf, mf_rhs, U, Uint);
  for (size_type i = 0; i < nb_dof_rhs; ++i) {
       Eint[i] = gmm::abs(Uint[i] - Vint[i]);
  }


  // computation of error on u.
  double errmax = 0.0, exactmax = 0.0;
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    if (ls_value(mf_rhs.point_of_dof(i)) >= 0.0) {
      errmax = std::max(errmax, Eint[i]);
      exactmax = std::max(exactmax, Vint[i]);
    }
    else Eint[i] = 0.0;
  cout << "Linfty error: (absolute) " << errmax << \
						"(relative) " <<100.0 * errmax / exactmax << "%" << endl;

  
  cout << "L2 error: " << getfem::asm_L2_dist(mim, mf, U, mf_rhs, Vint) <<
	"(relative) " << 100.0
    * getfem::asm_L2_dist(mim, mf, U, mf_rhs, Vint)
    / getfem::asm_L2_norm(mim, mf_rhs, Vint) << "%" << endl;
	
	
  cout << "H1 error: " << getfem::asm_H1_dist(mim, mf, U, mf_rhs, Vint) <<
	"(relative) " << 100.0
    * getfem::asm_H1_dist(mim, mf, U, mf_rhs, Vint)
    / getfem::asm_H1_norm(mim, mf_rhs, Vint) << "%" << endl;

/*
  // computation of error on multipliers.
cout << "OKKKK" << endl;

  gmm::resize(BA, nb_dof_mult, nb_dof_rhs); gmm::clear(BA);
  gmm::resize(KA, nb_dof_rhs, nb_dof_rhs);  gmm::clear(KA);
  gmm::resize(B, nb_dof_mult, nb_dof_mult); gmm::clear(B);
  
cout << "OKIDOC" << endl;
  
  asm_stabilization_mixed_term(BA, mimbounddown, mf_rhs, mf_mult, lsdown);

cout << "OKIDOC" << endl;

  getfem::asm_mass_matrix(B, mimbounddown, mf_mult, mf_mult);
  asm_stabilization_symm_term(KA, mimbounddown, mf_rhs, lsdown);

  scalar_type err_l2_mult =
    ( gmm::vect_sp(B, LAMBDA, LAMBDA) + gmm::vect_sp(KA, Vint, Vint)
      + 2 * gmm::vect_sp(BA, Vint, LAMBDA) ) / gmm::vect_sp(KA, Vint, Vint);


    
  cout << "L2 error on multipliers: "
       << gmm::sgn(err_l2_mult) * gmm::sqrt(gmm::abs(err_l2_mult)) * 100.0
       << "%" << endl;
  // cout << "L2^2 max on multipliers: " << gmm::vect_sp(KA, Vint, Vint);
  // cout << "  LAMBDA^2: " << gmm::vect_sp(B, LAMBDA, LAMBDA);
  // cout << "  Double produit: " <<  2*gmm::vect_sp(BA, Vint, LAMBDA)<<endl;



*/

  // exporting solution in vtk format.
  {
    getfem::vtk_export exp("xfem_dirichlet.vtk", (2==1));
    exp.exporting(mf); 
    exp.write_point_data(mf, U, "solution");
    //cout << "export done, you can view the data file with (for example)\n"
    //  "mayavi -d xfem_dirichlet.vtk -f WarpScalar -m BandedSurfaceMap "
    //  "-m Outline\n";
  }
  // exporting exact solution in vtk format.
  {
    getfem::vtk_export exp("xfem_dirichlet_exact.vtk", (2==1));
    exp.exporting(mf); 
    exp.write_point_data(mf, Uex, "solution");

    //cout << "export done, you can view the data file with (for example)\n"
    //  "mayavi -d xfem_dirichlet.vtk -f WarpScalar -m BandedSurfaceMap "
    //  "-m Outline\n";
  }
  // exporting error in vtk format.
  {
    getfem::vtk_export exp("xfem_dirichlet_error.vtk", (2==1));
    exp.exporting(mf_rhs); 
    exp.write_point_data(mf_rhs, Eint, "error");
    //cout << "export done, you can view the data file with (for example)\n"
    //  "mayavi -d xfem_dirichlet_error.vtk -f WarpScalar -m BandedSurfaceMap "
    //  "-m Outline\n";
  }
/*  // exporting multipliers in vtk format.
  {
    getfem::vtk_export exp("xfem_dirichlet_mult.vtk", (2==1));
    exp.exporting(mf_mult); 
    exp.write_point_data(mf_mult, LAMBDA, "multipliers");
    //cout << "export done, you can view the data file with (for example)\n"
    //  "mayavi -d xfem_dirichlet_mult.vtk -f WarpScalar -m BandedSurfaceMap "
    //  "-m Outline\n";
  }
*/

  lsmf.write_to_file("xfem_dirichlet_ls.mf", true);
  gmm::vecsave("xfem_dirichlet_ls.U", ls.values());


/*
  unsigned nrefine = mf.linked_mesh().convex_index().card() < 200 ? 32 : 4;
  if (0) {
    //cout << "saving the slice solution for matlab\n";
    getfem::stored_mesh_slice sl, sl0,sll;
    
    
    getfem::mesh_slicer slicer(mf.linked_mesh());
    getfem::slicer_build_stored_mesh_slice sbuild(sl);
    getfem::mesh_slice_cv_dof_data<plain_vector> mfU(mf,U);
    getfem::slicer_isovalues iso(mfU, 0.0, 0);
    getfem::slicer_build_stored_mesh_slice sbuild0(sl0);
    
    slicer.push_back_action(sbuild);  // full slice in sl
    slicer.push_back_action(iso);     // extract isosurface 0
    slicer.push_back_action(sbuild0); // store it into sl0
    slicer.exec(nrefine, mf.convex_index());
    
    getfem::mesh_slicer slicer2(mf.linked_mesh());
    getfem::mesh_slice_cv_dof_data<plain_vector> 
      mfL(ls.get_mesh_fem(), ls.values());
    getfem::slicer_isovalues iso2(mfL, 0.0, 0);
    getfem::slicer_build_stored_mesh_slice sbuildl(sll);
    slicer2.push_back_action(iso2);     // extract isosurface 0
    slicer2.push_back_action(sbuildl); // store it into sl0
    slicer2.exec(nrefine, mf.convex_index());
    
    sl.write_to_file("xfem_dirichlet.sl", true);
    sl0.write_to_file("xfem_dirichlet.sl0", true);
    sll.write_to_file("xfem_dirichlet.sll", true);
    plain_vector UU(sl.nb_points()), LL(sll.nb_points()); 
    sl.interpolate(mf, U, UU);
    gmm::vecsave("xfem_dirichlet.slU", UU);
    // gmm::scale(LAMBDA, 0.005);
    sll.interpolate(mf_mult, LAMBDA, LL);
    gmm::vecsave("xfem_dirichlet.slL", LL);
  }
*/
  return 0; 
}
