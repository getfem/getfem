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
 * Solve the Signorini problem with static Coulomb friction.
 * 
 * Research program.
 */


#include "getfem/getfem_assembling.h"
#include "getfem/getfem_export.h"
#include "getfem/getfem_import.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_Coulomb_friction.h"
#include "getfem/getfem_derivatives.h"
#include "gmm/gmm.h"
#include <fstream>

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;


/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::short_type;
using bgeot::dim_type;
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. 
 * default types of getfem_model_solvers.h
 */
typedef gmm::wsvector<scalar_type> sparse_vector;
typedef gmm::row_matrix<sparse_vector> sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

typedef getfem::model_state<sparse_matrix, sparse_matrix,
			    plain_vector> MODEL_STATE;

/*
 * structure for the friction problem
 */
struct friction_problem {

  enum { DIRICHLET_BOUNDARY, CONTACT_BOUNDARY, NEUMANN_BOUNDARY };
  getfem::mesh mesh;         /* the mesh */
  getfem::mesh_im  mim;      /* integration methods.                         */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the friction solution     */
  getfem::mesh_fem mf_l;     /* mesh_fem for the multipliers.                */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_vm;    /* mesh_fem used for the VonMises stress        */

  getfem::mesh mesh_refu;    /* mesh used for a reference solution.          */
  getfem::mesh_im  mim_refu; /* integration method for error computation.    */
  getfem::mesh_fem mf_refu;  /* mesh_fem used for a reference solution.      */

  scalar_type lambda, mu;    /* Lamé coefficients.                           */
  scalar_type rho, PG;       /* density, and gravity                         */
  scalar_type friction_coef; /* friction coefficient.                        */
  scalar_type gamma;       /* Stabilization parameter Barbosa-Hughes method. */

  scalar_type residual;      /* max residual for the iterative solvers       */
  
  size_type N, noisy, method, population, Dirichlet, Neumann;
  size_type contact_condition;
  std::vector<size_type> blocked_dofs;
  sparse_matrix Brigid;

  scalar_type r;
  scalar_type overlap, subdomsize;
  scalar_type Dirichlet_ratio, Neumann_intensity;
  bool dxexport, mlexport, compare_with_ref;

  plain_vector UREF;

  std::string datafilename, reffilename;
  bgeot::md_param PARAM;

  void solve(void);
  void init(void);
  friction_problem(void)
    : mim(mesh), mf_u(mesh), mf_l(mesh), mf_rhs(mesh),
      mf_vm(mesh), mim_refu(mesh_refu), mf_refu(mesh_refu) {}
};


namespace getfem {

  class position_vector : public getfem::nonlinear_elem_term {
    unsigned N;
    bgeot::base_vector coeff;
    bgeot::multi_index sizes_;
  public:
    position_vector(unsigned NN) : N(NN)
    { sizes_.resize(1); sizes_[0] = short_type(N); }
    const bgeot::multi_index &sizes(size_type) const {  return sizes_; }
    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t)
    { for (size_type i = 0; i < N; ++i) t[i] = ctx.xreal()[i]; }
  };
  
  // Compute the matrix for the momentums
  template<class MAT>
  void asm_momentum_matrix
  (const MAT &B_, const mesh_im &mim, const mesh_fem &mf,
   const mesh_region &rg = mesh_region::all_convexes()) {
    // "non linéaires".
    MAT &B = const_cast<MAT &>(B_);

    position_vector nterm(mf.linked_mesh().dim());
    unsigned line = 0;

    std::vector<scalar_type> V(mf.nb_dof());

    for (unsigned k = 1; k <= mf.get_qdim(); ++k) {
      std::stringstream msg;
      msg << "V(#1)+=comp(vBase(#1))(:," << k << ");";
      generic_assembly assem(msg.str());
      assem.push_mi(mim);
      assem.push_mf(mf);
      assem.push_vec(V);
      assem.assembly(rg);
      gmm::copy(gmm::row_vector(V),
		gmm::sub_matrix(B, gmm::sub_interval(line++, 1),
				gmm::sub_interval(0, gmm::mat_ncols(B))));
    }
    
    for (unsigned k1 = mf.get_qdim(); k1 >= 1; --k1)
      for (unsigned k2 = mf.get_qdim(); k2 >= k1+1; --k2) {
	std::stringstream msg;
	msg << "t=comp(vBase(#1).NonLin(#1));"
	    << "V(#1) += t(:," << k1 << "," << k2 << ")"
	    << "- t(:," << k2 << "," << k1 << ");";
	generic_assembly assem(msg.str());
	assem.push_mi(mim);
	assem.push_mf(mf);
	assem.push_vec(V);
	assem.push_nonlinear_term(&nterm);
	assem.assembly(rg);
	gmm::copy(gmm::row_vector(V),
		  gmm::sub_matrix(B, gmm::sub_interval(line++, 1),
				  gmm::sub_interval(0, gmm::mat_ncols(B))));
      } 
    
  }








}


/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void friction_problem::init(void) {
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string FEM_TYPE_L  = PARAM.string_value("FEM_TYPE_L",
					       "FEM name for multipliers");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "FEM_TYPE_L="  << FEM_TYPE_L << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  std::string meshname
    (PARAM.string_value("MESHNAME", "Nom du fichier de maillage"));

  /* First step : build the mesh */
  bgeot::pgeometric_trans pgt = bgeot::geometric_trans_descriptor(MESH_TYPE);
  N = pgt->dim();
  if (meshname.compare(0,5, "splx:")==0) {
    std::vector<size_type> nsubdiv(N);
    std::fill(nsubdiv.begin(),nsubdiv.end(),
	      PARAM.int_value("NX", "Number of space steps "));
    getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			      PARAM.int_value("MESH_NOISED") != 0);
  
    bgeot::base_matrix M(N,N);
    for (size_type i=0; i < N; ++i) {
      static const char *t[] = {"LX","LY","LZ"};
      M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
    }
    mesh.transformation(M);
  }
  else getfem::import_mesh(meshname, mesh);
  mesh.optimize_structure();

  double h = 0;
  for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i)
    h += mesh.convex_radius_estimate(i);
  h /= double(mesh.convex_index().card());
  cout << "mean h = " << 2*h << endl;
  
  for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i)
    if (mesh.convex_quality_estimate(i) < 0.1)
      GMM_WARNING1("Bad quality of convex " << i << " : "
		   << mesh.convex_quality_estimate(i));

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL");
  if (residual == 0.) residual = 1e-10;

  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  rho = PARAM.real_value("RHO", "Density");
  PG = PARAM.real_value("PG", "Gravity constant");
  friction_coef = PARAM.real_value("FRICTION_COEF", "Friction coefficient");
  Dirichlet = PARAM.int_value("DIRICHLET","Dirichlet condition or not");
  Dirichlet_ratio = PARAM.real_value("DIRICHLET_RATIO",
				     "parameter for Dirichlet condition");
  cout << "Dirichlet_ratio = " << Dirichlet_ratio << endl;
  Neumann = PARAM.int_value("NEUMANN","Neumann condition or not");
  Neumann_intensity = PARAM.real_value("NEUMANN_INTENSITY",
				       "Intensity of the applied force");
  dxexport = (PARAM.int_value("DX_EXPORT", "Exporting on OpenDX format") != 0);
  mlexport = (PARAM.int_value("ML_EXPORT", "Exporting on Matlab format") != 0);

  r = PARAM.real_value("R", "augmentation parameter");
  method = PARAM.int_value("METHOD", "solve method");
  subdomsize = PARAM.real_value("SUBDOMSIZE");
  overlap = PARAM.real_value("OVERLAP");
  contact_condition = PARAM.int_value("CONTACT_CONDITION",
				      "type of contact condition");
  if (contact_condition == 2)
    gamma = PARAM.real_value("GAMMA", "Stabilization parameter");
  noisy = PARAM.int_value("NOISY", "verbosity of iterative methods");
  mf_u.set_qdim(dim_type(N));
  mf_l.set_qdim(dim_type(N));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  getfem::pfem pf_l = getfem::fem_descriptor(FEM_TYPE_L);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_l.set_finite_element(mesh.convex_index(), pf_l);
  mf_vm.set_classical_discontinuous_finite_element(1);
  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u);
  }
  else
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));

  /* set boundary conditions */
  base_node Pmin(N), Pmax(N);
  mesh.bounding_box(Pmin, Pmax);
  base_node center = (Pmin + Pmax )/2.0;
  std::cout << "Reperage des bord de contact et Dirichlet\n";  
  for (dal::bv_visitor cv(mesh.convex_index()); !cv.finished(); ++cv) {
    size_type nf = mesh.structure_of_convex(cv)->nb_faces();
    for (short_type f = 0; f < nf; f++) {
      if (!mesh.is_convex_having_neighbour(cv, f)) {
	base_small_vector un = mesh.normal_of_face_of_convex(cv, f);
	un /= gmm::vect_norm2(un);	
	base_node pt = mesh.points_of_face_of_convex(cv,f)[0];
	if (un[N-1] < -0.000001 && (N != 3 || (gmm::vect_dist2(pt, center)
			   > .99*sqrt(25. + 15*15) && pt[N-1] < 20.1)))
	  mesh.region(CONTACT_BOUNDARY).add(cv, f); 
	if ((un[N-1] > 0.1 && (Dirichlet == 1))
	    || (un[0] < -0.5 && (Dirichlet == 2)))
	  mesh.region(DIRICHLET_BOUNDARY).add(cv, f);
	if (un[N-1] > 0.1 && (Neumann == 1))
	  mesh.region(NEUMANN_BOUNDARY).add(cv, f);
      }
    }
  }
  if (Dirichlet == 3) { // kill rigid motions.
    
    if (N == 3) {
      
      gmm::row_matrix<gmm::wsvector<scalar_type> >
	BBR(N+(N*(N-1))/2, mf_u.nb_dof()), BR(N-1+(N*(N-1))/2, mf_u.nb_dof());
      getfem::asm_momentum_matrix(BBR, mim, mf_u);
      gmm::sub_interval SUBN(0, N-1), SUBU(0, mf_u.nb_dof());
      gmm::copy(gmm::sub_matrix(BBR, SUBN, SUBU),
		gmm::sub_matrix(BR, SUBN, SUBU));
      gmm::copy(gmm::sub_matrix(BBR, gmm::sub_interval(N, (N*(N-1))/2), SUBU),
		gmm::sub_matrix(BR, gmm::sub_interval(N-1, (N*(N-1))/2), SUBU));
      gmm::resize(Brigid, N-1+(N*(N-1))/2, mf_u.nb_dof());
      gmm::copy(BR, Brigid);
    }
    else {
      
      base_node top = center; top[N-1] = Pmax[N-1];
      GMM_ASSERT1(!mf_u.is_reduced(), "To be adapted");
      for (size_type k = 0; k < N-1; ++k) {
	scalar_type dist_to_center = gmm::vect_dist2(Pmin, center);
	scalar_type dist_to_top = gmm::vect_dist2(Pmin, top);
	size_type indmid = size_type(-1), indtop = size_type(-1);
	for (size_type i = 0; i < mf_u.nb_dof(); i+=N) {
	  
	  if (gmm::abs(mf_u.point_of_basic_dof(i)[k]) < 1e-7) {
	    if (gmm::vect_dist2(mf_u.point_of_basic_dof(i), center)<dist_to_center) {
	      indmid = i;
	      dist_to_center = gmm::vect_dist2(mf_u.point_of_basic_dof(i), center); 
	    }
	    if (gmm::vect_dist2(mf_u.point_of_basic_dof(i), top) < dist_to_top) {
	      indtop = i;
	      dist_to_top = gmm::vect_dist2(mf_u.point_of_basic_dof(i), top); 
	    }
	  }
	}
	GMM_ASSERT1(indmid != size_type(-1) && indtop != size_type(-1),
		    "No dof found to kill rigid motions : " << indmid
		    << " : " << indtop);
	cout << "Blocking dof " << indmid+k << " : "
	     << mf_u.point_of_basic_dof(indmid+k) << endl;
	cout << "Blocking dof " << indtop+k << " : "
	     << mf_u.point_of_basic_dof(indtop+k) << endl;
	blocked_dofs.push_back(indmid+k); 
	blocked_dofs.push_back(indtop+k);
      }
    }
  }

  // loading an eventual reference solution

  reffilename = PARAM.string_value("REFSOLFILENAME");
  compare_with_ref = (reffilename.size() > 0);
  if (compare_with_ref) {
    cout << "Comparison with " << reffilename << endl;
    mesh_refu.read_from_file(reffilename + ".meshfem");
    mf_refu.read_from_file(reffilename + ".meshfem");
    mim_refu.set_integration_method(mesh_refu.convex_index(), ppi);
    gmm::resize(UREF, mf_refu.nb_dof());
    gmm::vecload(reffilename + ".U", UREF);
  }
}

template <typename VEC1, typename VEC2>
void calcul_von_mises(const getfem::mesh_fem &mf_u, const VEC1 &U,
		      const getfem::mesh_fem &mf_vm, VEC2 &VM,
		      scalar_type mu=1) {
  assert(mf_vm.get_qdim() == 1); 
  unsigned N = mf_u.linked_mesh().dim(); assert(N == mf_u.get_qdim());
  std::vector<scalar_type> DU(mf_vm.nb_dof() * N * N);

  getfem::compute_gradient(mf_u, mf_vm, U, DU);
  
  gmm::resize(VM, mf_vm.nb_dof());
  scalar_type vm_min, vm_max;
  for (size_type i=0; i < mf_vm.nb_dof(); ++i) {
    VM[i] = 0;
    scalar_type sdiag = 0.;
    for (unsigned j=0; j < N; ++j) {
      sdiag += DU[i*N*N + j*N + j];
      for (unsigned k=0; k < N; ++k) {
	scalar_type e = .5*(DU[i*N*N + j*N + k] + DU[i*N*N + k*N + j]);
	VM[i] += e*e;	
      }      
    }
    VM[i] -= 1./N * sdiag * sdiag;
    vm_min = (i == 0 ? VM[0] : std::min(vm_min, VM[i]));
    vm_max = (i == 0 ? VM[0] : std::max(vm_max, VM[i]));
    assert(VM[i] > -1e-6);
  }
  cout << "Von Mises : min=" << 4*mu*mu*vm_min << ", max="
       << 4*mu*mu*vm_max << "\n";
  gmm::scale(VM, 4*mu*mu);
}

/**************************************************************************/
/*  Function determining the situation of a solution.                     */
/**************************************************************************/

template <class Matrix, class Vector1, class Vector2>
void situation_of(const Vector1 &U, const Vector2 &gap,
		  const Matrix &BN, const Matrix &BT,
		  std::vector<int> &situation) {
  size_type nbc = gmm::mat_nrows(BN), N1 = gmm::mat_nrows(BT) / nbc;
  plain_vector RLN(nbc), RLT(nbc*N1);
  gmm::mult(BN, U, gmm::scaled(gap, -1.0), RLN);
  
  gmm::mult(BT, U, RLT);
  for (size_type j = 0; j < nbc; ++j) {
    if (RLN[j] < -1E-10) situation[j] = 0;
    else if (gmm::vect_norm2
	     (gmm::sub_vector(RLT, gmm::sub_interval(j*N1, N1))) < 1E-10)
      situation[j] = 1;
    else if (N1 > 1 || RLT[j] < 0) situation[j] = 2;
    else situation[j] = 3;
  }
}

/**************************************************************************/
/*  Structure for the Newton Additive Schwarz.                            */
/**************************************************************************/

void build_vB(std::vector<sparse_matrix> &vB,
	      getfem::mesh_fem &mf_u, dal::bit_vector dofgammac,
	      dal::bit_vector dofgammad, double subdomsize, double overlap) {
  //  size_type hsize = F.size(), nb_dof = mf_u.nb_dof();
  // size_type nbc = coulomb.nb_contact_nodes();
  // size_type nbc_coarse = 0;
  // if (usecoarse) nbc_coarse = coulomb_coarse.nb_contact_nodes();
  
  size_type N = mf_u.linked_mesh().dim();

  GMM_ASSERT1(!mf_u.is_reduced(), "To be adapted");

  std::vector<base_node> pts;
  size_type i;
  for (i = 0; i < mf_u.nb_dof(); ++i)
    pts.push_back(mf_u.point_of_basic_dof(i));

  //   std::vector<size_type> links(nbc * N);
  //   std::vector<size_type> links_coarse(nbc_coarse * N);
    
  for (dal::bv_visitor j(dofgammad); !j.finished(); ++j, ++i)
    pts.push_back(mf_u.point_of_basic_dof(j));

  for (dal::bv_visitor j(dofgammac); !j.finished(); ++j) {
    if ((j % N) == 0) { pts.push_back(mf_u.point_of_basic_dof(j)); ++i; }
  }

  for (dal::bv_visitor j(dofgammac); !j.finished(); ++j) {
    if ((j % N) != 0) { pts.push_back(mf_u.point_of_basic_dof(j)); ++i; }
  }

  //   for (size_type i = 0; i < nbc; ++i) { // Ne marche pas dans tous cas ...
  //     size_type j 
  //       = gmm::vect_const_begin(gmm::mat_row(coulomb.get_BN(), i)).index();
  //     pts[nb_dof + i] = mf_u.point_of_basic_dof(j);
  //     //    links[i] = j;
  //     for (size_type k = 0; k < N-1;  ++k) {
  //       pts[nb_dof + nbc + i*(N-1) + k] = mf_u.point_of_basic_dof(j);
  //       // links[nbc + i*(N-1) + k] = j+k-N+1;
  //     }
  //   }
  
  //   if (usecoarse) {
  //     for (size_type i = 0; i < nbc_coarse; ++i) {
  //       size_type j=gmm::vect_const_begin
  //                  (gmm::mat_row(coulomb_coarse.get_BN(),i)).index();
  //       links_coarse[i] = j;
  //       for (size_type k = 0; k < N-1;  ++k)
  // 	links_coarse[nbc_coarse + i*(N-1) + k] = j+k-N+1;
  //     }
  //   }

  // cout << "links = " << links << endl;
  // cout << "links_coarse = " << links_coarse << endl;
  // cout << "nb_dof = " << nb_dof << endl;
  // cout << "nb_dof_coarse = " << mf_coarse.nb_dof() << endl;
 
  gmm::rudimentary_regular_decomposition(pts, subdomsize, overlap, vB);
  size_type nsd = vB.size();

  cout << "Number of subdomains : " << nsd << endl;
 
  //   if (usecoarse) {
  //     vB.resize(nsd+1);
  //     cout << "interpolation coarse mesh\n";
  //     size_type nb_dof_coarse = mf_coarse.nb_dof();
  //     gmm::resize(vB[nsd], nb_dof, nb_dof_coarse);
  //     getfem::interpolation(mf_coarse, mf_u, vB[nsd], true);
  //     gmm::resize(vB[nsd], hsize, nb_dof_coarse+nbc_coarse*N);
  
  //     gmm::unsorted_sub_index si1(links), si2(links_coarse);
  //     gmm::sub_interval si3(nb_dof, nbc*N), si4(nb_dof_coarse,nbc_coarse*N);
  
  //     sparse_matrix Maux(nbc*N, nbc_coarse*N);
  //     gmm::copy(gmm::sub_matrix(vB[nsd], si1, si2),Maux);
  
  //     // cout << "Maux = " << Maux << endl;
  //     gmm::copy(Maux, gmm::sub_matrix(vB[nsd], si3, si4));
  
  //     // cout << "vB[nsd] = " << vB[nsd] << endl;
  
  //     // gmm::copy(gmm::sub_matrix(vB[nsd], si1, si2),Maux);
  //     // cout << "Maux = " << Maux << endl;
  
  //     ++nsd;
  //   }
  cout << "There is " << nsd << " sub-domains\n";
}

struct Coulomb_NewtonAS_struct
  : public gmm::NewtonAS_struct<sparse_matrix, sparse_matrix> {
  
  std::vector<sparse_matrix> vB;
  getfem::mdbrick_abstract<MODEL_STATE> *coulomb;
  MODEL_STATE MS;

  size_type size(void) { return coulomb->nb_dof(); }
  const std::vector<sparse_matrix> &get_vB() { return vB; }
    
  // very inefficient, to be done again.

  void compute_tangent_matrix(sparse_matrix &M, Vector &x) {
    gmm::copy(x, MS.state());
    coulomb->compute_tangent_matrix(MS);
    gmm::copy(MS.tangent_matrix(), M);
  }
  void compute_sub_tangent_matrix(sparse_matrix &Mloc, Vector &x,
				  size_type i) {
    gmm::copy(x, MS.state());
    coulomb->compute_tangent_matrix(MS);
    
    sparse_matrix aux(gmm::mat_ncols(vB[i]),
		      gmm::mat_ncols(MS.tangent_matrix()));
    gmm::mult(gmm::transposed(vB[i]), MS.tangent_matrix(), aux);
    gmm::mult(aux, vB[i], Mloc);
  }
  void compute_sub_F(Vector &fi, Vector &x, size_type i) {
    gmm::copy(x, MS.state());
    coulomb->compute_residual(MS);
    gmm::mult(gmm::transposed(vB[i]), MS.residual(), fi);
  }
  void compute_F(Vector &f, Vector &x) {
    gmm::copy(x, MS.state());
    coulomb->compute_residual(MS);
    gmm::copy(MS.residual(), f);
  }
  Coulomb_NewtonAS_struct(getfem::mdbrick_abstract<MODEL_STATE> &C) : coulomb(&C),MS(C) {}
  
};


/**************************************************************************/
/*  Stiffness matrix for Barbosa-Hughes method and linear elasticity.     */
/**************************************************************************/
/* attention, il manque le \gamma pour le moment.                         */

namespace getfem {

  template<class MAT, class VECT>
  void asm_stiffmatrix_for_stabilization_on_linear_elasticity
  (const MAT &RM_, const mesh_im &mim, const mesh_fem &mf,
   const mesh_fem &mf_data, const VECT &LAMBDA, const VECT &MU,
   const mesh_region &rg) { // à simplifier, faire des réductions par termes
    // "non linéaires".
    MAT &RM = const_cast<MAT &>(RM_);
    GMM_ASSERT1(mf_data.get_qdim() == 1,
		"invalid data mesh fem (Qdim=1 required)");
    
    GMM_ASSERT1(mf.get_qdim() == mf.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");
    
    generic_assembly
      assem("lambda=data$1(#2); mu=data$2(#2);"
	    "t=comp(Normal().vGrad(#1).Normal().Normal().vGrad(#1)"
	    ".Normal().Base(#2).Base(#2));"
	    "M(#1,#1)+= sym(t(i,:,i,j,j,k,:,k,l,l,p,q).mu(p).mu(q)*4"
	    "+ t(i,:,i,j,j,k,:,l,l,k,p,q).mu(p).lambda(q)*2"
	    "+ t(i,:,j,j,i,k,:,k,l,l,p,q).lambda(p).mu(q)*2"
	    "+ t(i,:,j,j,i,k,:,l,l,k,p,q).lambda(p).lambda(q))");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mf_data);
    assem.push_data(LAMBDA);
    assem.push_data(MU);
    assem.push_mat(RM);
    assem.assembly(rg);
  }
  
  template<class MAT, class VECT>
  void asm_mixed_matrix_for_stabilization_on_linear_elasticity
  (const MAT &RM_, const mesh_im &mim, const mesh_fem &mf,
   const mesh_fem &mf_l, const mesh_fem &mf_data, const VECT &LAMBDA,
   const VECT &MU, const mesh_region &rg) {
    MAT &RM = const_cast<MAT &>(RM_);

    GMM_ASSERT1(mf_data.get_qdim() == 1,
		"invalid data mesh fem (Qdim=1 required)");
    GMM_ASSERT1(mf.get_qdim() == mf.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");
    
    generic_assembly
      assem("lambda=data$1(#3); mu=data$2(#3);"
	    "t=comp(vBase(#2).Normal().Normal().vGrad(#1).Normal().Base(#3));"
	    "M(#2,#1)+= t(:,k,k,i,:,i,j,j,p).mu(p)*2"
	    "+ t(:,k,k,i,:,j,j,i,p).lambda(p)");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mf_l);
    assem.push_mf(mf_data);
    assem.push_data(LAMBDA);
    assem.push_data(MU);
    assem.push_mat(RM);
    assem.assembly(rg);
  }
  
  // A mettre dans getfem_modeling.h (et à changer en additional matrices).

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_additional_matrix : public mdbrick_abstract<MODEL_STATE>  {
  public :
    
    TYPEDEF_MODEL_STATE_TYPES;
    
  protected :
    
    mdbrick_abstract<MODEL_STATE> &sub_problem;
    T_MATRIX M;
    size_type num_fem1, num_fem2;
    
    virtual void proper_update(void) { }
    
  public :
    
    template <typename MAT> void set_matrix(const MAT &MM)
    { gmm::resize(M, gmm::mat_nrows(MM),gmm::mat_ncols(MM)); gmm::copy(MM,M); }
    
    const T_MATRIX &get_matrix() { return M; }
    
    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      if (gmm::mat_nrows(M) > 0) {
	const mesh_fem &mf_u1 = *(this->mesh_fems[num_fem1]);
	size_type i1 = this->mesh_fem_positions[num_fem1], nbd1=mf_u1.nb_dof();
	const mesh_fem &mf_u2 = *(this->mesh_fems[num_fem2]);
	size_type i2 = this->mesh_fem_positions[num_fem2], nbd2=mf_u2.nb_dof();
	gmm::sub_interval SUBI(i0+i1, nbd1);
	gmm::sub_interval SUBJ(i0+i2, nbd2);
	gmm::add(M, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
      }
    }
    
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				     size_type) {
      if (gmm::mat_nrows(M) > 0) {
	const mesh_fem &mf_u1 = *(this->mesh_fems[num_fem1]);
	size_type i1 = this->mesh_fem_positions[num_fem1], nbd1=mf_u1.nb_dof();
	const mesh_fem &mf_u2 = *(this->mesh_fems[num_fem2]);
	size_type i2 = this->mesh_fem_positions[num_fem2], nbd2=mf_u2.nb_dof();
	gmm::sub_interval SUBI(i0+i1, nbd1);
	gmm::sub_interval SUBJ(i0+i2, nbd2);
	gmm::mult_add(M, gmm::sub_vector(MS.state(), SUBJ),
		      gmm::sub_vector(MS.residual(), SUBI));
      }
    }
    
    mdbrick_additional_matrix(mdbrick_abstract<MODEL_STATE> &problem,
			      size_type numfem1 = 0, size_type numfem2 = 0)
      : sub_problem(problem), num_fem1(numfem1), num_fem2(numfem2) {
      this->add_sub_brick(sub_problem);
      this->force_update();
    }
    
  };
  
}


/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

void friction_problem::solve(void) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  N = mesh.dim();
  cout << "Number of dof for u: " << mf_u.nb_dof() << endl;

  // Linearized elasticity brick.
  getfem::mdbrick_isotropic_linearized_elasticity<MODEL_STATE> ELAS(mim, mf_u, lambda,mu);

  getfem::mdbrick_additional_matrix<MODEL_STATE> AUGMENTATION(ELAS);

  // Defining the volumic source term brick.
  plain_vector F(nb_dof_rhs * N);
  plain_vector f(N);
  if (Dirichlet == 2) f[0] = -rho*PG; else f[N-1] = -rho*PG;
  for (size_type i = 0; i < nb_dof_rhs; ++i)
      gmm::copy(f,gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  getfem::mdbrick_source_term<MODEL_STATE> VOL_F(AUGMENTATION, mf_rhs, F);
  
  // Defining the applied force source term brick
  gmm::clear(f);
  if (Neumann == 1) f[N-1] = Neumann_intensity;
  for (size_type i = 0; i < nb_dof_rhs; ++i)
      gmm::copy(f,gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  getfem::mdbrick_source_term<MODEL_STATE> NEUMANN_F(VOL_F, mf_rhs, F, NEUMANN_BOUNDARY);
 
  GMM_ASSERT1(!mf_rhs.is_reduced(), "To be adapted");
  // Dirichlet condition brick.
  gmm::clear(F);
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    F[i*N+N-1] = Dirichlet_ratio * mf_rhs.point_of_basic_dof(i)[N-1];

  getfem::mdbrick_constraint<MODEL_STATE> RIGID_MOTIONS(NEUMANN_F);
  if (Dirichlet == 3) {
    if (blocked_dofs.size() > 0) {
      sparse_matrix CB(blocked_dofs.size(), mf_u.nb_dof());
      for (size_type i = 0; i < blocked_dofs.size(); ++i)
	CB(i, blocked_dofs[i]) = scalar_type(1);
      RIGID_MOTIONS.set_constraints(CB, plain_vector(blocked_dofs.size()));
    }
    else if (gmm::mat_nrows(Brigid) > 0)
      RIGID_MOTIONS.set_constraints(Brigid,
				    plain_vector(gmm::mat_nrows(Brigid)));
  }

  getfem::mdbrick_Dirichlet<MODEL_STATE> DIRICHLET(RIGID_MOTIONS, DIRICHLET_BOUNDARY);
  DIRICHLET.rhs().set(mf_rhs, F);
//   if (method == 2)
//     DIRICHLET.set_constraints_type(getfem::AUGMENTED_CONSTRAINTS);
//   else
//     DIRICHLET.set_constraints_type(getfem::ELIMINATED_CONSTRAINTS);

    
  // contact condition for Lagrange elements
  GMM_ASSERT1(!mf_u.is_reduced(), "To be adapted");
  dal::bit_vector cn, dn;
  size_type nbc;
  if (contact_condition == 0) {
    cn = mf_u.basic_dof_on_region(CONTACT_BOUNDARY); // cout << "cn = " << cn << endl;
    dn = mf_u.basic_dof_on_region(DIRICHLET_BOUNDARY); // cout << "dn = " << dn << endl;
  }
  else {
    cn = mf_l.basic_dof_on_region(CONTACT_BOUNDARY);
    dn = mf_l.basic_dof_on_region(DIRICHLET_BOUNDARY);
  }
  cn.setminus(dn); // not to be done for P0 ...
  nbc = cn.card()/N;

  cout << "Nbc = " << nbc << endl;

  sparse_matrix BN(nbc, mf_u.nb_dof());
  sparse_matrix MMB(nbc, nbc);
  sparse_matrix BT((N-1)*nbc, mf_u.nb_dof());
  plain_vector gap(nbc);
  sparse_matrix AUG_M;

  switch (contact_condition) {
  case 0:
    {
      size_type jj = 0;
      for (dal::bv_visitor i(cn); !i.finished(); ++i)
	if (i % N == 0) {
	  BN(jj, i+N-1) = -1.;
	  gap[jj] = mf_u.point_of_basic_dof(i)[N-1];
	  for (size_type k = 0; k < N-1; ++k) BT((N-1)*jj+k, i+k) = 1.;
	  ++jj;
	}
      sparse_matrix BB(mf_u.nb_dof(), mf_u.nb_dof());
      getfem::asm_mass_matrix(BB, mim, mf_u, mf_u, CONTACT_BOUNDARY);
      std::vector<size_type> ind;
      for (dal::bv_visitor i(cn); !i.finished(); ++i) 
	if ((i%N) == N-1) ind.push_back(i);
      gmm::sub_index SUBI(ind);
      gmm::copy(gmm::sub_matrix(BB, SUBI, SUBI), MMB);
    }
    break;
  case 1: case 2:
    {
      sparse_matrix BB(mf_l.nb_dof(), mf_u.nb_dof());
      std::vector<std::vector<size_type> > ind(N);
      for (dal::bv_visitor i(cn); !i.finished(); ++i) 
	ind[i%N].push_back(i);
      GMM_ASSERT1(ind[0].size() == nbc, "Internal error");
      gmm::sub_index SUBI(ind[N-1]);
      gmm::sub_interval SUBJ(0, mf_u.nb_dof());
      getfem::asm_mass_matrix(BB, mim, mf_l, mf_u, CONTACT_BOUNDARY);
      gmm::copy(gmm::scaled(gmm::sub_matrix(BB, SUBI, SUBJ), -1.0), BN);
      
      for (size_type k = 0; k < N-1; ++k) {
	gmm::sub_index SUBIk(ind[k]);
	gmm::copy(gmm::sub_matrix(BB, SUBIk, SUBJ),
		  gmm::sub_matrix(BT, gmm::sub_slice(k, nbc, N-1), SUBJ));
      }

      plain_vector pregap(mf_u.nb_dof());
      for (size_type i=0; i < mf_u.nb_dof(); ++i)
        if (i % N == N-1) pregap[i] = -mf_u.point_of_basic_dof(i)[N-1];
      gmm::mult(BN, pregap, gap);
      // gmm::clear(gap);

      if (contact_condition == 2) {
	sparse_matrix AUG_MM(mf_l.nb_dof(), mf_l.nb_dof());
	getfem::asm_mass_matrix(AUG_MM, mim, mf_l, mf_l, CONTACT_BOUNDARY);
	gmm::scale(AUG_MM, gamma);
	gmm::resize(AUG_M, nbc, nbc);
	gmm::copy(gmm::sub_matrix(AUG_MM, SUBI, SUBI), AUG_M);

	sparse_matrix AUG_C(mf_u.nb_dof(), mf_u.nb_dof());
	plain_vector LAMBDA(nb_dof_rhs, lambda), MU(nb_dof_rhs, mu);
	getfem::asm_stiffmatrix_for_stabilization_on_linear_elasticity
	  (AUG_C, mim, mf_u, mf_rhs, LAMBDA, MU, CONTACT_BOUNDARY);
	gmm::scale(AUG_C, gamma);
	AUGMENTATION.set_matrix(gmm::scaled(AUG_C, -1.0));
	
	sparse_matrix AUG_D(mf_l.nb_dof(), mf_u.nb_dof());
	getfem::asm_mixed_matrix_for_stabilization_on_linear_elasticity
	  (AUG_D, mim, mf_u, mf_l, mf_rhs, LAMBDA, MU, CONTACT_BOUNDARY);
	gmm::scale(AUG_D, gamma);
	gmm::add(gmm::scaled(gmm::sub_matrix(AUG_D, SUBI, SUBJ), 1.0), BN);
      }

    }
    break;
  default: GMM_ASSERT1(false, "Unknown contact condition option");
  }

  getfem::mdbrick_Coulomb_friction<MODEL_STATE>
    FRICTION(DIRICHLET, BN, gap, friction_coef, BT);
  FRICTION.set_r(r);
  if (contact_condition == 2) FRICTION.set_augmented_matrix(AUG_M);
  
  cout << "Total number of variables: " << FRICTION.nb_dof() << endl;
  MODEL_STATE MS(FRICTION);
 
  switch (method) {
  case 0 :
    {
      MS.adapt_sizes(FRICTION);
      if (Dirichlet == 3) gmm::fill(MS.state(), scalar_type(-1));
      gmm::iteration iter(residual, int(noisy), 40000);
      
      getfem::default_newton_line_search ls;
      getfem::model_problem<MODEL_STATE> mdpb(MS,FRICTION,ls);
      
      int l = gmm::get_warning_level();
      gmm::set_warning_level(1);
      getfem::useful_types<MODEL_STATE>::plsolver_type plinsolv
	= getfem::default_linear_solver(FRICTION);
      if (N == 3)
	plinsolv.reset
	  (new getfem::linear_solver_gmres_preconditioned_ilu<sparse_matrix,
	   plain_vector>());
      getfem::classical_Newton(mdpb, iter, *plinsolv);
      gmm::set_warning_level(l);
    }
    break;
  case 1 :
    {
      GMM_ASSERT1(!mf_u.is_reduced(), "To be adapted");
      gmm::iteration iter(residual, int(noisy));
      iter.set_maxiter(1000000);
      Coulomb_NewtonAS_struct NS(FRICTION);
      build_vB(NS.vB, mf_u, mf_u.basic_dof_on_region(CONTACT_BOUNDARY),
	       mf_u.basic_dof_on_region(DIRICHLET_BOUNDARY), subdomsize, overlap);
      gmm::fill_random(MS.state());
      gmm::Newton_additive_Schwarz(NS, MS.state(), iter,
				   gmm::identity_matrix(),
				   gmm::using_superlu(), gmm::using_gmres());
    }
    break;
  }

  std::vector<scalar_type> U(mf_u.nb_dof());
  gmm::copy(ELAS.get_solution(MS), U);
  
  std::vector<int> situation1(nbc);
  situation_of(gmm::sub_vector(MS.state(),
			       gmm::sub_interval(0, mf_u.nb_dof())),
	       gap, BN, BT, situation1);
  cout << "situation of solution : ";
  for (size_type j = 0; j < nbc; ++j) cout << situation1[j];
  cout << endl;
  cout << "Final residual : " <<  MS.reduced_residual_norm() << endl;
  cout << "Norm of solution : " << gmm::vect_norm2(U) << endl;

  
  plain_vector LN1(nbc), LT1(nbc*(N-1));
  gmm::copy(FRICTION.get_LN(MS), LN1);
  gmm::copy(FRICTION.get_LT(MS), LT1);

  if (contact_condition == 0) { // to be done also for the friction ...
    gmm::iteration itercg(1e-12, 1);
    plain_vector LL(nbc);
    gmm::cg(MMB, LL, LN1, gmm::identity_matrix(), itercg);
    gmm::copy(LL, LN1);
  }

  plain_vector BIGLN(mf_l.nb_dof());
  size_type j = 0;
  for (dal::bv_visitor i(cn); !i.finished(); ++i)
    if (i % N == N-1) BIGLN[i] = LN1[j++];

  cout << "contact stress : " << LN1 << endl;
  cout << "friction stress : " << LT1 << endl;

  if (mlexport) {
    getfem::mesh_fem mf_printed_vm(mesh);
    mf_printed_vm.set_finite_element(mesh.convex_index(),
				     getfem::classical_discontinuous_fem
				     (mesh.trans_of_convex(0), 3));
    mf_u.write_to_file(datafilename + ".meshfem", true);
    mf_l.write_to_file(datafilename + ".meshfem_l", true);
    mf_printed_vm.write_to_file(datafilename + ".meshfem_vm", true);
    gmm::vecsave(datafilename + ".U", U);
    
    plain_vector VM(mf_printed_vm.nb_dof());
    getfem::interpolation_von_mises(mf_u, mf_printed_vm, U, VM);
    gmm::vecsave(datafilename + ".VM", VM);

    getfem::slicer_boundary a0(mf_l.linked_mesh()); //,
    //			       getfem::slicer_none::static_instance(),
    //		       CONTACT_BOUNDARY);
    getfem::stored_mesh_slice sl;
    getfem::slicer_build_stored_mesh_slice a1(sl);
    getfem::mesh_slicer slicer(mf_l.linked_mesh());
    slicer.push_back_action(a0);
    slicer.push_back_action(a1);
    unsigned nrefine = 4;
    slicer.exec(nrefine, mf_l.convex_index());
    sl.write_to_file(datafilename + ".sl", true);
    
    gmm::vecsave(datafilename + ".BIGLN", BIGLN);
    plain_vector LLN(sl.nb_points()*N);
    sl.interpolate(mf_l, BIGLN, LLN);
    gmm::vecsave(datafilename + ".LN", LLN);
  }

  if (compare_with_ref) { // error with reference solution
    plain_vector UR(mf_refu.nb_dof());
    getfem::interpolation(mf_u, mf_refu, U, UR, true);
    gmm::add(gmm::scaled(UREF, -1.0), UR);

    getfem::mesh_fem mf_refl(mesh_refu);
    mf_refl.read_from_file(reffilename + ".meshfem_l");
    plain_vector BIGRLN(mf_refl.nb_dof()), BIGREFLN(mf_refl.nb_dof());
    gmm::vecload(reffilename + ".BIGLN", BIGREFLN);
    getfem::interpolation(mf_l, mf_refl, BIGLN, BIGRLN, true);
    gmm::add(gmm::scaled(BIGREFLN, -1.0), BIGRLN);

    scalar_type l2gc = getfem::asm_L2_norm(mim_refu, mf_refl, BIGRLN,
					   CONTACT_BOUNDARY);
    cout << "L2 error on contact pressure = " << l2gc << endl;
    scalar_type l2 = getfem::asm_L2_norm(mim_refu, mf_refu, UR);
    scalar_type h1 = getfem::asm_H1_norm(mim_refu, mf_refu, UR);
    cout << "L2 error = " << l2 << endl
	 << "H1 error = " << h1 << endl
	 << "Linfty error = " << gmm::vect_norminf(UR) << endl;

    l2gc = getfem::asm_L2_norm(mim_refu, mf_refl, BIGREFLN,
			       CONTACT_BOUNDARY);
    l2 = getfem::asm_L2_norm(mim_refu, mf_refu, UREF);
    h1 = getfem::asm_H1_norm(mim_refu, mf_refu, UREF);
    cout << "L2 contact pressure ref = " << l2gc << endl
	 << "L2 ref = " << l2 << endl
	 << "H1 ref = " << h1 << endl
	 << "Linfty ref = " << gmm::vect_norminf(UREF) << endl;

    gmm::vecsave(reffilename + "_error.U", UR);
  }

  if (dxexport) {
    getfem::dx_export exp(datafilename + ".dx", false);
    getfem::stored_mesh_slice sl;
    if (N <= 2)
      sl.build(mesh, getfem::slicer_none(),4);
    else
      sl.build(mesh, getfem::slicer_boundary(mesh),4);
    exp.exporting(sl,true);
    exp.exporting_mesh_edges();
    exp.write_point_data(mf_u, U);
    exp.serie_add_object("deformationsteps");
    std::vector<scalar_type> VM;
    calcul_von_mises(mf_u, U, mf_vm, VM, mu);
    exp.write_point_data(mf_vm, VM); 
    exp.serie_add_object("vonmisessteps");    
  }
}

  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  // gmm::set_warning_level(1);

  friction_problem p;
  p.PARAM.read_command_line(argc, argv);
  p.init();
  p.solve();
  
  return 0; 
}
