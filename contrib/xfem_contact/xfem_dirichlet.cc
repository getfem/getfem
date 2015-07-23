/*===========================================================================

 Copyright (C) 2002-2015 Yves Renard, Julien Pommier.

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

/*
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
#include "getfem/bgeot_mesh_structure.h"
#include "getfem/getfem_config.h"

#if GETFEM_HAVE_METIS_OLD_API
extern "C" void METIS_PartGraphKway(int *, int *, int *, int *, int *, int *,
                                    int *, int *, int *, int *, int *);
extern "C" void METIS_PartGraphRecursive(int *, int *, int *, int *, int *, int *,
                                         int *, int *, int *, int *, int *);
extern "C" void METIS_mCPartGraphKway(int *, int *, int *, int *, int *, int *, int *,
                                      int *, int *, float *, int *, int *, int *);
extern "C" void METIS_mCPartGraphRecursive(int *, int *, int *, int *, int *, int *, int *,
                                           int *, int *, int *, int *, int *);
#elif GETFEM_HAVE_METIS
# include <metis.h>
#endif

using std::endl; using std::cout;


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

/* 
 * Exact solution 
 */
double Radius;
int u_version;
double u_alpha = 1.5;
double u_B = 20.0;
double u_n = 8.0;
double dtheta = M_PI/36;
double sol5_a = 0.075, sol5_b = 30.0;

double u_exact(const base_node &p) {
  double R = Radius, r=gmm::vect_norm2(p), T=atan2(p[1], p[0])+dtheta;
  switch (u_version) {
  case 0: {
    double sum = std::accumulate(p.begin(), p.end(), double(0));
    
    return 5.0 * sin(sum) * (r*r - Radius*Radius);
  }
  case 1: {
    double A=u_alpha, n=u_n;
    return (R*R - r*r *(1+A*(1.0 + sin(n*T))));
  }
  case 2: {
    double A=u_alpha, n=u_n, B=u_B;
    return (R*R - r*r *(1+A*(1.0 + sin(n*T)))) * cos(B*r);      
  }
  case 3: {
    double A=u_alpha, n=u_n;
    return 5*(R*R*R*R - r*r*r*r*(1+A*(1.0 + sin(n*T))));      
  }
  case 4:{
    return 5.0 * (r*r*r - Radius*Radius*Radius);
  }
  case 5:{
    double a = sol5_a, b = sol5_b;
    T = atan2(p[1], p[0]);
    return 5.0 * r*r*(r - Radius*(1 - a*sin(b*T)));
  }
  case 6:{
    double sum = std::accumulate(p.begin(), p.end(), double(0));
    
    return 5.0 * sin(sum) * (r*r*r - Radius*Radius*Radius);
  }
  case 7:{
    double sum = std::accumulate(p.begin(), p.end(), double(0));
    
    return 5.0 * sin(sum) * (r*r*r*r - Radius*Radius*Radius*Radius);
  }
  case 8:{
    double rho=gmm::sqrt(p[0]*p[0]+ p[1]*p[1]+ p[2]*p[2]); 
    
    return 5.0 * ( Radius*Radius* Radius-rho*rho*rho);
  }
  }
  GMM_ASSERT1(false, "Invalid exact solution");
}

double g_exact(const base_node &p) {
  // value of the normal derivative. Due to the pb idem than the opposite of
  // norm of the gradient.
  double R = Radius, r=gmm::vect_norm2(p);
  switch (u_version) {
  case 0: {
    double sum=std::accumulate(p.begin(), p.end(), double(0)), norm_sqr = r*r;
    if (norm_sqr < 1e-10) norm_sqr = 1e-10;
    return 5.0 * (sum * cos(sum) * (norm_sqr - R*R)
		  + 2.0 * norm_sqr * sin(sum)) / sqrt(norm_sqr);
  }
  case 1: {
    double A=u_alpha, T=atan2(p[1], p[0])+dtheta, n=u_n;
    return - sqrt(r*r*pow(2.0*sin(T)+2.0*sin(T)*A+2.0*sin(T)*A*sin(n*T)+cos(T)*A*cos(n*T)*n,2.0)+r*r*pow(-2.0*cos(T)-2.0*cos(T)*A-2.0*cos(T)*A*sin(n*T)+sin(T)*A*cos(n*T)*n,2.0));
  } 
  case 2: {
    double A=u_alpha, T=atan2(p[1], p[0])+dtheta, n=u_n,B=u_B;
    if (gmm::abs(r) < 1e-10) r = 1e-10;
    return -(4.0*r*cos(B*r)+8.0*r*cos(B*r)*A+8.0*r*cos(B*r)*A*A
	     +2.0*sin(B*r)*B*R*R-2.0*sin(B*r)*B*r*r
	     +r*A*A*pow(cos(n*T),2.0)*n*n*cos(B*r)+8.0*r*cos(B*r)*A*A*sin(n*T)
	     -4.0*sin(B*r)*B*r*r*A*sin(n*T)
	     +2.0*sin(B*r)*B*r*r*A*A*pow(cos(n*T),2.0)
	     -4.0*r*cos(B*r)*A*A*pow(cos(n*T),2.0)+8.0*r*cos(B*r)*A*sin(n*T)
	     -4.0*sin(B*r)*B*r*r*A*A*sin(n*T)-4.0*sin(B*r)*B*r*r*A*A
	     -4.0*sin(B*r)*B*r*r*A+2.0*sin(B*r)*B*R*R*A
	     +2.0*sin(B*r)*B*R*R*A*sin(n*T))
      / sqrt(A*A*pow(cos(n*T),2.0)*n*n +4.0+8.0*A+8.0*A*sin(n*T)
	     +8.0*A*A+8.0*A*A*sin(n*T)-4.0*A*A*pow(cos(n*T),2.0));
  }
  case 3: {
    double A=u_alpha, T=atan2(p[1], p[0])+dtheta, n=u_n;
    return -5*r*r*r*sqrt(16.0 + 32.0*A*A + 32.0*A*sin(n*T) + 32.0*A
			 + 32.0*A*A*sin(n*T) - 16*gmm::sqr(A*cos(n*T))
			 + gmm::sqr(A*cos(n*T)*n));
  }   
  case 4: {
    r=gmm::vect_norm2(p);
    return 15*r*r;
  }
  case 5:{
    double a = sol5_a, b = sol5_b;
    double T = atan2(p[1], p[0]);
    return 5.0*r*sqrt( gmm::sqr(3*r-2*Radius*(1-a*sin(b*T)))
		       + gmm::sqr(R*cos(b*T)*b*a));
  }  
  case 6: {
    double sum=std::accumulate(p.begin(), p.end(), double(0)), T=atan2(p[1], p[0]);
    return 5*cos(sum)*(cos(T)+sin(T))*(r*r*r-Radius*Radius*Radius)+15*r*r*sin(sum);
    //  return 5.0*sqrt( gmm::sqr(5*cos(sum)*(cos(T)+sin(T))*(r*r*r-Radius*Radius*Radius)+15*r*r*sin(sum))
    //	     + gmm::sqr(5*cos(sum)*(cos(T)-sin(T))*(r*r*r-Radius*Radius*Radius))) ;
  }  
  case 7: {
    double sum=std::accumulate(p.begin(), p.end(), double(0)), T=atan2(p[1], p[0]);
    return 5*cos(sum)*(cos(T)+sin(T))*(r*r*r*r-Radius*Radius*Radius*Radius)+20*r*r*r*sin(sum);
  }  
  case 8: {
    double rho=gmm::sqrt(p[0]*p[0]+ p[1]*p[1]+ p[2]*p[2]); 
    return -15*rho*rho; ;
  }
  }
  return 0;
}

double rhs(const base_node &p) {
  double R = Radius, r=gmm::vect_norm2(p);
  switch (u_version) {
  case 0: {
    double sum = std::accumulate(p.begin(), p.end(), double(0));
    double norm_sqr = r*r, N = double(gmm::vect_size(p));
    return 5.0 * (N * sin(sum) * (norm_sqr - R*R-2.0)
		  - 4.0 * sum * cos(sum));
  }
  case 1: {
    double A=u_alpha, T=atan2(p[1], p[0])+dtheta, n=u_n;
    return -(-4.0-4.0*A-4.0*A*sin(n*T)+A*sin(n*T)*n*n);
  }
  case 2: {
    double A=u_alpha, T=atan2(p[1], p[0])+dtheta, n=u_n, B=u_B;
    if (gmm::abs(r) < 1e-10) r = 1e-10;
    return (4.0*r*cos(B*r)*A*sin(n*T)-5.0*sin(B*r)*B*r*r*A+4.0*r*cos(B*r)
	    +4.0*r*cos(B*r)*A+sin(B*r)*B*R*R-5.0*sin(B*r)*B*r*r
	    -5.0*sin(B*r)*B*r*r*A*sin(n*T)-r*r*r*cos(B*r)*B*B
	    -r*A*sin(n*T)*n*n*cos(B*r)+r*cos(B*r)*B*B*R*R
	    -r*r*r*cos(B*r)*B*B*A-r*r*r*cos(B*r)*B*B*A*sin(n*T))/r;
  }
  case 3: {
    double A=u_alpha, T=atan2(p[1], p[0])+dtheta, n=u_n;
    return -5*r*r*(-16.0-16.0*A*sin(n*T)-16.0*A+A*sin(n*T)*n*n);
  }
  case 4: {
    r=gmm::vect_norm2(p);
    return -45.0*r;
  }
  case 5:{
    double a = sol5_a, b = sol5_b;
    double T = atan2(p[1], p[0]);
    return 5.0*(4.0*Radius-9.0*r-4.0*Radius*sin(b*T)*a+Radius*sin(b*T)*b*b*a);
  }  
  case 6: {
    double sum = std::accumulate(p.begin(), p.end(), double(0));
    return 5.0 * (2*sin(sum)*(r*r*r-Radius*Radius*Radius)-6*r*sum*cos(sum)-9*r*sin(sum));
  }
  case 7: {
    double sum = std::accumulate(p.begin(), p.end(), double(0));
    return 5.0 * (2*sin(sum)*(r*r*r*r-Radius*Radius*Radius*Radius)-8*r*r*sum*cos(sum)-16*r*r*sin(sum));
  }
  case 8: {
    double rho=gmm::sqrt(p[0]*p[0]+ p[1]*p[1]+ p[2]*p[2]); 
    return 60*rho;
  }
  }
  return 0;
}

double ls_value(const base_node &p) {
  double R = Radius, r=gmm::vect_norm2(p);
  switch (u_version) {
  case 0: case 4: case 6: case 7:  return gmm::vect_norm2_sqr(p) - R*R;
  case 1: case 2: {
    double A=u_alpha, T=atan2(p[1],p[0])+dtheta, n=u_n;
    return -(R*R - r*r*(1+A*(1.0 + sin(n*T)))) / 15.0;
  }
  case 3: {
    double A=u_alpha, T=atan2(p[1], p[0])+dtheta, n=u_n;
    return -(R*R*R*R - r*r*r*r*(1+A*(1.0 + sin(n*T)))) / 4.0;      
  } 
  case 5:{
    double a = sol5_a, b = sol5_b;
    double T = atan2(p[1], p[0]);
    return r*r*(r - Radius*(1 - a*sin(b*T))) / (Radius*sqrt(Radius+Radius*Radius*b*b*a*a));
  } 
  case 8:{
    double rho=gmm::sqrt(p[0]*p[0]+ p[1]*p[1]+ p[2]*p[2]); 
    return (rho*rho*rho - Radius*Radius*Radius);
  }
    //  case 6:{
    //return r*r*r-R*R*R;
    //}
  }
  return 0;
}

/*
 * Test procedure
 */

void test_mim(getfem::mesh_im_level_set &mim, getfem::mesh_fem &mf_rhs,
	      bool bound) {
  if (!u_version) {
    unsigned N =  unsigned(mim.linked_mesh().dim());
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
  std::vector<scalar_type> U;
  size_type N;
  base_matrix gradU;
  bgeot::base_vector coeff;
  bgeot::multi_index sizes_;
public:
  level_set_unit_normal(const getfem::mesh_fem &mf_, const VECT1 &U_) 
    : mf(mf_), U(mf_.nb_basic_dof()), N(mf_.linked_mesh().dim()),
      gradU(1, N) {
    sizes_.resize(1); sizes_[0] = short_type(N);
    mf.extend_vector(U_, U);
  }
  const bgeot::multi_index &sizes(size_type) const {  return sizes_; }
  virtual void compute(getfem::fem_interpolation_context& ctx,
		       bgeot::base_tensor &t) {
    size_type cv = ctx.convex_num();
    coeff.resize(mf.nb_basic_dof_of_element(cv));
    gmm::copy
      (gmm::sub_vector(U,gmm::sub_index(mf.ind_basic_dof_of_element(cv))),
       coeff);
    ctx.pf()->interpolation_grad(ctx, coeff, gradU, 1);
    scalar_type norm = gmm::vect_norm2(gmm::mat_row(gradU, 0));
    for (size_type i = 0; i < N; ++i) t[i] = gradU(0, i) / norm;
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

  getfem::generic_assembly assem("t=comp(Base(#2).Grad(#1).NonLin(#3));"
				 "M(#2, #1)+= t(:,:,i,i)");
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
    assem("t=comp(Grad(#1).NonLin(#2).Grad(#1).NonLin(#2));"
	  "M(#1, #1)+= sym(t(:,i,i,:,j,j))");
  assem.push_mi(mim);
  assem.push_mf(mf);
  assem.push_mf(ls.get_mesh_fem());
  assem.push_mat(RM);
  assem.push_nonlinear_term(&nterm);
  assem.assembly(rg);
}



/**************************************************************/
/* assembling patch vector                                     */
/**************************************************************/

template<class VEC>
void asm_patch_vector     
(const VEC &RM_, const getfem::mesh_im &mim, const getfem::mesh_fem &mf_mult,
 const getfem::mesh_region &rg = getfem::mesh_region::all_convexes()) {
  VEC &RM = const_cast<VEC &>(RM_);
    
  getfem::generic_assembly assem("t=comp(Base(#1)); V(#1)+= t(:);");
  assem.push_mi(mim);
  assem.push_mf(mf_mult);
  assem.push_vec(RM);
  // assem.set("RM(#1)+=comp(Base(#1))()");
  assem.assembly(rg);

}
/**************************************************************/
/* assembling patch matrix                                     */
/**************************************************************/

template<class MAT>
void asm_stabilization_patch_term
(const MAT &RM_, const getfem::mesh &mesh, const getfem::mesh_im &mimbounddown, const getfem::mesh_fem &mf_mult,
 scalar_type ratio_size, scalar_type h ){
  MAT &M1 = const_cast<MAT &>(RM_);
  
  /****************************************************/
  /*        " select patch "                          */
  /****************************************************/
  
  
  
  // assemby patch vector
  const getfem::mesh_fem &mf_P0 = getfem::classical_mesh_fem(mesh, 0);
  size_type nbe = mf_P0.nb_dof();
  int ne = 0;
  double size_of_crack = 0;
  plain_vector Patch_Vector(nbe);
  asm_patch_vector(Patch_Vector, mimbounddown, mf_P0);
  // cout<<"patch_vectot="<< Patch_Vector<<endl;
  dal::bit_vector  Patch_element_list, Patch_dof_ind;
  for (size_type i = 0; i < nbe; ++i) {
    if (Patch_Vector[i] != scalar_type(0)){
      size_type cv = mf_P0.first_convex_of_basic_dof(i);
      Patch_element_list.add(cv);
      Patch_dof_ind.add(i);
      ne++;
      size_of_crack=size_of_crack + Patch_Vector[i];
    }
  }
  // cout<<"Path_element_list="<< Patch_element_list <<endl;
  //cout<<"Path_dof_ind="<< Patch_dof_ind <<endl;
  cout<<"number of element in patch="<< ne <<endl;
  std::vector<int> xadj(ne+1), adjncy, numelt(ne), part(ne);
  std::vector<int> vwgt(ne), indelt(mesh.convex_index().last_true()+1);
  std::vector<double> vwgtt(ne);
  int j = 0, k = 0;
  for (dal::bv_visitor ic(Patch_element_list); !ic.finished(); ++ic, j++) {
    numelt[j] = int(ic);
    indelt[ic] = j;
  }
  j = 0;
  for (dal::bv_visitor ic(Patch_element_list); !ic.finished(); ++ic, j++) {
    size_type ind_dof_patch = mf_P0.ind_basic_dof_of_element(ic)[0];
    vwgt[indelt[ic]] = int(1000000*Patch_Vector[ind_dof_patch]);
    vwgtt[indelt[ic]] = Patch_Vector[ind_dof_patch];
    xadj[j] = k;
    bgeot::mesh_structure::ind_set s;
    mesh.neighbours_of_convex(ic, s);
    for (bgeot::mesh_structure::ind_set::iterator it = s.begin(); it != s.end(); ++it) {
      if (Patch_element_list.is_in(*it)) { adjncy.push_back(indelt[*it]); ++k; }
    }
  }
  
  xadj[j] = k;
  // cout<<"xadj="<<xadj<<endl;
  //cout<<"adjncy="<<adjncy<<endl;
  //cout<<"vwgt="<<vwgt<<endl;
  
  cout<<"ratio size beween mesh and coarse mesh= "<< ratio_size <<endl;

  int nparts = 1;
#if GETFEM_HAVE_METIS
  nparts = int(size_of_crack/(ratio_size*h));
# ifdef GETFEM_HAVE_METIS_OLD_API
  std::vector<int> adjwgt(k); // actually Metis would also accept NULL instead of an empty array
  int wgtflag = 2, numflag = 0, edgecut;
  // float ubvec[1] = {1.03f};
  int options[5] = {0,0,0,0,0};
  //METIS_mCPartGraphKway(&ne, &ncon, &(xadj[0]), &(adjncy[0]), &(vwgt[0]), &(adjwgt[0]), &wgtflag,
  //		    &numflag, &nparts, &(ubvec[0]),  options, &edgecut, &(part[0]));
  //METIS_mCPartGraphRecursive(&ne, &ncon, &(xadj[0]), &(adjncy[0]), &(vwgt[0]), &(adjwgt[0]), &wgtflag,
  //			 &numflag, &nparts,  options, &edgecut, &(part[0]));
  //METIS_PartGraphKway(&ne, &(xadj[0]), &(adjncy[0]), &(vwgt[0]), &(adjwgt[0]), &wgtflag,
  //	  &numflag, &nparts, options, &edgecut, &(part[0]));
  METIS_PartGraphRecursive(&ne, &(xadj[0]), &(adjncy[0]), &(vwgt[0]), &(adjwgt[0]), &wgtflag,
                           &numflag, &nparts, options, &edgecut, &(part[0]));
# else
  int ncon = 1, edgecut;
  int options[METIS_NOPTIONS] = { 0 };
  METIS_SetDefaultOptions(options);
  METIS_PartGraphRecursive(&ne, &ncon, &(xadj[0]), &(adjncy[0]), &(vwgt[0]), 0, 0,
                           &nparts, 0, 0, options, &edgecut, &(part[0]));
# endif
  //cout<<"size_of_mesh="<<h<<endl;
  //cout<<"size_of_crack="<< size_of_crack <<endl;
  cout<<"nb_partition="<<nparts<<endl;
  // cout<<"partition="<<part<<endl;
  //cout<<"edgecut="<<edgecut<<endl;
#else
  GMM_ASSERT1(false, "METIS not linked");
#endif


  /**************************************************************/
  /*        Assembly matrices                                   */
  /**************************************************************/
  
  
  std::vector<double> size_patch(nparts);
  size_type nb_dof_mult=mf_mult.nb_dof();
  sparse_matrix M0(nb_dof_mult, nbe);
  getfem::asm_mass_matrix(M0, mimbounddown, mf_mult, mf_P0);
  
  for (size_type i=0; i < size_type(ne); i++) {
    size_patch[part[i]]= size_patch[part[i]] + vwgtt[i];	  
  }
  
  //cout<<"size_patch="<<size_patch<<endl;
  
  sparse_row_matrix MAT_aux(nparts, nb_dof_mult);
  for (size_type r=0; r < nbe; r++) {
    size_type cv = mf_P0.first_convex_of_basic_dof(r);
    gmm::add(gmm::mat_col(M0, r), gmm::mat_row(MAT_aux, part[indelt[cv]]));
  }
  
  sparse_row_matrix MAT_proj(nbe, nb_dof_mult);
  
  for (size_type r=0; r < nbe; r++) {
    size_type cv = mf_P0.first_convex_of_basic_dof(r);
    int p=part[indelt[cv]];
    gmm::copy(gmm::scaled(gmm::mat_row(MAT_aux, p), 1./size_patch[p]),
	      gmm::mat_row(MAT_proj, r));
  }
  
  gmm::mult(M0, MAT_proj, M1);
  
}





/*
 * Elementary extrapolation matrices
 */

void compute_mass_matrix_extra_element
(base_matrix &M, const getfem::mesh_im &mim, const getfem::mesh_fem &mf,
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
  // cout << "M = " << M << endl;
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
  
  // Load the mesh
  getfem::mesh mesh;
  // std::string MESH_FILE = PARAM.string_value("MESH_FILE", "Mesh file");
  //getfem::import_mesh(MESH_FILE, mesh);
  //unsigned N = mesh.dim();

  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type");
  bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Nomber of space steps "));
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			    PARAM.int_value("MESH_NOISED") != 0);
  //  base_small_vector tt(N); tt[1] = -0.5;
  //  mesh.translation(tt); 
  // center the mesh in (0, 0).
  base_node Pmin(N), Pmax(N);
  mesh.bounding_box(Pmin, Pmax);
  Pmin += Pmax; Pmin /= -2.0;
  // Pmin[N-1] = -Pmax[N-1];
  mesh.translation(Pmin);
  cout<<"Creating mesh done"<<endl; 
  scalar_type h = 2. * mesh.minimal_convex_radius_estimate();
  cout << "h = " << h << endl;
  // Level set definition
  unsigned lsdeg = unsigned(PARAM.int_value("LEVEL_SET_DEGREE", "level set degree"));
  bool simplify_level_set = 
    (PARAM.int_value("SIMPLIFY_LEVEL_SET",
		     "simplification or not of the level sets") != 0);
  Radius = PARAM.real_value("RADIUS", "Domain radius");
  getfem::level_set ls(mesh, dim_type(lsdeg));
  getfem::level_set lsup(mesh, dim_type(lsdeg), true), lsdown(mesh, dim_type(lsdeg), true);
  const getfem::mesh_fem &lsmf = ls.get_mesh_fem();
  for (unsigned i = 0; i < lsmf.nb_dof(); ++i) {
    lsup.values()[i] = lsdown.values()[i] = ls.values()[i]
      = ls_value(lsmf.point_of_basic_dof(i));
    if(N==2){
      lsdown.values(1)[i] = lsmf.point_of_basic_dof(i)[1];
      lsup.values(1)[i] = -lsmf.point_of_basic_dof(i)[1];
    }else{
      lsdown.values(2)[i] = lsmf.point_of_basic_dof(i)[2];
      lsup.values(2)[i] = -lsmf.point_of_basic_dof(i)[2];
    }
  }
  
  if (simplify_level_set) {
    scalar_type simplify_rate = std::min(0.03, 0.05 * sqrt(h));
    cout << "Simplification of level sets with rate: " <<
      simplify_rate << endl;
    ls.simplify(simplify_rate);
    lsup.simplify(simplify_rate);
    lsdown.simplify(simplify_rate); 
  }
  
  getfem::mesh_level_set mls(mesh), mlsup(mesh), mlsdown(mesh);
  mls.add_level_set(ls);
  mls.adapt();
  mlsup.add_level_set(lsup);
  mlsup.adapt();
  mlsdown.add_level_set(lsdown);
  mlsdown.adapt();
  
  getfem::mesh mcut;
  mls.global_cut_mesh(mcut);
  mcut.write_to_file("cut.mesh");
  
  // Integration method on the domain
  std::string IM = PARAM.string_value("IM", "Mesh file");
  std::string IMS = PARAM.string_value("IM_SIMPLEX", "Mesh file");
  int intins = getfem::mesh_im_level_set::INTEGRATE_INSIDE;
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
  getfem::mesh_im_level_set mimbounddown(mlsdown, intbound,
					 getfem::int_method_descriptor(IMS));
  mimbounddown.set_integration_method(mesh.convex_index(),
				      getfem::int_method_descriptor(IM));
  mimbounddown.adapt();
  getfem::mesh_im_level_set mimboundup(mlsup, intbound,
				       getfem::int_method_descriptor(IMS));
  mimboundup.set_integration_method(mesh.convex_index(),
				    getfem::int_method_descriptor(IM));
  mimboundup.adapt();
  
  // Finite element method for the unknown
  getfem::mesh_fem pre_mf(mesh);
  std::string FEM = PARAM.string_value("FEM", "finite element method");
  pre_mf.set_finite_element(mesh.convex_index(),
			    getfem::fem_descriptor(FEM));
  getfem::partial_mesh_fem mf(pre_mf);

  dal::bit_vector kept_dof;
  std::set<size_type> cols;
//   sparse_matrix BRB(pre_mf.nb_dof(), pre_mf.nb_dof());
//   getfem::asm_mass_matrix(BRB, mim, pre_mf);
//   cout << "Selecting dofs" << endl;
//   gmm::range_basis(BRB, cols);
//   for (std::set<size_type>::iterator it = cols.begin();
//        it != cols.end(); ++it)
//     kept_dof.add(*it);
  kept_dof = getfem::select_dofs_from_im(pre_mf, mim);
  cout << "Dof Selected : " << kept_dof.card() << " / "
       << pre_mf.nb_dof() << endl;  
  dal::bit_vector rejected_elt;
  for (dal::bv_visitor cv(mim.convex_index()); !cv.finished(); ++cv)
    if (mim.int_method_of_element(cv) == getfem::im_none())
      rejected_elt.add(cv);
  mf.adapt(kept_dof, rejected_elt);
  size_type nb_dof = mf.nb_dof();
  
  // Finite element method for the rhs
  getfem::mesh_fem mf_rhs(mesh);
  std::string FEMR = PARAM.string_value("FEM_RHS", "finite element method");
  mf_rhs.set_finite_element(mesh.convex_index(),
			    getfem::fem_descriptor(FEMR));
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  cout << "nb_dof_rhs = " << nb_dof_rhs << endl;
  
  // A P0 finite element
  const getfem::mesh_fem &mf_P0 = getfem::classical_mesh_fem(mesh, 0);
  
  // Finite element method for the multipliers
  getfem::mesh_fem mf_mult(mesh);
  std::string FEMM = PARAM.string_value("FEM_MULT", "fem for multipliers");
  mf_mult.set_finite_element(mesh.convex_index(),
				 getfem::fem_descriptor(FEMM));
  // getfem::partial_mesh_fem mf_mult(pre_mf_mult);
  cout << "NX="                      << PARAM.int_value("NX", "Nomber of space steps ")     << "\n"; 
  cout << "MESH_TYPE="               << MESH_TYPE                 << "\n";
  cout << "FEM_TYPE="                << FEM                       << "\n";
  cout << "FEM_MULT="                << FEMM                      << "\n";

  int stabilized_dirichlet =
    int(PARAM.int_value("STABILIZED_DIRICHLET", "Stabilized version of "
			"Dirichlet condition or not"));
  cout << "stabilized_dirichlet ="<< stabilized_dirichlet << "\n";

  // Range_basis call
  cols.clear();
  sparse_matrix BRBB;
  if (stabilized_dirichlet) {
    gmm::resize(BRBB, mf_mult.nb_dof(), mf_mult.nb_dof());
    getfem::asm_mass_matrix(BRBB, mimbounddown, mf_mult, mf_mult);
    // cout << BRBB << endl;
  } else {
    gmm::resize(BRBB, nb_dof, mf_mult.nb_dof());
    getfem::asm_mass_matrix(BRBB, mimbounddown, mf, mf_mult);
  }
  cout << "Selecting dofs for the multiplier" << endl;
  // cout << "nb_dof = " << mf.nb_dof() << endl;
  //cout << "nb_dof_mult = " << mf_mult.nb_dof() << endl;
  gmm::range_basis(BRBB, cols, 1e-12);
  mf_mult.reduce_to_basic_dof(cols);
  // penser à l'optimisation sur les mailles ...
  
  
  // kept_dof_mult = select_dofs_from_im(pre_mf_mult, mimbounddown, N-1);
  // mf_mult.adapt(kept_dof_mult, rejected_elt);
  size_type nb_dof_mult = mf_mult.nb_dof();
  cout << "nb_dof_mult = " << nb_dof_mult << endl;

  // Mass matrix on the boundary
  sparse_matrix B2(mf_rhs.nb_dof(), nb_dof);
  getfem::asm_mass_matrix(B2, mimboundup, mf_rhs, mf);
  sparse_matrix B(nb_dof_mult, nb_dof);
  getfem::asm_mass_matrix(B, mimbounddown, mf_mult, mf);

  scalar_type dir_gamma0(0);
  sparse_matrix MA(nb_dof_mult, nb_dof_mult), KA(nb_dof, nb_dof);
  sparse_matrix BA(nb_dof_mult, nb_dof);
  sparse_row_matrix M1(nb_dof_mult, nb_dof_mult);
  if (stabilized_dirichlet > 0) {
    
    sparse_row_matrix E1(nb_dof, nb_dof);
    
    if (stabilized_dirichlet == 3) {
      
      std::string datafilename;
      datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
      mesh.write_to_file(datafilename + ".mesh");
      cout<<"save mesh done"<<endl;
      scalar_type tpa = PARAM.real_value("TPA", "Type of patch assembly");

      if(tpa){
	
	scalar_type ratio_size = PARAM.real_value("RATIO_GR_MESH", "ratio size between mesh and patches");
	//       cout<<"ratio size beween mesh and coarse mesh= "<< ratio_size <<endl;
	
	asm_stabilization_patch_term(M1, mesh, mimbounddown, mf_mult, ratio_size, h);
       
	//	cout<<'patch matrix='<<M1<<endl; 
      }else{
      
      
      /****************************************************/
      /*        " select patch "                          */
      /****************************************************/
      


      // assemby patch vector
      size_type nbe = mf_P0.nb_dof();
      int ne = 0;
      double size_of_crack = 0;
      plain_vector Patch_Vector(nbe);
      asm_patch_vector(Patch_Vector, mimbounddown, mf_P0);
      // cout<<"patch_vectot="<< Patch_Vector<<endl;
      dal::bit_vector  Patch_element_list, Patch_dof_ind;
      for (size_type i = 0; i < nbe; ++i) {
	if (Patch_Vector[i] != scalar_type(0)){
	  size_type cv = mf_P0.first_convex_of_basic_dof(i);
	  Patch_element_list.add(cv);
	  Patch_dof_ind.add(i);
	  ne++;
	  size_of_crack=size_of_crack + Patch_Vector[i];
	}
      }
      // cout<<"Path_element_list="<< Patch_element_list <<endl;
      //cout<<"Path_dof_ind="<< Patch_dof_ind <<endl;
      cout<<"number of element in patch="<< ne <<endl;
      std::vector<int> xadj(ne+1), adjncy, numelt(ne), part(ne);
      std::vector<int> vwgt(ne), indelt(mesh.convex_index().last_true()+1);
      std::vector<double> vwgtt(ne);
      int j = 0, k = 0;
      for (dal::bv_visitor ic(Patch_element_list); !ic.finished(); ++ic, j++) {
	numelt[j] = int(ic);
	indelt[ic] = j;
      }
      j = 0;
      for (dal::bv_visitor ic(Patch_element_list); !ic.finished(); ++ic, j++) {
	size_type ind_dof_patch = mf_P0.ind_basic_dof_of_element(ic)[0];
	vwgt[indelt[ic]] = int(1000000*Patch_Vector[ind_dof_patch]);
	vwgtt[indelt[ic]] = Patch_Vector[ind_dof_patch];
	xadj[j] = k;
	bgeot::mesh_structure::ind_set s;
	mesh.neighbours_of_convex(ic, s);
	for (bgeot::mesh_structure::ind_set::iterator it = s.begin(); it != s.end(); ++it) {
	  if (Patch_element_list.is_in(*it)) { adjncy.push_back(indelt[*it]); ++k; }
	}
      }

      xadj[j] = k;
      // cout<<"xadj="<<xadj<<endl;
      //cout<<"adjncy="<<adjncy<<endl;
      //cout<<"vwgt="<<vwgt<<endl;

      scalar_type ratio_size = PARAM.real_value("RATIO_GR_MESH", "ratio size between mesh and patches");
      cout<<"ratio size beween mesh and coarse mesh= "<< ratio_size <<endl;

      int nparts = 1;
#if GETFEM_HAVE_METIS
      nparts = int(size_of_crack/(ratio_size*h));
# ifdef GETFEM_HAVE_METIS_OLD_API
      std::vector<int> adjwgt(k); // actually Metis would also accept NULL instead of an empty array
      int wgtflag = 2, numflag = 0, edgecut;
      int options[5] = {0,0,0,0,0};
      // float ubvec[1] = {1.03f};
      //METIS_mCPartGraphKway(&ne, &ncon, &(xadj[0]), &(adjncy[0]), &(vwgt[0]), &(adjwgt[0]), &wgtflag,
      //		    &numflag, &nparts, &(ubvec[0]),  options, &edgecut, &(part[0]));
      //METIS_mCPartGraphRecursive(&ne, &ncon, &(xadj[0]), &(adjncy[0]), &(vwgt[0]), &(adjwgt[0]), &wgtflag,
      //			 &numflag, &nparts,  options, &edgecut, &(part[0]));
      //METIS_PartGraphKway(&ne, &(xadj[0]), &(adjncy[0]), &(vwgt[0]), &(adjwgt[0]), &wgtflag,
      //	  &numflag, &nparts, options, &edgecut, &(part[0]));
      METIS_PartGraphRecursive(&ne, &(xadj[0]), &(adjncy[0]), &(vwgt[0]), &(adjwgt[0]), &wgtflag,
                               &numflag, &nparts, options, &edgecut, &(part[0]));
# else
      int ncon = 1, edgecut;
      int options[METIS_NOPTIONS] = { 0 };
      METIS_SetDefaultOptions(options);
      METIS_PartGraphRecursive(&ne, &ncon, &(xadj[0]), &(adjncy[0]), &(vwgt[0]), 0, 0,
                               &nparts, 0, 0, options, &edgecut, &(part[0]));
# endif
      //cout<<"size_of_mesh="<<h<<endl;
      //cout<<"size_of_crack="<< size_of_crack <<endl;
      cout<<"nb_partition="<<nparts<<endl;
      // cout<<"partition="<<part<<endl;
      //cout<<"edgecut="<<edgecut<<endl;
#else
  GMM_ASSERT1(false, "METIS not linked");
#endif


      /**************************************************************/
      /*        Assembly matrices                                   */
      /**************************************************************/


      std::vector<double> size_patch(nparts);

      sparse_matrix M0(nb_dof_mult, nbe);
      getfem::asm_mass_matrix(M0, mimbounddown, mf_mult, mf_P0);

      for (size_type i=0; i < size_type(ne); i++) {
	size_patch[part[i]]= size_patch[part[i]] + vwgtt[i];	  
      }

      //cout<<"size_patch="<<size_patch<<endl;

      sparse_row_matrix MAT_aux(nparts, nb_dof_mult);
      for (size_type r=0; r < nbe; r++) {
	size_type cv = mf_P0.first_convex_of_basic_dof(r);
	gmm::add(gmm::mat_col(M0, r), gmm::mat_row(MAT_aux, part[indelt[cv]]));
      }
      
      sparse_row_matrix MAT_proj(nbe, nb_dof_mult);

      for (size_type r=0; r < nbe; r++) {
	size_type cv = mf_P0.first_convex_of_basic_dof(r);
	int p=part[indelt[cv]];
	gmm::copy(gmm::scaled(gmm::mat_row(MAT_aux, p), 1./size_patch[p]),
		  gmm::mat_row(MAT_proj, r));
      }
     
      gmm::mult(M0, MAT_proj, M1);


      }


  
//       sparse_matrix MAT_proj(nbe, nb_dof_mult);

//       for (int r=0; r<nbe; r++){
// 	for (int jj=0; jj<nb_dof_mult; jj++){
// 	  size_type cv = mf_P0.first_convex_of_basic_dof(r);
// 	  int p=part[indelt[cv]];
// 	  for (int i=0; i<ne; i++){
// 	    if ( part[i]== p) {
// 	      size_type ind_dof_patch= mf_P0.ind_basic_dof_of_element(numelt[i])[0];
// 	      MAT_proj(r,jj) += M0(jj,ind_dof_patch);
// 	    }
// 	  }
//  	  MAT_proj(r,jj) /=  size_patch[p];
// 	}
//       }
//       gmm::mult(M0, MAT_proj, M1);

      // cout<<"M0="<<M0<<endl;
      //cout<<"MAT_proj="<<MAT_proj<<endl;
      // cout<<"M1="<<M1<<endl;

    }//end stabilized_dirichlet == 3
    
    



    if (stabilized_dirichlet == 2) {
      
      /***************************************************************/ 
      /*           Computation of the extrapolation operator.        */
      /***************************************************************/

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
	size_type cv = mf_P0.first_convex_of_basic_dof(i);
	ratios[cv] = gmm::abs(MC1(i,i)) / gmm::abs(MC2(i,i));
	if (ratios[cv] > 0 && ratios[cv] < min_ratio) elt_black_list.add(cv);
      }
      
	
      sparse_matrix EO(mf.nb_basic_dof(), mf.nb_basic_dof());
      sparse_row_matrix T1(nb_dof, nb_dof);
      sparse_row_matrix EX(mf.nb_basic_dof(), mf.nb_basic_dof());
      asm_mass_matrix(EO, uncutmim, pre_mf);
//       sparse_matrix MM1(nbe, nb_dof);
//       getfem::asm_mass_matrix(MM1, mim, mf_P0, mf);

      for (size_type i = 0; i < mf.nb_basic_dof(); ++i) {
	bool found = false;
	
// 	gmm::linalg_traits<gmm::linalg_traits<sparse_matrix>::sub_col_type >
// 	  ::iterator it = gmm::vect_begin(gmm::mat_col(MM1, i)),
// 	             ite = gmm::vect_end(gmm::mat_col(MM1, i));
// 	for ( ;  it != ite; ++it)
// 	  if (!elt_black_list.is_in(it.index())) found = true;

	getfem::mesh::ind_cv_ct ct = mf.convex_to_basic_dof(i);
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
	compute_mass_matrix_extra_element(Mloc, uncutmim, pre_mf, i, cv2);
	for (size_type ii = 0; ii < gmm::mat_nrows(Mloc); ++ii) 
	  for (size_type jj = 0; jj < gmm::mat_ncols(Mloc); ++jj)
	    EX(mf.ind_basic_dof_of_element(i)[ii],
	       mf.ind_basic_dof_of_element(cv2)[jj])
	      += Mloc(ii, jj);
      }

      gmm::copy(gmm::identity_matrix(), E1);
      gmm::copy(gmm::identity_matrix(), T1);
      for (dal::bv_visitor i(dof_black_list); !i.finished(); ++i)
	gmm::mult(gmm::transposed(mf.extension_matrix()), gmm::mat_row(EX, i),
		  gmm::mat_row(T1, i));

      plain_vector BE(mf.nb_basic_dof()), BS(mf.nb_basic_dof()), BBS(nb_dof);
      for (dal::bv_visitor i(dof_black_list); !i.finished(); ++i) {
	BE[i] = scalar_type(1);
	// TODO: store LU decomp.
	double rcond; 
	gmm::SuperLU_solve(EO, BS, BE, rcond);
	gmm::mult(gmm::transposed(mf.extension_matrix()), BS, BBS);
	gmm::mult(gmm::transposed(T1), BBS, gmm::mat_row(E1, i));
	BE[i] = scalar_type(0);
      }
      gmm::clean(E1, 1e-13);
      

      cout << "Extrapolation operator computed" << endl;

    }

    dir_gamma0 = PARAM.real_value("DIRICHLET_GAMMA0",
				  "Stabilization parameter for "
				  "Dirichlet condition");
   
    if (stabilized_dirichlet == 3) {
      getfem::asm_mass_matrix(MA, mimbounddown, mf_mult);
      gmm::scale(M1, dir_gamma0 );
      gmm::scale(MA, -dir_gamma0 );
      gmm::add(M1, MA);

    }else{
      getfem::asm_mass_matrix(MA, mimbounddown, mf_mult);
      asm_stabilization_mixed_term(BA, mimbounddown, mf, mf_mult, lsdown);
      asm_stabilization_symm_term(KA, mimbounddown, mf, lsdown);
      gmm::scale(MA, -dir_gamma0 * h);
      gmm::scale(BA, -dir_gamma0 * h);
      gmm::scale(KA, -dir_gamma0 * h);
    }
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
  test_mim(mimbounddown, mf_rhs, true);

  /***************************************************************/ 
  /*          New brick system                                   */
  /***************************************************************/
  getfem::model model;
  
  model.add_fem_variable("u", mf);
  model.add_fem_variable("Lambda", mf_mult);
  getfem::add_Laplacian_brick(model, mim, "u");
  
  
  if (stabilized_dirichlet > 0){
    getfem::add_explicit_matrix(model, "Lambda", "Lambda", MA);
    if (stabilized_dirichlet != 3) {
      getfem::add_explicit_matrix(model, "u", "u", KA);
    }
  }
  
  //Volumic source term
  plain_vector F(nb_dof_rhs);
  getfem::interpolation_function(mf_rhs, F, rhs);
  model.add_initialized_fem_data("VolumicData", mf_rhs, F);
  getfem::add_source_term_brick(model, mim, "u", "VolumicData");
  
  // Neumann condition
  getfem::interpolation_function(mf_rhs, F, g_exact);
  plain_vector R(nb_dof);
  gmm::mult(gmm::transposed(B2), F, R);
  getfem::add_explicit_rhs(model, "u", R);
  
  // Dirichlet condition
  
  getfem::add_constraint_with_multipliers(model, "u", "Lambda", B, plain_vector(nb_dof_mult));
 
  //Solving the problem
  cout<< "Stabilized_parameter="<< dir_gamma0 <<endl;
  gmm::iteration iter(1e-9, 1, 40000);

  //getfem::standard_solve(model, iter);


  getfem::default_newton_line_search lsear;
  getfem::standard_solve(model, iter, getfem::rselect_linear_solver(model,"mumps"), lsear);


  plain_vector U(nb_dof);
  gmm::copy(model.real_variable("u"), U);
  plain_vector LAMBDA(nb_dof_mult);
  gmm::copy(model.real_variable("Lambda"), LAMBDA);
  
  // interpolation of the solution on mf_rhs
  GMM_ASSERT1(!mf_rhs.is_reduced(), "To be adapted");
  plain_vector Uint(nb_dof_rhs), Vint(nb_dof_rhs), Eint(nb_dof_rhs),lambda_int(nb_dof_rhs);
  getfem::interpolation(mf, mf_rhs, U, Uint);
  //getfem::interpolation(mf_mult, mf_rhs, LAMBDA, lambda_int);
  for (size_type i = 0; i < nb_dof_rhs; ++i) {
    Vint[i] = u_exact(mf_rhs.point_of_basic_dof(i));
    Eint[i] = gmm::abs(Uint[i] - Vint[i]);
  }
  /***************************************************************/ 
  /*          computation of error on u.                         */
  /***************************************************************/

  double errmax = 0.0, exactmax = 0.0;
  for (size_type i = 0; i < nb_dof_rhs; ++i)
    if (ls_value(mf_rhs.point_of_basic_dof(i)) <= 0.0) {
      errmax = std::max(errmax, Eint[i]);
      exactmax = std::max(exactmax, Vint[i]);
    }
    else Eint[i] = 0.0;

  mf_rhs.write_to_file("xfem_dirichlet.mf", true);
  mf_mult.write_to_file("xfem_dirichlet.mfmult", true);
  gmm::vecsave("xfem_dirichlet_exact.U", Vint);
  gmm::vecsave("xfem_dirichlet.U", Uint);
  gmm::vecsave("xfem_dirichlet.L", LAMBDA);
  gmm::vecsave("xfem_dirichlet.map_error", Eint);

  // cout << "Linfty error: " << 100.0 * errmax / exactmax << "%" << endl;
  cout << "L2 error: " << 100.0
    * getfem::asm_L2_dist(mim, mf, U, mf_rhs, Vint)
    / getfem::asm_L2_norm(mim, mf_rhs, Vint) << "%" << endl;
  cout << "H1 error: " << 100.0
    * getfem::asm_H1_dist(mim, mf, U, mf_rhs, Vint)
    / getfem::asm_H1_norm(mim, mf_rhs, Vint) << "%" << endl;

  /*********************************************************************/  
  /*         computation of error on multipliers.                      */
  /*********************************************************************/
  gmm::resize(BA, nb_dof_mult, nb_dof_rhs); gmm::clear(BA);
  gmm::resize(KA, nb_dof_rhs, nb_dof_rhs);  gmm::clear(KA);
  gmm::resize(B, nb_dof_mult, nb_dof_mult); gmm::clear(B);
  asm_stabilization_mixed_term(BA, mimbounddown, mf_rhs, mf_mult, lsdown);
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
  if(1){
    
    /*********************************************************************/  
    /*         exporting solution in vtk format.                         */
    /*********************************************************************/
    {
      getfem::vtk_export exp("xfem_dirichlet.vtk", (2==1));
    exp.exporting(mf); 
    exp.write_point_data(mf, U, "solution");
    cout << "export done, you can view the data file with (for example)\n"
      "mayavi -d xfem_dirichlet.vtk -f WarpScalar -m BandedSurfaceMap "
      "-m Outline\n";
    }
    // exporting error in vtk format.
    {
      getfem::vtk_export exp("xfem_dirichlet_error.vtk", (2==1));
      exp.exporting(mf_rhs); 
      exp.write_point_data(mf_rhs, Eint, "error");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d xfem_dirichlet_error.vtk -f WarpScalar -m BandedSurfaceMap "
	"-m Outline\n";
    }
    // exporting multipliers in vtk format.
    {
      getfem::vtk_export exp("xfem_dirichlet_mult.vtk", (2==1));
      exp.exporting(mf_mult); 
      exp.write_point_data(mf_mult, LAMBDA, "multipliers");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d xfem_dirichlet_mult.vtk -f WarpScalar -m BandedSurfaceMap "
	"-m Outline\n";
    }
    
    lsmf.write_to_file("xfem_dirichlet_ls.mf", true);
    
    gmm::vecsave("xfem_dirichlet_ls.U", ls.values());
    
    //save solution
    
    /*********************************************************************/  
    /*         Save the solution.                                        */
    /*********************************************************************/
    mf.write_to_file("xfem_dirichlet.mfsol", true);
    gmm::vecsave("xfem_dirichlet.Usol", U);
    
    mf_mult.write_to_file("xfem_dirichlet.mf_mult", true);
    gmm::vecsave("xfem_dirichlet.Lsol", LAMBDA);
    cout << "saving solution done"<<endl;
    
    unsigned nrefine = mf.linked_mesh().convex_index().card() < 200 ? 32 : 4;
    /*********************************************************************/  
    /*         Creationg Slice                                           */
    /*********************************************************************/
    
    if (N==2) {
      cout << "saving the slice solution for matlab\n";
      getfem::stored_mesh_slice sl, sl0,sll;
      
      
      getfem::mesh_slicer slicer(mf.linked_mesh());
      getfem::slicer_build_stored_mesh_slice sbuild(sl);
      getfem::mesh_slice_cv_dof_data<plain_vector> mfU(mf, U);
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
      slicer2.push_back_action(sbuildl); // store it into sll
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
    }else{
      cout << "saving the slice solution for matlab\n";
      // Create slice of the sphere to plot the Solution in the half sphere
      // getfem::slicer_boundary exbond(mf.linked_mesh());//extract boundary
      getfem::stored_mesh_slice  sl, sll, sl0, slU, slsph;
      
    
      
      getfem::mesh_slicer slicer(mf.linked_mesh());
      getfem::slicer_build_stored_mesh_slice sbuild(sl);
      getfem::mesh_slice_cv_dof_data<plain_vector> mfU(mf, U);
      getfem::slicer_isovalues iso(mfU, 0.0, 0);
      getfem::slicer_build_stored_mesh_slice sbuild0(sl0);
      
      slicer.push_back_action(sbuild);  // full slice in sl
      slicer.push_back_action(iso);     // extract isosurface 0
      slicer.push_back_action(sbuild0); // store it into sl0
      slicer.exec(nrefine, mf.convex_index());
      


      getfem::mesh_slicer slicer2(ls.get_mesh_fem().linked_mesh());
      getfem::mesh_slice_cv_dof_data<plain_vector> 
	mfL(ls.get_mesh_fem(), ls.values());
      getfem::slicer_isovalues  iso2(mfL, 0.0, 0);
      //getfem::slicer_half_space hs(base_node(0,0,0), base_node(0,1,0),-1);
      // getfem::slicer_intersect  sect(iso2,hs);
      getfem::slicer_build_stored_mesh_slice sbuildlu(slU);
      slicer2.push_back_action(sbuild);                 // Full slice in slfulll
      slicer2.push_back_action(iso2);                   // extract isosurface 0
      // slicer2.push_back_action(hs);                  // cut with half space
      // slicer2.push_back_action(exbond);              // extract boundary
      slicer2.push_back_action(sbuildlu);               // store it into slU
      slicer2.exec(nrefine, ls.get_mesh_fem().convex_index());
      
      sl.write_to_file("xfem_dirichlet.sl", true);
      sl0.write_to_file("xfem_dirichlet.sl0", true);
      slU.write_to_file("xfem_dirichlet.slU", true);
      plain_vector  UU(slU.nb_points());
      
      slU.interpolate(mf, U, UU);
      gmm::vecsave("xfem_dirichlet.sl_U", UU);
      
   
      
      // Create slice of the sphere to plot the multiplier at the boundary
      
      getfem::mesh_slicer slicer3(mf.linked_mesh());
      //getfem::mesh_slice_cv_dof_data<plain_vector> 
      //  mfL(ls.get_mesh_fem(), ls.values());
      getfem::slicer_isovalues  iso3(mfL, 0.0, 0);
      getfem::slicer_half_space hs(base_node(0,0,0), base_node(0,1,0),-1);
      //    getfem::slicer_intersect  sect3(iso3,hs);
      getfem::slicer_build_stored_mesh_slice sbuildl(sll);
      slicer3.push_back_action(sbuild);            // Full slice in sl
      slicer3.push_back_action(iso3);              // extract isosurface 0
      slicer3.push_back_action(hs);                // cut with half space
      slicer3.push_back_action(sbuildl);           // store it into sll
      slicer3.exec(nrefine, mf.convex_index());
      
      sll.write_to_file("xfem_dirichlet.sll", true);
      plain_vector  LL(sll.nb_points()), UUb(sll.nb_points());
      
      sll.interpolate(mf, U, UUb);
      gmm::vecsave("xfem_dirichlet.slUb", UUb);
      
      gmm::scale(LAMBDA, 0.005);
      sll.interpolate(mf_mult, LAMBDA, LL);
      gmm::vecsave("xfem_dirichlet.slL", LL);
      
      getfem::slicer_boundary exbond1(mf.linked_mesh());//extract boundary
      getfem::mesh_slicer slicer4(mf.linked_mesh());
      //getfem::slicer_build_stored_mesh_slice sbuildfull(sl);
      getfem::slicer_sphere sph(base_node(0,0), Radius, -1 );
      // getfem::slicer_half_space hs1(base_node(0,0,0), base_node(0,1,0),-1);
      //getfem::slicer_intersect  sect4(sph,hs1);
      getfem::slicer_build_stored_mesh_slice sbuildsph(slsph);
      
      slicer4.push_back_action(sbuild);  // Full slice in sl
      slicer4.push_back_action(sph);         // extract sphere 
      //slicer4.push_back_action(hs1);         // cut with half space
      slicer4.push_back_action(exbond1);      // extract boundary
      slicer4.push_back_action(sbuildsph);   // store it into slsph
      
      slicer4.exec(nrefine, mf.convex_index());
      slsph.write_to_file("xfem_dirichlet.slsph", true);
    
      plain_vector UUs(slsph.nb_points()); 
      slsph.interpolate(mf, U, UUs);
      gmm::vecsave("xfem_dirichlet.slUsph", UUs);
    }
    
    /********************************************/
    /*exacte solution                           */
    /********************************************/
    plain_vector UE(nb_dof_rhs);
    plain_vector UEE(nb_dof);
    for (size_type i = 0; i < nb_dof_rhs; ++i) {
      UE[i] = u_exact(mf_rhs.point_of_basic_dof(i));
    }
    mf.write_to_file("xfem_dirichlet.mfE", true);
    getfem::interpolation(mf_rhs, mf, UE, UEE);
    gmm::vecsave("xfem_dirichlet_exact.UE", UEE);
    getfem::stored_mesh_slice sle, sl0e,slle;
    getfem::mesh_slicer slicere(mf.linked_mesh());
    getfem::slicer_build_stored_mesh_slice sbuilde(sle);
    getfem::mesh_slice_cv_dof_data<plain_vector> mfUe(mf, UEE);
    getfem::slicer_isovalues isoe(mfUe, 0.0, 0);
    getfem::slicer_build_stored_mesh_slice sbuild0e(sl0e);
    
    slicere.push_back_action(sbuilde);  // full slice in sle
    slicere.push_back_action(isoe);     // extract isosurface 0
    slicere.push_back_action(sbuild0e); // store it into sl0e
    slicere.exec(nrefine, mf.convex_index());
    
    
    
    getfem::mesh_slicer slicer2e(mf.linked_mesh());
    getfem::mesh_slice_cv_dof_data<plain_vector> 
      mfLe(ls.get_mesh_fem(), ls.values());
    getfem::slicer_isovalues iso2e(mfLe, 0.0, 0);
    getfem::slicer_build_stored_mesh_slice sbuildle(slle);
    slicer2e.push_back_action(iso2e);     // extract isosurface 0
    slicer2e.push_back_action(sbuildle); // store it into sl0e
    slicer2e.exec(nrefine, mf.convex_index());
    
    sle.write_to_file("xfem_dirichlet.sle", true);
    sl0e.write_to_file("xfem_dirichlet.sl0e", true);
    slle.write_to_file("xfem_dirichlet.slle", true);
    plain_vector UUE(sle.nb_points());
    sle.interpolate(mf, UEE, UUE);
    gmm::vecsave("xfem_dirichlet.slUE", UUE);
  }
  
  return 0; 
}
