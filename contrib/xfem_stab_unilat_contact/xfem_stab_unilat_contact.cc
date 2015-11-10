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

/**
 * Goal : stabilization of unilateral contact with Xfem.
 *
 * Research program.
 */

#include "gmm/gmm.h"
#include "getfem/getfem_assembling.h" /* import assembly methods      */
#include "getfem/getfem_export.h"     /* export functions             */
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_mesh_fem_product.h"
#include "getfem/getfem_mesh_fem_global_function.h"
#include "getfem/getfem_partial_mesh_fem.h"
#include "getfem/getfem_contact_and_friction_nodal.h"
#include "getfem/getfem_import.h"
#include "getfem/getfem_mesh_fem_sum.h"
#include "gmm/gmm_inoutput.h"

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

using std::endl; using std::cout; using std::cerr;

/* some GetFEM++ types that we will be using */
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

/******************************************************************************/
/* level set unit normal *                                                    */
/******************************************************************************/
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

/******************************************************************************/
/* h                                                    */
/******************************************************************************/
template<typename VECT1> class nonlin_h
  : public getfem::nonlinear_elem_term {
  const getfem::mesh &m;
  size_type N;
  bgeot::multi_index sizes_;
  size_type cv_old;
  scalar_type h;
public:
  nonlin_h(const getfem::mesh &m_) 
    : m(m_), N(m.dim())
      {
    sizes_.resize(1); sizes_[0] = 1;
    cv_old = size_type(-1); h = 0.;
  }
  const bgeot::multi_index &sizes(size_type) const {  return sizes_; }
  virtual void compute(getfem::fem_interpolation_context& ctx,
		       bgeot::base_tensor &t) {
    size_type cv = ctx.convex_num();
    if (cv != cv_old) {
      h = m.convex_radius_estimate(cv);
      cv_old = cv;
    }
    t[0] = 1./h;
  }
};




/******************************************************************************/
/* level set unit tangent                                                     */
/******************************************************************************/
template<typename VECT1> class level_set_unit_tang 
  : public getfem::nonlinear_elem_term {
  const getfem::mesh_fem &mf;
  std::vector<scalar_type> U;
  size_type N;
  base_matrix gradU;
  bgeot::base_vector coeff;
  bgeot::multi_index sizes_;
public:
  level_set_unit_tang(const getfem::mesh_fem &mf_, const VECT1 &U_) 
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
    GMM_ASSERT1(N == 2, "Sorry, to be done for N != 2");
    t[0] = gradU(0, 1) / norm; // this expression is correct only for the case of 2 dimensional
    t[1]=-gradU(0, 0) / norm; // this expression is correct only for the case of 2 dimensional
  }
};




/**********************************************************************************/
/* asembling stabilised mixed term                                                */
/**********************************************************************************/

template<class MAT>
void asm_stabilization_mixed_term
(const MAT &RM_, const getfem::mesh_im &mim, const getfem::mesh_fem &mf, 
 const getfem::mesh_fem &mf_mult, getfem::level_set &ls,
 scalar_type lambda, scalar_type mu,
 const getfem::mesh_region &rg = getfem::mesh_region::all_convexes()) {
  MAT &RM = const_cast<MAT &>(RM_);
  
  level_set_unit_normal<std::vector<scalar_type> >
    nterm(ls.get_mesh_fem(), ls.values());


  plain_vector LAMBDA(1, lambda), MU(1, mu);
  
  getfem::generic_assembly assem("lambda=data$1(1); mu=data$2(1);"
                                 "t=comp(Base(#2).NonLin(#3).vGrad(#1).NonLin(#3));"
				 "M(#2, #1)+= t(:,i,:,j,j,i).lambda(1)"
				 "+t(:,i,:,i,j,j).mu(1)*2");
  assem.push_mi(mim);
  assem.push_mf(mf);
  assem.push_mf(mf_mult);
  assem.push_mf(ls.get_mesh_fem());
  assem.push_data(LAMBDA);
  assem.push_data(MU);
  assem.push_mat(RM);
  assem.push_nonlinear_term(&nterm);
  
  
  gmm::clear(RM); 
  ls.set_shift(1e-7);
  assem.assembly(rg);
  ls.set_shift(-1e-7);
  assem.assembly(rg);
  ls.set_shift(0.);
}

/**********************************************************************************/
/* asembling stabilised tangent mixed term                                                */
/**********************************************************************************/

template<class MAT>
void asm_stabilization_tang_mixed_term
(const MAT &RM_, const getfem::mesh_im &mim, const getfem::mesh_fem &mf, 
 const getfem::mesh_fem &mf_mult, getfem::level_set &ls,
 scalar_type lambda, scalar_type mu,
 const getfem::mesh_region &rg = getfem::mesh_region::all_convexes()) {
  MAT &RM = const_cast<MAT &>(RM_);
  
  level_set_unit_normal<std::vector<scalar_type> >
    nterm(ls.get_mesh_fem(), ls.values());

  level_set_unit_tang<std::vector<scalar_type> >
    taterm(ls.get_mesh_fem(), ls.values());

  plain_vector LAMBDA(1, lambda), MU(1, mu);
  
  getfem::generic_assembly assem("lambda=data$1(1); mu=data$2(1);"
                                 "t=comp(Base(#2).NonLin$2(#3).vGrad(#1).NonLin$1(#3));"
				 "M(#2, #1)+= t(:,i,:,i,j,j).mu(1)"
				 "+t(:,i,:,j,i,j).mu(1)");
  
  assem.push_mi(mim);
  assem.push_mf(mf);
  assem.push_mf(mf_mult);
  assem.push_mf(ls.get_mesh_fem());
  assem.push_data(LAMBDA);
  assem.push_data(MU);
  assem.push_mat(RM);
  assem.push_nonlinear_term(&nterm);
  assem.push_nonlinear_term(&taterm);
  

  gmm::clear(RM); 
  ls.set_shift(1e-7);
  assem.assembly(rg);
  ls.set_shift(-1e-7);
  assem.assembly(rg);
  ls.set_shift(0.);
}

/**************************************************************************************/
/* asembling masse matrix for  mixed term                                             */
/**************************************************************************************/

template<class MAT>
void asm_mass_matrix_mixed_term
(const MAT &RM_, const getfem::mesh_im &mim, const getfem::mesh_fem &mf, 
 const getfem::mesh_fem &mf_mult, getfem::level_set &ls,
 const getfem::mesh_region &rg = getfem::mesh_region::all_convexes()) {
  MAT &RM = const_cast<MAT &>(RM_);
  
  level_set_unit_normal<std::vector<scalar_type> >
    nterm(ls.get_mesh_fem(), ls.values());
  
  getfem::generic_assembly assem("t=comp(Base(#2).vBase(#1).NonLin(#3));"
				 "M(#2,#1)+= t(:,:,i,i)");
  assem.push_mi(mim);
  assem.push_mf(mf);
  assem.push_mf(mf_mult);
  assem.push_mf(ls.get_mesh_fem());
  assem.push_mat(RM);
  assem.push_nonlinear_term(&nterm);
  
  
  gmm::clear(RM);
  ls.set_shift(1e-7);
  assem.assembly(rg);
  gmm::scale(RM, scalar_type(-1));
  ls.set_shift(-1e-7);
  assem.assembly(rg);
  ls.set_shift(0.);
  // cout << "RM = " << RM << endl;

}

/**************************************************************************************/
/* asembling masse matrix for tangent mixed term                                      */
/**************************************************************************************/

template<class MAT>
void asm_mass_matrix_tang_mixed_term
(const MAT &RM_, const getfem::mesh_im &mim, const getfem::mesh_fem &mf, 
 const getfem::mesh_fem &mf_mult, getfem::level_set &ls,
 const getfem::mesh_region &rg = getfem::mesh_region::all_convexes()) {
  MAT &RM = const_cast<MAT &>(RM_);
  
  level_set_unit_tang<std::vector<scalar_type> >
    taterm(ls.get_mesh_fem(), ls.values());
  
  getfem::generic_assembly assem("t=comp(Base(#2).vBase(#1).NonLin(#3));"
				 "M(#2,#1)+= t(:,:,i,i)");
  assem.push_mi(mim);
  assem.push_mf(mf);
  assem.push_mf(mf_mult);
  assem.push_mf(ls.get_mesh_fem());
  assem.push_mat(RM);
  assem.push_nonlinear_term(&taterm);
  
  
  gmm::clear(RM);
  ls.set_shift(1e-7);
  assem.assembly(rg);
  gmm::scale(RM, scalar_type(-1));
  ls.set_shift(-1e-7);
  assem.assembly(rg);
  ls.set_shift(0.);
  //cout << "RM = " << RM << endl;

}

/******************************************************************************************/
/* assembling stabilised symetric term                                                     */
/******************************************************************************************/

template<class MAT>
void asm_stabilization_symm_term
(const MAT &RM_, const getfem::mesh_im &mim, const getfem::mesh_fem &mf, 
 getfem::level_set &ls, scalar_type lambda, scalar_type mu,
 const getfem::mesh_region &rg = getfem::mesh_region::all_convexes()) {
  MAT &RM = const_cast<MAT &>(RM_);
  
  level_set_unit_normal<std::vector<scalar_type> >
    nterm(ls.get_mesh_fem(), ls.values());
  
  plain_vector LAMBDA(1, lambda), MU(1, mu);

  getfem::generic_assembly
    assem("lambda=data$1(1); mu=data$2(1);"
          "t=comp(NonLin(#2).vGrad(#1).NonLin(#2).NonLin(#2).vGrad(#1).NonLin(#2));"
	  "M(#1, #1)+= sym(t(i,:,j,j,i,k,:,l,l,k).lambda(1).lambda(1)"
          "+t(i,:,j,j,i,k,:,k,l,l).lambda(1).mu(1)*2"
	  "+t(i,:,i,j,j,k,:,l,l,k).lambda(1).mu(1)*2"
          "+t(i,:,i,j,j,k,:,k,l,l).mu(1).mu(1)*4)");
  assem.push_mi(mim);
  assem.push_mf(mf);
  assem.push_mf(ls.get_mesh_fem());
  assem.push_data(LAMBDA);
  assem.push_data(MU);
  assem.push_mat(RM);
  assem.push_nonlinear_term(&nterm);

  gmm::clear(RM);
  ls.set_shift(1e-7);
  assem.assembly(rg);
  ls.set_shift(-1e-7);
  assem.assembly(rg);
  ls.set_shift(0.);
}

/******************************************************************************************/
/* asembling stabilised symetric tangent term                                             */
/******************************************************************************************/

template<class MAT>
void asm_stabilization_symm_tang_term
(const MAT &RM_, const getfem::mesh_im &mim, const getfem::mesh_fem &mf, 
 getfem::level_set &ls, scalar_type lambda, scalar_type mu,
 const getfem::mesh_region &rg = getfem::mesh_region::all_convexes()) {
  MAT &RM = const_cast<MAT &>(RM_);
  
  level_set_unit_normal<std::vector<scalar_type> >
    nterm(ls.get_mesh_fem(), ls.values());

 level_set_unit_tang<std::vector<scalar_type> >
    taterm(ls.get_mesh_fem(), ls.values());
  
  plain_vector LAMBDA(1, lambda), MU(1, mu);

  getfem::generic_assembly
    assem("lambda=data$1(1); mu=data$2(1);"
          "t=comp(NonLin$2(#2).vGrad(#1).NonLin$1(#2).NonLin$2(#2).vGrad(#1).NonLin$1(#2));"
	  "M(#1, #1)+= sym(t(i,:,i,j,j,k,:,k,l,l).mu(1).mu(1)"
	  "+t(i,:,i,j,j,k,:,l,k,l).mu(1).mu(1)"
	  "+t(i,:,j,i,j,k,:,k,l,l).mu(1).mu(1)"
          "+t(i,:,j,i,j,k,:,l,k,l).mu(1).mu(1))");
  assem.push_mi(mim);
  assem.push_mf(mf);
  assem.push_mf(ls.get_mesh_fem());
  assem.push_data(LAMBDA);
  assem.push_data(MU);
  assem.push_mat(RM);
  assem.push_nonlinear_term(&nterm);
  assem.push_nonlinear_term(&taterm);


  gmm::clear(RM);
  ls.set_shift(1e-7);
  assem.assembly(rg);
  ls.set_shift(-1e-7);
  assem.assembly(rg);
  ls.set_shift(0.);
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
 scalar_type ratio_size, scalar_type
#ifdef GETFEM_HAVE_METIS
h
#endif
){
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
   cout<<"Path_element_list="<< Patch_element_list <<endl;
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
#ifdef GETFEM_HAVE_METIS
  nparts = int(size_of_crack/(ratio_size*h));
# ifdef GETFEM_HAVE_METIS_OLD_API
  std::vector<int> adjwgt(k); // actually Metis would also accept NULL instead of an empty array
  int wgtflag = 2, numflag = 0, edgecut;
  int  options[5] = {0,0,0,0,0};
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
  cout<<"size_of_crack="<< size_of_crack <<endl;
  cout<<"nb_partition="<<nparts<<endl;
  cout<<"partition="<<part<<endl;
  cout<<"edgecut="<<edgecut<<endl;
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


/**************************************************************************************/
/* asembling masse matrix for  inf-sup condition                                            */
/**************************************************************************************/

template<class MAT>
void asm_mass_matrix_for_inf_sup
(const MAT &RM_, const getfem::mesh_im &mim, const getfem::mesh_fem &mf_mult,
 getfem::level_set &ls, const getfem::mesh_region &rg = getfem::mesh_region::all_convexes()) {
  MAT &RM = const_cast<MAT &>(RM_);
  
//  level_set_unit_normal<std::vector<scalar_type> >
//    nterm(ls.get_mesh_fem(), ls.values());
  nonlin_h<std::vector<scalar_type> >
    nlinh(ls.get_mesh_fem().linked_mesh());
  
  
  getfem::generic_assembly assem("t=comp(NonLin(#2).Base(#1).Base(#1));"
				 "M(#1,#1)+= t(i,:,:)");
  assem.push_mi(mim);
  assem.push_mf(mf_mult);
  assem.push_mf(ls.get_mesh_fem());
  assem.push_mat(RM);
 // assem.push_nonlinear_term(&nterm);
  assem.push_nonlinear_term(&nlinh);
  
  gmm::clear(RM);
  assem.assembly(rg);

}



/**************************************************************************/
/*structure of contact problem                                 */
/**************************************************************************/

struct unilateral_contact_problem {
  
  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN1_BOUNDARY_NUM = 1, NEUMANN2_BOUNDARY_NUM=2, NEUMANN3_BOUNDARY_NUM=3, NEUMANN4_BOUNDARY_NUM=4};
  getfem::mesh mesh;  /* the mesh */
  getfem::level_set ls;      /* The two level sets defining the crack.       */
  getfem::mesh_level_set mls, mls_bound;       /* two meshs level set                    */
  getfem::mesh_im_level_set mim, mimbound;    /* the integration methods.              */
  getfem::mesh_fem mf_pre_u, mf_dir;
  getfem::mesh_fem_level_set mfls_u;
  getfem::mesh_fem_global_function mf_sing_u;
  
  getfem::mesh_fem mf_partition_of_unity;
  getfem::mesh_fem_product mf_product;
  getfem::mesh_fem_sum mf_u_sum;
  
  getfem::mesh_fem mf_contt; /* mesh_fem for the contact multiplier     */
  getfem::mesh_fem_level_set mfls_cont;   /* mesh_fem for the multiplier contact enriched with H.   */
  getfem::mesh_fem_sum mf_cont_sum;
  
  
  base_small_vector cracktip;
  
  
  getfem::mesh_fem& mf_u() { return mf_u_sum; }
  getfem::mesh_fem& mf_cont() { return mf_contt; }
  getfem::mesh_fem& mf_pre_uu() { return mf_pre_u; }
  
  
  scalar_type mu, lambda;    /* Lame coeff                   */
  
  int dgr;                   /* Order of enrichement for u */ 
  scalar_type stabilized_problem;         
  scalar_type  friction_coeff;
  scalar_type tresca_threshold;
  scalar_type ratio_size;
  
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  
  
  scalar_type residual;      /* max residual for the iterative solvers      */
  bool contact_only;
  bool Tresca_version;
  bool Test_of_coer, strmesh, rangeP_P;
  scalar_type cutoff_radius, cutoff_radius1, cutoff_radius0, enr_area_radius;
  
  size_type cutoff_func;
  
  typedef enum { NO_ENRICHMENT=0, 
		 FIXED_ZONE=1, 
		 GLOBAL_WITH_CUTOFF=2 } enrichment_option_enum;
  enrichment_option_enum enrichment_option;
  
  scalar_type h, cont_gamma0, R; // mesh parameter
  
  
  std::string datafilename;
  
  int reference_test;
  
  bgeot::md_param PARAM;
  
  bool solve(plain_vector &U, plain_vector &LAMBDA, plain_vector &LAMBDA_T);
  void init(void);
  unilateral_contact_problem(void) : ls(mesh, 1, true), mls(mesh), mls_bound(mesh), mim(mls),
		                     mf_pre_u(mesh), mf_dir(mesh),
				     mfls_u(mls, mf_pre_u),
				     mf_sing_u(mesh),
				     mf_partition_of_unity(mesh),
				     mf_product(mf_partition_of_unity, mf_sing_u),
				     mf_u_sum(mesh), mf_contt(mesh), mfls_cont(mls, mf_contt), 				      
				     mf_cont_sum(mesh),
				     /*mf_pe(mesh),*/ mf_rhs(mesh)
  {}

};//end structure of unilateral contact problem          


std::string name_of_dof(getfem::pdof_description dof) {
  char s[200];
  sprintf(s, "UnknownDof[%p]", (void*)dof);
  for (dim_type d = 0; d < 4; ++d) {
    if (dof == getfem::lagrange_dof(d)) {
      sprintf(s, "Lagrange[%d]", d); goto found;
    }
    if (dof == getfem::normal_derivative_dof(d)) {
      sprintf(s, "D_n[%d]", d); goto found;
    }
    if (dof == getfem::global_dof(d)) {
      sprintf(s, "GlobalDof[%d]", d);
    }
    if (dof == getfem::mean_value_dof(d)) {
      sprintf(s, "MeanValue[%d]", d);
    }
    if (getfem::dof_xfem_index(dof) != 0) {
      sprintf(s, "Xfem[idx:%d]", int(dof_xfem_index(dof)));
    }
    
    for (dim_type r = 0; r < d; ++r) {
      if (dof == getfem::derivative_dof(d, r)) {
	sprintf(s, "D_%c[%d]", "xyzuvw"[r], d); goto found;
      }
      for (dim_type t = 0; t < d; ++t) {
	if (dof == getfem::second_derivative_dof(d, r, t)) {
	  sprintf(s, "D2%c%c[%d]", "xyzuvw"[r], "xyzuvw"[t], d); 
	  goto found;
	}
      }
    }
  }
 found:
  return s;
}

/****************************************************************************************************/
/*Initialisation of unilateral contact problem                                                      */
/****************************************************************************************************/
void  unilateral_contact_problem::init(void) {
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string FEM_TYPE_CONT  = PARAM.string_value("FEM_TYPE_cont","FEM name mult contact viriable");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
						       "Name of simplex integration method");
  std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION");
  
  enrichment_option = enrichment_option_enum(PARAM.int_value("ENRICHMENT_OPTION",
							     "Enrichment option"));
  cout << "MESH_TYPE="               << MESH_TYPE                 << "\n";
  cout << "FEM_TYPE="                << FEM_TYPE                  << "\n";
  cout << "FEM_TYPE_CONT="           << FEM_TYPE_CONT             << "\n";
  cout << "INTEGRATION="             << INTEGRATION               << "\n";
  cout << " SIMPLEX_INTEGRATION="    << SIMPLEX_INTEGRATION        << "\n";
  cout << "SINGULAR_INTEGRATION="    << SINGULAR_INTEGRATION      << "\n";

  
  dgr = int(PARAM.int_value("dgr", "Enrichement order ofu"));
  
  
  mu    = PARAM.real_value("Mu", "Lame coefficient mu"); 
  lambda =PARAM.real_value("Lambda", "Lame coefficient lambda");
  
  
  /* First step : build the mesh */
  strmesh=PARAM.int_value("strmesh", "Structured mesh or not");
  
  bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  std::vector<size_type> nsubdiv(N);
  if (strmesh) {
   std::fill(nsubdiv.begin(),nsubdiv.end(),
	      PARAM.int_value("NX", "Nomber of space steps "));
    getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			      PARAM.int_value("MESH_NOISED") != 0);
    base_small_vector tt(N); tt[1] = -0.5;
    mesh.translation(tt); 
    cout<<"Creting mesh done"<<endl;     
    
    cracktip.resize(2); // Cordinate of cracktip
    cracktip[0] = 0.5;
    cracktip[1] = 0.;
  }else{
    std::string MESH_FILE = PARAM.string_value("MESH_FILE","Mesh file");
    mesh.read_from_file(MESH_FILE);
    cout<<"Import non-conforming mesh"<<":"<<MESH_FILE <<endl;     
  }
  scalar_type refinement_radius;
  refinement_radius
    = PARAM.real_value("REFINEMENT_RADIUS", "Refinement Radius");
  size_type refinement_process;
  refinement_process
    = PARAM.int_value("REFINEMENT_PROCESS", "Refinement process");
  
  if (refinement_radius > 0) {
    for (size_type ref = 0; ref < refinement_process; ++ref){
      dal::bit_vector conv_to_refine;
      for(dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i){
	for(size_type j=0; j < 3; ++j)
	  if(gmm::vect_dist2(mesh.points_of_convex(i)[j],cracktip)
	     < refinement_radius )
	    conv_to_refine.add(i);
      }
      mesh.Bank_refine(conv_to_refine);
      
      refinement_radius = refinement_radius/3.;
      cout <<"refining process step " << ref << " ... refining "
	   << conv_to_refine.size() <<" convexes..." << endl;
    }
    cout << "refinement process completed." << endl ;
  }
  
  mesh.write_to_file("toto.mesh");
  
  h = mesh.maximal_convex_radius_estimate();
  // h = 2. * mesh.minimal_convex_radius_estimate();

  cout << "h = " << h << endl;
  
  cont_gamma0 = PARAM.real_value("CONTACT_GAMMA0",
				 "Stabilization parameter for "
				 "contact condition");
  R = PARAM.real_value( "R", "augmented parameter");	  
  
  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL");
  if (residual == 0.) residual = 1e-10;
  enr_area_radius = PARAM.real_value("RADIUS_ENR_AREA",
				     "radius of the enrichment area");
  
  reference_test = int(PARAM.int_value("REFERENCE_TEST", "Reference test")); 
  
  
  cutoff_func = PARAM.int_value("CUTOFF_FUNC", "cutoff function");
  
  cutoff_radius = PARAM.real_value("CUTOFF", "Cutoff");
  cutoff_radius1 = PARAM.real_value("CUTOFF1", "Cutoff1");
  cutoff_radius0 = PARAM.real_value("CUTOFF0", "Cutoff0");
  mf_u().set_qdim(dim_type(N));
  

  friction_coeff= PARAM.real_value("FRICTION_COEFF", "Friction_coeff");
  tresca_threshold= PARAM.real_value("TRESCA_THRESHOLD", "Tresca_threshold");
  ratio_size= PARAM.real_value("RATIO_SIZE", "ratio size beween mesh and coarse mesh");
  
  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pfem pf_mult_cont = 
    getfem::fem_descriptor(FEM_TYPE_CONT);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);
  getfem::pintegration_method simp_ppi = 
    getfem::int_method_descriptor(SIMPLEX_INTEGRATION);
  getfem::pintegration_method sing_ppi = (SINGULAR_INTEGRATION.size() ?
					  getfem::int_method_descriptor(SINGULAR_INTEGRATION) : 0);
  
  mim.set_integration_method(mesh.convex_index(), ppi);
  mls.add_level_set(ls);
  
  mim.set_simplex_im(simp_ppi, sing_ppi);
  mf_pre_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_dir.set_finite_element(mesh.convex_index(), pf_u);
  mf_dir.set_qdim(dim_type(N));
  mf_partition_of_unity.set_classical_finite_element(1);
  
  // Integration method on the boudary
  
  mls_bound.add_level_set(ls);
  int intbound = getfem::mesh_im_level_set::INTEGRATE_BOUNDARY;
  mimbound.init_with_mls(mls_bound, intbound, simp_ppi, sing_ppi);
  mimbound.set_integration_method(mesh.convex_index(), ppi);
  
  
  
  
  contact_only = (PARAM.int_value("CONTACT_ONLY"," contact_only or not.") != 0);
  Tresca_version=(PARAM.int_value("TRESCA_VERSION"," tresca version or not.") != 0);
  cout<< "Tresca Version="<< Tresca_version<< endl;
  stabilized_problem = PARAM.real_value("STABILIZED_PROBLEM"," stabilized_problem or not.");
  Test_of_coer= (PARAM.int_value("STABILIZED_PROBLEM"," stabilized_problem or not.") != 0);
  rangeP_P =
    (PARAM.int_value("rangeP_P"," The creteria of choise in range bases is performed by comparison to p  or not.") != 0);
  
  mf_contt.set_finite_element(mesh.convex_index(), pf_mult_cont);
  
  
  
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
    
    base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);
    
    if (un[0]  > 0.5) mesh.region(NEUMANN1_BOUNDARY_NUM).add(i.cv(), i.f());
    if (un[1]  > 0.5) mesh.region(NEUMANN2_BOUNDARY_NUM).add(i.cv(), i.f());
    if (un[0]  < -0.5) mesh.region(NEUMANN3_BOUNDARY_NUM).add(i.cv(), i.f());
    if (un[1]  < -0.5) mesh.region(NEUMANN4_BOUNDARY_NUM).add(i.cv(), i.f());
  }
  
}//end intialisation


/****************************************************************************************************/
/*Defining level set function                                                                       */
/****************************************************************************************************/

base_small_vector ls_function(const base_node P, int num = 0) {
  scalar_type x = P[0], y = P[1], x0=0, y0=0, phi, ll;
  base_small_vector res(2), cracktip(2), t(2), xx(2);
  cracktip[0]=0.5; cracktip[1]= 0.;
  switch (num) {
    
  case 0: {
    xx[0]= P[0]; xx[1]=P[1];
    t[0]= cracktip[0] - x0;
    t[1]= cracktip[1] - y0;
    t= t * (1./ gmm::vect_norm2(t));
    ll= gmm::sqr(t[0]* t[0] + t[1]*t[1]);
    phi= t[0]*y - t[1]*x + x0*cracktip[1] - y0*cracktip[0];
    res[0]= phi/ll;
    res[1]= gmm::vect_sp(xx-cracktip, t);
  } break;
    
  case 1: {
    res[0] = y;
    res[1] = -0.5 + x;
  } break;
  case 2: {
    res[0] = x - 0.25;
    res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
  } break;
  default: assert(0);
  }
  return res;
}//end ls_function



/****************************************************************************************************/
/*Inf-Sup condition                                                                                 */
/****************************************************************************************************/

struct matrix_G {
  
  const sparse_matrix &B;
  const sparse_matrix &S;
  mutable plain_vector W1, W2;
  
  gmm::SuperLU_factor<scalar_type> SLUF;

  matrix_G(const sparse_matrix &BB, const sparse_matrix &SS)
    : B(BB), S(SS), W1(gmm::mat_nrows(SS)), W2(gmm::mat_nrows(SS)) {
    SLUF.build_with(SS);
  }
  
};



template <typename vector1, typename vector2>
void mult(const matrix_G &G, const vector1 &X, vector2 &Y) {
  gmm::mult(gmm::transposed(G.B), X, G.W1);
  // gmm::iteration it(1E-6, 0);
  // gmm::cg(G.S, G.W2, G.W1,  gmm::identity_matrix(), it);
  G.SLUF.solve(G.W2, G.W1);
  gmm::mult(G.B, G.W2, Y);
}



template <typename vector1, typename vector2>
void mult(const matrix_G &G, const vector1 &X, const vector2 &b, vector2 &Y)
{ mult(G, X, Y); gmm::add(b, Y); }

scalar_type smallest_eigen_value(const sparse_matrix &B,
				 const sparse_matrix &M,
				 const sparse_matrix &S) {
  
  size_type n = gmm::mat_nrows(M);
  scalar_type lambda;
  plain_vector V(n), W(n), V2(n);
  gmm::fill_random(V2);
  matrix_G G(B, S);
  
  do {
    gmm::copy(V2, V);
    gmm::scale(V, 1./gmm::vect_norm2(V));
    gmm::mult(M, V, W);
    
    gmm::iteration it(1E-3, 0);
    gmm::cg(G, V2, W,  gmm::identity_matrix(), it);    
    lambda = gmm::vect_norm2(V2);
    
    //  compute the Rayleigh quotient
    //     mult(G, V2, W);
    //     scalar_type lambda2 = gmm::vect_sp(V2, W);
    //     gmm::mult(M, V2, W);
    //     lambda2 /= gmm::vect_sp(V2, W);
    //     cout << "lambda2 = " << sqrt(lambda2) << endl;
    
    cout << "lambda = " << sqrt(1./lambda) << endl;
    cout << "residu = " << gmm::vect_dist2(V2, gmm::scaled(V, lambda)) << endl;
    
  } while (gmm::vect_dist2(V2, gmm::scaled(V, lambda)) > 1E-3);
  
  return sqrt(1./lambda);

}



/************************************************************************************************/
/*     solv contact_problem                                                                     */
/************************************************************************************************/

bool  unilateral_contact_problem::solve(plain_vector &U, plain_vector &LAMBDA, plain_vector &LAMBDA_T ) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();
  ls.reinit();  
  
  for (size_type d = 0; d < ls.get_mesh_fem().nb_basic_dof(); ++d) {
    ls.values(0)[d] = ls_function(ls.get_mesh_fem().point_of_basic_dof(d), 0)[0];
    ls.values(1)[d] = ls_function(ls.get_mesh_fem().point_of_basic_dof(d), 0)[1];
  }
  ls.touch();
  
  mls.adapt();
  mls_bound.adapt();
  mim.adapt();
  mimbound.adapt();
  mfls_u.adapt();
  mfls_cont.adapt();
  
  bool load_global_fun =  0;
  
  
  cout << "Setting up the singular functions for the enrichment\n";
  cout << dgr << endl;
  
  std::vector<getfem::pglobal_function> vfunc(dgr*4);
  //std::vector<getfem::pglobal_function> vfunc_u(6);
  //std::vector<getfem::pglobal_function> vfunc_p(dgrp*2);
  if (!load_global_fun) {
    std::cout << "Using default singular functions\n";
    for (unsigned i = 0; i < vfunc.size(); ++i){
      /* use the singularity */
      
      getfem::abstract_xy_function *s = 
	new getfem::crack_singular_xy_function(i);
      
      if (enrichment_option != FIXED_ZONE ) {
	/* use the product of the singularity function
	   with a cutoff */
	getfem::abstract_xy_function *c = 
	  new getfem::cutoff_xy_function(int(cutoff_func),
					 cutoff_radius, 
					 cutoff_radius1,cutoff_radius0);
	s  = new getfem::product_of_xy_functions(*s, *c);
       }
      vfunc[i]=getfem::global_function_on_level_set(ls, *s);           
    }
  }
  
  mf_sing_u.set_functions(vfunc);
  
  switch (enrichment_option) {
    
  case FIXED_ZONE :
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
      if (enriched_dofs.card() < 3)
	GMM_WARNING0("There is " << enriched_dofs.card() <<
		     " enriched dofs for the crack tip");
      mf_product.set_enrichment(enriched_dofs);
      mf_u_sum.set_mesh_fems(mf_product, mfls_u);
      //      mf_product_cont.set_enrichment(enriched_dofs);
      mf_cont_sum.set_mesh_fems(mf_contt);
    }
    break;
    
    
  case GLOBAL_WITH_CUTOFF :
    {
      if(cutoff_func == 0)
	cout<<"Using exponential Cutoff..."<<endl;
      else
	cout<<"Using Polynomial Cutoff..."<<endl;
      mf_u_sum.set_mesh_fems(mf_sing_u, mfls_u);
      mf_cont_sum.set_mesh_fems(mf_contt);
    } break;
    
  case NO_ENRICHMENT:
    {
      mf_u_sum.set_mesh_fems(mfls_u);
      mf_cont_sum.set_mesh_fems(mf_contt);
    } break;
    
  }// end switch
  
  
  /*..................................................................................*/
  /* Range_basis call                                                                 */
  /*..................................................................................*/
  
  std::set<size_type> cols;
  cols.clear(); 
  if (rangeP_P){
    sparse_matrix BRBB(mf_cont().nb_dof(), mf_cont().nb_dof());
    asm_mass_matrix(BRBB, mimbound, mf_cont(), mf_cont());
    cout << "Selecting dofs for the multiplier PP" << endl;
    cout << "nb_dof_mult = " << mf_cont().nb_dof() << endl;
    gmm::range_basis(BRBB, cols);
    mf_cont().reduce_to_basic_dof(cols);
  }else{
    sparse_matrix BRBB(mf_pre_uu().nb_dof(), mf_cont().nb_dof());
    mf_pre_uu().set_qdim(1);
    asm_mass_matrix(BRBB, mimbound, mf_pre_uu(), mf_cont());
    cout << "Selecting dofs for the multiplier UP" << endl;
    cout << "nb_dof_mult = " << mf_cont().nb_dof() << endl;
    gmm::range_basis(BRBB, cols);
    mf_cont().reduce_to_basic_dof(cols);
    mf_pre_uu().set_qdim(dim_type(N));
}
  
  
  
  
  size_type nb_dof = mf_u().nb_dof();
  cout << "nb_dof = " << nb_dof << endl;
  size_type nb_dof_cont = mf_cont().nb_dof();
  cout << "nb_dof_cont.... = " << nb_dof_cont << endl;
  cout << "nb_basic_dof_cont.... = " <<mf_cont().nb_basic_dof() << endl;
  
  U.resize(nb_dof);
  LAMBDA.resize(nb_dof_cont);
  LAMBDA_T.resize(nb_dof_cont);
  
  // Find the dofs on the upper right and lower right corners
  cout << "Find the dofs on the upper right and lower right corners" << endl;
  
  scalar_type d1 = 1.0, d2 = 1.0;
  size_type icorner1 = size_type(-1), icorner2 = size_type(-1);
  base_node corner1 = base_node(1.0, -0.5);
  base_node corner2 = base_node(1.0, 0.5);
  GMM_ASSERT1(!(mf_u().is_reduced()), "To be adapted for reduced fems");
  for (size_type i = 0; i < mf_u().nb_basic_dof(); i+=N) {
    scalar_type dd1 = gmm::vect_dist2(mf_u().point_of_basic_dof(i), corner1);
    if (dd1 < d1) { icorner1 = i; d1 = dd1; }
    scalar_type dd2 = gmm::vect_dist2(mf_u().point_of_basic_dof(i), corner2);
    if (dd2 < d2) { icorner2 = i; d2 = dd2; }
  }//end for
  
  GMM_ASSERT1(((d1 < 1E-8) && (d2 < 1E-8)),
	      "Upper right or lower right corners not found d1 = "
	      << d1 << " d2 = " << d2);
  
  // Assembling mixed term for contact problem
  
  getfem::CONTACT_B_MATRIX BN(nb_dof_cont, nb_dof);
  getfem::CONTACT_B_MATRIX BT(nb_dof_cont, nb_dof);
  asm_mass_matrix_mixed_term(BN, mimbound, mf_u(), mf_cont(), ls);
  if (!contact_only) {
    asm_mass_matrix_tang_mixed_term(BT, mimbound, mf_u(), mf_cont(), ls);
  }
  
  getfem::CONTACT_B_MATRIX CA(nb_dof_cont, nb_dof);
  // Assembling stabilized mixed term for contact problem
  //cont_gamma0=cont_gamma/h;
  // if (!contact_only){
 getfem::CONTACT_B_MATRIX CAT(nb_dof_cont, nb_dof);
 // }
 cout<<"Cont_gamma0="<< cont_gamma0 <<endl;
  if (stabilized_problem==1) {
    cout<<"stabilized_problem="<<stabilized_problem<<endl;
        cout<< "Assembling stabilized mixed term for contact problem"<<endl;
    asm_stabilization_mixed_term(CA, mimbound, mf_u(), mf_cont(), ls, lambda, mu);
    gmm::scale(CA, -cont_gamma0 * h);
    if (!contact_only){
      cout<< "Assembling stabilized tangent mixed term for friction problem"<<endl;
      asm_stabilization_tang_mixed_term(CAT, mimbound, mf_u(), mf_cont(), ls, lambda, mu);
      gmm::scale(CAT, -cont_gamma0 * h);
    }
  }
  
  // Assembling symetrique term for stabilized contact problem
  
  sparse_matrix KA(nb_dof, nb_dof);
  // if (!contact_only){
  sparse_matrix KAT(nb_dof, nb_dof);
    //}
  if (stabilized_problem==1) {
    cout<<"Assembling symetrique term for stabilized contact problem"<<endl;
    asm_stabilization_symm_term(KA, mimbound, mf_u(), ls, lambda, mu);
    gmm::scale(KA, -cont_gamma0 * h);
    if (!contact_only){
      cout<<"Assembling symetrique tangent term for stabilized friction problem"<<endl;
      asm_stabilization_symm_tang_term(KAT, mimbound, mf_u(), ls, lambda, mu);
      gmm::scale(KAT, -cont_gamma0 * h);
    }
  }
  
  //Assembling mass term for stabilized problem
  
  getfem::CONTACT_B_MATRIX MA(nb_dof_cont, nb_dof_cont);
  getfem::CONTACT_B_MATRIX MAT(nb_dof_cont, nb_dof_cont);
  getfem::CONTACT_B_MATRIX M1(nb_dof_cont, nb_dof_cont);
  getfem::CONTACT_B_MATRIX M1T(nb_dof_cont, nb_dof_cont);
  if (stabilized_problem!=0) {
    cout<<"Assembling mass term for stabilized problem"<<endl;
    getfem::asm_mass_matrix(MA, mimbound, mf_cont());
    if (!contact_only){
      cout<<"Assembling mass tangent term for stabilized problem"<<endl;
      getfem::asm_mass_matrix(MAT, mimbound, mf_cont());
     
    }
    if (stabilized_problem==1) {
      gmm::scale(MA, cont_gamma0 * h);
      //cout << "MA = " << MA << endl;
      if (!contact_only){ 
	gmm::scale(MAT, cont_gamma0 * h);
      }
    }
    if (stabilized_problem==2) {
      cout<<"Patch Stabilization"<<endl;
      //gmm::scale(MA, cont_gamma0);
      asm_stabilization_patch_term(M1, mesh, mimbound, mf_cont(), ratio_size, 2.*h);
      gmm::scale(M1, -cont_gamma0);
      gmm::scale(MA, cont_gamma0);
      gmm::add(M1, MA);
      // cout << "MA = " << MA << endl;
      if (!contact_only){ 
	//gmm::scale(MAT, cont_gamma0);
	asm_stabilization_patch_term(M1T, mesh, mimbound, mf_cont(), ratio_size, 2.*h);
	gmm::scale(M1T, -cont_gamma0);
	gmm::scale(MAT, cont_gamma0);
	gmm::add(M1T, MAT);
      }
    }
  }
  
  
  getfem::model model;
  
  // Main unknown of the problem.
  model.add_fem_variable("u", mf_u());
  model.add_fem_variable("Lambda", mf_cont()); // Adding normal contact variable
  if (!contact_only){
    model.add_fem_variable("Lambda_t", mf_cont()); // Adding friction contact variable
  }
  
  // Linearized elasticity brick.
  model.add_initialized_scalar_data("lambda", lambda);
  model.add_initialized_scalar_data("mu", mu);
  if (!contact_only){
    cout<<"Friction coeff="<< friction_coeff<<endl;
    model.add_initialized_scalar_data("Friction_coeff", friction_coeff);
    if (Tresca_version){
      cout<<"tresca_threshold="<<tresca_threshold <<endl;
      model.add_initialized_scalar_data("Tresca_threshold", tresca_threshold);
    }
    
    
  }
  model.add_initialized_scalar_data("augmentation_parameter", R/h);
  //model.add_initialized_scalar_data("augmentation_parameter", R);
  getfem::add_isotropic_linearized_elasticity_brick
    (model, mim, "u", "lambda", "mu");
  //  model.add_initialized_scalar_data("augmentation_parameter", R/h);
  // model.add_initialized_scalar_data
  //("augmentation_parameter", mu * (3*lambda + 2*mu) / (h*(lambda + mu)) );  // r ~= Young modulus
  
  
    
  if (stabilized_problem!=0) {
    if (stabilized_problem==1){
      // cout << "KA = " << KA << endl;
      cout << "gamma = " << cont_gamma0 * h << endl;
      getfem::add_explicit_matrix(model, "u", "u", KA);
      // Defining the contact condition.
      gmm::add(CA, BN);
      if (!contact_only){
	getfem::add_explicit_matrix(model, "u", "u", KAT);
	gmm::add(CAT, BT); 
      }

    }

    if (contact_only){
      getfem::add_Hughes_stab_basic_contact_brick(model, "u", "Lambda",
						  "augmentation_parameter", BN, MA);
    }else{
      if (Tresca_version){
	getfem::add_Hughes_stab_basic_contact_brick
	  (model, "u", "Lambda", "Lambda_t", "augmentation_parameter",
	   BN, BT, MA, MAT, "Friction_coeff","","",1, Tresca_version, "Tresca_threshold");
      }else{
	getfem::add_Hughes_stab_basic_contact_brick
	  (model, "u", "Lambda", "Lambda_t", "augmentation_parameter",
	   BN, BT, MA, MAT, "Friction_coeff");
      }
    }

  }else{
    if (contact_only){		
      getfem::add_basic_contact_brick(model, "u", "Lambda",
				      "augmentation_parameter", BN);
    }else{
      if (Tresca_version){
     	getfem::add_basic_contact_brick
	  (model, "u", "Lambda", "Lambda_t",
	   "augmentation_parameter", BN, BT, "Friction_coeff","","",1,Tresca_version,"Tresca_threshold");
      }else{
	getfem::add_basic_contact_brick
	  (model, "u", "Lambda", "Lambda_t",
	   "augmentation_parameter", BN, BT, "Friction_coeff");
      }
      
    }
    
  }
  

 // if (stabilized_problem) {
//     // cout << "KA = " << KA << endl;
//     cout << "gamma = " << cont_gamma0 * h << endl;
//     getfem::add_explicit_matrix(model, "u", "u", KA);
//     // Defining the contact condition.
//     gmm::add(CA, BN); 
//     if (contact_only){
//       getfem::add_Hughes_stab_basic_contact_brick(model, "u", "Lambda",
// 						  "augmentation_parameter",
// 						  BN, MA);
//     }else { 
//       getfem::add_explicit_matrix(model, "u", "u", KAT);
//       // Defining the contact condition.
//       gmm::add(CAT, BT); 
//       getfem::add_Hughes_stab_basic_contact_brick
// 	(model, "u", "Lambda", "Lambda_t", "augmentation_parameter",
// 	 BN, BT, MA, MAT, "friction_coeff");
//     }
//   }  else {		
//     getfem::add_basic_contact_brick(model, "u", "Lambda",
// 				    "augmentation_parameter", BN);
    
//     if (!contact_only){
//       getfem::add_basic_contact_brick
// 	(model, "u", "Lambda", "Lambda_t",
// 	 "augmentation_parameter", BN, BT, "friction_coeff","","",0);
      
//     }
//   }
  
  
  // Defining the Neumann condition right hand side.
  
  std::vector<scalar_type> F(nb_dof_rhs * N * N);
  std::vector<scalar_type> FF(nb_dof_rhs * N * N);
  
  // Neumann condition brick.
  
  for(size_type i = 0; i < F.size(); i=i+N*N) {
    base_node pt = mf_rhs.point_of_basic_dof(i / (N*N));
    for(size_type j = 0; j < N; ++j){ F[i+j+j*N] =sin(2*M_PI*pt[1])*0.4;}
    
    // scalar_type coeff = pt[1] > 0. ? -1 : 1;
    // for(size_type j = 0; j < N; ++j){ F[i+j+j*N] =  coeff*0.1; FF[i+j+j*N] = (pt[0]-0.5)*2.0;}
    
  }
  
  cout <<"Applied  Neumann condition"<< endl;
  model.add_initialized_fem_data("NeumannData", mf_rhs, F);
  model.add_initialized_fem_data("NeumannData1", mf_rhs, FF);
  //getfem::add_normal_source_term_brick
  //      (model, mim, "u", "NeumannData1", NEUMANN4_BOUNDARY_NUM);
  getfem::add_normal_source_term_brick
    (model, mim, "u", "NeumannData", NEUMANN3_BOUNDARY_NUM);
  //    getfem::add_normal_source_term_brick
  //      (model, mim, "u", "NeumannData1", NEUMANN2_BOUNDARY_NUM);
  getfem::add_normal_source_term_brick
    (model, mim, "u", "NeumannData", NEUMANN1_BOUNDARY_NUM);
  
  
  std::vector<scalar_type> FFF(mf_rhs.nb_dof()*N, 0.0);
  for(size_type i = 0; i < FFF.size(); i+=N){ 
    base_node pt = mf_rhs.point_of_basic_dof(i/N);
    FFF[i+1]=3.5*cos(pt[0]*2.*M_PI) * pt[1]*pt[0]*(1-pt[0]) ;
  }
  model.add_initialized_fem_data("VolumicForce", mf_rhs, FFF);
  getfem::add_source_term_brick(model, mim, "u",  "VolumicForce");
  
  
  //  std::vector<scalar_type> FFF(mf_rhs.nb_dof()*N);
  //   for(size_type i = 0; i <FFF.size(); i++){ FFF[i]=0;}
  //   model.add_initialized_fem_data("DirichletData", mf_rhs, FFF);
  //   getfem::add_Dirichlet_condition_with_multipliers
  //     (model, mim, "u", mf_dir,  NEUMANN4_BOUNDARY_NUM, "DirichletData");
  
  cout<<"kill rigid motion"<<endl;
  GMM_ASSERT1(N==2, "To be corrected for 2D computation");
  sparse_matrix BB(3, mf_u().nb_dof());
  BB(0, icorner1) = 1.0;
  BB(1, icorner1+1) = 1.0;
  BB(2, icorner2) = 1.0;
  size_type size(3);
  std::vector<scalar_type> LRH(size);
  model.add_fixed_size_variable("dir", size);
  getfem::add_constraint_with_multipliers(model, "u", "dir", BB, LRH);

  //verified coercivit condition
  
  

  
  
  
  
  //classical inf-sup condition
  
  //   if (PARAM.int_value("INF_SUP_COMP")) {
  
  //      cout << "Sparse matrices computation for the test of inf-sup condition"
  // 	 << endl;
  
  //   //  sparse_matrix MUU(nb_dof, nb_dof);
  //     sparse_matrix MPP(nb_dof_cont, nb_dof_cont);
  //     sparse_matrix BBB(nb_dof_cont, nb_dof);
  //     sparse_matrix MTILD(nb_dof_cont, nb_dof_cont);
  //     plain_vector VVV(nb_dof),VVV1(nb_dof),VVV0(nb_dof), dd1(nb_dof_cont), dd2(nb_dof_cont), YYY(nb_dof_cont), ZZZ(nb_dof_cont);
  //     scalar_type lllambda;
//     asm_mass_matrix_for_inf_sup(MTILD, mimbound, mf_cont(), ls); // int_gamma_c 1/h psi_i psi_j
//     asm_mass_matrix_mixed_term(BBB, mimbound, mf_u(), mf_cont(), ls);
//    // getfem::asm_mass_matrix(MUU, mimbound, mf_u());
//     getfem::asm_mass_matrix(MPP, mimbound, mf_cont());
//     gmm::fill_random(VVV);
//   //  gmm::scaled(VVV,gmm::vect_norm2(VVV));//x_n/norm(x_n)
// do{
//     gmm::mult(BBB, VVV, dd1);
  
//    gmm::lu_solve(MPP, YYY,dd1);
  
//    gmm::mult(MTILD, YYY, dd2);
//    gmm::lu_solve(gmm::transposed(MPP), ZZZ,dd2);

//    gmm::mult(gmm::transposed(BBB), ZZZ, VVV1);
//    lllambda = gmm::vect_sp(VVV1,VVV);//gmm::vect_norm2(VVV);
//    gmm::copy(VVV,VVV0);
//    gmm::copy(VVV1,VVV);
//    cout << "lambda = " << lllambda << endl;
// }while (gmm::vect_dist2(VVV1, gmm::scaled(VVV0, lllambda )) > 1E-3);

// //cout << "lambda = " << sqrt(1./lllambda) << endl;
// cout << "residu = " << gmm::vect_dist2(VVV1, gmm::scaled(VVV, lllambda)) << endl;
// //     gmm::copy(gmm::lu_inverse(Bis), Biss);
// //     gmm::copy(BN, Mis);
// //     gmm::mult(Biss, Mis, BissMis);
// //     gmm::mult(Miss, BissMis, BBissMis);
// //     gmm::mult(gmm::transposed(BissMis), BBissMis, Siss);


// //     cout << "Inf-sup condition test" << endl;
// //     scalar_type lllambda = smallest_eigen_value(Siss, Siss, Sis);
// //     cout << "The inf-sup test gives " << lllambda << endl;

//   }


  
  // Generic solve.

  gmm::iteration iter(residual, 1, 2000);
  cout << "Solving..." << endl;
  iter.init();
  //getfem::standard_solve(model, iter);

  getfem::simplest_newton_line_search lse; // (size_t(-1), 6.0/5.0, 0.2, 3.0/5.0);
  //getfem::standard_solve(model, iter,
  //	getfem::default_linear_solver<getfem::model_real_sparse_matrix,
  //			 getfem::model_real_plain_vector>(model) , lse);

  getfem::standard_solve(model, iter, getfem::rselect_linear_solver(model,"mumps"), lse);

  gmm::Harwell_Boeing_save("toto.hb", model.real_tangent_matrix());


  //Test of coercivity
  if (PARAM.int_value("Test_of_coer")) {
    cout << "Compute the smallest eigenvalue" << endl;
    size_type  nb_up;
    gmm::sub_interval interval; 
    nb_up= gmm::nnz(model.real_variable("u"));
    interval = model.interval_of_variable("u");
    cout<< "nb_up=" << nb_up << endl;
    cout <<"nbdof_in_u="<<  gmm::nnz(model.real_variable("u"))<<endl;
    cout <<"nbdof_in_p="<<  gmm::nnz(model.real_variable("Lambda"))<<endl;
    sparse_matrix KKK(gmm::nnz(model.real_variable("u")), gmm::nnz(model.real_variable("u"))), 
		       KKKK(gmm::nnz(model.real_variable("u")), gmm::nnz(model.real_variable("u")));


    gmm::copy(gmm::sub_matrix(model.real_tangent_matrix(),interval, interval) ,KKK); 
    plain_vector VVV(gmm::nnz(model.real_variable("u"))), VVV1(gmm::nnz(model.real_variable("u"))) ;
    gmm::fill_random(VVV);   
    scalar_type lllambda(0), lllambda_aux(1), lllambda_max(0);
    do{
      lllambda_aux= lllambda;
      // gmm::copy(VVV,VVV1);
      gmm::copy(gmm::scaled(VVV, (1./gmm::vect_norm2(VVV))),VVV1);
      // cout<<"norm="<<gmm::vect_norm2(VVV1)<<endl;
      gmm::mult(KKK,VVV1,VVV);
      lllambda=gmm::vect_sp(VVV,VVV1); 
      // cout<<"lllambda="<< lllambda<<endl;
    }while (gmm::abs(lllambda-lllambda_aux ) > 1E-12);
    lllambda_max=lllambda;
    cout<<"lllambda_max="<< lllambda_max<<endl;


    gmm::copy(gmm::identity_matrix(),KKKK);
    gmm::scale(KKKK,-2.*lllambda_max);
    gmm::add(KKKK,KKK);
    //  cout<<KKKK<<endl;



    gmm::fill_random(VVV);   
    lllambda=0;
    lllambda_aux=1;
    do{
      lllambda_aux= lllambda;
      gmm::copy(gmm::scaled(VVV,(1./gmm::vect_norm2(VVV))),VVV1);
      gmm::mult(KKK,VVV1,VVV);
      lllambda=gmm::vect_sp(VVV,VVV1); 
      // cout<<"lllambda="<< lllambda<<endl;
    }while (gmm::abs(lllambda-lllambda_aux ) > 1E-10);
    cout<<"lllambda_min="<<  lllambda+2.*lllambda_max <<endl;
      }
  




  gmm::resize(U, mf_u().nb_dof());
  gmm::copy(model.real_variable("u"), U);
  gmm::resize(LAMBDA, mf_cont().nb_dof());
  gmm::copy(model.real_variable("Lambda"), LAMBDA);
  
  if (!contact_only) {
    gmm::resize(LAMBDA_T, mf_cont().nb_dof());
    gmm::copy(model.real_variable("Lambda_t"), LAMBDA_T);


  }
  
  return (iter.converged());
  
  
}//end solv






/* ************************************************************************************/
/* Main program                                                                       */
/**************************************************************************************/

int main(int argc, char *argv[]) {
  
  
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.
  
  //getfem::getfem_mesh_level_set_noisy();
  
  
  unilateral_contact_problem p;
  p.PARAM.read_command_line(argc, argv);
  p.init();
  p.mesh.write_to_file(p.datafilename + ".mesh");
  
  
  plain_vector U, Lambda, Lambda_t;
  // if (!p.contact_only) {
   //}
  // if (p.contact_only) {
  //  if (!p.solve(U, Lambda)) GMM_ASSERT1(false,"Solve has failed");
  //}else{
  if (!p.solve(U, Lambda, Lambda_t)) GMM_ASSERT1(false,"Friction Solve has failed")
      //   }
  bgeot::md_param PARAM1;
  PARAM1.read_command_line(argc, argv);
  std::string INTEGRATION1 = PARAM1.string_value("INTEGRATION",
						 "Name of integration method");
  std::string SIMPLEX_INTEGRATION1 = PARAM1.string_value("SIMPLEX_INTEGRATION",
							 "Name of simplex integration method");
  std::string SINGULAR_INTEGRATION1 = PARAM1.string_value("SINGULAR_INTEGRATION","Name of singular integration method");
  
  getfem::mesh mcut;
  p.mls.global_cut_mesh(mcut);
  unsigned Q = p.mf_u().get_qdim();
  
  cout << p.mf_u().nb_basic_dof() << endl;
  cout << p.mf_cont().nb_basic_dof() << endl;
   cout << p.mf_cont().nb_dof() << endl;
  //construction of mesh refined
  
  // getfem::stored_mesh_slice slr;
  getfem::mesh mcut_refined;
  
  unsigned NX = unsigned(p.PARAM.int_value("NX")), nn;
  if (NX < 6) nn = 24;
  else if (NX < 12) nn = 6;
  else if (NX < 30) nn = 3;
  else nn = 3;
  
  /* choose an adequate slice refinement based on the distance
     to the crack tip */
  std::vector<bgeot::short_type> nrefine(mcut.convex_index().last_true()+1);
  for (dal::bv_visitor cv(mcut.convex_index()); !cv.finished(); ++cv) {
    scalar_type dmin=0, d;
    base_node Pmin,Pp;
    for (unsigned i=0; i < mcut.nb_points_of_convex(cv); ++i) {
      Pp = mcut.points_of_convex(cv)[i];
      d = gmm::vect_norm2(ls_function(Pp));
      if (d < dmin || i == 0) { dmin = d; Pmin = Pp; }
    }
    
    if (dmin < 1e-5)
      nrefine[cv] = short_type(nn*8);
    else if (dmin < .1) 
      nrefine[cv] = short_type(nn*2);
    else nrefine[cv] = short_type(nn);
    // if (dmin < .01)
    //  cout << "cv: "<< cv << ", dmin = " << dmin << "Pmin=" << Pmin
    //       << " " << nrefine[cv] << "\n";
  }
  
  {
    getfem::mesh_slicer slicer(mcut); 
    getfem::slicer_build_mesh bmesh(mcut_refined);
    slicer.push_back_action(bmesh);
    slicer.exec(nrefine, getfem::mesh_region::all_convexes());
  }
  
  if (1) { 
    
    getfem::mesh_fem mf(mcut, dim_type(Q));
    mf.set_classical_discontinuous_finite_element(2, 1E-5);
    plain_vector VV(mf.nb_dof());
    getfem::interpolation(p.mf_u(), mf, U, VV);
    
    getfem::mesh_fem mf_cont(mcut, 1);
    mf_cont.set_classical_discontinuous_finite_element(2, 1E-5);
    plain_vector PPP(mf_cont.nb_dof());
    getfem::interpolation(p.mf_cont(), mf_cont, Lambda, PPP);
    
    //if (!p.contact_only) {
     
    //}
    getfem::level_set ls(mcut, 1, true);
    for (size_type d = 0; d < ls.get_mesh_fem().nb_basic_dof(); ++d) {
      ls.values(0)[d] =ls_function(ls.get_mesh_fem().point_of_basic_dof(d), 0)[0];
      ls.values(1)[d] =ls_function(ls.get_mesh_fem().point_of_basic_dof(d), 0)[1];
    }
    ls.touch();
    
    const getfem::mesh_fem &lsmf = ls.get_mesh_fem();
    lsmf.write_to_file(p.datafilename + "_ls.mf", true);
    gmm::vecsave(p.datafilename + "_ls.U", ls.values());
    
    
    
    
    cout << "saving the slice solution for matlab\n";
    getfem::stored_mesh_slice sl, sl0,sll;
    
    
    getfem::mesh_slicer slicer1(mf.linked_mesh());
    getfem::slicer_build_stored_mesh_slice sbuild(sl);
    //   getfem::mesh_slice_cv_dof_data<plain_vector> mfU(mf,VV);
    //   getfem::slicer_isovalues iso(mfU, 0.0, 0);
    getfem::slicer_build_stored_mesh_slice sbuild0(sl0);
    
    slicer1.push_back_action(sbuild);  // full slice in sl
    //  slicer1.push_back_action(iso);     // extract isosurface 0
    slicer1.push_back_action(sbuild0); // store it into sl0
    slicer1.exec(nrefine, mf.convex_index());
    
    getfem::mesh_slicer slicer2(mf.linked_mesh());
    getfem::mesh_slice_cv_dof_data<plain_vector> 
      mfL(ls.get_mesh_fem(), ls.values());
    getfem::slicer_isovalues iso2(mfL, 0.0, 0);
    getfem::slicer_build_stored_mesh_slice sbuildl(sll);
    slicer2.push_back_action(iso2);     // extract isosurface 0
    slicer2.push_back_action(sbuildl); // store it into sl0
    slicer2.exec(nrefine, mf.convex_index());

    if (p.reference_test) {
    sl.write_to_file(p.datafilename + "_ref.sl", true);
    sl0.write_to_file(p.datafilename + "_ref.sl0", true);
    sll.write_to_file(p.datafilename + "_ref.sll", true);
    plain_vector LL(sll.nb_points()); 
    //  gmm::scale(PP, 0.005);
    sll.interpolate(mf_cont, PPP, LL);
    gmm::vecsave(p.datafilename + "_ref.slL", LL);

    if (!p.contact_only) {
      getfem::mesh_fem mf_tang_cont(mcut, 1);
      mf_tang_cont.set_classical_discontinuous_finite_element(2, 1E-5);
      plain_vector PPPT(mf_tang_cont.nb_dof());
      getfem::interpolation(p.mf_cont(), mf_tang_cont, Lambda_t, PPPT);


      plain_vector LLT(sll.nb_points()); 
      sll.interpolate(mf_tang_cont, PPPT, LLT);
      gmm::vecsave(p.datafilename + "_ref.slLT", LLT);
    }
    //    plain_vector UU(sl.nb_points()), LL(sll.nb_points()); 
    //    sl.interpolate(mf, V, UU);
    //    gmm::vecsave("xfem_stab_unilat_contact.slU", UU);
    }else{
      sl.write_to_file(p.datafilename + ".sl", true);
      sl0.write_to_file(p.datafilename + ".sl0", true);
      sll.write_to_file(p.datafilename + ".sll", true);
      plain_vector LL(sll.nb_points()); 
           //  gmm::scale(PP, 0.005);
      sll.interpolate(mf_cont, PPP, LL);
      gmm::vecsave(p.datafilename + ".slL", LL);

      if (!p.contact_only) {
	getfem::mesh_fem mf_tang_cont(mcut, 1);
	mf_tang_cont.set_classical_discontinuous_finite_element(2, 1E-5);
	plain_vector PPPT(mf_tang_cont.nb_dof());
	getfem::interpolation(p.mf_cont(), mf_tang_cont, Lambda_t, PPPT);
	
	plain_vector LLT(sll.nb_points()); 
	sll.interpolate(mf_tang_cont, PPPT, LLT);
	gmm::vecsave(p.datafilename + ".slLT", LLT);
      }
      //    plain_vector UU(sl.nb_points()), LL(sll.nb_points()); 
      //    sl.interpolate(mf, V, UU);
      //    gmm::vecsave("xfem_stab_unilat_contact.slU", UU);
    }
  } 	  
  

  
  
  //save desplacement
  if (p.reference_test) {
    cout << "Saving the reference desplacement" << endl;
    getfem::mesh_fem mf(mcut_refined, dim_type(Q));
    mf.set_classical_discontinuous_finite_element(2, 1E-5);
    plain_vector V(mf.nb_dof());
    getfem::interpolation(p.mf_u(), mf, U, V);
    
    mf.write_to_file(p.datafilename + ".meshfem_refined", true);
    gmm::vecsave(p.datafilename + ".U_refined", V);
    
    //Save normal contact multiplier
    cout << "Saving the reference contact multiplier" << endl;
    
    getfem::mesh_fem mf_cont(mcut_refined, 1);
    mf_cont.set_classical_discontinuous_finite_element(2, 1E-5);
    plain_vector PP(mf_cont.nb_dof());
    
    getfem::interpolation(p.mf_cont(), mf_cont, Lambda, PP);
    
    mf_cont.write_to_file(p.datafilename + ".cont_meshfem_refined", true);
    gmm::vecsave(p.datafilename + ".cont_refined", PP);
    p.mf_cont().write_to_file(p.datafilename + ".contt_meshfem", true);
    gmm::vecsave(p.datafilename + ".contt", Lambda);
    
    //Save tangent contact multiplier
    cout << "Saving the reference tangant contact multiplier" << endl;
    
    getfem::mesh_fem mf_tang_cont(mcut_refined, 1);
    mf_tang_cont.set_classical_discontinuous_finite_element(2, 1E-5);
    plain_vector PPT(mf_tang_cont.nb_dof());
    
    getfem::interpolation(p.mf_cont(), mf_tang_cont, Lambda_t, PPT);
    
    mf_tang_cont.write_to_file(p.datafilename + ".tang_cont_meshfem_refined", true);
    gmm::vecsave(p.datafilename + ".tang_cont_refined", PPT);
    p.mf_cont().write_to_file(p.datafilename + ".tang_contt_meshfem", true);
    gmm::vecsave(p.datafilename + ".tang_contt", Lambda_t);
    
  }else{
    cout << "Saving the solution" << endl;
    cout << "Saving desplacement" << endl;
    getfem::mesh_fem mf(mcut, dim_type(Q));
    mf.set_classical_discontinuous_finite_element(2, 1E-5);
    plain_vector V(mf.nb_dof());
    getfem::interpolation(p.mf_u(), mf, U, V);
    
    mf.write_to_file(p.datafilename + ".meshfem", true);
    gmm::vecsave(p.datafilename + ".U", V);
    
    //Save normal contact multiplier
    cout << "Saving normal contact multiplier" << endl;
    
    getfem::mesh_fem mf_cont(mcut, 1);
    mf_cont.set_classical_discontinuous_finite_element(2, 1E-5);
    plain_vector PP(mf_cont.nb_dof());
    
    getfem::interpolation(p.mf_cont(), mf_cont, Lambda, PP);
    
    mf_cont.write_to_file(p.datafilename + ".cont_meshfem", true);
    gmm::vecsave(p.datafilename + ".cont", PP);
    p.mf_cont().write_to_file(p.datafilename + ".contt_meshfem", true);
    gmm::vecsave(p.datafilename + ".contt", Lambda); 
    
    //Save tangent contact multiplier
    cout<<"Saving tangent contact multiplier"<< endl;
    
    getfem::mesh_fem mf_tang_cont(mcut, 1);
    mf_tang_cont.set_classical_discontinuous_finite_element(2, 1E-5);
    plain_vector PPT(mf_tang_cont.nb_dof());
    
    getfem::interpolation(p.mf_cont(), mf_tang_cont, Lambda_t, PPT);
    
    mf_tang_cont.write_to_file(p.datafilename + ".tang_cont_meshfem", true);
    gmm::vecsave(p.datafilename + ".tang_cont", PPT);
    p.mf_cont().write_to_file(p.datafilename + ".tang_contt_meshfem", true);
    gmm::vecsave(p.datafilename + ".tang_contt", Lambda_t); 
  }
  
  


 cout << "plot Von mises" << endl;


  getfem::mesh_im mim_refined(mcut); 
  mim_refined.set_integration_method(getfem::int_method_descriptor
				     ("IM_TRIANGLE(6)"));

  mcut.write_to_file(p.datafilename + ".meshvm");

  getfem::mesh_fem mf_refined(mcut,dim_type(Q));
  mf_refined.set_classical_discontinuous_finite_element(4, 0.001);
  plain_vector W(mf_refined.nb_dof());
  
  getfem::interpolation(p.mf_u(), mf_refined, U, W);

  mf_refined.write_to_file(p.datafilename + ".meshfemuvm", true);
  gmm::vecsave(p.datafilename + ".Uvm", W);
 
  getfem::mesh_fem mf_vm(mcut,  1);
  mf_vm.set_classical_discontinuous_finite_element(4, 0.001);
  plain_vector Vm(mf_vm.nb_dof());

  
  cout << "compute Von_mises" << endl;
  getfem::interpolation_von_mises(mf_refined, mf_vm, W, Vm);

  gmm::vecsave(p.datafilename + ".VM",Vm);
  mf_vm.write_to_file(p.datafilename + ".meshfemvm", true);  






  //Compute error 
  
  if(p.PARAM.int_value("ERROR_TO_REF_SOL") == 1){
    cout << "Computing error with respect to a reference solution..." << endl;
    std::string REFERENCE_MF = p.datafilename + ".meshfem_refined";
    std::string REFERENCE_U = p.datafilename + ".U_refined";
    std::string REFERENCE_MF_cont = p.datafilename + ".cont_meshfem_refined";
    std::string REFERENCE_cont = p.datafilename + ".cont_refined";
    std::string REFERENCE_MF_tang_cont = p.datafilename + ".tang_cont_meshfem_refined";
    std::string REFERENCE_tang_cont = p.datafilename + ".tang_cont_refined";
    
    
    cout << "Load reference displacement from "
	 << REFERENCE_MF << " and " << REFERENCE_U << "\n";
    getfem::mesh ref_m; 
    ref_m.read_from_file(REFERENCE_MF);
    getfem::mesh_fem ref_mf(ref_m); 
    ref_mf.read_from_file(REFERENCE_MF);
    plain_vector ref_U(ref_mf.nb_dof());
    gmm::vecload(REFERENCE_U, ref_U);
    
    cout << "Load reference pressure from "
	 << REFERENCE_MF_cont << " and " << REFERENCE_cont << "\n";
    getfem::mesh ref_m_cont; 
    ref_m_cont.read_from_file(REFERENCE_MF_cont);
    getfem::mesh_fem ref_mf_cont(ref_m_cont); 
    ref_mf_cont.read_from_file(REFERENCE_MF_cont);
    plain_vector ref_cont(ref_mf_cont.nb_dof());
    gmm::vecload(REFERENCE_cont, ref_cont);

   
    
   

    // if (!p.contact_only) {
    
    cout << "Load reference tangent pressure from "
	 << REFERENCE_MF_tang_cont << " and " << REFERENCE_tang_cont << "\n";
    getfem::mesh ref_m_tang_cont; 
    ref_m_tang_cont.read_from_file(REFERENCE_MF_tang_cont);
    getfem::mesh_fem ref_mf_tang_cont(ref_m_tang_cont);
    ref_mf_tang_cont.read_from_file(REFERENCE_MF_tang_cont);
    
    plain_vector ref_tang_cont(ref_mf_tang_cont.nb_dof());
    gmm::vecload(REFERENCE_tang_cont, ref_tang_cont);
    // }
    
    getfem::level_set ls(ref_m_cont, 1, true);
    
    
    
    for (size_type d = 0; d < ls.get_mesh_fem().nb_basic_dof(); ++d) {
      ls.values(0)[d] =ls_function(ls.get_mesh_fem().point_of_basic_dof(d), 0)[0];
      ls.values(1)[d] =ls_function(ls.get_mesh_fem().point_of_basic_dof(d), 0)[1];
    }
    ls.touch();
    getfem::mesh_level_set mls(ref_m);
    getfem::mesh_level_set mls_cont(ref_m_cont);
    getfem::mesh_level_set mls_tang_cont(ref_m_tang_cont);
    mls.add_level_set(ls);
    mls_cont.add_level_set(ls);
    mls_tang_cont.add_level_set(ls);
    mls.adapt();
    mls_cont.adapt();
    mls_tang_cont.adapt();

    //int intbound = getfem::mesh_im_level_set::INTEGRATE_BOUNDARY;
    //mimbound.init_with_mls(mls, intbound, simp_ppi, sing_ppi);
    //mimbound.set_integration_method(mesh.convex_index(), ppi);
    
    getfem::mesh_im_level_set mimm(mls);
    mimm.set_integration_method(ref_m.convex_index(),
				getfem::int_method_descriptor(INTEGRATION1));
    mimm.set_simplex_im(getfem::int_method_descriptor(SIMPLEX_INTEGRATION1), 
			getfem::int_method_descriptor(SINGULAR_INTEGRATION1));
    
    // Integration methods on the boudary
    int intboundd = getfem::mesh_im_level_set::INTEGRATE_BOUNDARY;
    getfem::mesh_im_level_set mimboundd(mls_cont, intboundd,
					getfem::int_method_descriptor(SIMPLEX_INTEGRATION1));
    mimboundd.set_integration_method(ref_m_cont.convex_index(),
				     getfem::int_method_descriptor(INTEGRATION1));
    mimboundd.adapt();



    getfem::mesh_im_level_set mimbounddt(mls_tang_cont, intboundd,
					 getfem::int_method_descriptor(SIMPLEX_INTEGRATION1));
    if (!p.contact_only) {
     
      mimbounddt.set_integration_method( ref_m_tang_cont.convex_index(),
					getfem::int_method_descriptor(INTEGRATION1));
      mimbounddt.adapt();
    }



    //Interpolation of  U on a reference mesh 
    getfem::mesh_im ref_mim(ref_m);
    getfem::pintegration_method ppi = 
      getfem::int_method_descriptor("IM_TRIANGLE(6)");
    ref_mim.set_integration_method(ref_m.convex_index(), ppi);
    
    plain_vector interp_U(ref_mf.nb_dof());
    getfem::interpolation(p.mf_u(), ref_mf, U, interp_U);
    
    
    plain_vector interp_U_error(ref_mf.nb_dof());
    gmm::add(interp_U, gmm::scaled(ref_U, -1.), interp_U_error);
    gmm::vecsave(p.datafilename+".U_map_error", interp_U_error);
    
    cout << "To ref L2 ERROR on U:"
	 << getfem::asm_L2_dist(mimm, ref_mf, interp_U,
				ref_mf, ref_U) << endl;
    
    cout << "To ref H1 ERROR on U:"
	 << getfem::asm_H1_dist(mimm, ref_mf, interp_U,
				ref_mf, ref_U) << endl;
    gmm::add(gmm::scaled(interp_U, -1.), ref_U);
    gmm::vecsave(p.datafilename + ".diff_ref", ref_U);
    
    //Interpolation of  Lambda on a reference mesh 
    
    plain_vector interp_cont(ref_mf_cont.nb_dof());
    getfem::interpolation(p.mf_cont(), ref_mf_cont, Lambda, interp_cont);
    
    plain_vector interp_cont_error(ref_mf_cont.nb_dof());
    gmm::add(interp_cont, gmm::scaled(ref_cont, -1.), interp_cont_error);
    gmm::vecsave(p.datafilename+".cont_map_error", interp_cont_error);
    
    cout << "To ref L2 ERROR on P:"
	 << getfem::asm_L2_dist(mimboundd, ref_mf_cont, interp_cont,
				ref_mf_cont, ref_cont) << endl;
    //Interpolation of  Lambda_t on a reference mesh 
    plain_vector interp_tang_cont(ref_mf_tang_cont.nb_dof());
    if (!p.contact_only) {
     
      getfem::interpolation(p.mf_cont(), ref_mf_tang_cont, Lambda_t, interp_tang_cont);
      
      plain_vector interp_tang_cont_error(ref_mf_tang_cont.nb_dof());
      gmm::add(interp_tang_cont, gmm::scaled(ref_tang_cont, -1.), interp_tang_cont_error);
      gmm::vecsave(p.datafilename+".tang_cont_map_error", interp_tang_cont_error);
      
      cout << "To ref L2 ERROR on P tangent:"
	   << getfem::asm_L2_dist(mimbounddt, ref_mf_tang_cont, interp_tang_cont,
				  ref_mf_tang_cont, ref_tang_cont) << endl;
      
    }
  }
  
  return 0; 
  
  
}


