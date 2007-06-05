

#include "getfem/getfem_modeling.h"
#include "getfem/getfem_inter_element.h"

using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrices. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;


class mass_matrix_jump : public getfem::compute_on_inter_element {

 protected :

  sparse_matrix &M;
  
  virtual void compute_on_gauss_point
    (getfem::fem_interpolation_context ctx1, getfem::pfem pf1,
     getfem::fem_interpolation_context ctx2, getfem::pfem pf2,
     getfem::papprox_integration pai1) {
    
    scalar_type w = pai1->integration_coefficients()[ctx1.ii()];
    
    size_type cv1 = ctx1.convex_num();
    size_type cv2 = ctx2.convex_num();
    getfem::base_tensor t1, t2;
    pf1->real_base_value(ctx1, t1);
    pf2->real_base_value(ctx2, t2);

    for (size_type i = 0; i < pf1->nb_dof(cv1); ++i)
      for (size_type j = 0; j < pf1->nb_dof(cv1); ++j)
	M(mf.ind_dof_of_element(cv1)[i],
	  mf.ind_dof_of_element(cv1)[j]) += w * t1[i] * t1[j];
       
    for (size_type i = 0; i < pf2->nb_dof(cv2); ++i)
      for (size_type j = 0; j < pf2->nb_dof(cv2); ++j)
	M(mf.ind_dof_of_element(cv2)[i],
	  mf.ind_dof_of_element(cv2)[j]) += w * t2[i] * t2[j];
    
    for (size_type i = 0; i < pf1->nb_dof(cv1); ++i)
      for (size_type j = 0; j < pf2->nb_dof(cv2); ++j) {
	M(mf.ind_dof_of_element(cv1)[i],
	  mf.ind_dof_of_element(cv2)[j]) -= w * t1[i] * t2[j];
	M(mf.ind_dof_of_element(cv2)[j],
	  mf.ind_dof_of_element(cv1)[i]) -= w * t1[i] * t2[j];
      }
  }

 public :
  
  mass_matrix_jump(sparse_matrix &MM, getfem::mesh_im &mmim,
		   getfem::mesh_fem &mmf)
    : compute_on_inter_element(mmim, mmf), M(MM) {
    GMM_ASSERT1(mf.get_qdim() > 1,
	       "Vectorial elements not taken into account ... to be done");
  }

};



int main(void) {

  getfem::mesh m;

  getfem::mesh_fem mf(m);

  getfem::mesh_im mim(m);

  sparse_matrix M;
  mass_matrix_jump MMJ(M, mim, mf);


  getfem::mesh_region r_border;


  for (getfem::mr_visitor v(r_border); !v.finished(); ++v)
    MMJ.compute_on_face(v.cv(), v.f());





}
