/*===========================================================================
 
 Copyright (C) 2007-2012 Yves Renard, Julien Pommier.
 
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

#include "getfem/getfem_models.h"
#include "getfem/getfem_import.h"
#include "getfem/getfem_inter_element.h"
#include "getfem/getfem_export.h"
#include "gmm/gmm_superlu_interface.h"

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

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
  scalar_type coeff_;
  
  virtual void compute_on_gauss_point
    (getfem::fem_interpolation_context ctx1, getfem::pfem pf1,
     getfem::fem_interpolation_context ctx2, getfem::pfem pf2,
     getfem::papprox_integration pai1) {
    
    scalar_type w = pai1->integration_coefficients()[ctx1.ii()];
    base_small_vector up(mf.linked_mesh().dim());
    const base_matrix& B = ctx1.B();
    gmm::mult(B, pgt1->normals()[f1], up);
    scalar_type norm = gmm::vect_norm2(up);
    scalar_type J = ctx1.J() * norm;

    size_type cv1 = ctx1.convex_num();
    size_type cv2 = ctx2.convex_num();
    getfem::base_tensor t1, t2;
    pf1->real_base_value(ctx1, t1);
    pf2->real_base_value(ctx2, t2);

    for (size_type i = 0; i < pf1->nb_dof(cv1); ++i)
      for (size_type j = 0; j < pf1->nb_dof(cv1); ++j)
	M(mf.ind_basic_dof_of_element(cv1)[i],
	  mf.ind_basic_dof_of_element(cv1)[j])
	  += coeff_ * J * w * t1[i] * t1[j];
       
    for (size_type i = 0; i < pf2->nb_dof(cv2); ++i)
      for (size_type j = 0; j < pf2->nb_dof(cv2); ++j)
	M(mf.ind_basic_dof_of_element(cv2)[i],
	  mf.ind_basic_dof_of_element(cv2)[j])
	  += coeff_ * J * w * t2[i] * t2[j];
    
    for (size_type i = 0; i < pf1->nb_dof(cv1); ++i)
      for (size_type j = 0; j < pf2->nb_dof(cv2); ++j) {
	M(mf.ind_basic_dof_of_element(cv1)[i],
	  mf.ind_basic_dof_of_element(cv2)[j])
	  -= coeff_ * J * w * t1[i] * t2[j];
	M(mf.ind_basic_dof_of_element(cv2)[j],
	  mf.ind_basic_dof_of_element(cv1)[i])
	  -= coeff_ * J * w * t1[i] * t2[j];
      }
  }

 public :
  
  mass_matrix_jump(sparse_matrix &MM, const getfem::mesh_im &mmim,
		   const getfem::mesh_fem &mmf, scalar_type coeff = 1)
    : compute_on_inter_element(mmim, mmf), M(MM), coeff_(coeff) {
    GMM_ASSERT1(mf.get_qdim() <= 1,
	       "Vectorial elements not taken into account ... to be done");
  }

};

/*! post_proc

  A short post_proc class for simple post-processing operations.
*/
class post_proc {

protected:
  const getfem::mesh_fem *_mf; //!< Mesh fem structure for postprocessing.

public:
  /*!
    \brief A simple constructor. Just linking the mesh fem structure.
    \param mf the mesh fem structure for post processing.
  */
  post_proc(const getfem::mesh_fem &mf){ _mf = &mf; }; 
    
  /*!
    \brief Return the difference of potential V(x2)-V(x1).
    \param X the vector to represent.
    \param x1 the coordinates of the first point.
    \param x2 the coordinates of the second point.
    \param DDP the value of the difference of potential computed.
  */
  template <class typeX, class Scalar>
  void emf(const typeX &X,
	   const getfem::base_node &x1,
	   const getfem::base_node &x2,
	   Scalar &DDP) {

    getfem::base_node vecdir(x2-x1);
    getfem::slicer_half_space a0(x1, vecdir, 0);
    getfem::slicer_half_space a1(x1, getfem::base_node
				 (vecdir[1], -vecdir[0]), 0);
    getfem::slicer_half_space a2(x2, vecdir, 0);
    // We use here two parallel 1D slicers and one perpendicular to the
    // others.
    getfem::stored_mesh_slice sl, sl2;
    getfem::slicer_build_stored_mesh_slice a3(sl), a4(sl2);
    getfem::mesh_slicer slicer((this->_mf)->linked_mesh()),
      slicer2((this->_mf)->linked_mesh());
    slicer.push_back_action(a2);
    slicer.push_back_action(a1);
    slicer.push_back_action(a3);
    slicer.exec(1);
    
    typeX Xinterp(sl.nb_points(), 0.);
    sl.interpolate(*(this->_mf), X, Xinterp); // Data are interpolated on the
					      // 0D slice.
    slicer2.push_back_action(a0);
    slicer2.push_back_action(a1);
    slicer2.push_back_action(a4);
    slicer2.exec(1);
    
    typeX Xinterp2(sl2.nb_points(), 0.);
    sl2.interpolate(*(this->_mf), X, Xinterp2); // Data are interpolated on
						// the 0D slice.
    DDP = Xinterp[0] - Xinterp2[0];
  }

};

int main(void) {

  getfem::mesh m;
  // Import mesh from gmsh.
  getfem::import_mesh("gmshv2:square.msh", m);
  // 2d problem.
  getfem::maybe_remove_last_dimension(m);

  // Mesh fem variables for the solution and the data.
  getfem::mesh_fem mf(m), mfd(m);
  mf.set_finite_element(m.convex_index(),
			getfem::fem_descriptor("FEM_PK(2,2)"));
  mfd.set_finite_element(m.convex_index(),
			 getfem::fem_descriptor("FEM_PK_DISCONTINUOUS(2,0)"));

  // Dissocation between the upper and the lower region.
  int LOWER_REGION = 1, UPPER_REGION = 3,
    LOWER_BOUNDARY = 10, UPPER_BOUNDARY = 30;
  getfem::mesh_region rg = m.region(UPPER_REGION);
  for (getfem::mr_visitor i(rg); !i.finished(); ++i) {
    mf.set_dof_partition(i.cv(), UPPER_REGION);
  }
  rg = m.region(LOWER_REGION);
  for (getfem::mr_visitor i(rg); !i.finished(); ++i) {
    mf.set_dof_partition(i.cv(), LOWER_REGION);
  }

  // Extraction of the boundary separating the upper and lower region.
  getfem::mesh_region r_border;
  getfem::outer_faces_of_mesh(m, r_border);
  
  getfem::outer_faces_of_mesh(m, m.region(UPPER_REGION).index(), 
			      m.region(UPPER_BOUNDARY));
  getfem::outer_faces_of_mesh(m, m.region(LOWER_REGION).index(), 
			      m.region(LOWER_BOUNDARY));
  for (getfem::mr_visitor v(r_border); !v.finished(); ++v) {
    m.region(UPPER_BOUNDARY).sup(v.cv(), v.f());
  }
  
  // Selection of the integration method.
  getfem::mesh_im mim(m);
  mim.set_integration_method(m.convex_index(), 
			     getfem::int_method_descriptor("IM_TRIANGLE(5)"));

  // Creation of the matrix for the problem assembly.
  sparse_matrix M(mf.nb_dof(), mf.nb_dof());
  // Creation of the RHS and of the solution vectors.
  plain_vector B(mf.nb_dof(), 0.), X(mf.nb_dof(), 0.);
  // Creation of a vector containing the medium coefficients.
  plain_vector A(mfd.nb_dof(), 0.);
  
  // Characteristic constants for the different regions.
  double slow(1.), supper(1.);
  rg = m.region(1);
  GMM_ASSERT1(!mfd.is_reduced(), "To be adapted");
  for (getfem::mr_visitor i(rg); !i.finished(); ++i) {
    A[ mfd.ind_basic_dof_of_element( i.cv() )[0] ] = slow;
  }
  rg = m.region(3);
  for (getfem::mr_visitor i(rg); !i.finished(); ++i) {
    A[ mfd.ind_basic_dof_of_element( i.cv() )[0] ] = supper;
  }

  // Assembly of the stiffness matrix for the laplacian.
  getfem::asm_stiffness_matrix_for_laplacian(M, mim, mf, mfd, A);

  // Assembly of the jump on in the interface.  
  GMM_ASSERT1(!mf.is_reduced(), "To be adapted");
  mass_matrix_jump MMJ(M, mim, mf, 2.);
  getfem::mesh_region interface = m.region(UPPER_BOUNDARY);
  for (getfem::mr_visitor v(interface); !v.finished(); ++v)
    MMJ.compute_on_face(v.cv(), v.f());

  // Set the Dirichlet constraints.
  for (unsigned int j = 0; j < mf.nb_dof(); j++) {
    if ( mf.point_of_basic_dof(j)[1] == .5 ) X[j] = 1.;
  }
  getfem::assembling_Dirichlet_condition(M, B, mf, 1010, X);

  // Solution of the linear system.
  double condest;
  gmm::SuperLU_solve(M, X, B, condest);

  // Post-processing of the solution.
  post_proc sortie(mf);
  // The emf is computed between the points x1 and x2.
  getfem::base_node x1(0., -0.0001), x2(0., .0001);
  scalar_type ddp;
  sortie.emf(X, x1, x2, ddp);
  std::cout << "Jump of the potential V on the interface :" << std::endl;
  std::cout << "V" << x2 << " - " << "V" << x1
       << " = " << ddp << std::endl;
}
