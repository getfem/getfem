/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2004-2015 Yves Renard
 
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


#include "getfem/getfem_linearized_plates.h"
#include "getfem/getfem_contact_and_friction_common.h" // for vectorize_base_tensor


namespace getfem {


  // For the moment, projection onto RT0 element
  // and with no storage.

  class elementary_projection_transformation
    : public virtual_elementary_transformation {

  public:

    virtual void give_transformation(const mesh_fem &mf, size_type cv,
                                     base_matrix &M) const{


      
      // Obtaining the fem descriptors
      pfem pf1 = mf.fem_of_element(cv);
      size_type N = pf1->dim();
      size_type qmult =  N / pf1->target_dim();

      bool simplex = false;
      if (pf1->ref_convex(cv) == bgeot::simplex_of_reference(dim_type(N))) {
        simplex = true;
      } else if (pf1->ref_convex(cv)
                 == bgeot::parallelepiped_of_reference(dim_type(N))) {
        simplex = false;
      } else {
        GMM_ASSERT1(false, "Cannot adapt the method for such an element.");
      }

      GMM_ASSERT1(pf1->is_equivalent(), "For tau-equivalent fem only."); // Indeed no, for the moment ...

      std::stringstream fem_desc;
      fem_desc << "FEM_RT0" << (simplex ? "":"Q") << "(" << N << ")";
      pfem pf2 = fem_descriptor(fem_desc.str());

      // Obtaining a convenient integration method
      size_type degree = pf1->estimated_degree() *  pf2->estimated_degree();
      bgeot::pgeometric_trans pgt = mf.linked_mesh().trans_of_convex(cv);
      papprox_integration pim
        = classical_approx_im(pgt, dim_type(degree))->approx_method();

      // Computation of mass matrices
      size_type ndof1 = pf1->nb_dof(cv) * qmult;
      size_type ndof2 = pf2->nb_dof(0);
      base_matrix M1(ndof1, ndof1), M2(ndof2, ndof2), B(ndof1, ndof2);
      base_matrix aux0(ndof1, ndof1), aux1(ndof1, ndof2), aux2(ndof1, ndof2);
      base_matrix aux3(ndof2, ndof2);

      
      base_matrix G;
      bgeot::vectors_to_base_matrix(G, mf.linked_mesh().points_of_convex(cv));
      fem_interpolation_context ctx1(pgt, pf1, base_node(N), G, cv);
      fem_interpolation_context ctx2(pgt, pf2, base_node(N), G, cv);

      base_tensor t1, t2;
      base_matrix tv1, tv2;
        
      for (size_type i = 0; i < pim->nb_points_on_convex(); ++i) {

        scalar_type coeff = pim->coeff(i); // Mult by ctx.J() not useful here
        ctx1.set_xref(pim->point(i));
        ctx2.set_xref(pim->point(i));    
        pf1->real_base_value(ctx1, t1);
        vectorize_base_tensor(t1, tv1, ndof1, pf1->target_dim(), N);
        pf2->real_base_value(ctx2, t2);
        vectorize_base_tensor(t2, tv2, ndof2, pf2->target_dim(), N);

        if (N == 2) { // H curl ... ! and for triangles ?
          gmm::clear(tv2);
          tv2(0,0) = scalar_type(1);
          tv2(1,1) = scalar_type(1);
          tv2(2,0) = pim->point(i)[1];
          tv2(3,1) = pim->point(i)[0];
        }
       
        gmm::mult(tv1, gmm::transposed(tv1), aux0);
        gmm::add(gmm::scaled(aux0, coeff), M1);
        gmm::mult(tv2, gmm::transposed(tv2), aux3);
        gmm::add(gmm::scaled(aux3, coeff), M2);
        gmm::mult(tv1, gmm::transposed(tv2), aux1);
        gmm::add(gmm::scaled(aux1, coeff), B);
      }
      
      
      // Computation of M
      gmm::lu_inverse(M1);
      gmm::lu_inverse(M2);
      gmm::mult(M1, B, aux1);
      gmm::mult(aux1, M2, aux2);
      GMM_ASSERT1(gmm::mat_nrows(M) == ndof1,
                  "Element not convenient for projection");
      gmm::mult(aux2, gmm::transposed(B), M);
      gmm::clean(M, 1E-15);
      // cout << "M = " << M << endl;
    }
  };



  void add_RT0_projection(model &md) {

     pelementary_transformation p
      = new elementary_projection_transformation();

    md.add_elementary_transformation("RT0_projection", p);
  }







}  /* end of namespace getfem.                                             */

