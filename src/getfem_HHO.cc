/*===========================================================================

 Copyright (C) 2019-2019 Yves Renard

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


#include "getfem/getfem_HHO.h"


namespace getfem {

  THREAD_SAFE_STATIC bgeot::geotrans_precomp_pool HHO_pgp_pool;
  THREAD_SAFE_STATIC fem_precomp_pool HHO_pfp_pool;

  // To be optimized:
  // - The fact that (when pf2->target_dim() = 1) the
  //   problem can be solved componentwise can be more exploited in
  //   avoiding the computation of the whole matrix M2.
  // - The vectorization can be avoided in most cases
  

  class _HHO_reconstructed_gradient
    : public virtual_elementary_transformation {

  public:

    virtual void give_transformation(const mesh_fem &mf1, const mesh_fem &mf2,
                                     size_type cv, base_matrix &M) const {

      // The reconstructed Gradient "G" is described on mf2 and computed by
      // the formula on the element T :
      // \int_T G.w = \int_T Grad(v_T).w + \int_{dT}(v_{dT} - v_T).(w.n)
      // where "w" is the test function arbitrary in mf2, "v_T" is the field
      // inside the element whose gradient is to be reconstructed,
      // "v_{dT}" is the field on the boundary of T and "n" is the outward
      // unit normal.
      
      // Obtaining the fem descriptors
      pfem pf1 = mf1.fem_of_element(cv);
      pfem pf2 = mf2.fem_of_element(cv);
      pfem pfi = interior_fem_of_hho_method(pf1);

      size_type degree = std::max(pf1->estimated_degree(),
                                  pf2->estimated_degree());
      bgeot::pgeometric_trans pgt = mf1.linked_mesh().trans_of_convex(cv);
      papprox_integration pim
        = classical_approx_im(pgt, dim_type(2*degree))->approx_method();

      base_matrix G;
      bgeot::vectors_to_base_matrix(G, mf1.linked_mesh().points_of_convex(cv));

      bgeot::pgeotrans_precomp pgp
        = HHO_pgp_pool(pgt, pim->pintegration_points());
      pfem_precomp pfp1 = HHO_pfp_pool(pf1, pim->pintegration_points());
      pfem_precomp pfp2 = HHO_pfp_pool(pf2, pim->pintegration_points());
      pfem_precomp pfpi = HHO_pfp_pool(pfi, pim->pintegration_points());
      
      fem_interpolation_context ctx1(pgp, pfp1, 0, G, cv);
      fem_interpolation_context ctx2(pgp, pfp2, 0, G, cv);
      fem_interpolation_context ctxi(pgp, pfpi, 0, G, cv);

      size_type Q = mf1.get_qdim(), N = mf1.linked_mesh().dim();
      base_vector un(N);
      size_type qmult1 =  Q / pf1->target_dim();
      size_type ndof1 = pf1->nb_dof(cv) * qmult1;
      size_type qmult2 =  Q*N / pf2->target_dim();
      size_type ndof2 = pf2->nb_dof(cv) * qmult2;
      size_type qmulti =  Q / pfi->target_dim();
      size_type ndofi = pfi->nb_dof(cv) * qmulti;

      
      base_tensor t1, t2, ti, tv1;
      base_matrix tv2, tv1p, tvi;
      base_matrix M1(ndof2, ndof1), M2(ndof2, ndof2), M2inv(ndof2, ndof2);
      base_matrix aux2(ndof2, ndof2);

      // Integrals on the element : \int_T G.w (M2) and  \int_T Grad(v_T).w (M1)
      for (size_type ipt = 0; ipt < pim->nb_points_on_convex(); ++ipt) {
        ctx1.set_ii(ipt); ctx2.set_ii(ipt);
        scalar_type coeff = pim->coeff(ipt) * ctx1.J();
        
        ctx1.grad_base_value(t1);
        vectorize_grad_base_tensor(t1, tv1, ndof1, pf1->target_dim(), Q);

        ctx2.base_value(t2);
        vectorize_base_tensor(t2, tv2, ndof2, pf2->target_dim(), Q*N);
       
        gmm::mult(tv2, gmm::transposed(tv2), aux2);
        gmm::add(gmm::scaled(aux2, coeff), M2);

        for (size_type i = 0; i < ndof1; ++i) // To be optimized
          for (size_type j = 0; j < ndof2; ++j)
            for (size_type k = 0; k < Q*N; ++k)
              M1(j, i) += coeff * tv1.as_vector()[i+k*ndof1] * tv2(j, k);
      }

      // Integrals on the faces : \int_{dT}(v_{dT} - v_T).(w.n) (M1)
      for (short_type ifc = 0; ifc < pgt->structure()->nb_faces(); ++ifc) {
        ctx1.set_face_num(ifc); ctx2.set_face_num(ifc); ctxi.set_face_num(ifc);
        size_type first_ind = pim->ind_first_point_on_face(ifc);
        for (size_type ipt = 0; ipt < pim->nb_points_on_face(ifc); ++ipt) {
          ctx1.set_ii(first_ind+ipt);
          ctx2.set_ii(first_ind+ipt);
          ctxi.set_ii(first_ind+ipt);
          scalar_type coeff = pim->coeff(first_ind+ipt) * ctx1.J();
          
          ctx2.base_value(t2);
          vectorize_base_tensor(t2, tv2, ndof2, pf2->target_dim(), Q*N);

          ctx1.base_value(t1);
          vectorize_base_tensor(t1, tv1p, ndof1, pf1->target_dim(), Q);
          
          ctxi.base_value(ti);
          vectorize_base_tensor(ti, tvi, ndofi, pfi->target_dim(), Q);

          gmm::mult(ctx1.B(), pgt->normals()[ifc], un);

          for (size_type i = 0; i < ndof1; ++i) // To be optimized
            for (size_type j = 0; j < ndof2; ++j)
              for (size_type k1 = 0; k1 < Q; ++k1) {
                scalar_type b(0), a = coeff *
                  (tv1p(i, k1) - (i < ndofi ? tvi(i, k1) : 0.));
                for (size_type k2 = 0; k2 < N; ++k2)
                  b += a * tv2(j, k1 + k2*Q) * un[k2];
                M1(j, i) += b;
              }
        }
      }
      
      if (pf2->target_dim() == 1) {
        gmm::sub_slice I(0, ndof2/(N*Q), N*Q);
        gmm::lu_inverse(gmm::sub_matrix(M2, I, I));
        for (size_type i = 0; i < N*Q; ++i) {
          gmm::sub_slice I2(i, ndof2/(N*Q), N*Q);
          gmm::copy(gmm::sub_matrix(M2, I, I), gmm::sub_matrix(M2inv, I2, I2));
        }
      } else { gmm::copy(M2, M2inv); gmm::lu_inverse(M2inv); }
      
      
      gmm::mult(M2inv, M1, M);
      gmm::clean(M, gmm::vect_norminf(M.as_vector()) * 1E-13);
      // cout << "M = " << M << endl;
    }
  };

  void add_HHO_reconstructed_gradient(model &md, std::string name) {
    pelementary_transformation
      p = std::make_shared<_HHO_reconstructed_gradient>();
    md.add_elementary_transformation(name, p);
  }


  class _HHO_reconstructed_sym_gradient
    : public virtual_elementary_transformation {

  public:

    virtual void give_transformation(const mesh_fem &mf1, const mesh_fem &mf2,
                                     size_type cv, base_matrix &M) const {

      // The reconstructed symmetric Gradient "G" is described on mf2 and
      // computed by the formula on the element T :
      // \int_T G:w =   (1/2)*\int_T 0.5*Grad(v_T):(w+w^T)
      //              + (1/2)*\int_{dT}(v_{dT} - v_T).((w+w^T).n)
      // where "w" is the test function arbitrary in mf2, "v_T" is the field
      // inside the element whose gradient is to be reconstructed,
      // "v_{dT}" is the field on the boundary of T and "n" is the outward
      // unit normal.
      
      // Obtaining the fem descriptors
      pfem pf1 = mf1.fem_of_element(cv);
      pfem pf2 = mf2.fem_of_element(cv);
      pfem pfi = interior_fem_of_hho_method(pf1);

      size_type degree = std::max(pf1->estimated_degree(),
                                  pf2->estimated_degree());
      bgeot::pgeometric_trans pgt = mf1.linked_mesh().trans_of_convex(cv);
      papprox_integration pim
        = classical_approx_im(pgt, dim_type(2*degree))->approx_method();

      base_matrix G;
      bgeot::vectors_to_base_matrix(G, mf1.linked_mesh().points_of_convex(cv));

      bgeot::pgeotrans_precomp pgp
        = HHO_pgp_pool(pgt, pim->pintegration_points());
      pfem_precomp pfp1 = HHO_pfp_pool(pf1, pim->pintegration_points());
      pfem_precomp pfp2 = HHO_pfp_pool(pf2, pim->pintegration_points());
      pfem_precomp pfpi = HHO_pfp_pool(pfi, pim->pintegration_points());
      
      fem_interpolation_context ctx1(pgp, pfp1, 0, G, cv);
      fem_interpolation_context ctx2(pgp, pfp2, 0, G, cv);
      fem_interpolation_context ctxi(pgp, pfpi, 0, G, cv);

      size_type Q = mf1.get_qdim(), N = mf1.linked_mesh().dim();
      GMM_ASSERT1(Q == N, "This transformation works only for vector fields "
                  "having the same dimension as the domain");
      base_vector un(N);
      size_type qmult1 =  N / pf1->target_dim();
      size_type ndof1 = pf1->nb_dof(cv) * qmult1;
      size_type qmult2 =  N*N / pf2->target_dim();
      size_type ndof2 = pf2->nb_dof(cv) * qmult2;
      size_type qmulti =  N / pfi->target_dim();
      size_type ndofi = pfi->nb_dof(cv) * qmulti;

      
      base_tensor t1, t2, ti, tv1;
      base_matrix tv2, tv1p, tvi;
      base_matrix M1(ndof2, ndof1), M2(ndof2, ndof2), M2inv(ndof2, ndof2);
      base_matrix aux2(ndof2, ndof2);

      // Integrals on the element : \int_T G:w (M2)
      //                 and  (1/2)*\int_T 0.5*Grad(v_T):(w+w^T)
      for (size_type ipt = 0; ipt < pim->nb_points_on_convex(); ++ipt) {
        ctx1.set_ii(ipt); ctx2.set_ii(ipt);
        scalar_type coeff = pim->coeff(ipt) * ctx1.J();
        
        ctx1.grad_base_value(t1);
        vectorize_grad_base_tensor(t1, tv1, ndof1, pf1->target_dim(), N);

        ctx2.base_value(t2);
        vectorize_base_tensor(t2, tv2, ndof2, pf2->target_dim(), N*N);
       
        gmm::mult(tv2, gmm::transposed(tv2), aux2);
        gmm::add(gmm::scaled(aux2, coeff), M2);

        for (size_type i = 0; i < ndof1; ++i) // To be optimized
          for (size_type j = 0; j < ndof2; ++j)
            for (size_type k1 = 0; k1 < N; ++k1) 
              for (size_type k2 = 0; k2 < N; ++k2)
                M1(j, i) += coeff * tv1.as_vector()[i+(k1+k2*N)*ndof1]
                  * 0.5 * (tv2(j, k1+k2*N) + tv2(j, k2+k1*N));
      }

      // Integrals on the faces : (1/2)*\int_{dT}(v_{dT} - v_T).((w+w^T).n) (M1)
      for (short_type ifc = 0; ifc < pgt->structure()->nb_faces(); ++ifc) {
        ctx1.set_face_num(ifc); ctx2.set_face_num(ifc); ctxi.set_face_num(ifc);
        size_type first_ind = pim->ind_first_point_on_face(ifc);
        for (size_type ipt = 0; ipt < pim->nb_points_on_face(ifc); ++ipt) {
          ctx1.set_ii(first_ind+ipt);
          ctx2.set_ii(first_ind+ipt);
          ctxi.set_ii(first_ind+ipt);
          scalar_type coeff = pim->coeff(first_ind+ipt) * ctx1.J();
          
          ctx2.base_value(t2);
          vectorize_base_tensor(t2, tv2, ndof2, pf2->target_dim(), N*N);

          ctx1.base_value(t1);
          vectorize_base_tensor(t1, tv1p, ndof1, pf1->target_dim(), N);
          
          ctxi.base_value(ti);
          vectorize_base_tensor(ti, tvi, ndofi, pfi->target_dim(), N);

          gmm::mult(ctx1.B(), pgt->normals()[ifc], un);

          for (size_type i = 0; i < ndof1; ++i) // To be optimized
            for (size_type j = 0; j < ndof2; ++j)
              for (size_type k1 = 0; k1 < N; ++k1) {
                scalar_type b(0), a = coeff *
                  (tv1p(i, k1) - (i < ndofi ? tvi(i, k1) : 0.));
                for (size_type k2 = 0; k2 < N; ++k2)
                  b += a*0.5*(tv2(j, k1 + k2*N) + tv2(j, k2 + k1*N)) * un[k2];
                M1(j, i) += b;
              }
        }

      }
      
      if (pf2->target_dim() == 1) {
        gmm::sub_slice I(0, ndof2/(N*Q), N*Q);
        gmm::lu_inverse(gmm::sub_matrix(M2, I, I));
        for (size_type i = 0; i < N*Q; ++i) {
          gmm::sub_slice I2(i, ndof2/(N*Q), N*Q);
          gmm::copy(gmm::sub_matrix(M2, I, I), gmm::sub_matrix(M2inv, I2, I2));
        }
      } else 
        { gmm::copy(M2, M2inv); gmm::lu_inverse(M2inv); }
      
      gmm::mult(M2inv, M1, M);
      gmm::clean(M, gmm::vect_norminf(M.as_vector()) * 1E-13);
    }
  };

  void add_HHO_reconstructed_symmetrized_gradient(model &md, std::string name) {
    pelementary_transformation
      p = std::make_shared<_HHO_reconstructed_sym_gradient>();
    md.add_elementary_transformation(name, p);
  }

  

  class _HHO_reconstructed_value
    : public virtual_elementary_transformation {

  public:

    virtual void give_transformation(const mesh_fem &mf1, const mesh_fem &mf2,
                                     size_type cv, base_matrix &M) const {
      // The reconstructed variable "D" is described on mf2 and computed by
      // the formula on the element T :
      //   \int_T Grad(D).Grad(w) =   \int_T Grad(v_T).Grad(w)
      //                            + \int_{dT}(v_{dT} - v_T).(Grad(w).n)
      // with the constraint
      //   \int_T D = \int_T v_T
      // where "w" is the test function arbitrary in mf2, "v_T" is the field
      // inside the element whose gradient is to be reconstructed,
      // "v_{dT}" is the field on the boundary of T and "n" is the outward
      // unit normal.
      
      // Obtaining the fem descriptors
      pfem pf1 = mf1.fem_of_element(cv);
      pfem pf2 = mf2.fem_of_element(cv);
      pfem pfi = interior_fem_of_hho_method(pf1);

      size_type degree = std::max(pf1->estimated_degree(),
                                  pf2->estimated_degree());
      bgeot::pgeometric_trans pgt = mf1.linked_mesh().trans_of_convex(cv);
      papprox_integration pim
        = classical_approx_im(pgt, dim_type(2*degree))->approx_method();

      base_matrix G;
      bgeot::vectors_to_base_matrix(G, mf1.linked_mesh().points_of_convex(cv));

      bgeot::pgeotrans_precomp pgp
        = HHO_pgp_pool(pgt, pim->pintegration_points());
      pfem_precomp pfp1 = HHO_pfp_pool(pf1, pim->pintegration_points());
      pfem_precomp pfp2 = HHO_pfp_pool(pf2, pim->pintegration_points());
      pfem_precomp pfpi = HHO_pfp_pool(pfi, pim->pintegration_points());
      
      fem_interpolation_context ctx1(pgp, pfp1, 0, G, cv);
      fem_interpolation_context ctx2(pgp, pfp2, 0, G, cv);
      fem_interpolation_context ctxi(pgp, pfpi, 0, G, cv);

      size_type Q = mf1.get_qdim(), N = mf1.linked_mesh().dim();
      base_vector un(N);
      size_type qmult1 =  Q / pf1->target_dim();
      size_type ndof1 = pf1->nb_dof(cv) * qmult1;
      size_type qmult2 =  Q / pf2->target_dim();
      size_type ndof2 = pf2->nb_dof(cv) * qmult2;
      size_type qmulti =  Q / pfi->target_dim();
      size_type ndofi = pfi->nb_dof(cv) * qmulti;

      
      base_tensor t1, t2, ti, tv1, tv2, t1p, t2p;
      base_matrix tv1p, tv2p, tvi;
      base_matrix M1(ndof2, ndof1), M2(ndof2, ndof2), M2inv(ndof2, ndof2);
      base_matrix M3(Q, ndof1), M4(Q, ndof2);
      base_matrix aux1(ndof2, ndof1), aux2(ndof2, ndof2);

      // Integrals on the element : \int_T Grad(D).Grad(w) (M2)
      //                            \int_T Grad(v_T).Grad(w) (M1)
      //                            \int_T D (M4)  and \int_T v_T (M3)
      for (size_type ipt = 0; ipt < pim->nb_points_on_convex(); ++ipt) {
        ctx1.set_ii(ipt); ctx2.set_ii(ipt);
        scalar_type coeff = pim->coeff(ipt) * ctx1.J();
        
        ctx1.grad_base_value(t1);
        vectorize_grad_base_tensor(t1, tv1, ndof1, pf1->target_dim(), Q);

        ctx1.base_value(t1p);
        vectorize_base_tensor(t1p, tv1p, ndof1, pf1->target_dim(), Q);

        ctx2.grad_base_value(t2);
        vectorize_grad_base_tensor(t2, tv2, ndof2, pf2->target_dim(), Q);

        ctx2.base_value(t2p);
        vectorize_base_tensor(t2p, tv2p, ndof2, pf2->target_dim(), Q);

        for (size_type i = 0; i < ndof2; ++i) // To be optimized
          for (size_type j = 0; j < ndof2; ++j)
            for (size_type k = 0; k < Q*N; ++k)
              M2(j, i) += coeff * tv2.as_vector()[i+k*ndof2]
                                * tv2.as_vector()[j+k*ndof2];

        for (size_type i = 0; i < ndof2; ++i)
          for (size_type k = 0; k < Q; ++k)
            M4(k,  i) += coeff * tv2p(i, k);
              
        for (size_type i = 0; i < ndof1; ++i) // To be optimized
          for (size_type j = 0; j < ndof2; ++j)
            for (size_type k = 0; k < Q*N; ++k)
              M1(j, i) += coeff * tv1.as_vector()[i+k*ndof1]
                                * tv2.as_vector()[j+k*ndof2];

        for (size_type i = 0; i < ndof1; ++i)
          for (size_type k = 0; k < Q; ++k)
            M3(k,  i) += coeff * tv1p(i, k);

      }

      // Integrals on the faces : \int_{dT}(v_{dT} - v_T).(Grad(w).n) (M1)
      for (short_type ifc = 0; ifc < pgt->structure()->nb_faces(); ++ifc) {
        ctx1.set_face_num(ifc); ctx2.set_face_num(ifc); ctxi.set_face_num(ifc);
        size_type first_ind = pim->ind_first_point_on_face(ifc);
        for (size_type ipt = 0; ipt < pim->nb_points_on_face(ifc); ++ipt) {
          ctx1.set_ii(first_ind+ipt);
          ctx2.set_ii(first_ind+ipt);
          ctxi.set_ii(first_ind+ipt);
          scalar_type coeff = pim->coeff(first_ind+ipt) * ctx1.J();
          
          ctx2.grad_base_value(t2);
          vectorize_grad_base_tensor(t2, tv2, ndof2, pf2->target_dim(), Q);

          ctx1.base_value(t1);
          vectorize_base_tensor(t1, tv1p, ndof1, pf1->target_dim(), Q);
          
          ctxi.base_value(ti);
          vectorize_base_tensor(ti, tvi, ndofi, pfi->target_dim(), Q);

          gmm::mult(ctx1.B(), pgt->normals()[ifc], un);

          for (size_type i = 0; i < ndof1; ++i) // To be optimized
            for (size_type j = 0; j < ndof2; ++j)
              for (size_type k1 = 0; k1 < Q; ++k1) {
                scalar_type b(0), a = coeff *
                  (tv1p(i, k1) - (i < ndofi ? tvi(i, k1) : 0.));
                for (size_type k2 = 0; k2 < N; ++k2)
                  b += a * tv2.as_vector()[j+(k1 + k2*Q)*ndof2] * un[k2];
                M1(j, i) += b;
              }
        }

      }

      // Add the constraint with penalization
      gmm::mult(gmm::transposed(M4), M4, aux2);
      gmm::add (gmm::scaled(aux2, 1.E7), M2);
      gmm::mult(gmm::transposed(M4), M3, aux1);
      gmm::add (gmm::scaled(aux1, 1.E7), M1);
      
      if (pf2->target_dim() == 1 && Q > 1) {
        gmm::sub_slice I(0, ndof2/Q, Q);
        gmm::lu_inverse(gmm::sub_matrix(M2, I, I));
        for (size_type i = 0; i < Q; ++i) {
          gmm::sub_slice I2(i, ndof2/Q, Q);
          gmm::copy(gmm::sub_matrix(M2, I, I), gmm::sub_matrix(M2inv, I2, I2));
        }
      } else 
        { gmm::copy(M2, M2inv); gmm::lu_inverse(M2inv); }
      
      gmm::mult(M2inv, M1, M);
      gmm::clean(M, gmm::vect_norminf(M.as_vector()) * 1E-13);
    }
  };

  void add_HHO_reconstructed_value(model &md, std::string name) {
    pelementary_transformation
      p = std::make_shared<_HHO_reconstructed_value>();
    md.add_elementary_transformation(name, p);
  }


  class _HHO_reconstructed_sym_value
    : public virtual_elementary_transformation {

  public:

    virtual void give_transformation(const mesh_fem &mf1, const mesh_fem &mf2,
                                     size_type cv, base_matrix &M) const {
      // The reconstructed variable "D" is described on mf2 and computed by
      // the formula on the element T :
      //   \int_T Sym(Grad(D)).Grad(w) =   \int_T Sym(Grad(v_T)).Grad(w)
      //                            + \int_{dT}(v_{dT} - v_T).(Sym(Grad(w)).n)
      // with the constraints
      //   \int_T D = \int_T v_T
      //   \int_T Skew(Grad(D)) = 0.5\int_{dT}(n x v_{dT} - v_{dT} x n)
      // where "w" is the test function arbitrary in mf2, "v_T" is the field
      // inside the element whose gradient is to be reconstructed,
      // "v_{dT}" is the field on the boundary of T and "n" is the outward
      // unit normal.
      
      // Obtaining the fem descriptors
      pfem pf1 = mf1.fem_of_element(cv);
      pfem pf2 = mf2.fem_of_element(cv);
      pfem pfi = interior_fem_of_hho_method(pf1);

      size_type degree = std::max(pf1->estimated_degree(),
                                  pf2->estimated_degree());
      bgeot::pgeometric_trans pgt = mf1.linked_mesh().trans_of_convex(cv);
      papprox_integration pim
        = classical_approx_im(pgt, dim_type(2*degree))->approx_method();

      base_matrix G;
      bgeot::vectors_to_base_matrix(G, mf1.linked_mesh().points_of_convex(cv));

      bgeot::pgeotrans_precomp pgp
        = HHO_pgp_pool(pgt, pim->pintegration_points());
      pfem_precomp pfp1 = HHO_pfp_pool(pf1, pim->pintegration_points());
      pfem_precomp pfp2 = HHO_pfp_pool(pf2, pim->pintegration_points());
      pfem_precomp pfpi = HHO_pfp_pool(pfi, pim->pintegration_points());
      
      fem_interpolation_context ctx1(pgp, pfp1, 0, G, cv);
      fem_interpolation_context ctx2(pgp, pfp2, 0, G, cv);
      fem_interpolation_context ctxi(pgp, pfpi, 0, G, cv);

      size_type Q = mf1.get_qdim(), N = mf1.linked_mesh().dim();
      GMM_ASSERT1(Q == N, "This transformation works only for vector fields "
                  "having the same dimension as the domain");
      base_vector un(N);
      size_type qmult1 =  N / pf1->target_dim();
      size_type ndof1 = pf1->nb_dof(cv) * qmult1;
      size_type qmult2 =  N / pf2->target_dim();
      size_type ndof2 = pf2->nb_dof(cv) * qmult2;
      size_type qmulti =  N / pfi->target_dim();
      size_type ndofi = pfi->nb_dof(cv) * qmulti;

      
      base_tensor t1, t2, ti, tv1, tv2, t1p, t2p;
      base_matrix tv1p, tv2p, tvi;
      base_matrix M1(ndof2, ndof1), M2(ndof2, ndof2), M2inv(ndof2, ndof2);;
      base_matrix M3(N, ndof1), M4(N, ndof2);
      base_matrix M5(N*N, ndof1), M6(N*N, ndof2);
      base_matrix aux1(ndof2, ndof1), aux2(ndof2, ndof2);

      // Integrals on the element : \int_T Sym(Grad(D)).Grad(w) (M2)
      //                            \int_T Sym(Grad(v_T)).Grad(w) (M1)
      //                            \int_T D (M4)  and \int_T v_T (M3)
      //                            \int_T Skew(Grad(D)) (M6)
      for (size_type ipt = 0; ipt < pim->nb_points_on_convex(); ++ipt) {
        ctx1.set_ii(ipt); ctx2.set_ii(ipt);
        scalar_type coeff = pim->coeff(ipt) * ctx1.J();
        
        ctx1.grad_base_value(t1);
        vectorize_grad_base_tensor(t1, tv1, ndof1, pf1->target_dim(), N);

        ctx1.base_value(t1p);
        vectorize_base_tensor(t1p, tv1p, ndof1, pf1->target_dim(), N);

        ctx2.grad_base_value(t2);
        vectorize_grad_base_tensor(t2, tv2, ndof2, pf2->target_dim(), N);

        ctx2.base_value(t2p);
        vectorize_base_tensor(t2p, tv2p, ndof2, pf2->target_dim(), N);

        for (size_type i = 0; i < ndof2; ++i) // To be optimized
          for (size_type j = 0; j < ndof2; ++j)
            for (size_type k1 = 0; k1 < N; ++k1)
              for (size_type k2 = 0; k2 < N; ++k2)
                M2(j, i) += coeff * tv2.as_vector()[i+(k1+k2*N)*ndof2]
                  * 0.5 * (tv2.as_vector()[j+(k1+k2*N)*ndof2] +
                           tv2.as_vector()[j+(k2+k1*N)*ndof2]);

        for (size_type i = 0; i < ndof2; ++i)
          for (size_type k = 0; k < N; ++k)
            M4(k,  i) += coeff * tv2p(i, k);

        for (size_type i = 0; i < ndof2; ++i)
          for (size_type k1 = 0; k1 < N; ++k1)
            for (size_type k2 = 0; k2 < N; ++k2)
              M6(k1+k2*N, i) += 0.5*coeff*(tv2.as_vector()[i+(k1+k2*N)*ndof2] -
                                           tv2.as_vector()[i+(k2+k1*N)*ndof2]);
              
        for (size_type i = 0; i < ndof1; ++i) // To be optimized
          for (size_type j = 0; j < ndof2; ++j)
            for (size_type k1 = 0; k1 < N; ++k1)
              for (size_type k2 = 0; k2 < N; ++k2)
                M1(j, i) += coeff * tv1.as_vector()[i+(k1+k2*N)*ndof1]
                  * 0.5 * (tv2.as_vector()[j+(k1+k2*N)*ndof2] +
                           tv2.as_vector()[j+(k2+k1*N)*ndof2]);

        for (size_type i = 0; i < ndof1; ++i)
          for (size_type k = 0; k < N; ++k)
            M3(k,  i) += coeff * tv1p(i, k);

      }

      // Integrals on the faces : \int_{dT}(v_{dT} - v_T).(Sym(Grad(w)).n) (M1)
      //                          \int_{dT} Skew(n x v_{dT} - v_{dT} x n) (M5)
      for (short_type ifc = 0; ifc < pgt->structure()->nb_faces(); ++ifc) {
        ctx1.set_face_num(ifc); ctx2.set_face_num(ifc); ctxi.set_face_num(ifc);
        size_type first_ind = pim->ind_first_point_on_face(ifc);
        for (size_type ipt = 0; ipt < pim->nb_points_on_face(ifc); ++ipt) {
          ctx1.set_ii(first_ind+ipt);
          ctx2.set_ii(first_ind+ipt);
          ctxi.set_ii(first_ind+ipt);
          scalar_type coeff = pim->coeff(first_ind+ipt) * ctx1.J();
          
          ctx2.grad_base_value(t2);
          vectorize_grad_base_tensor(t2, tv2, ndof2, pf2->target_dim(), N);

          ctx1.base_value(t1);
          vectorize_base_tensor(t1, tv1p, ndof1, pf1->target_dim(), N);
          
          ctxi.base_value(ti);
          vectorize_base_tensor(ti, tvi, ndofi, pfi->target_dim(), N);

          gmm::mult(ctx1.B(), pgt->normals()[ifc], un);

          for (size_type i = 0; i < ndof1; ++i) // To be optimized
            for (size_type j = 0; j < ndof2; ++j)
              for (size_type k1 = 0; k1 < N; ++k1) {
                scalar_type b(0), a = coeff *
                  (tv1p(i, k1) - (i < ndofi ? tvi(i, k1) : 0.));
                for (size_type k2 = 0; k2 < N; ++k2)
                  b += a * 0.5 * (tv2.as_vector()[j+(k1 + k2*N)*ndof2] +
                                  tv2.as_vector()[j+(k2 + k1*N)*ndof2])* un[k2];
                M1(j, i) += b;
              }

          for (size_type i = 0; i < ndof1; ++i)
            for (size_type k1 = 0; k1 < N; ++k1)
              for (size_type k2 = 0; k2 < N; ++k2)
                M5(k1+k2*N, i) += 0.5 * coeff * (tv1p(i, k2) * un[k1] -
                                                 tv1p(i, k1) * un[k2]);
        }
      }

      // Add the constraint with penalization
      gmm::mult(gmm::transposed(M4), M4, aux2);
      gmm::add (gmm::scaled(aux2, 1.E7), M2);
      gmm::mult(gmm::transposed(M6), M6, aux2);
      gmm::add (gmm::scaled(aux2, 1.E7), M2);
      gmm::mult(gmm::transposed(M4), M3, aux1);
      gmm::add (gmm::scaled(aux1, 1.E7), M1);
      gmm::mult(gmm::transposed(M6), M5, aux1);
      gmm::add (gmm::scaled(aux1, 1.E7), M1);
      
      if (pf2->target_dim() == 1 && Q > 1) {
        gmm::sub_slice I(0, ndof2/Q, Q);
        gmm::lu_inverse(gmm::sub_matrix(M2, I, I));
        for (size_type i = 0; i < Q; ++i) {
          gmm::sub_slice I2(i, ndof2/Q, Q);
          gmm::copy(gmm::sub_matrix(M2, I, I), gmm::sub_matrix(M2inv, I2, I2));
        }
      } else 
        { gmm::copy(M2, M2inv); gmm::lu_inverse(M2inv); }

      gmm::mult(M2inv, M1, M);
      gmm::clean(M, gmm::vect_norminf(M.as_vector()) * 1E-13);
    }
  };

  void add_HHO_reconstructed_symmetrized_value(model &md, std::string name) {
    pelementary_transformation
      p = std::make_shared<_HHO_reconstructed_sym_value>();
    md.add_elementary_transformation(name, p);
  }


  class _HHO_stabilization
    : public virtual_elementary_transformation {

  public:

    virtual void give_transformation(const mesh_fem &mf1, const mesh_fem &mf2,
                                     size_type cv, base_matrix &M) const {
      // The reconstructed variable "S" is described on mf2 and computed by
      // S(v) = P_{\dT}(v_{dT} - D(v)  - P_T(v_T - D(v)))
      // where P__{\dT} et P_T are L2 projections on the boundary and on the
      // interior of T on the corresponding discrete spaces.
      // Note that P_{\dT}(v_{dT}) = v_{dT} and P_T(v_T) = v_T and D is
      // the reconstructed value on P^{k+1} given by the formula:
      //   \int_T Grad(D).Grad(w) =   \int_T Grad(v_T).Grad(w)
      //                            + \int_{dT}(v_{dT} - v_T).(Grad(w).n)
      // with the constraint
      //   \int_T D = \int_T v_T
      // where "w" is the test function arbitrary in mf2, "v_T" is the field
      // inside the element whose gradient is to be reconstructed,
      // "v_{dT}" is the field on the boundary of T and "n" is the outward
      // unit normal.
      // The implemented formula is
      // S(v) = v_{dT} - P_{\dT}D(v) - P_{\dT}(v_T) + P_{\dT}(P_T(D(v)))
      // by the mean of the projection matrix from P^{k+1} to the original space
      // and the projection matrix from interior space to the boundary space.
      // As it is built, S(v) is zero on interior dofs.
      
      GMM_ASSERT1(&mf1 == &mf2, "The HHO stabilization transformation is "
                  "only defined on the HHO space to itself");

      // Obtaining the fem descriptors
      pfem pf1 = mf1.fem_of_element(cv);
      short_type degree = pf1->estimated_degree();
      bgeot::pgeometric_trans pgt = mf1.linked_mesh().trans_of_convex(cv);
      pfem pf2 = classical_fem(pgt, short_type(degree + 1)); // Should be
                                         // changed for an interior PK method
      pfem pfi = interior_fem_of_hho_method(pf1);

      papprox_integration pim
        = classical_approx_im(pgt, dim_type(2*degree+2))->approx_method();

      base_matrix G;
      bgeot::vectors_to_base_matrix(G, mf1.linked_mesh().points_of_convex(cv));

      bgeot::pgeotrans_precomp pgp
        = HHO_pgp_pool(pgt, pim->pintegration_points());
      pfem_precomp pfp1 = HHO_pfp_pool(pf1, pim->pintegration_points());
      pfem_precomp pfp2 = HHO_pfp_pool(pf2, pim->pintegration_points());
      pfem_precomp pfpi = HHO_pfp_pool(pfi, pim->pintegration_points());
      
      fem_interpolation_context ctx1(pgp, pfp1, 0, G, cv);
      fem_interpolation_context ctx2(pgp, pfp2, 0, G, cv);
      fem_interpolation_context ctxi(pgp, pfpi, 0, G, cv);

      size_type Q = mf1.get_qdim(), N = mf1.linked_mesh().dim();
      base_vector un(N);
      size_type qmult1 =  Q / pf1->target_dim();
      size_type ndof1 = pf1->nb_dof(cv) * qmult1;
      size_type qmult2 =  Q / pf2->target_dim();
      size_type ndof2 = pf2->nb_dof(cv) * qmult2;
      size_type qmulti =  Q / pfi->target_dim();
      size_type ndofi = pfi->nb_dof(cv) * qmulti;

      
      base_tensor t1, t2, ti, tv1, tv2, t1p, t2p;
      base_matrix tv1p, tv2p, tvi;
      base_matrix M1(ndof2, ndof1), M2(ndof2, ndof2), M2inv(ndof2, ndof2);
      base_matrix M3(Q, ndof1), M4(Q, ndof2);
      base_matrix aux1(ndof2, ndof1), aux2(ndof2, ndof2);
      base_matrix M7(ndof1, ndof1), M7inv(ndof1, ndof1), M8(ndof1, ndof2);
      base_matrix M9(ndof1, ndof1), MD(ndof2, ndof1);

      // Integrals on the element : \int_T Grad(D).Grad(w) (M2)
      //                            \int_T Grad(v_T).Grad(w) (M1)
      //                            \int_T D (M4)  and \int_T v_T (M3)
      for (size_type ipt = 0; ipt < pim->nb_points_on_convex(); ++ipt) {
        ctx1.set_ii(ipt); ctx2.set_ii(ipt);
        scalar_type coeff = pim->coeff(ipt) * ctx1.J();
        
        ctx1.grad_base_value(t1);
        vectorize_grad_base_tensor(t1, tv1, ndof1, pf1->target_dim(), Q);

        ctx1.base_value(t1p);
        vectorize_base_tensor(t1p, tv1p, ndof1, pf1->target_dim(), Q);

        ctx2.grad_base_value(t2);
        vectorize_grad_base_tensor(t2, tv2, ndof2, pf2->target_dim(), Q);

        ctx2.base_value(t2p);
        vectorize_base_tensor(t2p, tv2p, ndof2, pf2->target_dim(), Q);

        for (size_type i = 0; i < ndof2; ++i) // To be optimized
          for (size_type j = 0; j < ndof2; ++j)
            for (size_type k = 0; k < Q*N; ++k)
              M2(j, i) += coeff * tv2.as_vector()[i+k*ndof2]
                                * tv2.as_vector()[j+k*ndof2];

        for (size_type i = 0; i < ndof2; ++i)
          for (size_type k = 0; k < Q; ++k)
            M4(k,  i) += coeff * tv2p(i, k);
              
        for (size_type i = 0; i < ndof1; ++i) // To be optimized
          for (size_type j = 0; j < ndof2; ++j)
            for (size_type k = 0; k < Q*N; ++k)
              M1(j, i) += coeff * tv1.as_vector()[i+k*ndof1]
                                * tv2.as_vector()[j+k*ndof2];

        for (size_type i = 0; i < ndof1; ++i)
          for (size_type k = 0; k < Q; ++k)
            M3(k,  i) += coeff * tv1p(i, k);

        for (size_type i = 0; i < ndof1; ++i) // To be optimized
          for (size_type j = 0; j < ndof1; ++j)
            for (size_type k = 0; k < Q; ++k)
              M7(i, j) += coeff * tv1p(i, k) * tv1p(j, k);
        
        for (size_type i = 0; i < ndof1; ++i) // To be optimized
          for (size_type j = 0; j < ndof2; ++j)
            for (size_type k = 0; k < Q; ++k)
              M8(i, j) += coeff * tv1p(i, k) * tv2p(j, k);

      }

      // Integrals on the faces : \int_{dT}(v_{dT} - v_T).(Grad(w).n) (M1)
      for (short_type ifc = 0; ifc < pgt->structure()->nb_faces(); ++ifc) {
        ctx1.set_face_num(ifc); ctx2.set_face_num(ifc); ctxi.set_face_num(ifc);
        size_type first_ind = pim->ind_first_point_on_face(ifc);
        for (size_type ipt = 0; ipt < pim->nb_points_on_face(ifc); ++ipt) {
          ctx1.set_ii(first_ind+ipt);
          ctx2.set_ii(first_ind+ipt);
          ctxi.set_ii(first_ind+ipt);
          scalar_type coeff = pim->coeff(first_ind+ipt) * ctx1.J();
          
          ctx2.grad_base_value(t2);
          vectorize_grad_base_tensor(t2, tv2, ndof2, pf2->target_dim(), Q);

          ctx2.base_value(t2p);
          vectorize_base_tensor(t2p, tv2p, ndof2, pf2->target_dim(), Q);

          ctx1.base_value(t1);
          vectorize_base_tensor(t1, tv1p, ndof1, pf1->target_dim(), Q);
          
          ctxi.base_value(ti);
          vectorize_base_tensor(ti, tvi, ndofi, pfi->target_dim(), Q);

          gmm::mult(ctx1.B(), pgt->normals()[ifc], un);

          for (size_type i = 0; i < ndof1; ++i) // To be optimized
            for (size_type j = 0; j < ndof2; ++j)
              for (size_type k1 = 0; k1 < Q; ++k1) {
                scalar_type b(0), a = coeff *
                  (tv1p(i, k1) - (i < ndofi ? tvi(i, k1) : 0.));
                for (size_type k2 = 0; k2 < N; ++k2)
                  b += a * tv2.as_vector()[j+(k1 + k2*Q)*ndof2] * un[k2];
                M1(j, i) += b;
              }

          for (size_type i = 0; i < ndof1; ++i) // To be optimized
            for (size_type j = 0; j < ndof1; ++j)
              for (size_type k = 0; k < Q; ++k)
                M7(i, j) += coeff * tv1p(i,k) * tv1p(j, k);

          for (size_type i = 0; i < ndof1; ++i) // To be optimized
            for (size_type j = 0; j < ndof2; ++j)
              for (size_type k = 0; k < Q; ++k)
                M8(i, j) += coeff * tv1p(i,k) * tv2p(j, k);

          for (size_type i = 0; i < ndof1; ++i) // To be optimized
            for (size_type j = 0; j < ndofi; ++j)
              for (size_type k = 0; k < Q; ++k)
                M9(i, j) += coeff * tv1p(i,k) * tvi(j, k); 
        }
      }

      // Add the constraint with penalization
      gmm::mult(gmm::transposed(M4), M4, aux2);
      gmm::add (gmm::scaled(aux2, 1.E7), M2);
      gmm::mult(gmm::transposed(M4), M3, aux1);
      gmm::add (gmm::scaled(aux1, 1.E7), M1);

      if (pf2->target_dim() == 1 && Q > 1) {
        gmm::sub_slice I(0, ndof2/Q, Q);
        gmm::lu_inverse(gmm::sub_matrix(M2, I, I));
        for (size_type i = 0; i < Q; ++i) {
          gmm::sub_slice I2(i, ndof2/Q, Q);
          gmm::copy(gmm::sub_matrix(M2, I, I), gmm::sub_matrix(M2inv, I2, I2));
        }
      } else 
        { gmm::copy(M2, M2inv); gmm::lu_inverse(M2inv); }
      
      if (pf1->target_dim() == 1 && Q > 1) {
        gmm::sub_slice I(0, ndof1/Q, Q);
        gmm::lu_inverse(gmm::sub_matrix(M7, I, I));
        for (size_type i = 0; i < Q; ++i) {
          gmm::sub_slice I2(i, ndof1/Q, Q);
          gmm::copy(gmm::sub_matrix(M7, I, I), gmm::sub_matrix(M7inv, I2, I2));
        }
      } else
        { gmm::copy(M7, M7inv); gmm::lu_inverse(M7inv); }
      
      gmm::mult(M2inv, M1, MD);
      gmm::clean(MD, gmm::vect_norminf(MD.as_vector()) * 1E-13);

      // S  = (I - inv(M7)*M9)(I - inv(M7)*M8*MD)
      base_matrix MPB(ndof1, ndof1);
      gmm::mult(M7inv, M9, MPB);
      gmm::copy(gmm::identity_matrix(), M9);
      gmm::add(gmm::scaled(MPB, scalar_type(-1)), M9);

      base_matrix MPC(ndof1, ndof1), MPD(ndof1, ndof1);
      gmm::mult(M8, MD, MPC);
      gmm::mult(M7inv, MPC, MPD);
      gmm::copy(gmm::identity_matrix(), M7);
      gmm::add(gmm::scaled(MPD, scalar_type(-1)), M7);

      gmm::mult(M9, M7, M);
      gmm::clean(M, 1E-13);
    }
  };

  void add_HHO_stabilization(model &md, std::string name) {
    pelementary_transformation
      p = std::make_shared<_HHO_stabilization>();
    md.add_elementary_transformation(name, p);
  }


  class _HHO_symmetrized_stabilization
    : public virtual_elementary_transformation {

  public:

    virtual void give_transformation(const mesh_fem &mf1, const mesh_fem &mf2,
                                     size_type cv, base_matrix &M) const {
      // The reconstructed variable "S" is described on mf2 and computed by
      // S(v) = P_{\dT}(v_{dT} - D(v)  - P_T(v_T - D(v)))
      // where P__{\dT} et P_T are L2 projections on the boundary and on the
      // interior of T on the corresponding discrete spaces.
      // Note that P_{\dT}(v_{dT}) = v_{dT} and P_T(v_T) = v_T and D is
      // the reconstructed value on P^{k+1} given by the formula:
      //   \int_T Sym(Grad(D)).Grad(w) =   \int_T Sym(Grad(v_T)).Grad(w)
      //                            + \int_{dT}(v_{dT} - v_T).(Sym(Grad(w)).n)
      // with the constraints
      //   \int_T D = \int_T v_T
      //   \int_T Skew(Grad(D)) = 0.5\int_{dT}(n x v_{dT} - v_{dT} x n)
      // where "w" is the test function arbitrary in mf2, "v_T" is the field
      // inside the element whose gradient is to be reconstructed,
      // "v_{dT}" is the field on the boundary of T and "n" is the outward
      // unit normal.
      // The implemented formula is
      // S(v) = v_{dT} - P_{\dT}D(v) - P_{\dT}(v_T) + P_{\dT}(P_T(D(v)))
      // by the mean of the projection matrix from P^{k+1} to the original space
      // and the projection matrix from interior space to the boundary space.
      // As it is built, S(v) is zero on interior dofs.
      
      GMM_ASSERT1(&mf1 == &mf2, "The HHO stabilization transformation is "
                  "only defined on the HHO space to itself");
      
      // Obtaining the fem descriptors
      pfem pf1 = mf1.fem_of_element(cv);
      short_type degree = pf1->estimated_degree();
      bgeot::pgeometric_trans pgt = mf1.linked_mesh().trans_of_convex(cv);
      pfem pf2 = classical_fem(pgt, short_type(degree + 1)); // Should be changed to an
                                                 // interior PK method
      pfem pfi = interior_fem_of_hho_method(pf1);

      papprox_integration pim
        = classical_approx_im(pgt, dim_type(2*degree+2))->approx_method();

      base_matrix G;
      bgeot::vectors_to_base_matrix(G, mf1.linked_mesh().points_of_convex(cv));

      bgeot::pgeotrans_precomp pgp
        = HHO_pgp_pool(pgt, pim->pintegration_points());
      pfem_precomp pfp1 = HHO_pfp_pool(pf1, pim->pintegration_points());
      pfem_precomp pfp2 = HHO_pfp_pool(pf2, pim->pintegration_points());
      pfem_precomp pfpi = HHO_pfp_pool(pfi, pim->pintegration_points());
      
      fem_interpolation_context ctx1(pgp, pfp1, 0, G, cv);
      fem_interpolation_context ctx2(pgp, pfp2, 0, G, cv);
      fem_interpolation_context ctxi(pgp, pfpi, 0, G, cv);

      size_type Q = mf1.get_qdim(), N = mf1.linked_mesh().dim();
      GMM_ASSERT1(Q == N, "This transformation works only for vector fields "
                  "having the same dimension as the domain");
      base_vector un(N);
      size_type qmult1 =  N / pf1->target_dim();
      size_type ndof1 = pf1->nb_dof(cv) * qmult1;
      size_type qmult2 =  N / pf2->target_dim();
      size_type ndof2 = pf2->nb_dof(cv) * qmult2;
      size_type qmulti =  N / pfi->target_dim();
      size_type ndofi = pfi->nb_dof(cv) * qmulti;

      
      base_tensor t1, t2, ti, tv1, tv2, t1p, t2p;
      base_matrix tv1p, tv2p, tvi;
      base_matrix M1(ndof2, ndof1), M2(ndof2, ndof2), M2inv(ndof2, ndof2);
      base_matrix M3(N, ndof1), M4(N, ndof2);
      base_matrix aux1(ndof2, ndof1), aux2(ndof2, ndof2);
      base_matrix M5(N*N, ndof1), M6(N*N, ndof2);
      base_matrix M7(ndof1, ndof1), M7inv(ndof1, ndof1), M8(ndof1, ndof2);
      base_matrix M9(ndof1, ndof1), MD(ndof2, ndof1);

      // Integrals on the element : \int_T Sym(Grad(D)).Grad(w) (M2)
      //                            \int_T Sym(Grad(v_T)).Grad(w) (M1)
      //                            \int_T D (M4)  and \int_T v_T (M3)
      //                            \int_T Skew(Grad(D)) (M6)
      for (size_type ipt = 0; ipt < pim->nb_points_on_convex(); ++ipt) {
        ctx1.set_ii(ipt); ctx2.set_ii(ipt);
        scalar_type coeff = pim->coeff(ipt) * ctx1.J();
        
        ctx1.grad_base_value(t1);
        vectorize_grad_base_tensor(t1, tv1, ndof1, pf1->target_dim(), N);

        ctx1.base_value(t1p);
        vectorize_base_tensor(t1p, tv1p, ndof1, pf1->target_dim(), N);

        ctx2.grad_base_value(t2);
        vectorize_grad_base_tensor(t2, tv2, ndof2, pf2->target_dim(), N);

        ctx2.base_value(t2p);
        vectorize_base_tensor(t2p, tv2p, ndof1, pf1->target_dim(), N);

        for (size_type i = 0; i < ndof2; ++i) // To be optimized
          for (size_type j = 0; j < ndof2; ++j)
           for (size_type k1 = 0; k1 < N; ++k1)
             for (size_type k2 = 0; k2 < N; ++k2)
               M2(j, i) += coeff * tv2.as_vector()[i+(k1+k2*N)*ndof2]
                 * 0.5 * (tv2.as_vector()[j+(k1+k2*N)*ndof2] +
                          tv2.as_vector()[j+(k2+k1*N)*ndof2]);
        
        for (size_type i = 0; i < ndof2; ++i)
          for (size_type k = 0; k < N; ++k)
            M4(k,  i) += coeff * tv2p(i, k);
              
        for (size_type i = 0; i < ndof2; ++i)
          for (size_type k1 = 0; k1 < N; ++k1)
            for (size_type k2 = 0; k2 < N; ++k2)
              M6(k1+k2*N, i) += 0.5*coeff*(tv2.as_vector()[i+(k1+k2*N)*ndof2] -
                                           tv2.as_vector()[i+(k2+k1*N)*ndof2]);

        for (size_type i = 0; i < ndof1; ++i) // To be optimized
          for (size_type j = 0; j < ndof2; ++j)
            for (size_type k1 = 0; k1 < N; ++k1)
              for (size_type k2 = 0; k2 < N; ++k2)
                M1(j, i) += coeff * tv1.as_vector()[i+(k1+k2*N)*ndof1]
                  * 0.5 * (tv2.as_vector()[j+(k1+k2*N)*ndof2] +
                           tv2.as_vector()[j+(k2+k1*N)*ndof2]);
        
        for (size_type i = 0; i < ndof1; ++i)
          for (size_type k = 0; k < N; ++k)
            M3(k,  i) += coeff * tv1p(i, k);

        for (size_type i = 0; i < ndof1; ++i) // To be optimized
          for (size_type j = 0; j < ndof1; ++j)
            for (size_type k = 0; k < N; ++k)
              M7(i, j) += coeff * tv1p(i,k) * tv1p(j, k);

        for (size_type i = 0; i < ndof1; ++i) // To be optimized
          for (size_type j = 0; j < ndof2; ++j)
            for (size_type k = 0; k < N; ++k)
              M8(i, j) += coeff * tv1p(i,k) * tv2p(j, k);

      }

      // Integrals on the faces : \int_{dT}(v_{dT} - v_T).(Grad(w).n) (M1)
      //                          \int_{dT} Skew(n x v_{dT} - v_{dT} x n) (M5)
      for (short_type ifc = 0; ifc < pgt->structure()->nb_faces(); ++ifc) {
        ctx1.set_face_num(ifc); ctx2.set_face_num(ifc); ctxi.set_face_num(ifc);
        size_type first_ind = pim->ind_first_point_on_face(ifc);
        for (size_type ipt = 0; ipt < pim->nb_points_on_face(ifc); ++ipt) {
          ctx1.set_ii(first_ind+ipt);
          ctx2.set_ii(first_ind+ipt);
          ctxi.set_ii(first_ind+ipt);
          scalar_type coeff = pim->coeff(first_ind+ipt) * ctx1.J();
          
          ctx2.grad_base_value(t2);
          vectorize_grad_base_tensor(t2, tv2, ndof2, pf2->target_dim(), N);

          ctx2.base_value(t2p);
          vectorize_base_tensor(t2p, tv2p, ndof2, pf2->target_dim(), N);

          ctx1.base_value(t1);
          vectorize_base_tensor(t1, tv1p, ndof1, pf1->target_dim(), N);
          
          ctxi.base_value(ti);
          vectorize_base_tensor(ti, tvi, ndofi, pfi->target_dim(), N);

          gmm::mult(ctx1.B(), pgt->normals()[ifc], un);

          for (size_type i = 0; i < ndof1; ++i) // To be optimized
            for (size_type j = 0; j < ndof2; ++j)
              for (size_type k1 = 0; k1 < N; ++k1) {
                scalar_type b(0), a = coeff *
                  (tv1p(i, k1) - (i < ndofi ? tvi(i, k1) : 0.));
                for (size_type k2 = 0; k2 < N; ++k2)
                  b += a * 0.5 * (tv2.as_vector()[j+(k1 + k2*N)*ndof2] +
                                  tv2.as_vector()[j+(k2 + k1*N)*ndof2])* un[k2];
                M1(j, i) += b;
              }

          for (size_type i = 0; i < ndof1; ++i)
            for (size_type k1 = 0; k1 < N; ++k1)
              for (size_type k2 = 0; k2 < N; ++k2)
                M5(k1+k2*N, i) += 0.5 * coeff * (tv1p(i, k2) * un[k1] -
                                                 tv1p(i, k1) * un[k2]);

          for (size_type i = 0; i < ndof1; ++i) // To be optimized
            for (size_type j = 0; j < ndof1; ++j)
              for (size_type k = 0; k < N; ++k)
                M7(i, j) += coeff * tv1p(i,k) * tv1p(j, k);

          for (size_type i = 0; i < ndof1; ++i) // To be optimized
            for (size_type j = 0; j < ndof2; ++j)
              for (size_type k = 0; k < N; ++k)
                M8(i, j) += coeff * tv1p(i,k) * tv2p(j, k);

          for (size_type i = 0; i < ndof1; ++i) // To be optimized
            for (size_type j = 0; j < ndofi; ++j)
              for (size_type k = 0; k < N; ++k)
                M9(i, j) += coeff * tv1p(i,k) * tvi(j, k);
        }
      }

      // Add the constraint with penalization
      gmm::mult(gmm::transposed(M4), M4, aux2);
      gmm::add (gmm::scaled(aux2, 1E7), M2);
      gmm::mult(gmm::transposed(M6), M6, aux2);
      gmm::add (gmm::scaled(aux2, 1E7), M2);
      gmm::mult(gmm::transposed(M4), M3, aux1);
      gmm::add (gmm::scaled(aux1, 1E7), M1);
      gmm::mult(gmm::transposed(M6), M5, aux1);
      gmm::add (gmm::scaled(aux1, 1E7), M1);

      if (pf2->target_dim() == 1 && Q > 1) {
        gmm::sub_slice I(0, ndof2/Q, Q);
        gmm::lu_inverse(gmm::sub_matrix(M2, I, I));
        for (size_type i = 0; i < Q; ++i) {
          gmm::sub_slice I2(i, ndof2/Q, Q);
          gmm::copy(gmm::sub_matrix(M2, I, I), gmm::sub_matrix(M2inv, I2, I2));
        }
      } else 
        { gmm::copy(M2, M2inv); gmm::lu_inverse(M2inv); }
      
      if (pf1->target_dim() == 1 && Q > 1) {
        gmm::sub_slice I(0, ndof1/Q, Q);
        gmm::lu_inverse(gmm::sub_matrix(M7, I, I));
        for (size_type i = 0; i < Q; ++i) {
          gmm::sub_slice I2(i, ndof1/Q, Q);
          gmm::copy(gmm::sub_matrix(M7, I, I), gmm::sub_matrix(M7inv, I2, I2));
        }
      } else
        { gmm::copy(M7, M7inv); gmm::lu_inverse(M7inv); }
      
      gmm::mult(M2inv, M1, MD);
      gmm::clean(MD, gmm::vect_norminf(MD.as_vector()) * 1E-13);
      
      // S  = (I - inv(M7)*M9)(I - inv(M7)*M8*MD)
      base_matrix MPB(ndof1, ndof1);
      gmm::mult(M7inv, M9, MPB);
      gmm::copy(gmm::identity_matrix(), M9);
      gmm::add(gmm::scaled(MPB, scalar_type(-1)), M9);

      base_matrix MPC(ndof1, ndof1), MPD(ndof1, ndof1);
      gmm::mult(M8, MD, MPC);
      gmm::mult(M7inv, MPC, MPD);
      gmm::copy(gmm::identity_matrix(), M7);
      gmm::add(gmm::scaled(MPD, scalar_type(-1)), M7);

      gmm::mult(M9, M7, M);
      gmm::clean(M, 1E-13);
    }
  };

  void add_HHO_symmetrized_stabilization(model &md, std::string name) {
    pelementary_transformation
      p = std::make_shared<_HHO_symmetrized_stabilization>();
    md.add_elementary_transformation(name, p);
  }





}  /* end of namespace getfem.                                             */

