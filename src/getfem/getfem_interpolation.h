/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2001-2015 Yves Renard, Julien Pommier
 
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
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/

/**@file getfem_interpolation.h
   @author Yves Renard    <Yves.Renard@insa-lyon.fr>
   @author Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date October 15, 2001.
   @brief Interpolation of fields from a mesh_fem onto another.
*/

#ifndef GETFEM_INTERPOLATION_H__
#define GETFEM_INTERPOLATION_H__

#include "getfem_mesh_fem.h"
#include "bgeot_geotrans_inv.h"
#include "dal_tree_sorted.h"
#include "getfem_im_data.h"
#include "getfem_torus.h"

namespace getfem {

  /* ********************************************************************* */
  /*                                                                       */
  /*        I. Distribution of a set of points on a mesh.                  */
  /*                                                                       */
  /* ********************************************************************* */

  class mesh_trans_inv : public bgeot::geotrans_inv {

  protected :
    typedef std::set<size_type>::const_iterator set_iterator;
    typedef std::map<size_type,size_type>::const_iterator map_iterator;

    const mesh &msh;
    std::vector<std::set<size_type> > pts_cvx;
    std::vector<base_node> ref_coords;
    std::map<size_type,size_type> ids;

  public :

    size_type nb_points_on_convex(size_type i) const
    { return pts_cvx[i].size(); }
    void points_on_convex(size_type i, std::vector<size_type> &itab) const;
    size_type point_on_convex(size_type cv, size_type i) const;
    const std::vector<base_node> &reference_coords(void) const { return ref_coords; }

    void add_point_with_id(base_node n, size_type id)
    { size_type ipt = add_point(n); ids[ipt] = id; }
    size_type id_of_point(size_type ipt) const;
    const mesh &linked_mesh(void) const { return msh; }

    /* extrapolation = 0 : Only the points inside the mesh are distributed.
     * extrapolation = 1 : Try to extrapolate the exterior points near the
     *                     boundary.
     * extrapolation = 2 : Extrapolate all the exterior points. Could be
     *                     expensive.
     *
     * if rg_source is provided only the corresponding part of the mesh is
     * taken into account and extrapolation is done with respect to the
     * boundary of the specified region. rg_source must contain only convexes.
     */
    void distribute(int extrapolation = 0,
                    mesh_region rg_source=mesh_region::all_convexes());
    mesh_trans_inv(const mesh &m, double EPS_ = 1E-12)
      : bgeot::geotrans_inv(EPS_), msh(m) {}
  };
  

  /* ********************************************************************* */
  /*                                                                       */
  /*        II. Interpolation of functions.                                */
  /*                                                                       */
  /* ********************************************************************* */


  template <typename VECT, typename F, typename M>
  inline void interpolation_function__(const mesh_fem &mf, VECT &V,
                                       F &f, const dal::bit_vector &dofs,
                                       const M &, gmm::abstract_null_type) {
    size_type Q = mf.get_qdim();
    GMM_ASSERT1(gmm::vect_size(V) == mf.nb_basic_dof() && Q == 1,
                "Dof vector has not the right size");
    for (dal::bv_visitor i(dofs); !i.finished(); ++i)
      V[i] = f(mf.point_of_basic_dof(i));
  }

  template <typename VECT, typename F, typename M>
  inline void interpolation_function__(const mesh_fem &mf, VECT &V,
                                       F &f, const dal::bit_vector &dofs,
                                       const M &v, gmm::abstract_vector) {
    size_type N = gmm::vect_size(v),  Q = mf.get_qdim();
    GMM_ASSERT1(gmm::vect_size(V) == mf.nb_basic_dof()*N/Q,
                "Dof vector has not the right size");
    for (dal::bv_visitor i(dofs); !i.finished(); ++i)
      if (i % Q == 0)
        gmm::copy(f(mf.point_of_basic_dof(i)),
                  gmm::sub_vector(V, gmm::sub_interval(i*N/Q, N)));
  }

  template <typename VECT, typename F, typename M>
  inline void interpolation_function__(const mesh_fem &mf, VECT &V,
                                       F &f, const dal::bit_vector &dofs,
                                       const M &mm, gmm::abstract_matrix) {
    // typedef typename gmm::linalg_traits<VECT>::value_type T;
    size_type Nr = gmm::mat_nrows(mm), Nc = gmm::mat_ncols(mm), N = Nr*Nc;
    size_type Q = mf.get_qdim();
    base_matrix m(Nr, Nc);
    GMM_ASSERT1(gmm::vect_size(V) == mf.nb_basic_dof()*N/Q,
                "Dof vector has not the right size");
    for (dal::bv_visitor i(dofs); !i.finished(); ++i)
      if (i % Q == 0) {
        gmm::copy(f(mf.point_of_basic_dof(i)), m);
        for (size_type j = 0; j < Nc; ++j)
          gmm::copy(gmm::mat_col(m, j),
                    gmm::sub_vector(V, gmm::sub_interval(i*N/Q+j*Nr, Nr)));
      }

  }

  template <typename VECT, typename F, typename M>
  inline void interpolation_function_(const mesh_fem &mf, VECT &V,
                                      F &f, const dal::bit_vector &dofs,
                                      const M &m) {
    interpolation_function__(mf, V, f, dofs, m,
                             typename gmm::linalg_traits<M>::linalg_type());
  }

#if GETFEM_PARA_LEVEL > 0
  template <typename T>
  void take_one_op(void *a, void *b, int *len, MPI_Datatype *) {
    T aa = *((T*)a);
    return aa ? aa : *((T*)b);
  }

  template <typename T>
  inline MPI_Op mpi_take_one_op(T) {
    static bool isinit = false;
    static MPI_Op op;
    if (!isinit) {
      MPI_Op_create(take_one_op<T>, true, &op);
      isinit = true;
    }
    return op;
  }
#endif

  // TODO : verify that rhs is a lagrange fem
  /**
     @brief interpolation of a function f on mf_target.
     - mf_target must be of lagrange type.
     - mf_target's qdim should be equal to the size of the return value of f,
       or equal to 1
     - V should have the right size
     CAUTION: with the parallized version (GETFEM_PARA_LEVEL >= 2) the
     resulting vector V is distributed.
  */
  template <typename VECT, typename F>
  void interpolation_function(mesh_fem &mf_target, const VECT &VV, F &f,
                         mesh_region rg=mesh_region::all_convexes()) {
    typedef typename gmm::linalg_traits<VECT>::value_type T;
    size_type qqdimt = gmm::vect_size(VV) / mf_target.nb_dof();
    std::vector<T> V(mf_target.nb_basic_dof()*qqdimt);
    mf_target.linked_mesh().intersect_with_mpi_region(rg);
    dal::bit_vector dofs = mf_target.basic_dof_on_region(rg);
    if (dofs.card() > 0)
      interpolation_function_(mf_target, V, f, dofs,
                              f(mf_target.point_of_basic_dof(dofs.first())));

    if (mf_target.is_reduced()) {
      for (size_type k = 0; k < qqdimt; ++k)
        gmm::mult(mf_target.reduction_matrix(),
                  gmm::sub_vector(V,
                                  gmm::sub_slice(k, mf_target.nb_basic_dof(),
                                                 qqdimt)),
                  gmm::sub_vector(const_cast<VECT &>(VV),
                                  gmm::sub_slice(k, mf_target.nb_dof(),
                                                 qqdimt)));
    }
    else
      gmm::copy(V, const_cast<VECT &>(VV));
  }

  /* ********************************************************************* */
  /*                                                                       */
  /*        III. Interpolation between two meshes.                         */
  /*                                                                       */
  /* ********************************************************************* */

  /* ------------------------------ Interface -----------------------------*/

  /**
     @brief interpolation/extrapolation of (mf_source, U) on mf_target.
     - mf_target must be of lagrange type.
     - mf_target's qdim should be equal to mf_source qdim, or equal to 1
     - U.size() >= mf_source.get_qdim()
     - V.size() >= (mf_target.nb_dof() / mf_target.get_qdim())
                   * mf_source.get_qdim()

     With extrapolation = 0 a strict interpolation is done, with extrapolation = 1
     an extrapolation of the exterior points near the boundary is done (if any)
     and with extrapolation = 2 all  exterior points are extrapolated (could be expensive).

     If both mesh_fem shared the same mesh object, a fast interpolation
     will be used.

     If rg_source and rg_target are provided the operation is restricted to
     these regions. rg_source must contain only convexes.
  */
  template<typename VECTU, typename VECTV>
  void interpolation(const mesh_fem &mf_source, const mesh_fem &mf_target,
                     const VECTU &U, VECTV &V, int extrapolation = 0,
                     double EPS = 1E-10,
                     mesh_region rg_source=mesh_region::all_convexes(),
                     mesh_region rg_target=mesh_region::all_convexes());

  /**
     @brief Build the interpolation matrix of mf_source on mf_target.
     the matrix M is
     such that (V = M*U) == interpolation(mf_source, mf_target, U, V).

     Useful for repeated interpolations.
     For performance reasons the matrix M is recommended to be either
     a row or a row and column matrix.

     If rg_source and rg_target are provided the operation is restricted to
     these regions. rg_source must contain only convexes.
   */
  template<typename MAT>
  void interpolation(const mesh_fem &mf_source, const mesh_fem &mf_target,
                     MAT &M, int extrapolation = 0, double EPS = 1E-10,
                     mesh_region rg_source=mesh_region::all_convexes(),
                     mesh_region rg_target=mesh_region::all_convexes());


  /* --------------------------- Implementation ---------------------------*/

  /*
     interpolation of a solution on same mesh.
     - &mf_target.linked_mesh() == &mf_source.linked_mesh()
     - mf_target must be of lagrange type.
     - mf_target's qdim should be equal to mf_source qdim, or equal to 1
     - U.size() >= mf_source.get_qdim()
     - V.size() >= (mf_target.nb_dof() / mf_target.get_qdim())
                   * mf_source.get_qdim()
  */
  template<typename VECTU, typename VECTV, typename MAT>
    void interpolation_same_mesh(const mesh_fem &mf_source,
                                 const mesh_fem &mf_target,
                                 const VECTU &UU, VECTV &VV,
                                 MAT &MM, int version) {

    typedef typename gmm::linalg_traits<VECTU>::value_type T;
    base_matrix G;
    dim_type qdim = mf_source.get_qdim();
    dim_type qqdim = dim_type(gmm::vect_size(UU)/mf_source.nb_dof());
    std::vector<T> val(qdim);
    std::vector<std::vector<T> > coeff;
    std::vector<size_type> dof_source;
    GMM_ASSERT1(qdim == mf_target.get_qdim() || mf_target.get_qdim() == 1,
                "Attempt to interpolate a field of dimension "
                << qdim << " on a mesh_fem whose Qdim is " << 
                int(mf_target.get_qdim()));
    size_type qmult = mf_source.get_qdim()/mf_target.get_qdim();
    size_type qqdimt = qqdim * mf_source.get_qdim()/mf_target.get_qdim();
    fem_precomp_pool fppool;
    std::vector<size_type> dof_t_passes(mf_target.nb_basic_dof());
    std::vector<T> U(mf_source.nb_basic_dof()*qqdim);
    std::vector<T> V(mf_target.nb_basic_dof()*qqdimt);
    gmm::row_matrix<gmm::rsvector<scalar_type> >
      M(mf_target.nb_basic_dof(), mf_source.nb_basic_dof());

    if (version == 0) mf_source.extend_vector(UU, U);

    /* we should sort convexes by their fem */
    for (dal::bv_visitor cv(mf_source.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt=mf_source.linked_mesh().trans_of_convex(cv);
      pfem pf_s = mf_source.fem_of_element(cv);
      if (!mf_target.convex_index().is_in(cv)) 
        continue;
      pfem pf_t = mf_target.fem_of_element(cv);
      size_type nbd_s = pf_s->nb_dof(cv);
      size_type nbd_t = pf_t->nb_dof(cv);
      mesh_fem::ind_dof_ct::const_iterator itdof;
      size_type cvnbdof = mf_source.nb_basic_dof_of_element(cv);

      bool discontinuous_source = false;
      for (size_type dof=0; dof < nbd_s; ++dof)
        if (!dof_linkable(pf_s->dof_types()[dof])) {
          discontinuous_source = true;
          break;
        }

      if (version == 0) {
        coeff.resize(qqdim);
        for (size_type qq=0; qq < qqdim; ++qq) {
          coeff[qq].resize(cvnbdof);
          itdof = mf_source.ind_basic_dof_of_element(cv).begin();
          for (size_type k = 0; k < cvnbdof; ++k, ++itdof) {
            coeff[qq][k] = U[(*itdof)*qqdim+qq];
          }
        }
      }
      if (pf_s->need_G()) 
        bgeot::vectors_to_base_matrix
          (G, mf_source.linked_mesh().points_of_convex(cv));

      GMM_ASSERT1(pf_t->target_dim() == 1,
                  "won't interpolate on a vector FEM... ");
      pfem_precomp pfp = fppool(pf_s, pf_t->node_tab(cv));
      fem_interpolation_context ctx(pgt,pfp,size_type(-1), G, cv,
                                    short_type(-1));
      itdof = mf_target.ind_basic_dof_of_element(cv).begin();
      if (version != 0) {
        const mesh_fem::ind_dof_ct &idct
          = mf_source.ind_basic_dof_of_element(cv);
        dof_source.assign(idct.begin(), idct.end());
      }
      for (size_type i = 0; i < nbd_t; ++i, itdof+=mf_target.get_qdim()) {
        size_type dof_t = *itdof*qmult;
        if (!discontinuous_source && dof_t_passes[*itdof] > 0) continue;
        dof_t_passes[*itdof] += 1;
        ctx.set_ii(i);
        if (version == 0) {
          for (size_type qq=0; qq < qqdim; ++qq) {
            pf_s->interpolation(ctx, coeff[qq], val, qdim);
            for (size_type k=0; k < qdim; ++k)
              V[(dof_t + k)*qqdim+qq] += val[k];
          }
        }
        else {
          base_matrix Mloc(qdim, mf_source.nb_basic_dof_of_element(cv));
          pf_s->interpolation(ctx, Mloc, qdim);
          for (size_type k=0; k < qdim; ++k) {
            for (size_type j=0; j < dof_source.size(); ++j) {
              M(dof_t + k, dof_source[j]) += Mloc(k, j);
            }
          }
        }
      }
    }

    // calculate averages for discontinuous source and continuous target
    for (size_type i = 0; i < mf_target.nb_basic_dof(); ++i) {
      size_type dof_t = i*qmult;
      scalar_type passes = scalar_type(dof_t_passes[i]);
      if (version == 0 && passes > scalar_type(0))
        for (size_type qq=0; qq < qqdim; ++qq)
          for (size_type k=0; k < qdim; ++k)
            V[(dof_t + k)*qqdim+qq] /= passes;
      else if (passes > scalar_type(0))
        for (size_type k=0; k < qdim; ++k)
          for (size_type j=0; j < dof_source.size(); ++j)
            gmm::scale(gmm::mat_row(M, dof_t + k), scalar_type(1)/passes);
    }
    
    if (version == 0)
      mf_target.reduce_vector(V, VV);
    else {
      if (mf_target.is_reduced())
        if (mf_source.is_reduced()) {
          gmm::row_matrix<gmm::rsvector<scalar_type> >
            MMM(mf_target.nb_dof(), mf_source.nb_basic_dof());
          gmm::mult(mf_target.reduction_matrix(), M, MMM);
          gmm::mult(MMM, mf_source.extension_matrix(), MM);
        }
        else
          gmm::mult(mf_target.reduction_matrix(), M, MM); 
      else
        if (mf_source.is_reduced())
          gmm::mult(M, mf_source.extension_matrix(), MM);
        else
          gmm::copy(M, MM);
    }
  }

  
  /*
     interpolation of a solution on another mesh.
     - mti contains the points where to interpolate.
     - the solution should be continuous.
   */
  template<typename VECTU, typename VECTV, typename MAT>
  void interpolation(const mesh_fem &mf_source,
                     mesh_trans_inv &mti,
                     const VECTU &UU, VECTV &V, MAT &MM,
                     int version, int extrapolation = 0,
                     dal::bit_vector *dof_untouched = 0,
                     mesh_region rg_source=mesh_region::all_convexes()) {

    typedef typename gmm::linalg_traits<VECTU>::value_type T;
    const mesh &msh(mf_source.linked_mesh());
    dim_type qdim_s = mf_source.get_qdim();
    dim_type qqdim = dim_type(gmm::vect_size(UU)/mf_source.nb_dof());

    std::vector<T> U(mf_source.nb_basic_dof()*qqdim);
    gmm::row_matrix<gmm::rsvector<scalar_type> >
      M(gmm::mat_nrows(MM), mf_source.nb_basic_dof());

    if (version == 0) mf_source.extend_vector(UU, U);

    mti.distribute(extrapolation, rg_source);
    std::vector<size_type> itab;    
    base_matrix G;

    /* interpolation */
    dal::bit_vector points_to_do; points_to_do.add(0, mti.nb_points());
    std::vector<T> val(qdim_s);
    std::vector<std::vector<T> > coeff;
    base_tensor Z;
    std::vector<size_type> dof_source;
    
    for (dal::bv_visitor cv(mf_source.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt = msh.trans_of_convex(cv);
      mti.points_on_convex(cv, itab);
      if (itab.size() == 0) continue;

      pfem pf_s = mf_source.fem_of_element(cv);
      if (pf_s->need_G()) 
        bgeot::vectors_to_base_matrix(G, msh.points_of_convex(cv));

      fem_interpolation_context ctx(pgt, pf_s, base_node(), G, cv,
                                    short_type(-1));
      if (version == 0) {
        coeff.resize(qqdim);
        size_type cvnbdof = mf_source.nb_basic_dof_of_element(cv);
        mesh_fem::ind_dof_ct::const_iterator itdof;
        for (size_type qq=0; qq < qqdim; ++qq) {
          coeff[qq].resize(cvnbdof);
          itdof = mf_source.ind_basic_dof_of_element(cv).begin();
          for (size_type k = 0; k < cvnbdof; ++k, ++itdof) {
            coeff[qq][k] = U[(*itdof)*qqdim+qq];
          }
        }
      }
      if (version != 0) {
        const mesh_fem::ind_dof_ct &idct
          = mf_source.ind_basic_dof_of_element(cv);
        dof_source.assign(idct.begin(), idct.end());
      }
      for (size_type i = 0; i < itab.size(); ++i) {
        size_type ipt = itab[i];
        if (points_to_do.is_in(ipt)) {
          points_to_do.sup(ipt);
          ctx.set_xref(mti.reference_coords()[ipt]);
          size_type dof_t = mti.id_of_point(ipt);
          size_type pos = dof_t * qdim_s;
          if (version == 0) {
            for (size_type qq=0; qq < qqdim; ++qq) {           
              pf_s->interpolation(ctx, coeff[qq], val, qdim_s);
              for (size_type k=0; k < qdim_s; ++k)
                V[(pos + k)*qqdim+qq] = val[k];
            }
            // Part to be improved if one wants in option to be able to
            // interpolate the gradient.
            //          if (PVGRAD) {
            // base_matrix grad(mdim, qdim);
            // pf_s->interpolation_grad(ctx,coeff,gmm::transposed(grad), qdim);
            // std::copy(grad.begin(), grad.end(), V.begin()+dof_t*qdim*mdim);
            // }
          } else {
            base_matrix Mloc(qdim_s, mf_source.nb_basic_dof_of_element(cv));
            pf_s->interpolation(ctx, Mloc, qdim_s);
            for (size_type k=0; k < qdim_s; ++k) {
              for (size_type j=0; j < gmm::mat_ncols(Mloc); ++j)
                M(pos+k, dof_source[j]) = Mloc(k,j);
                /* does not work with col matrices
                   gmm::clear(gmm::mat_row(M, pos+k));
                   gmm::copy(gmm::mat_row(Mloc, k),
                   gmm::sub_vector(gmm::mat_row(M, pos+k), isrc));
                */
            }
          }
        }
      }
    }
    if (points_to_do.card() != 0) {
      if (dof_untouched) {
        dof_untouched->clear();
        for (dal::bv_visitor ipt(points_to_do); !ipt.finished(); ++ipt)
          dof_untouched->add(mti.id_of_point(ipt));
      }
      else {
        dal::bit_vector dofs_to_do;
        for (dal::bv_visitor ipt(points_to_do); !ipt.finished(); ++ipt)
          dofs_to_do.add(mti.id_of_point(ipt));
        GMM_WARNING2("in interpolation (different meshes),"
                     << dofs_to_do.card() << " dof of target mesh_fem have "
                     << " been missed\nmissing dofs : " << dofs_to_do);
      }
    }

    if (version != 0) {
      if (mf_source.is_reduced())
        gmm::mult(M, mf_source.extension_matrix(), MM);
      else
        gmm::copy(M, MM);
    }
    
  }

  template<typename VECTU, typename VECTV>
  void interpolation(const mesh_fem &mf_source, mesh_trans_inv &mti,
                     const VECTU &U, VECTV &V, int extrapolation = 0,
                     dal::bit_vector *dof_untouched = 0,
                     mesh_region rg_source=mesh_region::all_convexes()) {
    base_matrix M;
    GMM_ASSERT1((gmm::vect_size(U) % mf_source.nb_dof()) == 0 &&
                gmm::vect_size(V)!=0, "Dimension of vector mismatch");
    interpolation(mf_source, mti, U, V, M, 0, extrapolation, dof_untouched, rg_source);
  }



  /*
     interpolation of a solution on another mesh.
     - mf_target must be of lagrange type.
     - the solution should be continuous..
   */
  template<typename VECTU, typename VECTV, typename MAT>
  void interpolation(const mesh_fem &mf_source, const mesh_fem &mf_target,
                     const VECTU &U, VECTV &VV, MAT &MM,
                     int version, int extrapolation,
                     double EPS,
                     mesh_region rg_source=mesh_region::all_convexes(),
                     mesh_region rg_target=mesh_region::all_convexes()) {

    //Check if it is a torus mesh_fem
    const torus_mesh_fem * pmf_torus = dynamic_cast<const torus_mesh_fem *>(&mf_target);
    if(pmf_torus != 0){
      interpolation_to_torus_mesh_fem(mf_source, *pmf_torus, 
                                      U, VV, MM, version, extrapolation, EPS, rg_source, rg_target);
      return;
    }

    typedef typename gmm::linalg_traits<VECTU>::value_type T;
    dim_type qqdim = dim_type(gmm::vect_size(U)/mf_source.nb_dof());
    size_type qqdimt = qqdim * mf_source.get_qdim()/mf_target.get_qdim();
    std::vector<T> V(mf_target.nb_basic_dof()*qqdimt);
    mf_target.extend_vector(VV,V);
    gmm::row_matrix<gmm::rsvector<scalar_type> >
      M(mf_target.nb_basic_dof(), mf_source.nb_dof());

    const mesh &msh(mf_source.linked_mesh());
    getfem::mesh_trans_inv mti(msh, EPS);
    size_type qdim_s = mf_source.get_qdim(), qdim_t = mf_target.get_qdim();
    GMM_ASSERT1(qdim_s == qdim_t || qdim_t == 1,
                "Attempt to interpolate a field of dimension "
                << qdim_s << " on a mesh_fem whose Qdim is " << qdim_t);
    
    /* test if the target mesh_fem is really of Lagrange type.         */
    for (dal::bv_visitor cv(mf_target.convex_index()); !cv.finished();++cv) {
      pfem pf_t = mf_target.fem_of_element(cv);
      GMM_ASSERT1(pf_t->target_dim() == 1 && pf_t->is_lagrange(),
                  "Target fem not convenient for interpolation");
    }
    /* initialisation of the mesh_trans_inv */
    if (rg_target.id() == mesh_region::all_convexes().id()) {
      size_type nbpts = mf_target.nb_basic_dof() / qdim_t;
      for (size_type i = 0; i < nbpts; ++i)
        mti.add_point(mf_target.point_of_basic_dof(i * qdim_t));
      interpolation(mf_source, mti, U, V, M, version, extrapolation);
    }
    else {
      for (dal::bv_visitor_c dof(mf_target.basic_dof_on_region(rg_target)); !dof.finished(); ++dof)
        if (dof % qdim_t == 0)
          mti.add_point_with_id(mf_target.point_of_basic_dof(dof), dof/qdim_t);
      interpolation(mf_source, mti, U, V, M, version, extrapolation, 0, rg_source);
    }

    if (version == 0)
      mf_target.reduce_vector(V, VV);
    else {
      if (mf_target.is_reduced())
        gmm::mult(mf_target.reduction_matrix(), M, MM); 
      else
        gmm::copy(M, MM);
    }

  }

  /*
     interpolation of a solution on another torus mesh.
     - the solution should be continuous.
   */
  template<typename VECTU, typename VECTV, typename MAT>
  void interpolation_to_torus_mesh_fem(const mesh_fem &mf_source, const torus_mesh_fem &mf_target,
                                       const VECTU &U, VECTV &VV, MAT &MM,
                                       int version, int extrapolation,
                                       double EPS,
                                       mesh_region rg_source=mesh_region::all_convexes(),
                                       mesh_region rg_target=mesh_region::all_convexes()) {

    typedef typename gmm::linalg_traits<VECTU>::value_type T;
    dim_type qqdim = dim_type(gmm::vect_size(U)/mf_source.nb_dof());
    size_type qqdimt = qqdim * mf_source.get_qdim()/mf_target.get_qdim();
    std::vector<T> V(mf_target.nb_basic_dof()*qqdimt);
    mf_target.extend_vector(VV,V);
    gmm::row_matrix<gmm::rsvector<scalar_type> >
      M(mf_target.nb_basic_dof(), mf_source.nb_dof());

    const mesh &msh(mf_source.linked_mesh());
    getfem::mesh_trans_inv mti(msh, EPS);
    size_type qdim_s = mf_source.get_qdim(), qdim_t = mf_target.get_qdim();
    GMM_ASSERT1(qdim_s == qdim_t || qdim_t == 1,
                "Attempt to interpolate a field of dimension "
                << qdim_s << " on a mesh_fem whose Qdim is " << qdim_t);
    
    /* test if the target mesh_fem is convenient for interpolation.         */
    for (dal::bv_visitor cv(mf_target.convex_index()); !cv.finished();++cv) {
      pfem pf_t = mf_target.fem_of_element(cv);

      GMM_ASSERT1(pf_t->target_dim() == 1 || 
                  (mf_target.get_qdim() == mf_target.linked_mesh().dim()),
                  "Target fem not convenient for interpolation");
    }
    /* initialisation of the mesh_trans_inv */
    if (rg_target.id() == mesh_region::all_convexes().id()) {
      size_type nbpts = mf_target.nb_basic_dof() / qdim_t;
      for (size_type i = 0; i < nbpts; ++i)
      {
        base_node p(msh.dim());
        for (size_type k=0; k < msh.dim(); ++k) p[k] = mf_target.point_of_basic_dof(i * qdim_t)[k];
        mti.add_point(p);
      }
      interpolation(mf_source, mti, U, V, M, version, extrapolation);
    }
    else {
      for (dal::bv_visitor_c dof(mf_target.basic_dof_on_region(rg_target)); !dof.finished(); ++dof)
        if (dof % qdim_t == 0)
        {
           base_node p(msh.dim());
          for (size_type k=0; k < msh.dim(); ++k) p[k] = mf_target.point_of_basic_dof(dof)[k];          
          mti.add_point_with_id(p, dof/qdim_t);
        }
      interpolation(mf_source, mti, U, V, M, version, extrapolation, 0, rg_source);
    }

    if (version == 0)
      mf_target.reduce_vector(V, VV);
    else {
      if (mf_target.is_reduced())
        gmm::mult(mf_target.reduction_matrix(), M, MM); 
      else
        gmm::copy(M, MM);
    }
  }



  template<typename VECTU, typename VECTV>
  void interpolation(const mesh_fem &mf_source, const mesh_fem &mf_target,
                     const VECTU &U, VECTV &V, int extrapolation,
                     double EPS,
                     mesh_region rg_source, mesh_region rg_target) {
    base_matrix M;
    GMM_ASSERT1((gmm::vect_size(U) % mf_source.nb_dof()) == 0
                && (gmm::vect_size(V) % mf_target.nb_dof()) == 0
                && gmm::vect_size(V) != 0, "Dimensions mismatch");
    if (&mf_source.linked_mesh() == &mf_target.linked_mesh() &&
        rg_source.id() == mesh_region::all_convexes().id() &&
        rg_target.id() == mesh_region::all_convexes().id())
      interpolation_same_mesh(mf_source, mf_target, U, V, M, 0);
    else 
      interpolation(mf_source, mf_target, U, V, M, 0, extrapolation, EPS,
                    rg_source, rg_target);
  }

  template<typename MAT>
  void interpolation(const mesh_fem &mf_source, const mesh_fem &mf_target,
                     MAT &M, int extrapolation, double EPS,
                     mesh_region rg_source, mesh_region rg_target) {
    GMM_ASSERT1(mf_source.nb_dof() == gmm::mat_ncols(M)
                && (gmm::mat_nrows(M) % mf_target.nb_dof()) == 0
                && gmm::mat_nrows(M) != 0, "Dimensions mismatch");
    std::vector<scalar_type> U, V;
    if (&mf_source.linked_mesh() == &mf_target.linked_mesh() &&
        rg_source.id() == mesh_region::all_convexes().id() &&
        rg_target.id() ==mesh_region::all_convexes().id())
      interpolation_same_mesh(mf_source, mf_target, U, V, M, 1);
    else 
      interpolation(mf_source, mf_target, U, V, M, 1, extrapolation, EPS,
                    rg_source, rg_target);
  }


  /**Interpolate mesh_fem data to im_data.
   The qdim of mesh_fem must be equal to im_data nb_tensor_elem.
   Both im_data and mesh_fem must reside in the same mesh. 
   Only convexes defined with both mesh_fem and im_data will be interpolated.
   Version 0 uses filtered region of im_data, while other versions use full region.
  */
  template <typename VECT>
  void interpolation_to_im_data(const mesh_fem &mf_source,
    const im_data &im_target,
    const VECT &nodal_data, VECT &int_pt_data, bool use_im_data_filtered = true){
    // typedef typename gmm::linalg_traits<const VECT>::value_type T;
    
    dim_type qdim = mf_source.get_qdim();
    GMM_ASSERT1(qdim == im_target.nb_tensor_elem(), 
                "Incompatible size of qdim for mesh_fem " << qdim
                << " and im_data " << im_target.nb_tensor_elem());
    GMM_ASSERT1(&mf_source.linked_mesh() == &im_target.linked_mesh_im().linked_mesh(),
                "mf_source and im_data do not share the same mesh.");

    base_matrix G;
    base_vector coeff;
    bgeot::base_tensor tensor_int_point(im_target.tensor_size());    
    fem_precomp_pool fppool;

    size_type nb_dof = mf_source.nb_dof();
    size_type nb_basic_dof = mf_source.nb_basic_dof();

    GMM_ASSERT1(nb_dof == gmm::vect_size(nodal_data), 
                "Provided nodal data size is " << gmm::vect_size(nodal_data)
                << " but expecting vector size of " << nb_dof);
    
    size_type size_im_data = (use_im_data_filtered)
                              ?im_target.nb_filtered_index() * im_target.nb_tensor_elem()
                              :im_target.nb_index() * im_target.nb_tensor_elem();
    GMM_ASSERT1(size_im_data == gmm::vect_size(int_pt_data), 
                "Provided im data size is " << gmm::vect_size(int_pt_data)
                << " but expecting vector size of " << size_im_data);
 
    VECT extended_nodal_data(nb_basic_dof);
    if (nb_dof == nb_basic_dof) gmm::copy(nodal_data, extended_nodal_data);
    else mf_source.extend_vector(nodal_data, extended_nodal_data);
    
    dal::bit_vector im_data_convex_index;
    if(use_im_data_filtered) im_data_convex_index = im_target.filtered_convex_index();
    else im_data_convex_index = im_target.linked_mesh_im().convex_index();

    for (dal::bv_visitor cv(im_data_convex_index); !cv.finished(); ++cv) {

      bgeot::pgeometric_trans pgt = mf_source.linked_mesh().trans_of_convex(cv);
      pfem pf_source = mf_source.fem_of_element(cv);
      if(pf_source == NULL) continue;

      mesh_fem::ind_dof_ct::const_iterator it_dof;
      size_type cv_nb_dof = mf_source.nb_basic_dof_of_element(cv);
      size_type nb_nodal_pt = cv_nb_dof / qdim;
      coeff.resize(cv_nb_dof);      
      getfem::slice_vector_on_basic_dof_of_element(mf_source, extended_nodal_data, cv, coeff);

      getfem::pintegration_method pim = im_target.linked_mesh_im().int_method_of_element(cv);
      if (pf_source->need_G()) 
        bgeot::vectors_to_base_matrix(G, pim->approx_method()->integration_points());

      pfem_precomp pfp = fppool(pf_source, &pim->approx_method()->integration_points());
      fem_interpolation_context ctx(pgt, pfp, size_type(-1), G, cv, short_type(-1));
      size_type index_int_pt(-1);
      
      if(use_im_data_filtered) index_int_pt = im_target.filtered_index_of_point(cv, 0);
      else index_int_pt = im_target.index_of_point(cv, 0);

      size_type nb_int_points = im_target.nb_points_of_element(cv);

      for (size_type i = 0; i < nb_int_points; ++i, index_int_pt+=qdim) {
        ctx.set_ii(i);        
        ctx.pf()->interpolation(ctx, coeff, tensor_int_point.as_vector(), qdim);
        im_target.set_tensor(int_pt_data, cv, i, tensor_int_point, use_im_data_filtered);
      }    
    }//end of convex loop
  }

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_INTERPOLATION_H__  */
