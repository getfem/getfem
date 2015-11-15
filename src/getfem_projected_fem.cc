/*===========================================================================

 Copyright (C) 2012-2015 Yves Renard, Konstantinos Poulios

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

#include "getfem/getfem_projected_fem.h"

namespace getfem {

  typedef bgeot::convex<base_node>::dref_convex_pt_ct dref_convex_pt_ct;
//  typedef bgeot::basic_mesh::ref_mesh_face_pt_ct ref_mesh_face_pt_ct;

  /* calculates the projection of a point on the face of a convex
   * Input:
   *   pgt : the geometric transformation of the convex
   *   G_cv: the nodes of the convex, stored in columns
   *   fc  : the face of the convex to project on
   *   pt  : the point to be projected
   * Output:
   *   proj_ref: the projected point in the reference element
  */
  void projection_on_convex_face
    (const bgeot::pgeometric_trans pgt, const base_matrix &G_cv,
     const short_type fc, const base_node &pt,
     base_node &proj_ref) {

    size_type N = gmm::mat_nrows(G_cv); // dimension of the target space
    size_type P = pgt->dim();           // dimension of the reference element space

    size_type nb_pts_cv = gmm::mat_ncols(G_cv);
    size_type nb_pts_fc = pgt->structure()->nb_points_of_face(fc);

    GMM_ASSERT1( N == pt.size(), "Dimensions mismatch");
    GMM_ASSERT1( nb_pts_cv == pgt->nb_points(), "Dimensions mismatch");

    bgeot::convex_ind_ct ind_pts_fc = pgt->structure()->ind_points_of_face(fc);

    base_matrix G_fc(N, nb_pts_fc);
    for (size_type i=0; i < nb_pts_fc; i++)
      gmm::copy(gmm::mat_col(G_cv,ind_pts_fc[i]),gmm::mat_col(G_fc,i));

    // Local base on reference face
    base_matrix base_ref_fc(P,P-1);
    {
      dref_convex_pt_ct dref_pts_fc = pgt->convex_ref()->dir_points_of_face(fc);
      GMM_ASSERT1( dref_pts_fc.size() == P, "Dimensions mismatch");
      for (size_type i = 0; i < P-1; ++i)
          gmm::copy(dref_pts_fc[i+1] - dref_pts_fc[0], gmm::mat_col(base_ref_fc,i));
    }

    proj_ref.resize(P);

    base_node proj(N); // the projected point in the real element
    base_node vres(P); // residual vector
    scalar_type res= 1.;

    // initial guess
    proj_ref = gmm::mean_value(pgt->convex_ref()->points_of_face(fc));

    base_vector val(nb_pts_fc);
    pgt->poly_vector_val(proj_ref, ind_pts_fc, val);
    gmm::mult(G_fc, val, proj);

    base_matrix K(N,P-1);

    base_matrix grad_fc(nb_pts_fc, P);
    base_matrix grad_fc1(nb_pts_fc, P-1);
    base_matrix B(N,P-1), BB(N,P), CS(P-1,P-1);

    scalar_type EPS = 10E-12;
    unsigned cnt = 50;
    while (res > EPS && --cnt) {
      // computation of the pseudo inverse matrix B at point proj_ref
      pgt->poly_vector_grad(proj_ref, ind_pts_fc, grad_fc);
      gmm::mult(grad_fc, base_ref_fc, grad_fc1);
      gmm::mult(G_fc, grad_fc1, K);
      gmm::mult(gmm::transposed(K), K, CS);
      gmm::lu_inverse(CS);
      gmm::mult(K, CS, B);
      gmm::mult(B, gmm::transposed(base_ref_fc), BB);

      // Projection onto the face of the convex
      gmm::mult_add(gmm::transposed(BB), pt - proj, proj_ref);
      pgt->poly_vector_val(proj_ref, ind_pts_fc, val);
      gmm::mult(G_fc, val, proj);

      gmm::mult(gmm::transposed(BB), pt - proj, vres);
      res = gmm::vect_norm2(vres);
    }
    GMM_ASSERT1( res <= EPS,
                "Iterative pojection on convex face did not converge");
  }


  /* calculates the normal at a specific point on the face of a convex
   * Input:
   *   pgt   : the geometric transformation of the convex
   *   G_cv  : the nodes of the convex, stored in columns
   *   fc    : the face of the convex to project on
   *   ref_pt: the point in the reference element
   * Output:
   *   normal: the surface normal in the real element corresponding at
   *           the location of ref_pt in the reference element
  */
  void normal_on_convex_face
    (const bgeot::pgeometric_trans pgt, const base_matrix &G_cv,
     const short_type fc, const base_node &ref_pt, base_node &normal) {

    size_type N = gmm::mat_nrows(G_cv); // dimension of the target space
    size_type P = pgt->dim();           // dimension of the reference element space

    size_type nb_pts_cv = gmm::mat_ncols(G_cv);
    size_type nb_pts_fc = pgt->structure()->nb_points_of_face(fc);

    GMM_ASSERT1( nb_pts_cv == pgt->nb_points(), "Dimensions mismatch");

    bgeot::convex_ind_ct ind_pts_fc = pgt->structure()->ind_points_of_face(fc);

    base_matrix G_fc(N, nb_pts_fc);
    for (size_type i=0; i < nb_pts_fc; i++)
      gmm::copy(gmm::mat_col(G_cv,ind_pts_fc[i]),gmm::mat_col(G_fc,i));

    // Local base on reference face
    base_matrix base_ref_fc(P,P-1);
    {
      dref_convex_pt_ct dref_pts_fc = pgt->convex_ref()->dir_points_of_face(fc);
      GMM_ASSERT1( dref_pts_fc.size() == P, "Dimensions mismatch");
      for (size_type i = 0; i < P-1; ++i)
          gmm::copy(dref_pts_fc[i+1] - dref_pts_fc[0], gmm::mat_col(base_ref_fc,i));
    }

    normal.resize(N);

    base_matrix K(N,P-1);
    { // calculate K at the final point
      base_matrix grad_fc(nb_pts_fc, P);
      base_matrix grad_fc1(nb_pts_fc, P-1);
      pgt->poly_vector_grad(ref_pt, ind_pts_fc, grad_fc);
      gmm::mult(grad_fc, base_ref_fc, grad_fc1);
      gmm::mult(G_fc, grad_fc1, K);
    }

    base_matrix KK(N,P);
    { // calculate KK
      base_matrix grad_cv(nb_pts_cv, P);
      pgt->poly_vector_grad(ref_pt, grad_cv);
      gmm::mult(G_cv, grad_cv, KK);
    }

    base_matrix bases_product(P-1, P);
    gmm::mult(gmm::transposed(K), KK, bases_product);

    for (size_type i = 0; i < P; ++i) {
      std::vector<size_type> ind(0);
      for (size_type j = 0; j < P; ++j)
        if (j != i ) ind.push_back(j);
      scalar_type det = gmm::lu_det(gmm::sub_matrix(bases_product,
                                                    gmm::sub_interval(0, P-1),
                                                    gmm::sub_index(ind)      ) );
      gmm::add(gmm::scaled(gmm::mat_col(KK, i), (i % 2) ? -det : +det ), normal);
    }

    // normalizing
    gmm::scale(normal, 1/gmm::vect_norm2(normal));

    // ensure that normal points outwards
    base_node cv_center(N), fc_center(N);
    for (size_type i=0; i < nb_pts_cv; i++)
      gmm::add(gmm::mat_col(G_cv,i), cv_center);
    for (size_type i=0; i < nb_pts_fc; i++)
      gmm::add(gmm::mat_col(G_fc,i), fc_center);
    gmm::scale(cv_center, scalar_type(1)/scalar_type(nb_pts_cv));
    gmm::scale(fc_center, scalar_type(1)/scalar_type(nb_pts_fc));
    if (gmm::vect_sp(normal, fc_center -cv_center) < 0)
      gmm::scale(normal, scalar_type(-1));
  }

  /* calculates the normal at a specific point of a convex in a higher
   * dimension space
   * Input:
   *   pgt   : the geometric transformation of the convex
   *   G_cv  : the nodes of the convex, stored in columns
   *   ref_pt: the point in the reference element
   * Output:
   *   normal: the surface normal in the real element corresponding at
   *           the location of ref_pt in the reference element
   *           (or one of the possible normals if the space dimension
   *           is more than one higher than the convex dimension)
  */
  void normal_on_convex
    (const bgeot::pgeometric_trans pgt, const base_matrix &G_cv,
     const base_node &ref_pt, base_node &normal) {

    size_type N = gmm::mat_nrows(G_cv); // dimension of the target space
    size_type P = pgt->dim();           // dimension of the reference element space

    GMM_ASSERT1( N == 2 || N == 3, "Normal on convexes calculation is supported "
                                   "only for space dimension equal to 2 or 3.");
    GMM_ASSERT1( P < N, "Normal on convex is defined only in a space of"
                        "higher dimension.");

    size_type nb_pts = gmm::mat_ncols(G_cv);
    base_matrix K(N,P);
    { // calculate K at the final point
      base_matrix grad_cv(nb_pts, P);
      pgt->poly_vector_grad(ref_pt, grad_cv);
      gmm::mult(G_cv, grad_cv, K);
    }

    gmm::resize(normal,N);
    if (P==1 && N == 2) {
      normal[0] = -K(1,0);
      normal[1] = K(0,0);
    }
    else if (P==1 && N == 3) {
      normal[0] = K(2,0)-K(1,0);
      normal[1] = K(0,0)-K(2,0);
      normal[2] = K(1,0)-K(0,0);
    }
    else if (P==2) {
      normal[0] = K(1,0)*K(2,1)-K(2,0)*K(1,1);
      normal[1] = K(2,0)*K(0,1)-K(0,0)*K(2,1);
      normal[2] = K(0,0)*K(1,1)-K(1,0)*K(0,1);
    }
    gmm::scale(normal, 1/gmm::vect_norm2(normal));
  }

  void projected_fem::build_kdtree(void) const {
    tree.clear();
    dal::bit_vector dofs=mf_source.basic_dof_on_region(rg_source);
    dofs.setminus(blocked_dofs);
    dim_type qdim=target_dim();
    for (dal::bv_visitor dof(dofs); !dof.finished(); ++dof)
        if (dof % qdim == 0)
            tree.add_point_with_id(mf_source.point_of_basic_dof(dof), dof);
  }

  bool projected_fem::find_a_projected_point(base_node pt, base_node &ptr_proj,
                                             size_type &cv_proj, short_type &fc_proj) const {

    bgeot::index_node_pair ipt;
    //scalar_type dist =
    tree.nearest_neighbor(ipt, pt);

    size_type cv_sel = size_type(-1);
    short_type fc_sel = short_type(-1);
    scalar_type is_in_sel(1e10);
    base_node proj_ref, proj_ref_sel;
    const getfem::mesh::ind_cv_ct cvs = mf_source.convex_to_basic_dof(ipt.i);
    for (size_type i=0; i < cvs.size(); ++i) {
      size_type cv = cvs[i];
      const bgeot::pgeometric_trans pgt = mf_source.linked_mesh().trans_of_convex(cv);
      if (rg_source.is_in(cv)) { // project on the convex
        bool gt_invertible;
        gic = bgeot::geotrans_inv_convex(mf_source.linked_mesh().convex(cv), pgt);
        gic.invert(pt, proj_ref, gt_invertible);
        if (gt_invertible) {
          scalar_type is_in = pgt->convex_ref()->is_in(proj_ref);
          if (is_in < is_in_sel) {
            is_in_sel = is_in;
            cv_sel = cv;
            fc_sel = short_type(-1);
            proj_ref_sel = proj_ref;
          }
        }
      }
      else { // project on convex faces
        mesh_region::face_bitset faces = rg_source.faces_of_convex(cv);
        if (faces.count() > 0) { // this should rarely be more than one face
          bgeot::vectors_to_base_matrix(G, mf_source.linked_mesh().points_of_convex(cv));
          short_type nbf = mf_source.linked_mesh().nb_faces_of_convex(cv);
          for (short_type f = 0; f < nbf; ++f) {
            if (faces.test(f)) {
              projection_on_convex_face(pgt, G, f, pt, proj_ref);
              scalar_type is_in = pgt->convex_ref()->is_in(proj_ref);
              if (is_in < is_in_sel) {
                is_in_sel = is_in;
                cv_sel = cv;
                fc_sel = f;
                proj_ref_sel = proj_ref;
              }
            }
          }
        }
      }
    }
    if (cv_sel != size_type(-1) && is_in_sel < 0.05) {  //FIXME
        cv_proj = cv_sel;
        fc_proj = fc_sel;
        ptr_proj = proj_ref_sel;
        return true;
    }
    return false;
  }

  void projected_fem::update_from_context(void) const {
    fictx_cv = size_type(-1);
    dim_ = dim_type(-1);

    dim_type N = mf_source.linked_mesh().dim();
    GMM_ASSERT1( N == mim_target.linked_mesh().dim(),
                 "Dimensions mismatch between the source and the target meshes");

    build_kdtree();

    elements.clear();
    ind_dof.resize(mf_source.nb_basic_dof());
    size_type max_dof = 0;
    if (rg_target.id() != mesh_region::all_convexes().id() &&
        rg_target.is_empty()) {
      dim_ = mim_target.linked_mesh().dim();
      return;
    }

    for (mr_visitor i(rg_target); !i.finished(); ++i) {
      size_type cv = i.cv(); // refers to the target mesh
      short_type f = i.f();  // refers to the target mesh

      dim_type dim__ = mim_target.linked_mesh().structure_of_convex(cv)->dim();
      if (dim_ == dim_type(-1)) {
        dim_ = dim__;
        if (i.is_face()) dim__ = dim_type(dim__ - 1);
        GMM_ASSERT1(dim__ < N, "The projection should take place in lower "
                              "dimensions than the mesh dimension. Otherwise "
                              "use the interpolated_fem object instead.");
      }
      else
        GMM_ASSERT1(dim_ == dim__,
                    "Convexes/faces of different dimension in the target mesh");

      pintegration_method pim = mim_target.int_method_of_element(cv);
      GMM_ASSERT1(pim->type() == IM_APPROX,
                  "You have to use approximated integration to project a fem");
      papprox_integration pai = pim->approx_method();
      bgeot::pgeometric_trans pgt = mim_target.linked_mesh().trans_of_convex(cv);
      bgeot::pgeotrans_precomp pgp =
        bgeot::geotrans_precomp(pgt, pai->pintegration_points(), 0);
      dal::bit_vector dofs;
      size_type last_cv = size_type(-1); // refers to the source mesh
      short_type last_f = short_type(-1); // refers to the source mesh
      size_type nb_pts = i.is_face() ? pai->nb_points_on_face(f) : pai->nb_points();
      size_type start_pt = i.is_face() ? pai->ind_first_point_on_face(f) : 0;
      elt_projection_data &e = elements[cv];
      base_node gpt(N);
      for (size_type k = 0; k < nb_pts; ++k) {
        pgp->transform(mim_target.linked_mesh().points_of_convex(cv),
                       start_pt + k, gpt);
        gausspt_projection_data &gppd = e.gausspt[start_pt + k];
        gppd.iflags = find_a_projected_point(gpt, gppd.ptref, gppd.cv, gppd.f) ? 1 : 0;
        if (gppd.iflags) {
          // calculate gppd.normal
          const bgeot::pgeometric_trans pgt_source = mf_source.linked_mesh().trans_of_convex(gppd.cv);
          bgeot::vectors_to_base_matrix(G, mf_source.linked_mesh().points_of_convex(gppd.cv));
          if (gppd.f != short_type(-1))
            normal_on_convex_face(pgt_source, G, gppd.f, gppd.ptref, gppd.normal);
          else
            normal_on_convex(pgt_source, G, gppd.ptref, gppd.normal);
          // calculate gppd.gap
          base_node ppt = pgt_source->transform(gppd.ptref, G);
          gppd.gap = gmm::vect_sp(gpt-ppt, gppd.normal);
        }

        if (gppd.iflags && (last_cv != gppd.cv || last_f != gppd.f)) {
          if (gppd.f == short_type(-1)) { // convex
            size_type nbdof = mf_source.nb_basic_dof_of_element(gppd.cv);
            for (size_type loc_dof = 0; loc_dof < nbdof; ++loc_dof) {
              size_type idof = mf_source.ind_basic_dof_of_element(gppd.cv)[loc_dof];
              if (!(blocked_dofs[idof])) dofs.add(idof);
            }
          }
          else { // convex face
            size_type nbdof = mf_source.nb_basic_dof_of_face_of_element(gppd.cv, gppd.f);
            for (size_type loc_dof = 0; loc_dof < nbdof; ++loc_dof) {
              size_type idof = mf_source.ind_basic_dof_of_face_of_element(gppd.cv, gppd.f)[loc_dof];
              if (!(blocked_dofs[idof])) dofs.add(idof);
            }
          }
          last_cv = gppd.cv;
          last_f = gppd.f;
        }
      }
      e.nb_dof = dofs.card();
      e.pim = pim;
      e.inddof.resize(dofs.card());
      max_dof = std::max(max_dof, dofs.card());
      size_type cnt = 0;
      for (dal::bv_visitor idof(dofs); !idof.finished(); ++idof)
        { e.inddof[cnt] = idof; ind_dof[idof] = cnt++; }
      for (size_type k = 0; k < nb_pts; ++k) {
        gausspt_projection_data &gppd = e.gausspt[start_pt + k];
        if (gppd.iflags) {
          if (gppd.f == short_type(-1)) { // convex
            size_type nbdof = mf_source.nb_basic_dof_of_element(gppd.cv);
            for (size_type loc_dof = 0; loc_dof < nbdof; ++loc_dof) {
              size_type idof = mf_source.ind_basic_dof_of_element(gppd.cv)[loc_dof];
              gppd.local_dof[loc_dof] = dofs.is_in(idof) ? ind_dof[idof]
                                                         : size_type(-1);
            }
          }
          else { // convex face
            size_type nbdof = mf_source.nb_basic_dof_of_face_of_element(gppd.cv, gppd.f);
            pfem pf = mf_source.fem_of_element(gppd.cv);
            bgeot::convex_ind_ct ind_pts_fc = pf->structure(gppd.cv)->ind_points_of_face(gppd.f);
            unsigned rdim = target_dim() / pf->target_dim();
            if (rdim == 1)
              for (size_type loc_dof = 0; loc_dof < nbdof; ++loc_dof) { // local dof with respect to the source convex face
                size_type idof = mf_source.ind_basic_dof_of_face_of_element(gppd.cv, gppd.f)[loc_dof];
                size_type loc_dof2 = ind_pts_fc[loc_dof]; // local dof with respect to the source convex
                gppd.local_dof[loc_dof2] = dofs.is_in(idof) ? ind_dof[idof]
                                                            : size_type(-1);
              }
            else
              for (size_type ii = 0; ii < nbdof/rdim; ++ii)
                for (size_type jj = 0; jj < rdim; ++jj) {
                  size_type loc_dof = ii*rdim + jj; // local dof with respect to the source convex face
                  size_type idof = mf_source.ind_basic_dof_of_face_of_element(gppd.cv, gppd.f)[loc_dof];
                  size_type loc_dof2 = ind_pts_fc[ii]*rdim + jj; // local dof with respect to the source convex
                  gppd.local_dof[loc_dof2] = dofs.is_in(idof) ? ind_dof[idof]
                                                              : size_type(-1);
                }
          }
        }
      }
    }
    /** setup global dofs, with dummy coordinates */
    base_node P(dim()); gmm::fill(P,1./20);
    node_tab_.resize(max_dof);
    std::fill(node_tab_.begin(), node_tab_.end(), P);
    pspt_valid = false;
    dof_types_.resize(max_dof);
    std::fill(dof_types_.begin(), dof_types_.end(),
              global_dof(dim()));

    /* ind_dof should be kept full of -1 ( real_base_value and
       grad_base_value expect that)
    */
    std::fill(ind_dof.begin(), ind_dof.end(), size_type(-1));
  }

  size_type projected_fem::nb_dof(size_type cv) const
  {
    context_check();
    GMM_ASSERT1(mim_target.linked_mesh().convex_index().is_in(cv),
                "Wrong convex number: " << cv);
    std::map<size_type,elt_projection_data>::const_iterator eit;
    eit = elements.find(cv);
    return (eit != elements.end()) ? eit->second.nb_dof : 0;
  }

  size_type projected_fem::index_of_global_dof(size_type cv, size_type i) const
  {
    std::map<size_type,elt_projection_data>::const_iterator eit;
    eit = elements.find(cv);
    GMM_ASSERT1(eit != elements.end(), "Wrong convex number: " << cv);
    return eit->second.inddof[i];
  }

  bgeot::pconvex_ref projected_fem::ref_convex(size_type cv) const
  { return mim_target.int_method_of_element(cv)->approx_method()->ref_convex(); }

  const bgeot::convex<base_node> &projected_fem::node_convex(size_type cv) const
  {
    GMM_ASSERT1(mim_target.linked_mesh().convex_index().is_in(cv),
                "Wrong convex number: " << cv);
    return *(bgeot::generic_dummy_convex_ref
             (dim(), nb_dof(cv),
              mim_target.linked_mesh().structure_of_convex(cv)->nb_faces()));
  }

  bgeot::pstored_point_tab projected_fem::node_tab(size_type)
    const {
    if (!pspt_valid)
      { pspt = bgeot::store_point_tab(node_tab_); pspt_valid = true; }
    return pspt;
  }

  void projected_fem::base_value(const base_node &, base_tensor &) const
  { GMM_ASSERT1(false, "No base values, real only element."); }
  void projected_fem::grad_base_value(const base_node &,
                                      base_tensor &) const
  { GMM_ASSERT1(false, "No grad values, real only element."); }
  void projected_fem::hess_base_value(const base_node &,
                                      base_tensor &) const
  { GMM_ASSERT1(false, "No hess values, real only element."); }

  inline void projected_fem::actualize_fictx(pfem pf, size_type cv,
                                             const base_node &ptr) const {
    if (fictx_cv != cv) {
      bgeot::vectors_to_base_matrix
        (G, mf_source.linked_mesh().points_of_convex(cv));
      fictx = fem_interpolation_context
        (mf_source.linked_mesh().trans_of_convex(cv), pf, base_node(), G, cv);
      fictx_cv = cv;
    }
    fictx.set_xref(ptr);
  }

  void projected_fem::real_base_value(const fem_interpolation_context& c,
                                      base_tensor &t, bool) const {
    std::map<size_type,elt_projection_data>::iterator eit;
    eit = elements.find(c.convex_num());
    if (eit == elements.end()) {
      mi2[1] = target_dim(); mi2[0] = short_type(0);
      t.adjust_sizes(mi2);
      std::fill(t.begin(), t.end(), scalar_type(0));
      return;
    }
//    GMM_ASSERT1(eit != elements.end(), "Wrong convex number: " << c.convex_num());
    elt_projection_data &e = eit->second;

    mi2[1] = target_dim(); mi2[0] = short_type(e.nb_dof);
    t.adjust_sizes(mi2);
    std::fill(t.begin(), t.end(), scalar_type(0));
    if (e.nb_dof == 0) return;

    std::map<size_type,gausspt_projection_data>::iterator git;
    git = e.gausspt.find(c.ii());
    if (c.have_pgp() &&
        (c.pgp()->get_ppoint_tab()
         == e.pim->approx_method()->pintegration_points()) &&
        git != e.gausspt.end()) {
      gausspt_projection_data &gppd = git->second;
      if (gppd.iflags & 1) {
        if (gppd.iflags & 2) {
          t = gppd.base_val;
          return;
        }
        size_type cv = gppd.cv;
        pfem pf = mf_source.fem_of_element(cv);
        actualize_fictx(pf, cv, gppd.ptref);
        pf->real_base_value(fictx, taux);
        unsigned rdim = target_dim() / pf->target_dim();
        std::map<size_type,size_type>::const_iterator ii;
        if (rdim == 1) // mdim == 0
          for (size_type i = 0; i < pf->nb_dof(cv); ++i) {
            ii = gppd.local_dof.find(i);
            if (ii != gppd.local_dof.end() && ii->second != size_type(-1))
              for (size_type j = 0; j < target_dim(); ++j)
                t(ii->second,j) = taux(i,j);
          }
        else // mdim == 1
          for (size_type i = 0; i < pf->nb_dof(cv); ++i)
            for (size_type j = 0; j < target_dim(); ++j) {
              ii = gppd.local_dof.find(i*rdim+j);
              if (ii != gppd.local_dof.end() && ii->second != size_type(-1))
                t(ii->second,j) = taux(i,0);
            }

        if (store_values) {
          gppd.base_val = t;
          gppd.iflags |= 2;
        }
      }
    }
    else {
      size_type cv;
      short_type f;
      if (find_a_projected_point(c.xreal(), ptref, cv, f)) {
        pfem pf = mf_source.fem_of_element(cv);
        actualize_fictx(pf, cv, ptref);
        pf->real_base_value(fictx, taux);

        for (size_type i = 0; i < e.nb_dof; ++i)
          ind_dof.at(e.inddof[i]) = i;

        unsigned rdim = target_dim() / pf->target_dim();
        if (rdim == 1) // mdim == 0
          for (size_type i = 0; i < pf->nb_dof(cv); ++i) {
            size_type ii = ind_dof.at(mf_source.ind_basic_dof_of_element(cv)[i]);
            if (ii != size_type(-1)) {
              for (size_type j = 0; j < target_dim(); ++j)
                t(ii,j) = taux(i,j);
            }
          }
        else // mdim == 1
          for (size_type i = 0; i < pf->nb_dof(cv); ++i)
            for (size_type j = 0; j < target_dim(); ++j) {
              size_type ij = ind_dof.at(mf_source.ind_basic_dof_of_element(cv)[i*rdim+j]);
              if (ij != size_type(-1))
                t(ij,j) = taux(i,0);
            }

        for (size_type i = 0; i < e.nb_dof; ++i)
          ind_dof[e.inddof[i]] = size_type(-1);
      }
    }

  }

  void projected_fem::real_grad_base_value(const fem_interpolation_context& c,
                                           base_tensor &t, bool) const {
    std::map<size_type,elt_projection_data>::iterator eit;
    eit = elements.find(c.convex_num());
    if (eit == elements.end()) {
    mi3[2] = mf_source.linked_mesh().dim(); mi3[1] = target_dim(); mi3[0] = short_type(0);
      t.adjust_sizes(mi2);
      std::fill(t.begin(), t.end(), scalar_type(0));
      return;
    }
//    GMM_ASSERT1(eit != elements.end(), "Wrong convex number: " << c.convex_num());
    elt_projection_data &e = eit->second;

    size_type N = mf_source.linked_mesh().dim();
    mi3[2] = short_type(N); mi3[1] = target_dim(); mi3[0] = short_type(e.nb_dof);
    t.adjust_sizes(mi3);
    std::fill(t.begin(), t.end(), scalar_type(0));
    if (e.nb_dof == 0) return;

    std::map<size_type,gausspt_projection_data>::iterator git;
    git = e.gausspt.find(c.ii());
    if (c.have_pgp() &&
        (c.pgp()->get_ppoint_tab()
         == e.pim->approx_method()->pintegration_points()) &&
        git != e.gausspt.end()) {
      gausspt_projection_data &gppd = git->second;
      if (gppd.iflags & 1) {
        if (gppd.iflags & 4) {
          t = gppd.grad_val;
          return;
        }
        size_type cv = gppd.cv;
        pfem pf = mf_source.fem_of_element(cv);
        actualize_fictx(pf, cv, gppd.ptref);
        pf->real_grad_base_value(fictx, taux);

        unsigned rdim = target_dim() / pf->target_dim();
        std::map<size_type,size_type>::const_iterator ii;
        if (rdim == 1) // mdim == 0
          for (size_type i = 0; i < pf->nb_dof(cv); ++i) {
            ii = gppd.local_dof.find(i);
            if (ii != gppd.local_dof.end() && ii->second != size_type(-1))
              for (size_type j = 0; j < target_dim(); ++j)
                for (size_type k = 0; k < N; ++k)
                  t(ii->second, j, k) = taux(i, j, k);
          }
        else // mdim == 1
          for (size_type i = 0; i < pf->nb_dof(cv); ++i)
            for (size_type j = 0; j < target_dim(); ++j) {
              ii = gppd.local_dof.find(i*rdim+j);
              if (ii != gppd.local_dof.end() && ii->second != size_type(-1))
                for (size_type k = 0; k < N; ++k)
                  t(ii->second, j, k) = taux(i, 0, k);
            }
        if (store_values) {
          gppd.grad_val = t;
          gppd.iflags |= 4;
        }
      }
    }
    else {
      size_type cv;
      short_type f;
      if (find_a_projected_point(c.xreal(), ptref, cv, f)) {
        pfem pf = mf_source.fem_of_element(cv);
        actualize_fictx(pf, cv, ptref);
        pf->real_grad_base_value(fictx, taux);
        for (size_type i = 0; i < e.nb_dof; ++i)
          ind_dof.at(e.inddof[i]) = i;

        unsigned rdim = target_dim() / pf->target_dim();
        if (rdim == 1) // mdim == 0
          for (size_type i = 0; i < pf->nb_dof(cv); ++i) {
            size_type ii = ind_dof.at(mf_source.ind_basic_dof_of_element(cv)[i]);
            if (ii != size_type(-1))
              for (size_type j = 0; j < target_dim(); ++j)
                for (size_type k = 0; k < N; ++k)
                  t(ii,j,k) = taux(i,j,k);
          }
        else // mdim == 1
          for (size_type i = 0; i < pf->nb_dof(cv); ++i)
            for (size_type j = 0; j < target_dim(); ++j) {
              size_type ij = ind_dof.at(mf_source.ind_basic_dof_of_element(cv)[i*rdim+j]);
              if (ij != size_type(-1))
                for (size_type k = 0; k < N; ++k)
                  t(ij,j,k) = taux(i,0,k);
            }

        for (size_type i = 0; i < e.nb_dof; ++i)
          ind_dof[e.inddof[i]] = size_type(-1);
      }
    }
  }

  void projected_fem::real_hess_base_value
  (const fem_interpolation_context&, base_tensor &, bool) const
  { GMM_ASSERT1(false, "Sorry, to be done."); }

  void projected_fem::projection_data(const fem_interpolation_context& c,
                                      base_node &normal, scalar_type &gap) const {
    std::map<size_type,elt_projection_data>::iterator eit;
    eit = elements.find(c.convex_num());

    if (eit != elements.end()) {
      elt_projection_data &e = eit->second;
      if (e.nb_dof == 0) { // return undefined normal vector and huge gap
        normal = base_node(c.N());
        gap = 1e12;
        return;
      }
      std::map<size_type,gausspt_projection_data>::iterator git;
      git = e.gausspt.find(c.ii());
      if (c.have_pgp() &&
          (c.pgp()->get_ppoint_tab()
           == e.pim->approx_method()->pintegration_points()) &&
          git != e.gausspt.end()) {
        gausspt_projection_data &gppd = git->second;
        if (gppd.iflags & 1) {
          normal = gppd.normal;
          gap = gppd.gap;
        }
        else { // return undefined normal vector and huge gap
          normal = base_node(c.N());
          gap = 1e12;
        }
        return;
      }
    }

    // new projection
    projection_data(c.xreal(), normal, gap);
  }

  void projected_fem::projection_data(const base_node& pt,
                                      base_node &normal, scalar_type &gap) const {
    size_type cv;
    short_type f;
    if (find_a_projected_point(pt, ptref, cv, f)) {
      const bgeot::pgeometric_trans pgt = mf_source.linked_mesh().trans_of_convex(cv);
      bgeot::vectors_to_base_matrix(G, mf_source.linked_mesh().points_of_convex(cv));
      if (f != short_type(-1))
        normal_on_convex_face(pgt, G, f, ptref, normal);
      else
        normal_on_convex(pgt, G, ptref, normal);
      base_node ppt = pgt->transform(ptref, G);
      gap = gmm::vect_sp(pt-ppt, normal);
    }
    else { // return undefined normal vector and huge gap
      normal = base_node(pt.size());
      gap = 1e12;
    }

  }

  dal::bit_vector projected_fem::projected_convexes() const {
    dal::bit_vector bv;
    std::map<size_type,elt_projection_data>::const_iterator eit;
    for (eit = elements.begin(); eit != elements.end(); ++eit) {
      std::map<size_type,gausspt_projection_data>::const_iterator git;
      for (git = eit->second.gausspt.begin(); git != eit->second.gausspt.end(); ++git) {
        if (git->second.iflags)
          bv.add(git->second.cv);
      }
    }
    return bv;
  }

  void projected_fem::gauss_pts_stats(unsigned &ming, unsigned &maxg,
                                      scalar_type &meang) const {
    std::vector<unsigned> v(mf_source.linked_mesh().convex_index().last_true()+1);
    std::map<size_type,elt_projection_data>::const_iterator eit;
    for (eit = elements.begin(); eit != elements.end(); ++eit) {
      std::map<size_type,gausspt_projection_data>::const_iterator git;
      for (git = eit->second.gausspt.begin(); git != eit->second.gausspt.end(); ++git) {
        if (git->second.iflags)
          v[git->second.cv]++;
      }
    }

    ming = 100000; maxg = 0; meang = 0;
    unsigned cntg = 0;
    for (dal::bv_visitor cv(mf_source.linked_mesh().convex_index());
         !cv.finished(); ++cv) {
      ming = std::min(ming, v[cv]);
      maxg = std::max(maxg, v[cv]);
      meang += v[cv];
      if (v[cv] > 0) ++cntg;
    }
    meang /= scalar_type(cntg);
  }

  size_type projected_fem::memsize() const {
    size_type sz = 0;
    sz += blocked_dofs.memsize();
    sz += sizeof(*this);
    sz += elements.size() * sizeof(elt_projection_data); // Wrong for std::map
    std::map<size_type,elt_projection_data>::const_iterator eit;
    for (eit = elements.begin(); eit != elements.end(); ++eit) {
      sz += eit->second.gausspt.size() * sizeof(gausspt_projection_data); // Wrong for std::map
      sz += eit->second.inddof.capacity() * sizeof(size_type);
      std::map<size_type,gausspt_projection_data>::const_iterator git;
      for (git = eit->second.gausspt.begin(); git != eit->second.gausspt.end(); ++git) {
        sz += git->second.local_dof.size() * sizeof(size_type); // Wrong for std::map
      }
    }
    return sz;
  }

  projected_fem::projected_fem(const mesh_fem &mf_source_,
                               const mesh_im &mim_target_,
                               size_type rg_source_,
                               size_type rg_target_,
                               dal::bit_vector blocked_dofs_, bool store_val)
    : mf_source(mf_source_), mim_target(mim_target_),
      rg_source(mf_source.linked_mesh().region(rg_source_)),
      rg_target(mim_target.linked_mesh().region(rg_target_)),
      store_values(store_val), blocked_dofs(blocked_dofs_), mi2(2), mi3(3) {
    this->add_dependency(mf_source);
    this->add_dependency(mim_target);
    is_pol = is_lag = false; es_degree = 5;
    is_equiv = real_element_defined = true;
    ntarget_dim = mf_source.get_qdim();

    update_from_context();
  }

  DAL_SIMPLE_KEY(special_projfem_key, pfem);

  pfem new_projected_fem(const mesh_fem &mf_source_, const mesh_im &mim_target_,
                         size_type rg_source_, size_type rg_target_,
                         dal::bit_vector blocked_dofs_, bool store_val) {
    pfem pf(new projected_fem(mf_source_, mim_target_, rg_source_, rg_target_,
			      blocked_dofs_, store_val));
    dal::pstatic_stored_object_key
      pk = std::make_shared<special_projfem_key>(pf);
    dal::add_stored_object(pk, pf);
    return pf;
  }


}  /* end of namespace getfem.                                            */

