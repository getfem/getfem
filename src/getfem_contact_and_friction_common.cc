/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2013-2013 Yves Renard, Konstantinos Poulios.

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

#include "getfem/getfem_contact_and_friction_common.h"
#include "getfem/getfem_generic_assembly.h"
#ifndef _WIN32
#include <unistd.h>
#endif

namespace getfem {

  bool boundary_has_fem_nodes(bool slave_flag, int nodes_mode) {
    return (slave_flag && nodes_mode) ||
           (!slave_flag && nodes_mode == 2);
  }

  void compute_normal(const fem_interpolation_context &ctx,
                      size_type face, bool in_reference_conf,
                      const model_real_plain_vector &coeff,
                      base_node &n0, base_node &n,
                      base_matrix &grad) {
      n0 = bgeot::compute_normal(ctx, face);
      if (in_reference_conf) {
        n = n0;
      } else {
        ctx.pf()->interpolation_grad(ctx, coeff, grad, dim_type(ctx.N()));
        gmm::add(gmm::identity_matrix(), grad);
        scalar_type J = gmm::lu_inverse(grad);
        if (J <= scalar_type(0)) GMM_WARNING1("Inverted element !" << J);
        gmm::mult(gmm::transposed(grad), n0, n);
        gmm::scale(n, gmm::sgn(J)); // Test
      }
  }

  void vectorize_base_tensor(const base_tensor &t, base_matrix &vt,
                             size_type ndof, size_type qdim, size_type N) {
    GMM_ASSERT1(qdim == N || qdim == 1, "mixed intrinsic vector and "
                "tensorised fem is not supported");
    gmm::resize(vt, ndof, N);
    ndof = (ndof*qdim)/N;
    if (qdim == 1) {
      gmm::clear(vt);
      base_tensor::const_iterator it = t.begin();
      for (size_type i = 0; i < ndof; ++i, ++it)
        for (size_type j = 0; j < N; ++j) vt(i*N+j, j) = *it;
    } else if (qdim == N) {
      gmm::copy(t.as_vector(), vt.as_vector());
    }
  }

  void vectorize_grad_base_tensor(const base_tensor &t, base_tensor &vt,
                                         size_type ndof, size_type qdim,
                                         size_type N) {
    GMM_ASSERT1(qdim == N || qdim == 1, "mixed intrinsic vector and "
                  "tensorised fem is not supported");
    vt.adjust_sizes(bgeot::multi_index(ndof, N, N));
    ndof = (ndof*qdim)/N;
    if (qdim == 1) {
      gmm::clear(vt.as_vector());
      base_tensor::const_iterator it = t.begin();
      for (size_type k = 0; k < N; ++k)
        for (size_type i = 0; i < ndof; ++i, ++it)
          for (size_type j = 0; j < N; ++j) vt(i*N+j, j, k) = *it;
    } else if (qdim == N) {
      gmm::copy(t.as_vector(), vt.as_vector());
    }
  }

  //=========================================================================
  //
  //  Structure which store the contact boundaries, rigid obstacles and
  //  computes the contact pairs in large sliding/large deformation.
  //
  //=========================================================================

  size_type multi_contact_frame::add_U
  (const model_real_plain_vector *U, const std::string &name,
   const model_real_plain_vector *w, const std::string &wname) {
    if (!U) return size_type(-1);
    size_type i = 0;
    for (; i < Us.size(); ++i) if (Us[i] == U) return i;
    Us.push_back(U);
    Ws.push_back(w);
    Unames.push_back(name);
    Wnames.push_back(wname);
    ext_Us.resize(Us.size());
    ext_Ws.resize(Us.size());
    return i;
  }

  size_type multi_contact_frame::add_lambda
  (const model_real_plain_vector *lambda, const std::string &name) {
    if (!lambda) return size_type(-1);
    size_type i = 0;
    for (; i < lambdas.size(); ++i) if (lambdas[i] == lambda) return i;
    lambdas.push_back(lambda);
    lambdanames.push_back(name);
    ext_lambdas.resize(lambdas.size());
    return i;
  }

  void multi_contact_frame::extend_vectors(void) {
    dal::bit_vector iU, ilambda;
    for (size_type i = 0; i < contact_boundaries.size(); ++i) {
      size_type ind_U = contact_boundaries[i].ind_U;
      if (!(iU[ind_U])) {
        const mesh_fem &mf = *(contact_boundaries[i].mfu);
        gmm::resize(ext_Us[ind_U], mf.nb_basic_dof());
        mf.extend_vector(*(Us[ind_U]), ext_Us[ind_U]);
        if (Ws[ind_U]) {
          gmm::resize(ext_Ws[ind_U], mf.nb_basic_dof());
          mf.extend_vector(*(Ws[ind_U]), ext_Ws[ind_U]);
        } else gmm::resize(ext_Ws[ind_U], 0);
        iU.add(ind_U);
      }
      size_type ind_lambda = contact_boundaries[i].ind_lambda;
      if (ind_lambda != size_type(-1) && !(ilambda[ind_lambda])) {
        const mesh_fem &mf = *(contact_boundaries[i].mflambda);
        gmm::resize(ext_lambdas[ind_lambda], mf.nb_basic_dof());
        mf.extend_vector(*(lambdas[ind_lambda]), ext_lambdas[ind_lambda]);
        ilambda.add(ind_lambda);
      }
    }
  }

  void multi_contact_frame::normal_cone_simplification(void) {
    if (nodes_mode) {
      scalar_type threshold = ::cos(cut_angle);
      for (size_type i = 0; i < boundary_points_info.size(); ++i) {
        normal_cone &nc = boundary_points_info[i].normals;
        if (nc.size() > 1) {
          base_small_vector n_mean = nc[0];
          for (size_type j = 1; j < nc.size(); ++j) n_mean += nc[j];
          scalar_type nn_mean = gmm::vect_norm2(n_mean);
          GMM_ASSERT1(nn_mean != scalar_type(0), "oupssss");
          if (nn_mean != scalar_type(0)) {
            gmm::scale(n_mean, scalar_type(1)/nn_mean);
            bool reduce = true;
            for (size_type j = 0; j < nc.size(); ++j)
              if (gmm::vect_sp(n_mean, nc[j]) < threshold)
                { reduce = false; break; }
            if (reduce) {
              boundary_points_info[i].normals = normal_cone(n_mean);
            }
          }
        }
      }
    }
  }

  bool multi_contact_frame::test_normal_cones_compatibility
  (const normal_cone &nc1, const normal_cone &nc2) {
    for (size_type i = 0; i < nc1.size(); ++i)
      for (size_type j = 0; j < nc2.size(); ++j)
        if (gmm::vect_sp(nc1[i], nc2[j]) < scalar_type(0))
          return true;
    return false;
  }

  bool multi_contact_frame::test_normal_cones_compatibility
  (const base_small_vector &n, const normal_cone &nc2) {
    for (size_type j = 0; j < nc2.size(); ++j)
      if (gmm::vect_sp(n, nc2[j]) < scalar_type(0))
        return true;
    return false;
  }

  bool multi_contact_frame::are_dof_linked(size_type ib1, size_type idof1,
                                           size_type ib2, size_type idof2) {
    const mesh_fem &mf1 = mfdisp_of_boundary(ib1);
    const mesh_fem &mf2 = mfdisp_of_boundary(ib2);
    if ( &(mf1.linked_mesh()) != &(mf2.linked_mesh())) return false;
    GMM_ASSERT1(!(mf1.is_reduced()) && !(mf2.is_reduced()),
                "Nodal strategy can only be applied for non reduced fems");
    const mesh::ind_cv_ct &ic1 = mf1.convex_to_basic_dof(idof1);
    const mesh::ind_cv_ct &ic2 = mf2.convex_to_basic_dof(idof2);
    bool lk = false;
    for (size_type i = 0; i < ic1.size(); ++i) aux_dof_cv.add(ic1[i]);
    for (size_type i = 0; i < ic2.size(); ++i)
      if (aux_dof_cv.is_in(ic2[i])) { lk = true; break; }
    for (size_type i = 0; i < ic1.size(); ++i) aux_dof_cv.sup(ic1[i]);
    return lk;
  }

  bool multi_contact_frame::is_dof_linked(size_type ib1, size_type idof1,
                                          size_type ib2, size_type cv) {
    const mesh_fem &mf1 = mfdisp_of_boundary(ib1);
    const mesh_fem &mf2 = mfdisp_of_boundary(ib2);
    if ( &(mf1.linked_mesh()) != &(mf2.linked_mesh())) return false;
    GMM_ASSERT1(!(mf1.is_reduced()) && !(mf2.is_reduced()),
                "Nodal strategy can only be applied for non reduced fems");
    const mesh::ind_cv_ct &ic1 = mf1.convex_to_basic_dof(idof1);
    for (size_type i = 0; i < ic1.size(); ++i)
      if (cv == ic1[i]) return true;
    return false;
  }

  void multi_contact_frame::add_potential_contact_face
  (size_type ip, size_type ib, size_type ie, short_type iff) {
    bool found = false;
    std::vector<face_info> &sfi = potential_pairs[ip];
    for (size_type k = 0; k < sfi.size(); ++k)
      if (sfi[k].ind_boundary == ib &&
          sfi[k].ind_element == ie &&
          sfi[k].ind_face == iff) found = true;

    if (!found) sfi.push_back(face_info(ib, ie, iff));
  }

  void multi_contact_frame::clear_aux_info(void) {
    boundary_points = std::vector<base_node>();
    boundary_points_info = std::vector<boundary_point>();
    element_boxes.clear();
    element_boxes_info = std::vector<influence_box>();
    potential_pairs = std::vector<std::vector<face_info> >();
  }

  multi_contact_frame::multi_contact_frame(size_type NN, scalar_type r_dist,
                                           bool dela, bool selfc,
                                           scalar_type cut_a,
                                           bool rayt, int nmode, bool refc)
    : N(NN), self_contact(selfc), ref_conf(refc), use_delaunay(dela),
      nodes_mode(nmode), raytrace(rayt), release_distance(r_dist),
      cut_angle(cut_a), EPS(1E-8), md(0), coordinates(N), pt_eval(N) {
    if (N > 0) coordinates[0] = "x";
    if (N > 1) coordinates[1] = "y";
    if (N > 2) coordinates[2] = "z";
    if (N > 3) coordinates[3] = "w";
    GMM_ASSERT1(N <= 4, "Complete the definition for contact in "
                  "dimension greater than 4");
  }

  multi_contact_frame::multi_contact_frame(const model &mdd, size_type NN,
                                           scalar_type r_dist,
                                           bool dela, bool selfc,
                                           scalar_type cut_a,
                                           bool rayt, int nmode, bool refc)
    : N(NN), self_contact(selfc), ref_conf(refc),
      use_delaunay(dela), nodes_mode(nmode), raytrace(rayt),
      release_distance(r_dist), cut_angle(cut_a), EPS(1E-8), md(&mdd),
      coordinates(N), pt_eval(N) {
    if (N > 0) coordinates[0] = "x";
    if (N > 1) coordinates[1] = "y";
    if (N > 2) coordinates[2] = "z";
    if (N > 3) coordinates[3] = "w";
    GMM_ASSERT1(N <= 4, "Complete the definition for contact in "
                  "dimension greater than 4");
  }

  size_type multi_contact_frame::add_obstacle(const std::string &obs) {
    size_type ind = obstacles.size();
    obstacles.push_back(obs);
    obstacles_velocities.push_back("");
#if GETFEM_HAVE_MUPARSER_MUPARSER_H || GETFEM_HAVE_MUPARSER_H

    mu::Parser mu;
    obstacles_parsers.push_back(mu);
    obstacles_parsers[ind].SetExpr(obstacles[ind]);
    for (size_type k = 0; k < N; ++k)
      obstacles_parsers[ind].DefineVar(coordinates[k], &pt_eval[k]);
#else
    GMM_ASSERT1(false, "You have to link muparser with getfem to deal "
                "with rigid body obstacles");
#endif
    return ind;
  }



  size_type multi_contact_frame::add_master_boundary
  (const mesh_im &mim, const mesh_fem *mfu,
   const model_real_plain_vector *U, size_type reg,
   const mesh_fem *mflambda, const model_real_plain_vector *lambda,
   const model_real_plain_vector *w,
   const std::string &vvarname,
   const std::string &mmultname, const std::string &wname) {
    GMM_ASSERT1(mfu->linked_mesh().dim() == N,
                "Mesh dimension is " << mfu->linked_mesh().dim()
                << "should be " << N << ".");
    GMM_ASSERT1(&(mfu->linked_mesh()) == &(mim.linked_mesh()),
                "Integration and finite element are not on the same mesh !");
    if (mflambda)
      GMM_ASSERT1(&(mflambda->linked_mesh()) == &(mim.linked_mesh()),
                  "Integration and finite element are not on the same mesh !");
    contact_boundary cb(reg, mfu, mim, add_U(U, vvarname, w, wname),
                        mflambda, add_lambda(lambda, mmultname));
    contact_boundaries.push_back(cb);
    return size_type(contact_boundaries.size() - 1);
  }

  size_type multi_contact_frame::add_slave_boundary
  (const mesh_im &mim, const mesh_fem *mfu,
   const model_real_plain_vector *U, size_type reg,
   const mesh_fem *mflambda, const model_real_plain_vector *lambda,
   const model_real_plain_vector *w,
   const std::string &vvarname,
   const std::string &mmultname, const std::string &wname) {
    size_type ind
      = add_master_boundary(mim, mfu, U, reg, mflambda, lambda, w,
                            vvarname, mmultname, wname);
    contact_boundaries[ind].slave = true;
    return ind;
  }


  size_type multi_contact_frame::add_master_boundary
  (const mesh_im &mim, size_type reg, const std::string &vvarname,
   const std::string &mmultname, const std::string &wname) {
    GMM_ASSERT1(md, "This multi contact frame object is not linked "
                "to a model");
    const mesh_fem *mfl(0);
    const model_real_plain_vector *l(0);
    if (mmultname.size()) {
      mfl = &(md->mesh_fem_of_variable(mmultname));
      l = &(md->real_variable(mmultname));
    }
    const model_real_plain_vector *w(0);
    if (wname.compare(vvarname) == 0) {
      GMM_ASSERT1(md->n_iter_of_variable(vvarname) > 1, "More than one "
                 "versions of the displacement variable were expected here");
      w = &(md->real_variable(vvarname,1));
    }
    else if (wname.size()) {
      GMM_ASSERT1(&(md->mesh_fem_of_variable(wname))
                 == &(md->mesh_fem_of_variable(vvarname)), "The previous "
                 "displacement should be defined on the same mesh_fem as the "
                 "current one");
      w = &(md->real_variable(wname));
    }
    return add_master_boundary(mim, &(md->mesh_fem_of_variable(vvarname)),
                               &(md->real_variable(vvarname)), reg, mfl, l, w,
                               vvarname, mmultname, wname);
  }

  size_type multi_contact_frame::add_slave_boundary
  (const mesh_im &mim, size_type reg, const std::string &vvarname,
   const std::string &mmultname, const std::string &wname) {
    GMM_ASSERT1(md, "This multi contact frame object is not linked "
                "to a model");
    const mesh_fem *mfl(0);
    const model_real_plain_vector *l(0);
    if (mmultname.size()) {
      mfl = &(md->mesh_fem_of_variable(mmultname));
      l = &(md->real_variable(mmultname));
    }
    const model_real_plain_vector *w(0);
    if (wname.compare(vvarname) == 0) {
      GMM_ASSERT1(md->n_iter_of_variable(vvarname) > 1, "More than one "
                 "versions of the displacement variable were expected here");
      w = &(md->real_variable(vvarname,1));
    }
    else if (wname.size()) {
      GMM_ASSERT1(&(md->mesh_fem_of_variable(wname))
                 == &(md->mesh_fem_of_variable(vvarname)), "The previous "
                 "displacement should be defined on the same mesh_fem as the "
                 "current one");
      w = &(md->real_variable(wname));
    }
    return add_slave_boundary(mim, &(md->mesh_fem_of_variable(vvarname)),
                              &(md->real_variable(vvarname)), reg, mfl, l, w,
                              vvarname, mmultname, wname);
  }


  void multi_contact_frame::compute_boundary_points(bool slave_only) {
    fem_precomp_pool fppool;
    base_matrix G;
    model_real_plain_vector coeff;

    for (size_type i = 0; i < contact_boundaries.size(); ++i)
      if (!slave_only || is_slave_boundary(i)) {
        size_type bnum = region_of_boundary(i);
        const mesh_fem &mfu = mfdisp_of_boundary(i);
        const mesh_im &mim = mim_of_boundary(i);
        const model_real_plain_vector &U = disp_of_boundary(i);
        const mesh &m = mfu.linked_mesh();
        bool on_fem_nodes =
          boundary_has_fem_nodes(is_slave_boundary(i), nodes_mode);

        base_node val(N), bmin(N), bmax(N);
        base_small_vector n0(N), n(N), n_mean(N);
        base_matrix grad(N,N);
        mesh_region region = m.region(bnum);
        GMM_ASSERT1(mfu.get_qdim() == N, "Wrong mesh_fem qdim");


        dal::bit_vector dof_already_interpolated;
        std::vector<size_type> dof_ind(mfu.nb_basic_dof());
        for (getfem::mr_visitor v(region,m); !v.finished(); ++v) {
          size_type cv = v.cv();
          bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
          pfem pf_s = mfu.fem_of_element(cv);

          if (!ref_conf)
            slice_vector_on_basic_dof_of_element(mfu, U, cv, coeff);
          bgeot::vectors_to_base_matrix
            (G, mfu.linked_mesh().points_of_convex(cv));

          pfem_precomp pfp(0);
          size_type nbptf(0);
          std::vector<size_type> indpt, indpfp;
          if (on_fem_nodes) {
            dim_type qqdim = mfu.get_qdim() / pf_s->target_dim();
            pfp = fppool(pf_s, pf_s->node_tab(cv));
            nbptf = pf_s->node_convex(cv).structure()->nb_points_of_face(v.f());
            indpt.resize(nbptf); indpfp.resize(nbptf);
            for (short_type ip = 0; ip < nbptf; ++ip) {
              indpt[ip] =
                mfu.ind_basic_dof_of_face_of_element(cv,v.f())[ip*qqdim];
              indpfp[ip] =
                pf_s->node_convex(cv).structure()->ind_points_of_face(v.f())[ip];
            }
          }
          else {
            pintegration_method pim = mim.int_method_of_element(cv);
            GMM_ASSERT1(pim, "Integration method should be defined");
            pfp = fppool(pf_s,&(pim->approx_method()->integration_points()));
            nbptf = pim->approx_method()->nb_points_on_face(v.f());
            indpt.resize(nbptf); indpfp.resize(nbptf);
            for (short_type ip = 0; ip < nbptf; ++ip)
              indpt[ip] = indpfp[ip] =
                pim->approx_method()->ind_first_point_on_face(v.f())+ip;
          }
          fem_interpolation_context ctx(pgt,pfp,size_type(-1),G,cv,v.f());

          for (short_type ip = 0; ip < nbptf; ++ip) {
            ctx.set_ii(indpfp[ip]);

            size_type ind = indpt[ip];
            if (!(on_fem_nodes && dof_already_interpolated[ind])) {
              if (!ref_conf) {
                pf_s->interpolation(ctx, coeff, val, dim_type(N));
                val += ctx.xreal();
              } else {
                val = ctx.xreal();
              }
              if (on_fem_nodes) dof_ind[ind] = boundary_points.size();

            }

            // unit normal vector computation
            compute_normal(ctx, v.f(), ref_conf, coeff, n0, n, grad);
            n /= gmm::vect_norm2(n);

            if (on_fem_nodes && dof_already_interpolated[ind]) {
              boundary_points_info[dof_ind[ind]].normals.add_normal(n);
            } else {
              boundary_points.push_back(val);
              boundary_points_info.push_back(boundary_point(ctx.xreal(), i, cv,
                                                            v.f(), ind, n));
            }

            if (on_fem_nodes) dof_already_interpolated.add(ind);
          }
        }
      }
  }

  void multi_contact_frame::compute_potential_contact_pairs_delaunay(void) {

    compute_boundary_points();
    normal_cone_simplification();
    potential_pairs = std::vector<std::vector<face_info> >();
    potential_pairs.resize(boundary_points.size());

    gmm::dense_matrix<size_type> simplexes;
    base_small_vector rr(N);
    // Necessary ?
    // for (size_type i = 0; i < boundary_points.size(); ++i) {
    //   gmm::fill_random(rr);
    //   boundary_points[i] += 1E-9*rr;
    // }
    getfem::delaunay(boundary_points, simplexes);

    // connectivity analysis
    for (size_type is = 0; is < gmm::mat_ncols(simplexes); ++is) {

      for (size_type i = 1; i <= N; ++i)
        for (size_type j = 0; j < i; ++j) {
          size_type ipt1 = simplexes(i, is), ipt2 = simplexes(j, is);
          boundary_point *pt_info1 = &(boundary_points_info[ipt1]);
          boundary_point *pt_info2 = &(boundary_points_info[ipt2]);
          size_type ib1 = pt_info1->ind_boundary;
          size_type ib2 = pt_info2->ind_boundary;
          bool sl1 = is_slave_boundary(ib1);
          bool sl2 = is_slave_boundary(ib2);
          if (!sl1 && sl2) { // The slave in first if any
            std::swap(ipt1, ipt2);
            std::swap(pt_info1, pt_info2);
            std::swap(ib1, ib2);
            std::swap(sl1, sl2);
          }
          size_type ir1 = region_of_boundary(ib1);
          size_type ir2 = region_of_boundary(ib2);
          const mesh_fem &mf1 = mfdisp_of_boundary(ib1);
          const mesh_fem &mf2 = mfdisp_of_boundary(ib2);

          // CRITERION 1 : The unit normal cone / vector are compatible
          //               and the two points are not in the same element.
          if (
              // slave-master case
              ((sl1 && !sl2)
               // master-master self-contact case
               || (self_contact && !sl1 && !sl2))
              // test of unit normal vectors or cones
              && test_normal_cones_compatibility(pt_info1->normals,
                                                 pt_info2->normals)
              // In case of self-contact, test if the two points share the
              // same element.
              && (sl1
                  || ((nodes_mode < 2)
                      && (( &(mf1.linked_mesh()) != &(mf2.linked_mesh()))
                          || (pt_info1->ind_element != pt_info2->ind_element)))
                  || ((nodes_mode == 2)
                      && !(are_dof_linked(ib1, pt_info1->ind_pt,
                                          ib2, pt_info2->ind_pt)))
                  )
              ) {

            // Store the potential contact pairs

            if (boundary_has_fem_nodes(sl2, nodes_mode)) {
              const mesh::ind_cv_ct &ic2
                = mf2.convex_to_basic_dof(pt_info2->ind_pt);
              for (size_type k = 0; k < ic2.size(); ++k) {
                mesh_region::face_bitset fbs
                  = mf2.linked_mesh().region(ir2).faces_of_convex(ic2[k]);
                short_type nbf = mf2.linked_mesh().nb_faces_of_convex(ic2[k]);
                for (short_type f = 0; f < nbf; ++f)
                  if (fbs.test(f))
                    add_potential_contact_face(ipt1,
                                               pt_info2->ind_boundary,
                                               ic2[k], f);
              }
            } else
              add_potential_contact_face(ipt1, pt_info2->ind_boundary,
                                         pt_info2->ind_element,
                                         pt_info2->ind_face);

            if (self_contact && !sl1 && !sl2) {
              if (boundary_has_fem_nodes(sl2, nodes_mode)) {
                const mesh::ind_cv_ct &ic1
                  = mf1.convex_to_basic_dof(pt_info1->ind_pt);
                for (size_type k = 0; k < ic1.size(); ++k) {
                  mesh_region::face_bitset fbs
                    = mf1.linked_mesh().region(ir1).faces_of_convex(ic1[k]);
                  short_type nbf = mf1.linked_mesh().nb_faces_of_convex(ic1[k]);
                  for (short_type f = 0; f < nbf; ++f)
                    if (fbs.test(f))
                      add_potential_contact_face(ipt2,
                                                 pt_info1->ind_boundary,
                                                 ic1[k], f);
                }
              } else
                add_potential_contact_face(ipt2, pt_info1->ind_boundary,
                                           pt_info1->ind_element,
                                           pt_info1->ind_face);
            }

          }

        }
    }
  }


  void multi_contact_frame::compute_influence_boxes(void) {
    fem_precomp_pool fppool;
    bool avert = false;
    base_matrix G;
    model_real_plain_vector coeff;

    for (size_type i = 0; i < contact_boundaries.size(); ++i)
      if (!is_slave_boundary(i)) {
        size_type bnum = region_of_boundary(i);
        const mesh_fem &mfu = mfdisp_of_boundary(i);
        const model_real_plain_vector &U = disp_of_boundary(i);
        const mesh &m = mfu.linked_mesh();

        base_node val(N), bmin(N), bmax(N);
        base_small_vector n0(N), n(N), n_mean(N);
        base_matrix grad(N,N);
        mesh_region region = m.region(bnum);
        GMM_ASSERT1(mfu.get_qdim() == N, "Wrong mesh_fem qdim");

        dal::bit_vector points_already_interpolated;
        std::vector<base_node> transformed_points(m.nb_max_points());
        for (getfem::mr_visitor v(region,m); !v.finished(); ++v) {
          size_type cv = v.cv();
          bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
          pfem pf_s = mfu.fem_of_element(cv);
          pfem_precomp pfp = fppool(pf_s, &(pgt->geometric_nodes()));
          if (!ref_conf)
            slice_vector_on_basic_dof_of_element(mfu, U, cv, coeff);
          bgeot::vectors_to_base_matrix
            (G, mfu.linked_mesh().points_of_convex(cv));
          fem_interpolation_context ctx(pgt,pfp,size_type(-1), G, cv,
                                        size_type(-1));

          size_type nb_pt_on_face = 0;
          dal::bit_vector points_on_face;
          bgeot::pconvex_structure cvs = pgt->structure();
          for (size_type k = 0; k < cvs->nb_points_of_face(v.f()); ++k)
            points_on_face.add(cvs->ind_points_of_face(v.f())[k]);

          gmm::clear(n_mean);
          size_type nbd_t = pgt->nb_points();
          for (short_type ip = 0; ip < nbd_t; ++ip) {
            size_type ind = m.ind_points_of_convex(cv)[ip];

            // computation of transformed vertex
            if (!(points_already_interpolated.is_in(ind))) {
              ctx.set_ii(ip);
              if (!ref_conf) {
                pf_s->interpolation(ctx, coeff, val, dim_type(N));
                val += ctx.xreal();
                transformed_points[ind] = val;
              } else {
                transformed_points[ind] = ctx.xreal();
              }
              points_already_interpolated.add(ind);
            } else {
              val = transformed_points[ind];
            }

            if (ip == 0) // computation of bounding box
              bmin = bmax = val;
            else {
              for (size_type k = 0; k < N; ++k) {
                bmin[k] = std::min(bmin[k], val[k]);
                bmax[k] = std::max(bmax[k], val[k]);
              }
            }

            // computation of unit normal vector if the vertex is on the face
            if (points_on_face[ip]) {
              compute_normal(ctx, v.f(), ref_conf, coeff, n0, n, grad);
              n /= gmm::vect_norm2(n);
              n_mean += n;
              ++nb_pt_on_face;
            }

          }

          // is nb_pt_on_face really necessary, is this possible to occur?
          GMM_ASSERT1(nb_pt_on_face,
                      "This element has no vertex on considered face !");

          // Computation of influence box :
          // offset of the bounding box relatively to the release distance
          scalar_type h = bmax[0] - bmin[0];
          for (size_type k = 1; k < N; ++k) h = std::max(h, bmax[k]-bmin[k]);
          if (h < release_distance/scalar_type(40) && !avert) {
            GMM_WARNING1("Found an element whose size is smaller than 1/40 "
                         "of the release distance. You should probably "
                         "adapt the release distance.");
            avert = true;
          }
          for (size_type k = 0; k < N; ++k)
            { bmin[k] -= release_distance; bmax[k] += release_distance; }

          // Store the influence box and additional information.
          element_boxes.add_box(bmin, bmax, element_boxes_info.size());
          n_mean /= gmm::vect_norm2(n_mean);
          element_boxes_info.push_back(influence_box(i, cv, v.f(), n_mean));
        }
      }
  }

  void multi_contact_frame::compute_potential_contact_pairs_influence_boxes(void) {
    compute_influence_boxes();
    compute_boundary_points(!self_contact); // vraiment necessaire ?
    normal_cone_simplification();
    potential_pairs = std::vector<std::vector<face_info> >();
    potential_pairs.resize(boundary_points.size());

    for (size_type ip = 0; ip < boundary_points.size(); ++ip) {

      bgeot::rtree::pbox_set bset;
      element_boxes.find_boxes_at_point(boundary_points[ip], bset);
      boundary_point *pt_info = &(boundary_points_info[ip]);
      const mesh_fem &mf1 = mfdisp_of_boundary(pt_info->ind_boundary);
      size_type ib1 = pt_info->ind_boundary;

      bgeot::rtree::pbox_set::iterator it = bset.begin();
      for (; it != bset.end(); ++it) {
        influence_box &ibx = element_boxes_info[(*it)->id];
        size_type ib2 = ibx.ind_boundary;
        const mesh_fem &mf2 = mfdisp_of_boundary(ib2);

        // CRITERION 1 : The unit normal cone / vector are compatible
        //               and the two points are not in the same element.
        if (
            test_normal_cones_compatibility(ibx.mean_normal,
                                            pt_info->normals)
            // In case of self-contact, test if the points and the face
            // share the same element.
            && (((nodes_mode < 2)
                 && (( &(mf1.linked_mesh()) != &(mf2.linked_mesh()))
                     || (pt_info->ind_element != ibx.ind_element)))
                || ((nodes_mode == 2)
                    && !(is_dof_linked(ib1, pt_info->ind_pt,
                                       ibx.ind_boundary, ibx.ind_element)))
                )
            ) {

          add_potential_contact_face(ip, ibx.ind_boundary, ibx.ind_element,
                                     ibx.ind_face);
        }
      }

    }
  }

  struct proj_pt_surf_cost_function_object {
    size_type N;
    scalar_type EPS;
    const base_node &x0, &x;
    fem_interpolation_context &ctx;
    const model_real_plain_vector &coeff;
    const std::vector<base_small_vector> &ti;
    bool ref_conf;
    mutable base_node dxy;
    mutable base_matrix grad, gradtot;

    scalar_type operator()(const base_small_vector& a) const {
      base_node xx = x0;
      for (size_type i= 0; i < N-1; ++i) xx += a[i] * ti[i];
      ctx.set_xref(xx);
      if (!ref_conf) {
        ctx.pf()->interpolation(ctx, coeff, dxy, dim_type(N));
        dxy += ctx.xreal() - x;
      } else
        dxy = ctx.xreal() - x;
      return gmm::vect_norm2(dxy)/scalar_type(2);
    }

    scalar_type operator()(const base_small_vector& a,
                           base_small_vector &grada) const {
      base_node xx = x0;
      for (size_type i = 0; i < N-1; ++i) xx += a[i] * ti[i];
      ctx.set_xref(xx);
      if (!ref_conf) {
        ctx.pf()->interpolation(ctx, coeff, dxy, dim_type(N));
        dxy += ctx.xreal() - x;
        ctx.pf()->interpolation_grad(ctx, coeff, grad, dim_type(N));
        gmm::add(gmm::identity_matrix(), grad);
        gmm::mult(grad, ctx.K(), gradtot);
      } else {
        dxy = ctx.xreal() - x;
        gmm::copy(ctx.K(), gradtot);
      }
      for (size_type i = 0; i < N-1; ++i)
        grada[i] = gmm::vect_sp(gradtot, ti[i], dxy);
      return gmm::vect_norm2(dxy)/scalar_type(2);
    }
    void operator()(const base_small_vector& a,
                    base_matrix &hessa) const {
      base_small_vector b = a;
      base_small_vector grada(N-1), gradb(N-1);
      (*this)(b, grada);
      for (size_type i = 0; i < N-1; ++i) {
        b[i] += EPS;
        (*this)(b, gradb);
        for (size_type j = 0; j < N-1; ++j)
          hessa(j, i) = (gradb[j] - grada[j])/EPS;
        b[i] -= EPS;
      }
    }

    proj_pt_surf_cost_function_object
    (const base_node &x00, const base_node &xx,
     fem_interpolation_context &ctxx,
     const model_real_plain_vector &coefff,
     const std::vector<base_small_vector> &tii,
     scalar_type EPSS, bool rc)
      : N(gmm::vect_size(x00)), EPS(EPSS), x0(x00), x(xx),
        ctx(ctxx), coeff(coefff), ti(tii), ref_conf(rc),
        dxy(N), grad(N,N), gradtot(N,N) {}

  };

  struct raytrace_pt_surf_cost_function_object {
    size_type N;
    const base_node &x0, &x;
    fem_interpolation_context &ctx;
    const model_real_plain_vector &coeff;
    const std::vector<base_small_vector> &ti;
    const std::vector<base_small_vector> &Ti;
    bool ref_conf;
    mutable base_node dxy;
    mutable base_matrix grad, gradtot;

    void operator()(const base_small_vector& a,
                    base_small_vector &res) const {
      base_node xx = x0;
      for (size_type i = 0; i < N-1; ++i) xx += a[i] * ti[i];
      ctx.set_xref(xx);
      if (!ref_conf) {
        ctx.pf()->interpolation(ctx, coeff, dxy, dim_type(N));
        dxy += ctx.xreal() - x;
      } else
        dxy = ctx.xreal() - x;
      for (size_type i = 0; i < N-1; ++i)
        res[i] = gmm::vect_sp(dxy, Ti[i]);
    }

    void operator()(const base_small_vector& a,
                    base_matrix &hessa) const {
      base_node xx = x0;
      for (size_type i = 0; i < N-1; ++i) xx += a[i] * ti[i];
      ctx.set_xref(xx);
      if (!ref_conf) {
        ctx.pf()->interpolation_grad(ctx, coeff, grad, dim_type(N));
        gmm::add(gmm::identity_matrix(), grad);
        gmm::mult(grad, ctx.K(), gradtot);
      } else {
        gmm::copy(ctx.K(), gradtot);
      }
      for (size_type i = 0; i < N-1; ++i)
        for (size_type j = 0; j < N-1; ++j)
          hessa(j, i) = gmm::vect_sp(gradtot, ti[i], Ti[j]);
    }


    raytrace_pt_surf_cost_function_object
    (const base_node &x00, const base_node &xx,
     fem_interpolation_context &ctxx,
     const model_real_plain_vector &coefff,
     const std::vector<base_small_vector> &tii,
     const std::vector<base_small_vector> &Tii,
     bool rc)
      : N(gmm::vect_size(x00)), x0(x00), x(xx),
        ctx(ctxx), coeff(coefff), ti(tii), Ti(Tii), ref_conf(rc),
        dxy(N), grad(N,N), gradtot(N,N) {}

  };

  // Ideas to improve efficiency :
  // - From an iteration to another, is it possible to simplify the
  //   computation ? For instance in testing the old contact pairs ...
  //   But how to detect new contact situations ?
  // - A pre-test before projection (for Delaunay) : if the distance to a
  //   node is greater than the release distance + h then give up.
  // - Case J3 of valid/invalid contact situations is not really taken into
  //   account. How to take it into account in a cheap way ?

  void multi_contact_frame::compute_contact_pairs(void) {
    base_matrix G, grad(N,N);
    model_real_plain_vector coeff;
    base_small_vector a(N-1), ny(N);
    base_node y(N);
    std::vector<base_small_vector> ti(N-1), Ti(N-1);
    size_type nbwarn(0);

    // double time = dal::uclock_sec();

    clear_aux_info();
    contact_pairs = std::vector<contact_pair>();

    if (!ref_conf) extend_vectors();

    bool only_slave(true), only_master(true);
    for (size_type i = 0; i < contact_boundaries.size(); ++i)
      if (is_slave_boundary(i)) only_master = false;
      else only_slave = false;

    if (only_master && !self_contact) {
      GMM_WARNING1("There is only master boundary and no self-contact to detect. Exiting");
      return;
    }

    if (only_slave) {
      compute_boundary_points();
      potential_pairs.resize(boundary_points.size());
    }
    else if (use_delaunay)
      compute_potential_contact_pairs_delaunay();
    else
      compute_potential_contact_pairs_influence_boxes();

    // cout << "Time for computing potential pairs: " << dal::uclock_sec() - time << endl; time = dal::uclock_sec();


    // Scan of potential pairs
    for (size_type ip = 0; ip < potential_pairs.size(); ++ip) {
      bool first_pair_found = false;
      const base_node &x = boundary_points[ip];
      boundary_point &bpinfo = boundary_points_info[ip];
      size_type ibx = bpinfo.ind_boundary;
      bool slx = is_slave_boundary(ibx);
      scalar_type d0 = 1E300, d1, d2;

      base_small_vector nx = bpinfo.normals[0];
      if (raytrace) {
        if (bpinfo.normals.size() > 1) { // take the mean normal vector
          for (size_type i = 1; i < bpinfo.normals.size(); ++i)
            gmm::add(bpinfo.normals[i], nx);
          scalar_type nnx = gmm::vect_norm2(nx);
          GMM_ASSERT1(nnx != scalar_type(0), "Invalid normal cone");
          gmm::scale(nx, scalar_type(1)/nnx);
        }
      }

      if (self_contact || slx) {
#if GETFEM_HAVE_MUPARSER_MUPARSER_H || GETFEM_HAVE_MUPARSER_H
        // Detect here the nearest rigid obstacle (taking into account
        // the release distance)
        size_type irigid_obstacle(-1);
        gmm::copy(x, pt_eval);
        for (size_type i = 0; i < obstacles.size(); ++i) {
          d1 = scalar_type(obstacles_parsers[i].Eval());
          if (gmm::abs(d1) < release_distance && d1 < d0) {

            for (size_type j=0; j < bpinfo.normals.size(); ++j) {
              gmm::add(gmm::scaled(bpinfo.normals[j], EPS), pt_eval);
              d2 =  scalar_type(obstacles_parsers[i].Eval());
              if (d2 < d1) { d0 = d1; irigid_obstacle = i; break; }
              gmm::copy(x, pt_eval);
            }
          }
        }

        if (irigid_obstacle != size_type(-1)) {

          gmm::copy(x, pt_eval);
          gmm::copy(x, y);
          size_type nit = 0, nb_fail = 0;
          scalar_type alpha(0), beta(0);
          d1 = d0;

          while (++nit < 50 && nb_fail < 3) {
            for (size_type k = 0; k < N; ++k) {
              pt_eval[k] += EPS;
              d2 = scalar_type(obstacles_parsers[irigid_obstacle].Eval());
              ny[k] = (d2 - d1) / EPS;
              pt_eval[k] -= EPS;
            }

            if (gmm::abs(d1) < 1E-13)
              break; // point already lies on the rigid obstacle surface

            // ajouter un test de divergence ...
            for (scalar_type lambda(1); lambda >= 1E-3; lambda /= scalar_type(2)) {
              if (raytrace) {
                alpha = beta - lambda * d1 / gmm::vect_sp(ny, nx);
                gmm::add(x, gmm::scaled(nx, alpha), pt_eval);
              } else {
                gmm::add(gmm::scaled(ny, -d1/gmm::vect_norm2_sqr(ny)), y, pt_eval);
              }
              d2 = scalar_type(obstacles_parsers[irigid_obstacle].Eval());
//               if (nit > 10)
//                 cout << "nit = " << nit << " lambda = " << lambda
//                      << " alpha = " << alpha << " d2 = " << d2
//                      << " d1  = " << d1 << endl;
              if (gmm::abs(d2) < gmm::abs(d1)) break;
            }
            if (raytrace &&
                gmm::abs(beta - d1 / gmm::vect_sp(ny, nx)) > scalar_type(500))
              nb_fail++;
            gmm::copy(pt_eval, y); beta = alpha; d1 = d2;
          }

          if (gmm::abs(d1) > 1E-8) {
            GMM_WARNING1("Projection/raytrace on rigid obstacle failed");
            continue;
          }

          // CRITERION 4 for rigid bodies : Apply the release distance
          if (gmm::vect_dist2(y, x) > release_distance)
            continue;

          gmm::copy(pt_eval, y);
          ny /= gmm::vect_norm2(ny);

          d0 = gmm::vect_dist2(y, x) * gmm::sgn(d0);
          contact_pair ct(x, nx, bpinfo, y, ny, irigid_obstacle, d0);

          contact_pairs.push_back(ct);
          first_pair_found = true;
        }
#else
        if (obstacles.size() > 0)
          GMM_WARNING1("Rigid obstacles are ignored. Recompile with "
                       "muParser to account for rigid obstacles");
#endif
      }

      // if (potential_pairs[ip].size())
      // cout << "number of potential pairs for point " << ip << " : " << potential_pairs[ip].size() << endl;
      for (size_type ipf = 0; ipf < potential_pairs[ip].size(); ++ipf) {
        // Point to surface projection. Principle :
        //  - One parametrizes first the face on the reference element by
        //    obtaining a point x_0 on that face and t_i, i=1..d-1 some
        //    orthonormals tangent vectors to the face.
        //  - Let y_0 be the point to be projected and y the searched
        //    projected point. Then one searches for the minimum of
        //    J = (1/2)|| y - x ||
        //    with
        //    y = \phi(x0 + a_i t_i)
        //    (with a summation on i), where \phi = I+u(\tau(x)), and \tau
        //    the geometric transformation between reference and real
        //    elements.
        //  - The gradient of J with respect to a_i is
        //    \partial_{a_j} J = (\phi(x0 + a_i t_i) - x)
        //                       . (\nabla \phi(x0 + a_i t_i) t_j
        //  - A Newton algorithm is applied.
        //  - If it fails, a BFGS is called.

        const face_info &fi = potential_pairs[ip][ipf];
        size_type ib = fi.ind_boundary;
        size_type cv = fi.ind_element;
        short_type iff = fi.ind_face;

        const mesh_fem &mfu = mfdisp_of_boundary(ib);
        const mesh &m = mfu.linked_mesh();
        pfem pf_s = mfu.fem_of_element(cv);
        bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);

        if (!ref_conf)
          slice_vector_on_basic_dof_of_element(mfu, disp_of_boundary(ib),
                                               cv, coeff);

        bgeot::vectors_to_base_matrix(G, m.points_of_convex(cv));

        const base_node &x0 = pf_s->ref_convex(cv)->points_of_face(iff)[0];
        fem_interpolation_context ctx(pgt, pf_s, x0, G, cv, iff);

        const base_small_vector &n0 = pf_s->ref_convex(cv)->normals()[iff];
        for (size_type k = 0; k < N-1; ++k) { // A basis for the face
          gmm::resize(ti[k], N);
          scalar_type norm(0);
          while(norm < 1E-5) {
            gmm::fill_random(ti[k]);
            ti[k] -= gmm::vect_sp(ti[k], n0) * n0;
            for (size_type l = 0; l < k; ++l)
              ti[k] -= gmm::vect_sp(ti[k], ti[l]) * ti[l];
            norm = gmm::vect_norm2(ti[k]);
          }
          ti[k] /= norm;
        }

        bool converged = false;
        scalar_type residual(0);


        if (raytrace) { // Raytrace search for y by a Newton algorithm

          base_small_vector res(N-1), res2(N-1), dir(N-1), b(N-1);

          base_matrix hessa(N-1, N-1);
          gmm::clear(a);

          for (size_type k = 0; k < N-1; ++k) {
            gmm::resize(Ti[k], N);
            scalar_type norm(0);
            while (norm < 1E-5) {
              gmm::fill_random(Ti[k]);
              Ti[k] -= gmm::vect_sp(Ti[k], nx) * nx;
              for (size_type l = 0; l < k; ++l)
                Ti[k] -= gmm::vect_sp(Ti[k], Ti[l]) * Ti[l];
              norm = gmm::vect_norm2(Ti[k]);
            }
            Ti[k] /= norm;
          }

          raytrace_pt_surf_cost_function_object pps(x0, x, ctx, coeff, ti, Ti,
                                                    ref_conf);

          pps(a, res);
          residual = gmm::vect_norm2(res);
          scalar_type residual2(0), det(0);
          bool exited = false;
          size_type nbfail = 0, niter = 0;
          for (;residual > 2E-12 && niter <= 30; ++niter) {

            for (size_type subiter(0);;) {
              pps(a, hessa);
              det = gmm::abs(gmm::lu_inverse(hessa, false));
              if (det > 1E-15) break;
              for (size_type i = 0; i < N-1; ++i)
                a[i] += gmm::random() * 1E-7;
              if (++subiter > 4) break;
            }
            if (det <= 1E-15) break;
            // Computation of the descent direction
            gmm::mult(hessa, gmm::scaled(res, scalar_type(-1)), dir);

            if (gmm::vect_norm2(dir) > scalar_type(10)) nbfail++;
            if (nbfail >= 4) break;

            // Line search
            scalar_type lambda(1);
            for (size_type j = 0; j < 5; ++j) {
              gmm::add(a, gmm::scaled(dir, lambda), b);
              pps(b, res2);
              residual2 = gmm::vect_norm2(res2);
              if (residual2 < residual) break;
              lambda /= ((j < 3) ? scalar_type(2) : scalar_type(5));
            }

            residual = residual2;
            gmm::copy(res2, res);
            gmm::copy(b, a);
            scalar_type dist_ref = gmm::vect_norm2(a);
//             if (niter == 15)
//               cout << "more than 15 iterations " << a
//                    << " dir " << dir << " nbfail : " << nbfail << endl;
            if (niter > 1 && dist_ref > 15) break;
            if (niter > 5 && dist_ref > 8) break;
            if ((niter > 1 && dist_ref > 7) || nbfail == 3) exited = true;
          }
          converged = (gmm::vect_norm2(res) < 2E-6);
          GMM_ASSERT1(!((exited && converged &&
                         pf_s->ref_convex(cv)->is_in(ctx.xref()) < 1E-6)),
                      "A non conformal case !! " << gmm::vect_norm2(res)
                      << " : " << nbfail << " : " << niter);

        } else { // Classical projection for y

          proj_pt_surf_cost_function_object pps(x0, x, ctx, coeff, ti,
                                                EPS, ref_conf);

          // Projection could be ameliorated by finding a starting point near
          // x (with respect to the integration method, for instance).

          // A specific (Quasi) Newton algorithm for computing the projection
          base_small_vector grada(N-1), dir(N-1), b(N-1);
          gmm::clear(a);
          base_matrix hessa(N-1, N-1);
          scalar_type det(0);

          scalar_type dist = pps(a, grada);
          for (size_type niter = 0;
               gmm::vect_norm2(grada) > 1E-12 && niter <= 50; ++niter) {

            for (size_type subiter(0);;) {
              pps(a, hessa);
              det = gmm::abs(gmm::lu_inverse(hessa, false));
              if (det > 1E-15) break;
              for (size_type i = 0; i < N-1; ++i)
                a[i] += gmm::random() * 1E-7;
              if (++subiter > 4) break;
            }
            if (det <= 1E-15) break;
            // Computation of the descent direction
            gmm::mult(hessa, gmm::scaled(grada, scalar_type(-1)), dir);

            // Line search
            for (scalar_type lambda(1);
                 lambda >= 1E-3; lambda /= scalar_type(2)) {
              gmm::add(a, gmm::scaled(dir, lambda), b);
              if (pps(b) < dist) break;
              gmm::add(a, gmm::scaled(dir, -lambda), b);
              if (pps(b) < dist) break;
            }
            gmm::copy(b, a);
            dist = pps(a, grada);
          }

          converged = (gmm::vect_norm2(grada) < 2E-6);

          if (!converged) { // Try with BFGS
            gmm::iteration iter(1E-12, 0 /* noisy*/, 100 /*maxiter*/);
            gmm::clear(a);
            gmm::bfgs(pps, pps, a, 10, iter, 0, 0.5);
            residual = gmm::abs(iter.get_res());
            converged = (residual < 2E-5);
          }
        }

        bool is_in = (pf_s->ref_convex(cv)->is_in(ctx.xref()) < 1E-6);

        if (is_in || (!converged && !raytrace)) {
          if (!ref_conf) {
            ctx.pf()->interpolation(ctx, coeff, y, dim_type(N));
            y += ctx.xreal();
          } else {
            y = ctx.xreal();
          }
        }

        // CRITERION 2 : The contact pair is eliminated when
        //               projection/raytrace do not converge.
        if (!converged) {
          if (!raytrace && nbwarn < 4) {
            GMM_WARNING3("Projection or raytrace algorithm did not converge "
                         "for point " << x << " residual " << residual
                         << " projection computed " << y);
            ++nbwarn;
          }
          continue;
        }

        // CRITERION 3 : The projected point is inside the element
        //               The test should be completed: If the point is outside
        //               the element, a rapid reprojection on the face
        //               (on the reference element, with a linear algorithm)
        //               can be applied and a test with a neigbhour element
        //               to decide if the point is in fact ok ...
        //               (to be done only if there is no projection on other
        //               element which coincides and with a test on the
        //               distance ... ?) To be specified (in this case,
        //               change xref).
        if (!is_in) continue;

        // CRITERION 4 : Apply the release distance
        scalar_type signed_dist = gmm::vect_dist2(y, x);
        if (signed_dist > release_distance) continue;

        // compute the unit normal vector at y and the signed distance.
        base_small_vector ny0(N);
        compute_normal(ctx, iff, ref_conf, coeff, ny0, ny, grad);
        // ny /= gmm::vect_norm2(ny); // Useful only if the unit normal is kept
        signed_dist *= gmm::sgn(gmm::vect_sp(x - y, ny));

        // CRITERION 5 : comparison with rigid obstacles
        // CRITERION 7 : smallest signed distance on contact pairs
        if (first_pair_found && contact_pairs.back().signed_dist < signed_dist)
            continue;

        // CRITERION 1 : again on found unit normal vector
        if (!(test_normal_cones_compatibility(ny, bpinfo.normals)))
            continue;

        // CRITERION 6 : for self-contact only : apply a test on
        //               unit normals in reference configuration.
        if (&m == &(mfdisp_of_boundary(ibx).linked_mesh())) {

          base_small_vector diff = bpinfo.ref_point - ctx.xreal();
          scalar_type ref_dist = gmm::vect_norm2(diff);

          if ( (ref_dist < scalar_type(4) * release_distance)
               && (gmm::vect_sp(diff, ny0) < - 0.01 * ref_dist) )
            continue;
        }

        contact_pair ct(x, nx, bpinfo, ctx.xref(), y, ny, fi, signed_dist);
        if (first_pair_found) {
          contact_pairs.back() = ct;
        } else {
          contact_pairs.push_back(ct);
          first_pair_found = true;
        }

      }
    }

    // cout << "Time for computing pairs: " << dal::uclock_sec() - time << endl; time = dal::uclock_sec();

    clear_aux_info();
  }

  //=========================================================================
  //
  //  Raytracing interpolate transformation for generic assembly
  //
  //=========================================================================
  
  class  raytracing_interpolate_transformation
    : public virtual_interpolate_transformation {

    // Structure describing a contact boundary
    struct contact_boundary {
      size_type region;            // Boundary number
      const getfem::mesh_fem *mfu; // F.e.m. for the displacement.
      std::string dispname;        // Variable name for the displacement
      mutable const model_real_plain_vector *U;      // Displacement
      mutable model_real_plain_vector U_unred; // Unreduced displacement
      bool slave;
 
      contact_boundary(void) {}
      contact_boundary(size_type r, const mesh_fem *mf, const std::string dn,
                       bool sl)
        : region(r), mfu(mf), dispname(dn), slave(sl) {}
    };

    struct face_box_info {     // Additional information for a face box
      size_type ind_boundary;  // Boundary number
      size_type ind_element;   // Element number
      short_type ind_face;     // Face number in element
      base_small_vector mean_normal;   // Mean outward normal unit vector
      face_box_info(void) {}
      face_box_info(size_type ib, size_type ie,
                    short_type iff, const base_small_vector &n)
        : ind_boundary(ib), ind_element(ie), ind_face(iff), mean_normal(n) {}
    };

    scalar_type release_distance;  // Limit distance beyond which the contact
                                   // will not be considered.
    
    std::vector<contact_boundary> contact_boundaries;
    typedef std::map<const mesh *, std::vector<size_type> > mesh_boundary_cor;
    mesh_boundary_cor boundary_for_mesh;
    
    struct obstacle {
      ga_function f;
      ga_function der_f;
      mutable base_vector x;
    };

    std::vector<obstacle> obstacles;
        
    mutable bgeot::rtree face_boxes;
    mutable std::vector<face_box_info> face_boxes_info;


    void compute_face_boxes(void) const { // called by init
      fem_precomp_pool fppool;
      base_matrix G;
      model_real_plain_vector coeff;
      face_boxes.clear();
      face_boxes_info.resize(0);

      for (size_type i = 0; i < contact_boundaries.size(); ++i) {
        const contact_boundary &cb =  contact_boundaries[i];
        if (! cb.slave) {
          size_type bnum = cb.region;
          const mesh_fem &mfu = *(cb.mfu);
          const model_real_plain_vector &U = *(cb.U);
          const mesh &m = mfu.linked_mesh();
          size_type N = m.dim();
          
          base_node val(N), bmin(N), bmax(N);
          base_small_vector n0_x(N), n_x(N), n0_y(N), n_y(N), n_mean(N);
          base_matrix grad(N,N);
          mesh_region region = m.region(bnum);
          GMM_ASSERT1(mfu.get_qdim() == N, "Wrong mesh_fem qdim");
          
          dal::bit_vector points_already_interpolated;
          std::vector<base_node> transformed_points(m.nb_max_points());
          for (getfem::mr_visitor v(region,m); !v.finished(); ++v) {
            size_type cv = v.cv();
            bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
            pfem pf_s = mfu.fem_of_element(cv);
            pfem_precomp pfp = fppool(pf_s, &(pgt->geometric_nodes()));
            slice_vector_on_basic_dof_of_element(mfu, U, cv, coeff);
            bgeot::vectors_to_base_matrix
              (G, mfu.linked_mesh().points_of_convex(cv));
            fem_interpolation_context ctx(pgt,pfp,size_type(-1), G, cv,
                                          size_type(-1));
            
            bgeot::pconvex_structure cvs = pgt->structure();
            size_type nb_pt_on_face = cvs->nb_points_of_face(v.f());
            GMM_ASSERT1(nb_pt_on_face >= 2, "This element has less than two "
                        "vertices on considered face !");
            gmm::clear(n_mean);
            
            for (size_type k = 0; k < nb_pt_on_face; ++k) {
              size_type ip = cvs->ind_points_of_face(v.f())[k];
              size_type ind = m.ind_points_of_convex(cv)[ip];
              
              // computation of transformed vertex
              if (!(points_already_interpolated.is_in(ind))) {
                ctx.set_ii(ip);
                pf_s->interpolation(ctx, coeff, val, dim_type(N));
                val += ctx.xreal();
                transformed_points[ind] = val;
                points_already_interpolated.add(ind);
              } else {
                val = transformed_points[ind];
              }
              
              if (k == 0) // computation of bounding box
                bmin = bmax = val;
              else {
                for (size_type l = 0; l < N; ++l) {
                  bmin[l] = std::min(bmin[l], val[l]);
                  bmax[l] = std::max(bmax[l], val[l]);
                }
              }
              
              // computation of unit normal vector if the vertex is on the face
              compute_normal(ctx, v.f(), false, coeff, n0_x, n_x, grad);
              n_x /= gmm::vect_norm2(n_x);
              n_mean += n_x;
            }
            
            // Security coefficient of 1.3 (for nonlinear transformations)
            scalar_type h = bmax[0] - bmin[0];
            for (size_type k = 1; k < N; ++k) h = std::max(h, bmax[k]-bmin[k]);
            for (size_type k = 0; k < N; ++k)
              { bmin[k] -= h * 0.15; bmax[k] += h * 0.15; }
            
            // Store the bounding box and additional information.
            face_boxes.add_box(bmin, bmax, face_boxes_info.size());
            n_mean /= gmm::vect_norm2(n_mean);
            face_boxes_info.push_back(face_box_info(i, cv, v.f(), n_mean));
          }
        }
      }
      
    }

  public:

    void add_rigid_obstacle(const model &md, const std::string &expr,
                            size_type N) {
      // cout << "adding rigid obstacle " << expr << endl;
      obstacles.push_back(obstacle());
      obstacles.back().f = ga_function(md, expr);
      gmm::resize(obstacles.back().x, N);
      obstacles.back().f.workspace().add_fixed_size_variable
        ("x", gmm::sub_interval(0, N), obstacles.back().x);
      obstacles.back().f.compile();
      obstacles.back().der_f = obstacles.back().f;
      obstacles.back().der_f.derivative("x");
    }

    void add_rigid_obstacle(const ga_workspace &workspace,
                            const std::string &expr, size_type N) {
      obstacles.push_back(obstacle());
      obstacles.back().f = ga_function(workspace, expr);
      gmm::resize(obstacles.back().x, N);
      obstacles.back().f.workspace().add_fixed_size_variable
        ("x", gmm::sub_interval(0, N), obstacles.back().x);
      obstacles.back().f.compile();
      obstacles.back().der_f = obstacles.back().f;
      obstacles.back().der_f.derivative("x");
    }

    void add_contact_boundary(const model &md, const mesh &m,
                              const std::string dispname,
                              size_type region, bool slave) {
      const mesh_fem *mf = 0;
      if (md.variable_group_exists(dispname)) {
        const std::vector<std::string> &t = md.variable_group(dispname);
        for (size_type i = 0; i < t.size(); ++i) {
          const mesh_fem *mf2 = md.pmesh_fem_of_variable(t[i]);
          if (mf2 && &(mf2->linked_mesh()) == &m)
            { mf = mf2; break; }
        }
      } else mf = md.pmesh_fem_of_variable(dispname);
      GMM_ASSERT1(mf, "Displacement should be a fem variable");
      contact_boundary cb(region, mf, dispname, slave);
      boundary_for_mesh[&(mf->linked_mesh())]
        .push_back(contact_boundaries.size());
      contact_boundaries.push_back(cb);
    }
    
    void add_contact_boundary(const ga_workspace &workspace, const mesh &m,
                              const std::string dispname,
                              size_type region, bool slave) {
      const mesh_fem *mf = 0;
      if (workspace.variable_group_exists(dispname)) {
        const std::vector<std::string> &t = workspace.variable_group(dispname);
        for (size_type i = 0; i < t.size(); ++i) {
          const mesh_fem *mf2 = workspace.associated_mf(t[i]);
          if (mf2 && &(mf2->linked_mesh()) == &m)
            { mf = mf2; break; }
        }
      } else mf = workspace.associated_mf(dispname);
      GMM_ASSERT1(mf, "Displacement should be a fem variable");
      contact_boundary cb(region, mf, dispname, slave);
      boundary_for_mesh[&(mf->linked_mesh())]
        .push_back(contact_boundaries.size());
      contact_boundaries.push_back(cb);
    }

    void extract_variables(const ga_workspace &workspace,
                           std::set<var_trans_pair> &vars,
                           bool ignore_data, const mesh &m_x,
                           const std::string &interpolate_name) const {
      
      bool expand_groups = !ignore_data;
      // const mesh_fem *mf = workspace.associated_mf(name);
      // GMM_ASSERT1(mf, "Internal error");
      // const mesh &m_x = mf->linked_mesh();

      mesh_boundary_cor::const_iterator it =  boundary_for_mesh.find(&m_x);
      GMM_ASSERT1(it != boundary_for_mesh.end(), "Raytracing interpolate "
                  "transformation: Mesh with no declared contact boundary");
      const std::vector<size_type> &boundaries_ind = it->second;
      for (size_type i = 0; i < boundaries_ind.size(); ++i) {
        const contact_boundary &cb =  contact_boundaries[boundaries_ind[i]];
        if (expand_groups && workspace.variable_group_exists(cb.dispname)) {
          const std::vector<std::string> &t=workspace.variable_group(cb.dispname);
          for (size_type j = 0; j < t.size(); ++j)
            vars.insert(var_trans_pair(t[j], ""));
        } else vars.insert(var_trans_pair(cb.dispname, ""));
      }

      for (size_type i = 0; i < boundaries_ind.size(); ++i) {
        const contact_boundary &cb =  contact_boundaries[boundaries_ind[i]];
        if (!(cb.slave)) {
          if (expand_groups && workspace.variable_group_exists(cb.dispname)) {
            const std::vector<std::string> &t
              = workspace.variable_group(cb.dispname);
          for (size_type j = 0; j < t.size(); ++j)
            vars.insert(var_trans_pair(t[j], interpolate_name));
          } else vars.insert(var_trans_pair(cb.dispname, interpolate_name));
        }
      }
    }

    void init(const ga_workspace &workspace) const {
      for (size_type i = 0; i < contact_boundaries.size(); ++i) {
        const contact_boundary &cb =  contact_boundaries[i];
        const mesh_fem &mfu = *(cb.mfu);
        if (mfu.is_reduced()) {
          gmm::resize(cb.U_unred, mfu.nb_basic_dof());
          mfu.extend_vector(workspace.value(cb.dispname), cb.U_unred);
          cb.U = &(cb.U_unred);
        } else {
          cb.U = &(workspace.value(cb.dispname));
        }
      }
      compute_face_boxes();
    };

    void finalize(void) const {
      face_boxes.clear();
      face_boxes_info = std::vector<face_box_info>();
      for (size_type i = 0; i < contact_boundaries.size(); ++i)
        contact_boundaries[i].U_unred = model_real_plain_vector();
    }

    int transform(const ga_workspace &/*workspace*/, const mesh &m_x,
                  fem_interpolation_context &ctx_x,
                  const base_small_vector &/*Normal*/,
                  const mesh **m_t,
                  size_type &cv, size_type &face_num, base_node &P_ref,
                  base_small_vector &N_y,
                  std::map<var_trans_pair, base_tensor> &derivatives,
                  bool compute_derivatives) const {
      size_type cv_x = ctx_x.convex_num();
      size_type face_x = ctx_x.face_num();
      GMM_ASSERT1(face_x != size_type(-1), "The contact transformation can "
                  "only be applied to a boundary");

      //
      // Find the right (slave) contact boundary
      //
      mesh_boundary_cor::const_iterator it =  boundary_for_mesh.find(&m_x);
      GMM_ASSERT1(it != boundary_for_mesh.end(),
                  "Mesh with no declared contact boundary");
      const std::vector<size_type> &boundaries_ind = it->second;
      size_type ib_x = size_type(-1);
      for (size_type i = 0; i < boundaries_ind.size(); ++i) {
        const contact_boundary &cb =  contact_boundaries[boundaries_ind[i]];
        if (m_x.region(cb.region).is_in(cv_x, face_x))
          { ib_x = boundaries_ind[i]; break; }
      }
      GMM_ASSERT1(ib_x != size_type(-1),
                  "No contact region found for this point");
      const contact_boundary &cb_x =  contact_boundaries[ib_x];
      const mesh_fem &mfu_x = *(cb_x.mfu); 
      pfem pfu_x = mfu_x.fem_of_element(cv_x);
      size_type N = mfu_x.linked_mesh().dim();
      GMM_ASSERT1(mfu_x.get_qdim() == N,
                  "Displacment field with wrong dimension");
      
      model_real_plain_vector coeff_x, coeff_y, stored_coeff_y;
      base_small_vector a(N-1), b(N-1), pt_x(N), pt_y(N), n_x(N);
      base_small_vector stored_pt_y(N), stored_n_y(N), stored_pt_y_ref(N);
      base_small_vector n0_x, n_y(N), n0_y, res(N-1), res2(N-1), dir(N-1);
      base_matrix G_x, G_y, grad(N,N), hessa(N-1, N-1);
      std::vector<base_small_vector> ti(N-1), Ti(N-1);
      scalar_type stored_signed_distance(0);
      std::string stored_dispname;
      scalar_type d0 = 1E300, d1, d2;
      const mesh *stored_m_y(0);
      size_type stored_cv_y(-1), stored_face_y(-1);
      fem_interpolation_context stored_ctx_y;

      //
      // Computation of the deformed point and unit normal vectors
      //
      slice_vector_on_basic_dof_of_element(mfu_x, *(cb_x.U), cv_x, coeff_x);
      bgeot::vectors_to_base_matrix(G_x, m_x.points_of_convex(cv_x));
      ctx_x.set_pf(pfu_x);
      pfu_x->interpolation(ctx_x, coeff_x, pt_x, dim_type(N));
      pt_x += ctx_x.xreal();
      compute_normal(ctx_x, face_x, false, coeff_x, n0_x, n_x, grad);
      n_x /= gmm::vect_norm2(n_x);

      //
      //  Determine the nearest rigid obstacle, taking into account
      //  the release distance.
      //

      bool first_pair_found = false;
      size_type irigid_obstacle(-1);
      for (size_type i = 0; i < obstacles.size(); ++i) {
        const obstacle &obs = obstacles[i];
        gmm::copy(pt_x, obs.x);
        const base_tensor &t = obs.f.eval();
        
        GMM_ASSERT1(t.size() == 1, "Obstacle level set function as to be "
                    "a scalar valued one");
        d1 = t[0];
        // cout << "d1 = " << d1 << endl;
        if (gmm::abs(d1) < release_distance && d1 < d0) {
          const base_tensor &t_der = obs.der_f.eval();
          // cout << "t_der.as_vector() = " << t_der.as_vector() << endl;
          if (gmm::vect_sp(t_der.as_vector(), n_x) < scalar_type(0)) 
            { d0 = d1; irigid_obstacle = i; gmm::copy(t_der.as_vector(),n_y); }
        }
      }
      
      if (irigid_obstacle != size_type(-1)) {
        // cout << "Testing obstacle " << irigid_obstacle << endl;
        const obstacle &obs = obstacles[irigid_obstacle];
        gmm::copy(pt_x, obs.x);
        gmm::copy(pt_x, pt_y);
        size_type nit = 0, nb_fail = 0;
        scalar_type alpha(0), beta(0);
        d1 = d0;
        
        while (gmm::abs(d1) > 1E-13 && ++nit < 50 && nb_fail < 3) {
          if (nit != 1) gmm::copy(obs.der_f.eval().as_vector(), n_y);

          for (scalar_type lambda(1); lambda >= 1E-3; lambda/=scalar_type(2)) {
            alpha = beta - lambda * d1 / gmm::vect_sp(n_y, n_x);
            gmm::add(pt_x, gmm::scaled(n_x, alpha), obs.x);
            d2 = obs.f.eval()[0];
            if (gmm::abs(d2) < gmm::abs(d1)) break;
          }
          if (gmm::abs(beta - d1 / gmm::vect_sp(n_y, n_x)) > scalar_type(500))
            nb_fail++;
          beta = alpha; d1 = d2;
        }
        gmm::copy(obs.x, pt_y);

        if (gmm::abs(d1) > 1E-8) {
           GMM_WARNING1("Raytrace on rigid obstacle failed");
        } // CRITERION 4 for rigid bodies : Apply the release distance
        else if (gmm::vect_dist2(pt_y, pt_x) <= release_distance) {
          n_y /= gmm::vect_norm2(n_y);
          d0 = gmm::vect_dist2(pt_y, pt_x) * gmm::sgn(d0);
          stored_pt_y = stored_pt_y_ref = pt_y; stored_n_y = n_y, 
          stored_signed_distance = d0;
          first_pair_found = true;
        } else irigid_obstacle = size_type(-1);
      }
      
      //
      // Determine the potential contact pairs with deformable bodies
      //
      bgeot::rtree::pbox_set bset;
      base_node bmin(pt_x), bmax(pt_x);
      for (size_type i = 0; i < N; ++i)
        { bmin[i] -= release_distance; bmax[i] += release_distance; }

      face_boxes.find_line_intersecting_boxes(pt_x, n_x, bmin, bmax, bset);

      //
      // Iteration on potential contact pairs and application
      // of selection criteria
      //
      bgeot::rtree::pbox_set::iterator it_cp = bset.begin();
      for (; it_cp != bset.end(); ++it_cp) {
        face_box_info &fbox_y = face_boxes_info[(*it_cp)->id];
        size_type ib_y = fbox_y.ind_boundary;
        const contact_boundary &cb_y =  contact_boundaries[ib_y];
        const mesh_fem &mfu_y = *(cb_y.mfu);
        const mesh &m_y = mfu_y.linked_mesh();
        size_type cv_y = fbox_y.ind_element;
        pfem pfu_y = mfu_y.fem_of_element(cv_y);
        size_type face_y = fbox_y.ind_face;
        bgeot::pgeometric_trans pgt_y= m_y.trans_of_convex(cv_y);

        // CRITERION 1 : The unit normal vector are compatible
        //               and the two points are not in the same element.
        if (gmm::vect_sp(fbox_y.mean_normal, n_x) >= scalar_type(0) ||
            (cv_x == cv_y && &m_x == &m_y))
          continue;

        //
        // Raytrace search for y by Newton's algorithm
        //

        bgeot::vectors_to_base_matrix(G_y, m_y.points_of_convex(cv_y));
        const base_node &Y0
          = pfu_y->ref_convex(cv_y)->points_of_face(short_type(face_y))[0];
        fem_interpolation_context ctx_y(pgt_y, pfu_y, Y0, G_y, cv_y, face_y);
        
        const base_small_vector &NY0
          = pfu_y->ref_convex(cv_y)->normals()[face_y];
        for (size_type k = 0; k < N-1; ++k) { // A basis for the face
          gmm::resize(ti[k], N);
          scalar_type norm(0);
          while(norm < 1E-5) {
            gmm::fill_random(ti[k]);
            ti[k] -= gmm::vect_sp(ti[k], NY0) * NY0;
            for (size_type l = 0; l < k; ++l)
              ti[k] -= gmm::vect_sp(ti[k], ti[l]) * ti[l];
            norm = gmm::vect_norm2(ti[k]);
          }
          ti[k] /= norm;
        }

        gmm::clear(a);
        
        for (size_type k = 0; k < N-1; ++k) {
          gmm::resize(Ti[k], N);
          scalar_type norm(0);
          while (norm < 1E-5) {
            gmm::fill_random(Ti[k]);
            Ti[k] -= gmm::vect_sp(Ti[k], n_x) * n_x;
            for (size_type l = 0; l < k; ++l)
              Ti[k] -= gmm::vect_sp(Ti[k], Ti[l]) * Ti[l];
            norm = gmm::vect_norm2(Ti[k]);
          }
          Ti[k] /= norm;
        }

        slice_vector_on_basic_dof_of_element(mfu_y, *(cb_y.U), cv_y, coeff_y);

        raytrace_pt_surf_cost_function_object pps(Y0, pt_x, ctx_y, coeff_y,
                                                  ti, Ti, false);
        pps(a, res);
        scalar_type residual = gmm::vect_norm2(res);
        scalar_type residual2(0), det(0);
        bool exited = false;
        size_type nbfail = 0, niter = 0;
        for (;residual > 2E-12 && niter <= 30; ++niter) {
          
          for (size_type subiter(0);;) {
            pps(a, hessa);
            det = gmm::abs(gmm::lu_inverse(hessa, false));
            if (det > 1E-15) break;
            for (size_type i = 0; i < N-1; ++i)
              a[i] += gmm::random() * 1E-7;
            if (++subiter > 4) break;
          }
          if (det <= 1E-15) break;
          // Computation of the descent direction
          gmm::mult(hessa, gmm::scaled(res, scalar_type(-1)), dir);
          
          if (gmm::vect_norm2(dir) > scalar_type(10)) nbfail++;
          if (nbfail >= 4) break;
          
          // Line search
          scalar_type lambda(1);
          for (size_type j = 0; j < 5; ++j) {
            gmm::add(a, gmm::scaled(dir, lambda), b);
            pps(b, res2);
            residual2 = gmm::vect_norm2(res2);
            if (residual2 < residual) break;
            lambda /= ((j < 3) ? scalar_type(2) : scalar_type(5));
          }
          
          residual = residual2;
          gmm::copy(res2, res);
          gmm::copy(b, a);
          scalar_type dist_ref = gmm::vect_norm2(a);
          
          if (niter > 1 && dist_ref > 15) break;
          if (niter > 5 && dist_ref > 8) break;
          if ((niter > 1 && dist_ref > 7) || nbfail == 3) exited = true;
        }
        bool converged = (gmm::vect_norm2(res) < 2E-6);
        bool is_in = (pfu_y->ref_convex(cv_y)->is_in(ctx_y.xref()) < 1E-6);
        GMM_ASSERT1(!(exited && converged && is_in),
                    "A non conformal case !! " << gmm::vect_norm2(res)
                    << " : " << nbfail << " : " << niter);
        
        if (is_in) {
          ctx_y.pf()->interpolation(ctx_y, coeff_y, pt_y, dim_type(N));
          pt_y += ctx_y.xreal();
        }

        // CRITERION 2 : The contact pair is eliminated when
        //               raytrace do not converge.
        if (!converged) continue;
        
        // CRITERION 3 : The raytraced point is inside the element
        if (!is_in) continue;

        // CRITERION 4 : Apply the release distance
        scalar_type signed_dist = gmm::vect_dist2(pt_y, pt_x);
        if (signed_dist > release_distance) continue;
        
        // compute the unit normal vector at y and the signed distance.
        compute_normal(ctx_y, face_y, false, coeff_y, n0_y, n_y, grad);
        n_y /= gmm::vect_norm2(n_y);
        signed_dist *= gmm::sgn(gmm::vect_sp(pt_x - pt_y, n_y));

        // CRITERION 5 : comparison with rigid obstacles
        // CRITERION 7 : smallest signed distance on contact pairs
        if (first_pair_found && stored_signed_distance < signed_dist)
          continue;

        // CRITERION 1 : again on found unit normal vector
        if (gmm::vect_sp(n_y, n_x) >= scalar_type(0)) continue;

        // CRITERION 6 : for self-contact only : apply a test on
        //               unit normals in reference configuration.
        if (&m_x == &m_y) {
          base_small_vector diff = ctx_x.xreal() - ctx_y.xreal();
          scalar_type ref_dist = gmm::vect_norm2(diff);
          if ( (ref_dist < scalar_type(4) * release_distance)
               && (gmm::vect_sp(diff, n0_y) < - 0.01 * ref_dist) )
            continue;
        }

        stored_pt_y = pt_y; stored_pt_y_ref = ctx_y.xref();
        stored_m_y = &m_y; stored_cv_y = cv_y; stored_face_y = face_y;
        stored_n_y = n_y;
        stored_ctx_y = ctx_y;
        stored_coeff_y = coeff_y;
        stored_signed_distance = signed_dist;
        stored_dispname = cb_y.dispname;
        first_pair_found = true;
        irigid_obstacle = size_type(-1);
      }

      int ret_type = 0;

      if (irigid_obstacle != size_type(-1)) {
        *m_t = 0; cv = face_num = size_type(-1);
        P_ref = stored_pt_y; N_y = stored_n_y;
        ret_type = 2;
      } else if (first_pair_found) {
        *m_t = stored_m_y; cv = stored_cv_y; face_num = stored_face_y;
        P_ref = stored_pt_y_ref;
        ret_type = 1;
      }

      // Note on derivatives of the transformation : for efficiency and
      // simplicity reasons, the derivative should be computed with
      // the value of corresponding test functions. This means that
      // for a transformation F(u) the conputed derivative is F'(u).Test_u
      // including the Test_u.
      if (compute_derivatives) {
        if (ret_type >= 1) {
          fem_interpolation_context &ctx_y = stored_ctx_y;
          size_type cv_y = 0;
          if (ret_type == 1) cv_y = ctx_y.convex_num();
          
          base_matrix I_nxny(N,N); // I - nx@ny/nx.ny
          gmm::copy(gmm::identity_matrix(), I_nxny);
          gmm::rank_one_update(I_nxny, n_x,
                               gmm::scaled(stored_n_y,scalar_type(-1)
                                           / gmm::vect_sp(n_x, stored_n_y)));
        
          // cout << "n_y = " << stored_n_y << " pt_x = " << pt_x << " pt_y = " << stored_pt_y << endl;

          // Computation of F_y
          base_matrix F_y(N,N), F_y_inv(N,N), M1(N, N), M2(N, N);
          pfem pfu_y = 0;
          if (ret_type == 1) {
            pfu_y = ctx_y.pf();
            pfu_y->interpolation_grad(ctx_y, stored_coeff_y, F_y, dim_type(N));
            gmm::add(gmm::identity_matrix(), F_y);
            gmm::copy(F_y, F_y_inv);
            gmm::lu_inverse(F_y_inv);
          } else {
            gmm::copy(gmm::identity_matrix(), F_y);
            gmm::copy(gmm::identity_matrix(), F_y_inv);
          }

          // Computation of F_x
          base_matrix F_x(N,N), F_x_inv(N,N);
          pfu_x->interpolation_grad(ctx_x, coeff_x, F_x, dim_type(N));
          gmm::add(gmm::identity_matrix(), F_x);
          gmm::copy(F_x, F_x_inv);
          gmm::lu_inverse(F_x_inv);
        

          base_tensor base_ux;
          base_matrix vbase_ux;
          ctx_x.base_value(base_ux);
          size_type qdim_ux = pfu_x->target_dim();
          size_type ndof_ux = pfu_x->nb_dof(cv_x) * N / qdim_ux;
          vectorize_base_tensor(base_ux, vbase_ux, ndof_ux, qdim_ux, N);
          
          base_tensor base_uy;
          base_matrix vbase_uy;
          size_type ndof_uy = 0;
          if (ret_type == 1) {
            ctx_y.base_value(base_uy);
            size_type qdim_uy = pfu_y->target_dim();
            ndof_uy = pfu_y->nb_dof(cv_y) * N / qdim_uy;
            vectorize_base_tensor(base_uy, vbase_uy, ndof_uy, qdim_uy, N);
          }
          
          base_tensor grad_base_ux, vgrad_base_ux;
          ctx_x.grad_base_value(grad_base_ux);
          vectorize_grad_base_tensor(grad_base_ux, vgrad_base_ux, ndof_ux,
                                     qdim_ux, N);

          // Derivative : F_y^{-1}*I_nxny*(Test_u(X)-Test_u(Y)+gDn_x[Test_u])
          //         with Dn_x[Test_u] =-(I-nx@nx)*F_x^{-T}*Grad_Test_u^{T}*n_x
          //         and I_nxny*(I - nx@nx) = I_nxny
          
          // F_y^{-1}*I_nxny*Test_u(X)
          gmm::mult(F_y_inv, I_nxny, M1);
          base_matrix der_x(ndof_ux, N);
          gmm::mult(vbase_ux, gmm::transposed(M1), der_x);
          
          //         for (size_type i = 0; i < ndof_ux; ++i)
          //           for (size_type j = 0; j < N; ++j)
          //             for (size_type k = 0; k < N; ++k) 
          //               der_x(i, j) += M1(j,k) * vbase_ux(i,k);

          // -F_y^{-1}*I_nxny*Test_u(Y)
          base_matrix der_y(ndof_uy, N);
          if (ret_type == 1) {
            gmm::mult(vbase_uy, gmm::transposed(M1), der_y);
            gmm::scale(der_y, scalar_type(-1));
          }
          
          //         for (size_type i = 0; i < ndof_uy; ++i)
          //           for (size_type j = 0; j < N; ++j)
          //             for (size_type k = 0; k < N; ++k) 
          //               der_y(i, j) -= M1(j,k) * vbase_uy(i,k);
          
          // F_y^{-1}*I_nxny*gDn_x[Test_u]
          gmm::mult(M1, gmm::transposed(F_x_inv), M2);
          for (size_type i = 0; i < ndof_ux; ++i)
            for (size_type j = 0; j < N; ++j)
              for (size_type k = 0; k < N; ++k)
                for (size_type l = 0; l < N; ++l)
                  der_x(i, j) -= M2(j, k) * vgrad_base_ux(i, l, k)
                    * n_x[l] * stored_signed_distance;

          for (std::map<var_trans_pair, base_tensor>::iterator itd
                 = derivatives.begin(); itd != derivatives.end(); ++itd) {
            if (cb_x.dispname.compare(itd->first.first) == 0 &&
                itd->first.second.size() == 0) {
              itd->second.adjust_sizes(ndof_ux, N);
              gmm::copy(der_x.as_vector(), itd->second.as_vector());
            } else if (ret_type == 1 &&
                       stored_dispname.compare(itd->first.first) == 0 &&
                       itd->first.second.size() != 0) {
              itd->second.adjust_sizes(ndof_uy, N);
              gmm::copy(der_y.as_vector(), itd->second.as_vector());
            } else itd->second.adjust_sizes(0, 0);
          }
        } else {
          for (std::map<var_trans_pair, base_tensor>::iterator itd
                 = derivatives.begin(); itd != derivatives.end(); ++itd)
            itd->second.adjust_sizes(0, 0);
        }
      }

      return ret_type;
    }
    
    raytracing_interpolate_transformation(scalar_type d)
      : release_distance(d) {}
  };

  void add_raytracing_transformation
  (model &md, const std::string &transname, scalar_type d) {
    pinterpolate_transformation p = new raytracing_interpolate_transformation(d);
    md.add_interpolate_transformation(transname, p);
  }

  void add_raytracing_transformation
  (ga_workspace &workspace, const std::string &transname, scalar_type d) {
    pinterpolate_transformation p = new raytracing_interpolate_transformation(d);
    workspace.add_interpolate_transformation(transname, p);
  }

  void add_master_contact_boundary_to_raytracing_transformation
  (model &md, const std::string &transname, const mesh &m,
   const std::string &dispname, size_type region) {
    raytracing_interpolate_transformation *p
      = dynamic_cast<raytracing_interpolate_transformation *>
      (const_cast<virtual_interpolate_transformation *>
       (&(*(md.interpolate_transformation(transname)))));
    p->add_contact_boundary(md, m, dispname, region, false);
  }

  void add_slave_contact_boundary_to_raytracing_transformation
  (model &md, const std::string &transname, const mesh &m,
   const std::string &dispname, size_type region) {
    raytracing_interpolate_transformation *p
      = dynamic_cast<raytracing_interpolate_transformation *>
      (const_cast<virtual_interpolate_transformation *>
       (&(*(md.interpolate_transformation(transname)))));
    p->add_contact_boundary(md, m, dispname, region, true);
  }

  void add_master_contact_boundary_to_raytracing_transformation
  (ga_workspace &workspace, const std::string &transname, const mesh &m,
   const std::string &dispname, size_type region) {
    raytracing_interpolate_transformation *p
      = dynamic_cast<raytracing_interpolate_transformation *>
      (const_cast<virtual_interpolate_transformation *>
       (&(*(workspace.interpolate_transformation(transname)))));
    p->add_contact_boundary(workspace, m, dispname, region, false);
  }

  void add_slave_contact_boundary_to_raytracing_transformation
  (ga_workspace &workspace, const std::string &transname, const mesh &m,
   const std::string &dispname, size_type region) {
    raytracing_interpolate_transformation *p
      = dynamic_cast<raytracing_interpolate_transformation *>
      (const_cast<virtual_interpolate_transformation *>
       (&(*(workspace.interpolate_transformation(transname)))));
    p->add_contact_boundary(workspace, m, dispname, region, true);
  }

  void add_rigid_obstacle_to_raytracing_transformation
  (model &md, const std::string &transname,
   const std::string &expr, size_type N) {
    raytracing_interpolate_transformation *p
      = dynamic_cast<raytracing_interpolate_transformation *>
      (const_cast<virtual_interpolate_transformation *>
       (&(*(md.interpolate_transformation(transname)))));
    p->add_rigid_obstacle(md, expr, N);
  }

  void add_rigid_obstacle_to_raytracing_transformation
  (ga_workspace &workspace, const std::string &transname,
   const std::string &expr, size_type N) {
    raytracing_interpolate_transformation *p
      = dynamic_cast<raytracing_interpolate_transformation *>
      (const_cast<virtual_interpolate_transformation *>
       (&(*(workspace.interpolate_transformation(transname)))));
    p->add_rigid_obstacle(workspace, expr, N);
  }



  //=========================================================================
  //
  //  Specific nonlinear operator of the high-level generic assembly langage
  //  dedicated to contact/friction
  //
  //=========================================================================

  // static void ga_init_scalar(bgeot::multi_index &mi) { mi.resize(0); }
  static void ga_init_vector(bgeot::multi_index &mi, size_type N)
  { mi.resize(1); mi[0] = N; }
  // static void ga_init_square_matrix(bgeot::multi_index &mi, size_type N)
  // { mi.resize(2); mi[0] = mi[1] = N; }
  

  // Transformed_unit_vector(Grad_u, n)  = (I+Grad_u)^{-T}n / ||(I+Grad_u)^{-T}n||
  struct Transformed_unit_vector : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 2 || args[0]->sizes().size() != 2
          || args[1]->size() != args[0]->sizes()[0]
          || args[0]->sizes()[0] != args[0]->sizes()[1]) return false;
      ga_init_vector(sizes, args[0]->sizes()[0]);
      return true;
    }
    
    // Value : (I+Grad_u)^{-T}n / ||(I+Grad_u)^{-T}n||
    void value(const arg_list &args, base_tensor &result) const {
      size_type N = args[0]->sizes()[0];
      base_matrix F(N, N);
      gmm::copy(args[0]->as_vector(), F.as_vector());
      gmm::add(gmm::identity_matrix(), F);
      gmm::lu_inverse(F);
      gmm::mult(gmm::transposed(F), args[1]->as_vector(), result.as_vector());
      gmm::scale(result.as_vector(),
                 scalar_type(1)/gmm::vect_norm2(result.as_vector()));
    }

    // Derivative / Grad_u: -(I - n@n)(I+Grad_u)^{-T}Test_Grad_u n
    // Implementation: A{ijk} = -G{ik}ndef{j}
    //                 with G = (I - n@n)(I+Grad_u)^{-T}
    //                 and ndef the transformed normal         
    // Derivative / n: ((I+Grad_u)^{-T}Test_n - ndef(ndef.Test_n))/||(I+Grad_u)^{-T}n||
    // Implementation: A{ij} = (F{ij} - ndef{i}ndef{j})/norm_ndef
    //                 with F = (I+Grad_u)^{-1}
    void derivative(const arg_list &args, size_type nder,
                    base_tensor &result) const {
      size_type N = args[0]->sizes()[0];
      base_matrix F(N, N), G(N, N);
      base_small_vector ndef(N), aux(N);
      gmm::copy(args[0]->as_vector(), F.as_vector());
      gmm::add(gmm::identity_matrix(), F);
      gmm::lu_inverse(F);
      gmm::mult(gmm::transposed(F), args[1]->as_vector(), ndef);
      scalar_type norm_ndef = gmm::vect_norm2(ndef);
      gmm::scale(ndef, scalar_type(1)/norm_ndef);
      gmm::copy(gmm::transposed(F), G);
      gmm::mult(F, ndef, aux);
      gmm::rank_one_update(G, gmm::scaled(ndef, scalar_type(-1)), aux);
      base_tensor::iterator it = result.begin();
      switch (nder) {
      case 1:
        for (size_type k = 0; k < N; ++k)
          for (size_type j = 0; j < N; ++j)
            for (size_type i = 0; i < N; ++i, ++it)
              *it = -G(i, k) * ndef[j];
        break;
      case 2:
        for (size_type j = 0; j < N; ++j)
          for (size_type i = 0; i < N; ++i, ++it)
            *it = (F(j,i) - ndef[i]*ndef[j])/norm_ndef;
        break;
      default: GMM_ASSERT1(false, "Internal error");
      }
      GMM_ASSERT1(it == result.end(), "Internal error");
    }
    
    // Second derivative : not implemented
    void second_derivative(const arg_list &, size_type, size_type,
                           base_tensor &) const {
      GMM_ASSERT1(false, "Sorry, second derivative not implemented");
    }
  };


  // Coulomb_friction_coupled_projection(lambda, n, Vs, g, f, r)
  struct Coulomb_friction_coupled_projection : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 6 || args[1]->size() != args[0]->size()
          || args[2]->size() != args[0]->size()
          || args[3]->size() != 1 || args[4]->size() > 3 || args[4]->size() == 0
          || args[5]->size() != 1 ) return false;
      ga_init_vector(sizes, args[0]->sizes()[0]);
      return true;
    }
    
    // Value : (lambda.n+rg)_- n - P_B(n, f(lambda.n+rg)_-)(lambda-r Vs)
    void value(const arg_list &args, base_tensor &result) const {
      const base_vector &lambda = *(args[0]);
      const base_vector &n = *(args[1]);
      const base_vector &Vs = *(args[2]);
      base_vector &F = result;
      scalar_type g = (*(args[3]))[0];
      const base_vector &f = *(args[4]);
      scalar_type r = (*(args[5]))[0];
      

      scalar_type nn = gmm::vect_norm2(n);
      scalar_type lambdan = gmm::vect_sp(lambda, n)/nn;
      scalar_type lambdan_aug = gmm::neg(lambdan + r * g);
      size_type s_f = gmm::vect_size(f);
      scalar_type tau = ((s_f >= 3) ? f[2] : scalar_type(0)) + f[0]*lambdan_aug;
      if (s_f >= 2) tau = std::min(tau, f[1]);
    
      if (tau > scalar_type(0)) {
        gmm::add(lambda, gmm::scaled(Vs, -r), F);
        scalar_type mu = gmm::vect_sp(F, n)/nn;
        gmm::add(gmm::scaled(n, -mu/nn), F);
        scalar_type norm = gmm::vect_norm2(F);
        if (norm > tau) gmm::scale(F, tau / norm);
      } else { gmm::clear(F); }
      
      gmm::add(gmm::scaled(n, -lambdan_aug/nn), F);
    }

    // Derivative / Grad_u: -(I - n@n)(I+Grad_u)^{-T}Test_Grad_u n
    // Implementation: A{ijk} = -G{kj}ndef{i}
    //                 with G = (I - n@n)(I+Grad_u)^{-T}
    //                 and ndef the transformed normal         
    // Derivative / n: ((I+Grad_u)^{-T}Test_n - ndef(ndef.Test_n))/||(I+Grad_u)^{-T}n||
    // Implementation: A{ij} = (F{ij} - ndef{i}ndef{j})/norm_ndef
    //                 with F = (I+Grad_u)^{-1}
    void derivative(const arg_list &args, size_type nder,
                    base_tensor &result) const { // Can be optimized ?
      size_type N = args[0]->size();
      const base_vector &lambda = *(args[0]);
      const base_vector &n = *(args[1]);
      const base_vector &Vs = *(args[2]);
      base_vector F(N), dg(N);
      base_matrix dVs(N,N), dn(N,N);
      scalar_type g = (*(args[3]))[0];
      const base_vector &f = *(args[4]);
      scalar_type r = (*(args[5]))[0];

      scalar_type nn = gmm::vect_norm2(n);
      scalar_type lambdan = gmm::vect_sp(lambda, n)/nn;
      scalar_type lambdan_aug = gmm::neg(lambdan + r * g);
      size_type s_f = gmm::vect_size(f);
      scalar_type tau = ((s_f >= 3) ? f[2] : scalar_type(0)) + f[0]*lambdan_aug;
      if (s_f >= 2) tau = std::min(tau, f[1]);
      scalar_type norm(0);
      
      if (tau > scalar_type(0)) {
        gmm::add(lambda, gmm::scaled(Vs, -r), F);
        scalar_type mu = gmm::vect_sp(F, n)/nn;
        gmm::add(gmm::scaled(n, -mu/nn), F);
        norm = gmm::vect_norm2(F);
        gmm::copy(gmm::identity_matrix(), dn);
        gmm::scale(dn, -mu/nn);
        gmm::rank_one_update(dn, gmm::scaled(n, mu/(nn*nn*nn)), n);
        gmm::rank_one_update(dn, gmm::scaled(n, scalar_type(-1)/(nn*nn)), F);
        gmm::copy(gmm::identity_matrix(), dVs);
        gmm::rank_one_update(dVs, n, gmm::scaled(n, scalar_type(-1)/(nn*nn)));
        
        if (norm > tau) {
          gmm::rank_one_update(dVs, F,
                               gmm::scaled(F, scalar_type(-1)/(norm*norm)));
          gmm::scale(dVs, tau / norm);
          gmm::copy(gmm::scaled(F, scalar_type(1)/norm), dg);
          gmm::rank_one_update(dn, gmm::scaled(F, mu/(norm*norm*nn)), F);
          gmm::scale(dn, tau / norm);
          gmm::scale(F, tau / norm);
        } // else gmm::clear(dg);
        
      } // else { gmm::clear(dg); gmm::clear(dVs); gmm::clear(F); gmm::clear(dn); }
      // At this stage, F = P_{B_T}, dVs = d_q P_{B_T}, dn = d_n P_{B_T}
      // and dg = d_tau P_{B_T}.

      
      base_tensor::iterator it = result.begin();
      switch (nder) {
      case 1: // Derivative with respect to lambda
        if (norm > tau && ((s_f <= 1) || tau < f[1]) && ((s_f <= 2) || tau > f[2]))
          gmm::rank_one_update(dVs, dg, gmm::scaled(n, -f[0]/nn));
        if (lambdan_aug > scalar_type(0))
          gmm::rank_one_update(dVs, n, gmm::scaled(n, scalar_type(1)/(nn*nn)));
        for (size_type j = 0; j < N; ++j)
          for (size_type i = 0; i < N; ++i, ++it)
            *it = dVs(i, j);
        break;
      case 2: // Derivative with respect to n
        if (norm > tau && ((s_f <= 1) || tau < f[1]) && ((s_f <= 2) || tau > f[2])) {
          gmm::rank_one_update(dn, dg, gmm::scaled(lambda, -f[0]/nn));
          gmm::rank_one_update(dn, dg, gmm::scaled(n, f[0]*lambdan/(nn*nn)));
        }
        if (lambdan_aug > scalar_type(0)) {
          gmm::rank_one_update(dn, gmm::scaled(n, scalar_type(1)/(nn*nn)), lambda);
          gmm::rank_one_update(dn,
                               gmm::scaled(n,(lambdan_aug-lambdan)/(nn*nn*nn)), n);
          for (size_type j = 0; j < N; ++j) dn(j,j) -= lambdan_aug/nn;
        }
        for (size_type j = 0; j < N; ++j)
          for (size_type i = 0; i < N; ++i, ++it)
            *it = dn(i, j);
        break;
      case 3:
        gmm::scale(dVs, -r);
        for (size_type j = 0; j < N; ++j)
          for (size_type i = 0; i < N; ++i, ++it)
            *it = dVs(i, j);
        break;
      case 4:
         if (norm > tau && ((s_f <= 1) || tau < f[1]) && ((s_f <= 2) || tau > f[2]))
           gmm::scale(dg, -f[0]*r);
         else
           gmm::clear(dg);
         if (lambdan_aug > scalar_type(0))
           gmm::add(gmm::scaled(n, r/nn), dg);
         for (size_type i = 0; i < N; ++i, ++it)
           *it = dg[i];
        break;
      case 5:
        if (norm > tau && ((s_f <= 1) || tau < f[1]) && ((s_f <= 2) || tau > f[2]))
          gmm::scale(dg, -f[0]*g);
        else
          gmm::clear(dg);
        gmm::mult_add(dVs, gmm::scaled(Vs, scalar_type(-1)), dg);
        if (lambdan_aug > scalar_type(0))
          gmm::add(gmm::scaled(n, g/nn), dg);
        for (size_type i = 0; i < N; ++i, ++it)
          *it = dg[i];
        break;
      case 6:
        base_small_vector dtau_df(s_f);
        if ((s_f <= 1) || tau < f[1]) dtau_df[0] = lambdan_aug;
        if (s_f >= 2 && tau == f[1]) dtau_df[1] = 1;
        if (s_f >= 3 && tau < f[1]) dtau_df[2] = 1;
        for (size_type j = 0; j < s_f; ++j)
          for (size_type i = 0; i < N; ++i, ++it)
            *it = dg[i] * dtau_df[j];
        break;
      }
      GMM_ASSERT1(it == result.end(), "Internal error");
    }
    
    // Second derivative : not implemented
    void second_derivative(const arg_list &, size_type, size_type,
                           base_tensor &) const {
      GMM_ASSERT1(false, "Sorry, second derivative not implemented");
    }
  };

  static bool init_predef_operators(void) {

    ga_predef_operator_tab &PREDEF_OPERATORS
      = dal::singleton<ga_predef_operator_tab>::instance();
    
    PREDEF_OPERATORS.add_method("Transformed_unit_vector",
                                new Transformed_unit_vector());
    PREDEF_OPERATORS.add_method("Coulomb_friction_coupled_projection",
                                new Coulomb_friction_coupled_projection());

    return true;
   }

  static bool predef_operators_initialized = init_predef_operators();


}  /* end of namespace getfem.                                             */
