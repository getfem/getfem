/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2011-2013 Yves Renard, Konstantinos Poulios.

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


namespace getfem {


  //=========================================================================
  //
  //  Structure which store the contact boundaries, rigid obstacles and
  //  computes the contact pairs in large sliding/large deformation
  //
  //=========================================================================

  size_type multi_contact_frame::add_U(const model_real_plain_vector &U) {
    size_type i = 0;
    for (; i < Us.size(); ++i) if (Us[i] == &U) return i;
    Us.push_back(&U);
    ext_Us.resize(Us.size());
    return i;
  }
  
  void multi_contact_frame::extend_vectors(void) {
    dal::bit_vector iU;
    for (size_type i = 0; i < contact_boundaries.size(); ++i) {
      size_type ind_U = contact_boundaries[i].ind_U;
      if (!(iU[ind_U])) {
        const mesh_fem &mf = *(contact_boundaries[i].mfu);
        gmm::resize(ext_Us[ind_U], mf.nb_basic_dof());
        mf.extend_vector(*(Us[ind_U]), ext_Us[ind_U]);
        iU.add(ind_U);
      }
    }
  }
  
  void multi_contact_frame::normal_cone_simplicication(void) {
    if (fem_nodes_mode) {
      scalar_type threshold = ::cos(cut_angle);
      for (size_type i = 0; i < boundary_points_info.size(); ++i) {
        normal_cone &nc = boundary_points_info[i].normals;
        if (nc.size() > 1) {
          base_small_vector n_mean = nc[0];
          for (size_type j = 1; j < nc.size(); ++j) n_mean += nc[j];
          n_mean /= scalar_type(nc.size());
          bool reduce = true;
          for (size_type j = 0; j < nc.size(); ++j)
            if (gmm::vect_sp(n_mean, nc[j]) < threshold)
              { reduce = false; break; }
          if (reduce) boundary_points_info[i].normals = normal_cone(n_mean);
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
    const mesh_fem &mf1 = mfu_of_boundary(ib1);
    const mesh_fem &mf2 = mfu_of_boundary(ib2);
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
    const mesh_fem &mf1 = mfu_of_boundary(ib1);
    const mesh_fem &mf2 = mfu_of_boundary(ib2);
    if ( &(mf1.linked_mesh()) != &(mf2.linked_mesh())) return false;
    GMM_ASSERT1(!(mf1.is_reduced()) && !(mf2.is_reduced()),
                "Nodal strategy can only be applied for non reduced fems");
    const mesh::ind_cv_ct &ic1 = mf1.convex_to_basic_dof(idof1);
    for (size_type i = 0; i < ic1.size(); ++i)
      if (cv == ic1[i]) return true;
    return false;
  }

  void multi_contact_frame::add_potential_contact_face
  (size_type ip, size_type ib, size_type ie, short_type i_f) {
    bool found = false;
    std::vector<face_info> &sfi = potential_pairs[ip];
    for (size_type k = 0; k < sfi.size(); ++k)
      if (sfi[k].ind_boundary == ib &&
          sfi[k].ind_element == ie &&
          sfi[k].ind_face == i_f) found = true;
    
    if (!found) sfi.push_back(face_info(ib, ie, i_f));
  }

  void multi_contact_frame::clear_aux_info(void) {
    boundary_points = std::vector<base_node>();
    boundary_points_info = std::vector<boundary_point>();
    element_boxes.clear();
    element_boxes_info = std::vector<influence_box>();
    potential_pairs = std::vector<std::vector<face_info> >();
  }

  multi_contact_frame::multi_contact_frame(size_type NN, scalar_type r_dist,
                                           int fem_nodes, bool dela,
                                           bool refc, bool selfc,
                                           scalar_type cut_a)
    : N(NN), self_contact(selfc), ref_conf(refc), use_delaunay(dela), 
      fem_nodes_mode(fem_nodes), release_distance(r_dist), cut_angle(cut_a),
      EPS(1E-8), md(0), coordinates(N), pt_eval(N) {
    if (N > 0) coordinates[0] = "x";
    if (N > 1) coordinates[1] = "y";
    if (N > 2) coordinates[2] = "z";
    if (N > 3) coordinates[3] = "w";
    GMM_ASSERT1(N <= 4, "Complete the definition for contact in "
                  "dimension greater than 4");
  }

  multi_contact_frame::multi_contact_frame(const model &mdd, size_type NN,
                                           scalar_type r_dist,
                                           int fem_nodes, bool dela,
                                           bool refc, bool selfc,
                                           scalar_type cut_a)
    : N(NN), self_contact(selfc), ref_conf(refc),
      use_delaunay(dela), fem_nodes_mode(fem_nodes), release_distance(r_dist),
      cut_angle(cut_a), EPS(1E-8), md(&mdd), coordinates(N), pt_eval(N) {
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

  size_type multi_contact_frame::add_slave_boundary
  (const mesh_im &mim, const mesh_fem &mfu,
   const model_real_plain_vector &U, size_type reg) {
    GMM_ASSERT1(mfu.linked_mesh().dim() == N,
                "Mesh dimension is " << mfu.linked_mesh().dim()
                << "should be " << N << ".");
    contact_boundary cb(reg, mfu, mim, add_U(U));
    size_type ind = contact_boundaries.size();
    contact_boundaries.push_back(cb);
    slave_boundaries.add(ind);
    return ind;
  }

  size_type multi_contact_frame::add_slave_boundary
  (const mesh_im &mim, const std::string &varname, size_type reg) {
    GMM_ASSERT1(md, "This multi contact frame object is not linked "
                "to a model");
    return add_slave_boundary(mim, md->mesh_fem_of_variable(varname),
                       md->real_variable(varname), reg);
  }

  size_type multi_contact_frame::add_master_boundary
  (const mesh_im &mim, const mesh_fem &mfu,
   const model_real_plain_vector &U, size_type reg) {
    GMM_ASSERT1(mfu.linked_mesh().dim() == N,
                "Mesh dimension is " << mfu.linked_mesh().dim()
                << "should be " << N << ".");
    contact_boundary cb(reg, mfu, mim, add_U(U));
    contact_boundaries.push_back(cb);
    return size_type(contact_boundaries.size() - 1);
  }

  size_type multi_contact_frame::add_master_boundary
  (const mesh_im &mim, const std::string &varname, size_type reg) {
    GMM_ASSERT1(md, "This multi contact frame object is not linked "
                "to a model");
    return add_master_boundary(mim, md->mesh_fem_of_variable(varname),
                       md->real_variable(varname), reg);
  }


  void multi_contact_frame::compute_influence_boxes(void) {
    fem_precomp_pool fppool;
    base_matrix G;
    model_real_plain_vector coeff;
    
    for (size_type i = 0; i < contact_boundaries.size(); ++i)
      if (!(slave_boundaries[i])) {
        size_type bnum = region_of_boundary(i);
        const mesh_fem &mfu = mfu_of_boundary(i);
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
          size_type nbd_t = pgt->nb_points();
          slice_vector_on_basic_dof_of_element(mfu, U, cv, coeff);
          bgeot::vectors_to_base_matrix
            (G, mfu.linked_mesh().points_of_convex(cv));
          
          pfem_precomp pfp = fppool(pf_s, &(pgt->geometric_nodes()));
          fem_interpolation_context ctx(pgt,pfp,size_type(-1), G, cv,
                                        size_type(-1));
          
          size_type nb_pt_on_face = 0;
          gmm::clear(n_mean);
          for (short_type ip = 0; ip < nbd_t; ++ip) {
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
            // computation of unit normal vector if the vertex is on the face
            bool is_on_face = false;
            bgeot::pconvex_structure cvs = pgt->structure();
            for (size_type k = 0; k < cvs->nb_points_of_face(v.f()); ++k)
              if (cvs->ind_points_of_face(v.f())[k] == ip) is_on_face = true;
            if (is_on_face) {
              ctx.set_ii(ip);
              n0 = bgeot::compute_normal(ctx, v.f());
              pf_s->interpolation_grad(ctx, coeff, grad, dim_type(N));
              gmm::add(gmm::identity_matrix(), grad);
              scalar_type J = gmm::lu_inverse(grad);
              if (J <= scalar_type(0))
                GMM_WARNING1("Inverted element !" << J);
              gmm::mult(gmm::transposed(grad), n0, n);
              n /= gmm::vect_norm2(n);
              n_mean += n;
              ++nb_pt_on_face;
            }
            
            if (ip == 0) // computation of bounding box
              bmin = bmax = val;
            else {
              for (size_type k = 0; k < N; ++k) {
                bmin[k] = std::min(bmin[k], val[k]);
                bmax[k] = std::max(bmax[k], val[k]);
              }
            }
          }
          
          GMM_ASSERT1(nb_pt_on_face,
                      "This element has not vertex on considered face !");
          
          // Computation of influence box :
          // offset of the bounding box relatively to the release distance
          scalar_type h = bmax[0] - bmin[0];
          for (size_type k = 1; k < N; ++k) h = std::max(h, bmax[k]-bmin[k]);
          if (h < release_distance/scalar_type(20))
            GMM_WARNING1("Found an element whose size is smaller than 1/20 "
                         "of the release distance. You should probably "
                         "adapt the release distance.");
          for (size_type k = 0; k < N; ++k)
            { bmin[k] -= release_distance; bmax[k] += release_distance; }
          
          // Store the influence box and additional information.
          element_boxes.add_box(bmin, bmax, element_boxes_info.size());
          n_mean /= gmm::vect_norm2(n_mean);
          element_boxes_info.push_back(influence_box(i, cv, v.f(), n_mean));
        }
      }
  }
  
  void multi_contact_frame::compute_boundary_points(bool slave_only) {
    fem_precomp_pool fppool;
    base_matrix G;
    model_real_plain_vector coeff;
    
    for (size_type i = 0; i < contact_boundaries.size(); ++i)
      if (!slave_only || slave_boundaries[i]) {
        size_type bnum = region_of_boundary(i);
        const mesh_fem &mfu = mfu_of_boundary(i);
        const mesh_im &mim = mim_of_boundary(i);
        const model_real_plain_vector &U = disp_of_boundary(i);
        const mesh &m = mfu.linked_mesh();
        bool on_fem_nodes = (fem_nodes_mode == 2 ||
                             (fem_nodes_mode == 1 && slave_boundaries[i]));
        
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
          pintegration_method pim = mim.int_method_of_element(cv);
          slice_vector_on_basic_dof_of_element(mfu, U, cv, coeff);
          bgeot::vectors_to_base_matrix
            (G, mfu.linked_mesh().points_of_convex(cv));
          dim_type qqdim = mfu.get_qdim() / pf_s->target_dim();;
          
          pfem_precomp pfp(0); size_type nbptf(0);
          if (on_fem_nodes) {
            pfp = fppool(pf_s, pf_s->node_tab(cv));
            nbptf =pf_s->node_convex(cv).structure()->nb_points_of_face(v.f());
          }
          else {
            pfp = fppool(pf_s,&(pim->approx_method()->integration_points()));
            nbptf = pim->approx_method()->nb_points_on_face(v.f());
          }
          fem_interpolation_context ctx(pgt,pfp,size_type(-1), G, cv, v.f());
          
          
          for (short_type ip = 0; ip < nbptf; ++ip) {
              
            size_type ind(0), inddof(0);
            if (on_fem_nodes) {
              inddof = mfu.ind_basic_dof_of_element(cv)[ip*qqdim];
              ind = pf_s->node_convex(cv).structure()
                ->ind_points_of_face(v.f())[ip];
            }
            else
              ind = pim->approx_method()->ind_first_point_on_face(v.f())+ip;
            ctx.set_ii(ind);
            
            if (!(on_fem_nodes && dof_already_interpolated[inddof])) {
              pf_s->interpolation(ctx, coeff, val, dim_type(N));
              val += ctx.xreal();            
              if (on_fem_nodes)  dof_ind[inddof] = boundary_points.size();
              
            }
            
            // unit normal vector computation
            n0 = bgeot::compute_normal(ctx, v.f());
            pf_s->interpolation_grad(ctx, coeff, grad, dim_type(N));
            gmm::add(gmm::identity_matrix(), grad);
            scalar_type J = gmm::lu_inverse(grad);
            if (J <= scalar_type(0)) GMM_WARNING1("Inverted element !" << J);
            gmm::mult(gmm::transposed(grad), n0, n);
            n /= gmm::vect_norm2(n);
            
            if (on_fem_nodes && dof_already_interpolated[inddof]) {
              boundary_points_info[dof_ind[inddof]].normals.add_normal(n);
            } else {
              boundary_points.push_back(val);
              boundary_points_info.push_back(boundary_point(ctx.xreal(), i, cv,
                                                            v.f(), inddof, n));
            }
            
            if (on_fem_nodes)  dof_already_interpolated.add(inddof);
          }
        } 
      } 
  }
  
  void multi_contact_frame::compute_potential_contact_pairs_delaunay(void) {
    
    compute_boundary_points();
    normal_cone_simplicication();
    potential_pairs = std::vector<std::vector<face_info> >();
    potential_pairs.resize(boundary_points.size());
    
    gmm::dense_matrix<size_type> simplexes;
    getfem::delaunay(boundary_points, simplexes);
    
    // connectivity analysis
    for (size_type is = 0; is < gmm::mat_ncols(simplexes); ++is) {
      // il faut tester ce qui se passe en cas de points répétés !!
      
      for (size_type i = 1; i <= N; ++i)
        for (size_type j = 0; j < i; ++j) {
          size_type ipt1 = simplexes(i, is), ipt2 = simplexes(j, is);
          boundary_point *pt_info1 = &(boundary_points_info[ipt1]);
          boundary_point *pt_info2 = &(boundary_points_info[ipt2]);
          size_type ib1 = pt_info1->ind_boundary;
          size_type ib2 = pt_info2->ind_boundary;
          bool sl1 = slave_boundaries[ib1];
          bool sl2 = slave_boundaries[ib2];
          if (!sl1 && sl2) { // The slave in first if any
            std::swap(ipt1, ipt2);
            std::swap(pt_info1, pt_info2);
            std::swap(sl1, sl2);
            std::swap(ib1, ib2);
          }
          const mesh_fem &mf1 = mfu_of_boundary(ib1);
          const mesh_fem &mf2 = mfu_of_boundary(ib2);
          
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
                  || ((fem_nodes_mode < 2)
                      && (( &(mf1.linked_mesh()) != &(mf2.linked_mesh()))
                          || (pt_info1->ind_element != pt_info2->ind_element)))
                  || ((fem_nodes_mode == 2)
                      && !(are_dof_linked(ib1, pt_info1->ind_dof,
                                          ib2, pt_info2->ind_dof)))
                  )
              ) {
            
            // Store the potential contact pairs
            
            if ((sl2 && fem_nodes_mode) || (!sl2 && fem_nodes_mode == 2)) {
              const mesh::ind_cv_ct &ic2
                = mf2.convex_to_basic_dof(pt_info2->ind_dof);
              for (size_type k = 0; k < ic2.size(); ++k) {
                size_type ir1 = region_of_boundary(ib1);
                mesh_region::face_bitset fbs
                  = mf2.linked_mesh().region(ir1).faces_of_convex(ic2[k]);
                for (short_type f = 0;
                     f < mf2.linked_mesh().nb_faces_of_convex(ic2[k]); ++f)
                  if (fbs.test(f))
                    add_potential_contact_face(ipt1,pt_info2->ind_boundary,
                                               ic2[k], f);
              }
            } else
              add_potential_contact_face(ipt1, pt_info2->ind_boundary,
                                         pt_info2->ind_element,
                                         pt_info2->ind_face);
            
            if (self_contact && !sl1 && !sl2) {
              if ((sl2 && fem_nodes_mode) || (!sl2 && fem_nodes_mode==2)) {
                const mesh::ind_cv_ct &ic1
                  = mf1.convex_to_basic_dof(pt_info1->ind_dof);
                for (size_type k = 0; k < ic1.size(); ++k) {
                  size_type ir2 = region_of_boundary(ib2);
                  mesh_region::face_bitset fbs
                    = mf1.linked_mesh().region(ir2).faces_of_convex(ic1[k]);
                  for (short_type f = 0;
                       f < mf1.linked_mesh().nb_faces_of_convex(ic1[k]);
                       ++f)
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
  
  void multi_contact_frame::compute_potential_contact_pairs_influence_boxes(void) {
    compute_influence_boxes();
    compute_boundary_points(!self_contact); // vraiment nécessaire ?
    normal_cone_simplicication();
    potential_pairs = std::vector<std::vector<face_info> >();
    potential_pairs.resize(boundary_points.size());
    
    for (size_type ip = 0; ip < boundary_points.size(); ++ip) {
      
      bgeot::rtree::pbox_set bset;
      element_boxes.find_boxes_at_point(boundary_points[ip], bset);
      boundary_point *pt_info = &(boundary_points_info[ip]);
      const mesh_fem &mf1 = mfu_of_boundary(pt_info->ind_boundary);
      
      bgeot::rtree::pbox_set::iterator it = bset.begin();
      for (; it != bset.end(); ++it) {
        // - faire un test de normale et de non élément courant
        influence_box &ibx = element_boxes_info[(*it)->id];
        size_type ib1 = ibx.ind_boundary;
        const mesh_fem &mf2 = mfu_of_boundary(ib1);
        
        // CRITERION 1 : The unit normal cone / vector are compatible
        //               and the two points are not in the same element.
        if (
            test_normal_cones_compatibility(ibx.mean_normal,
                                            pt_info->normals)
            // In case of self-contact, test if the points and the face
            // share the same element.
            && (((fem_nodes_mode < 2)
                 && (( &(mf1.linked_mesh()) != &(mf2.linked_mesh()))
                     || (pt_info->ind_element != ibx.ind_element)))
                || ((fem_nodes_mode == 2)
                    && !(is_dof_linked(ib1, pt_info->ind_dof,
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
    const base_node &x0, &y0;
    fem_interpolation_context &ctx;
    const model_real_plain_vector &coeff;
    const std::vector<base_small_vector> &ti;
    base_node val;
    base_matrix grad, gradtot;
    
    scalar_type operator()(const base_small_vector& a) {
      base_node x = x0;
      for (size_type i= 0; i < N-1; ++i) x += a[i] * ti[i];
      ctx.set_xref(x);
      ctx.pf()->interpolation(ctx, coeff, val, dim_type(N));
      base_node y = val + ctx.xreal();
      return gmm::vect_dist2(y, y0)/scalar_type(2);
    }
    void operator()(const base_small_vector& a, base_small_vector &grada) {
      base_node x = x0;
      for (size_type i = 0; i < N-1; ++i) x += a[i] * ti[i];
      ctx.set_xref(x);
      ctx.pf()->interpolation(ctx, coeff, val, dim_type(N));
      base_node dy = val + ctx.xreal() - y0;
      ctx.pf()->interpolation_grad(ctx, coeff, grad, dim_type(N));
      gmm::add(gmm::identity_matrix(), grad);
      gmm::mult(grad, ctx.K(), gradtot);
      for (size_type i = 0; i < N-1; ++i)
        grada[i] = gmm::vect_sp(gradtot, ti[i], dy);
    }
    
    proj_pt_surf_cost_function_object
    (const base_node &x00, const base_node &y00,
     fem_interpolation_context &ctxx,
     const model_real_plain_vector &coefff,
     const std::vector<base_small_vector> &tii)
      : N(gmm::vect_size(x00)), x0(x00), y0(y00), ctx(ctxx),
        coeff(coefff), ti(tii), val(N), grad(N,N), gradtot(N,N) {}
    
  };
  
  void multi_contact_frame::compute_contact_pairs(void) {
    base_matrix G, grad(N,N);
    model_real_plain_vector coeff;
    base_small_vector a(N-1), n(N);
   
    clear_aux_info();
    contact_pairs = std::vector<contact_pair>();
    
    extend_vectors();
    
    if (use_delaunay)
      compute_potential_contact_pairs_delaunay();
    else
      compute_potential_contact_pairs_influence_boxes();
    
    
    // Scan of potential pairs
    for (size_type ip = 0; ip < potential_pairs.size(); ++ip) {
      const base_node &y0 = boundary_points[ip];
      boundary_point &bpinfo = boundary_points_info[ip];
      // Detect here the nearest rigid obstacle (taking into account
      // the release distance)
      size_type irigid_obstacle(-1);
      scalar_type d0 = 1E300, d1, d2;
#if GETFEM_HAVE_MUPARSER_MUPARSER_H || GETFEM_HAVE_MUPARSER_H
      gmm::copy(y0, pt_eval);
      for (size_type i = 0; i < obstacles.size(); ++i) {
        d1 = scalar_type(obstacles_parsers[i].Eval());
        if (gmm::abs(d1) < release_distance && d1 < d0) {
          
          for (size_type j=0; j < bpinfo.normals.size(); ++j) {
            gmm::add(gmm::scaled(bpinfo.normals[j], EPS), pt_eval);
            d2 =  scalar_type(obstacles_parsers[i].Eval());
            if (d2 < d1) { d0 = d1; irigid_obstacle = i; break; }
            gmm::copy(y0, pt_eval);
          }
        }
      }
#else
      if (obstacles.size() > 0)
        GMM_WARNING1("Rigid obstacles are ignored. Recompile with "
                     "muParser to account for rigid obstacles");
#endif
      
      bool first_pair = true;
      for (size_type ipf = 0; ipf < potential_pairs[ip].size(); ++ipf) {
        // Point to surface projection. Principle :
        //  - One parametrizes first the face on the reference element by
        //    obtaining a point x_0 on that face and t_i, i=1..d-1 some
        //    orthonormals tangent vectors to the face.
        //  - Let y_0 be the point ot be projected and y the searched
        //    projected point. Then one searches for the minimum of
        //    J = (1/2)|| y - y0 ||
        //    with
        //    y = \phi(x0 + a_i t_i)
        //    (with a summation on i), where \phi = I+u(\tau(x)), and \tau 
        //    the geometric transformation between reference and real
        //    elements.
        //  - The gradient of J with respect to a_i is
        //    \partial_{a_j} J = (\phi(x0 + a_i t_i) - y0)
        //                       . (\nabla \phi(x0 + a_i t_i) t_j
        //  - BFGS is called.           
        
        const face_info &fi = potential_pairs[ip][ipf];
        size_type ib = fi.ind_boundary;
        size_type cv = fi.ind_element;
        short_type i_f = fi.ind_face;
        
        const mesh_fem &mfu = mfu_of_boundary(ib);
        const mesh &m = mfu.linked_mesh();
        pfem pf_s = mfu.fem_of_element(cv);
        bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
        
        
        slice_vector_on_basic_dof_of_element(mfu, disp_of_boundary(ib),
                                             cv, coeff);
        const base_node &x0 = pf_s->ref_convex(cv)->points_of_face(i_f)[0];
        bgeot::vectors_to_base_matrix(G, m.points_of_convex(cv));
        fem_interpolation_context ctx(pgt, pf_s, x0, G, cv, i_f);
        
        const base_small_vector &n0 = pf_s->ref_convex(cv)->normals()[i_f];
        std::vector<base_small_vector> ti(N-1);
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

        proj_pt_surf_cost_function_object pps(x0, y0, ctx, coeff, ti);
        
        gmm::iteration iter(2E-7, 0/* noisy*/, 200 /*maxiter*/);
        gmm::clear(a);
        gmm::bfgs(pps, pps, a, 10, iter);
        // Projection could be ameliorated by finding a starting point near
        // y0 (with respect to the integration method, for instance).
   

        // CRITERION 2 : The contact pair is eliminated when bfgs
        //               do not converge.
        if (!(iter.converged())) {
          GMM_WARNING2("Projection did not converge for point " << y0);
          continue;
        }
//         cout << "BFGS did converge : " << ctx.xreal() << endl;
//         cout << "xref : " << ctx.xref() << endl;
//         cout << "xref : " << a << endl;

          
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
        if (pf_s->ref_convex(cv)->is_in(ctx.xref()) > 1E-5) continue;
        
        base_node y(N);
        ctx.pf()->interpolation(ctx, coeff, y, dim_type(N));
        y += ctx.xreal();

        // CRITERION 4 : Apply the release distance
        if (gmm::vect_dist2(y, y0) > release_distance) continue;

        // compute the unit_normal and the signed distance ...
        base_small_vector n00 = bgeot::compute_normal(ctx, i_f);
        ctx.pf()->interpolation_grad(ctx, coeff, grad, dim_type(N));
        gmm::add(gmm::identity_matrix(), grad);
        scalar_type J = gmm::lu_inverse(grad);
        if (J <= scalar_type(0)) GMM_WARNING1("Inverted element !" << J);
        gmm::mult(gmm::transposed(grad), n00, n);
        n /= gmm::vect_norm2(n);

        scalar_type signed_dist = gmm::vect_sp(y - y0, n);
        // CRITERION 5 : comparison with rigid obstacles
        if (irigid_obstacle != size_type(-1) && signed_dist  > d0)
          continue;

        // CRITERION 6 : for self-contact only : apply a test on
        //               unit normals in reference configuration.
        size_type iby0 = bpinfo.ind_boundary;
        if (&m == &(mfu_of_boundary(iby0).linked_mesh())) {
          
          base_small_vector diff = bpinfo.ref_point - ctx.xreal();
          scalar_type ref_dist = gmm::vect_norm2(diff);
          
          if ( (ref_dist < scalar_type(4) * release_distance)
               && (gmm::vect_sp(diff, n00) < - 0.01 * ref_dist) )
            continue;
        }
        
        contact_pair ct(y0, bpinfo, ctx.xref(), y,  fi, signed_dist);
        if (first_pair) {
          contact_pairs.push_back(ct);
          first_pair = false;
        }
        else { 
          // CRITERION 7 : smallest signed distance on contact pairs
          if (contact_pairs.back().signed_dist > signed_dist)
            contact_pairs.back() = ct;
        }

      }
      if (first_pair && irigid_obstacle != size_type(-1)) {
        contact_pair ct(y0, bpinfo, irigid_obstacle, d0);
        
#if GETFEM_HAVE_MUPARSER_MUPARSER_H || GETFEM_HAVE_MUPARSER_H
      gmm::copy(y0, pt_eval);
      size_type nit = 0;
      while (gmm::abs(d0) > 1E-10 && ++nit < 1000) {
        for (size_type k = 0; k < N; ++k) {
          pt_eval[k] += EPS;
          d1 = scalar_type(obstacles_parsers[irigid_obstacle].Eval());
          n[k] = (d1 - d0) / EPS;
          pt_eval[k] -= EPS;
        }
        
        gmm::add(gmm::scaled(n, -d0 / gmm::vect_norm2_sqr(n)), pt_eval);
        // A simple line search could be added
        d0 = scalar_type(obstacles_parsers[irigid_obstacle].Eval());
      }
      GMM_ASSERT1(nit < 1000, "Projection on rigid obstacle did not converge");

      ct.master_point.resize(N);
      gmm::copy(pt_eval, ct.master_point);
      
#endif
        contact_pairs.push_back(ct);
      }
    }
      
      
    clear_aux_info();
  }



}  /* end of namespace getfem.                                             */
