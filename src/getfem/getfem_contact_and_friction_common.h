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
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/

/** @file getfem_contact_and_friction_common.h
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @author Konstantinos Poulios <logari81@googlemail.com>
    @date November, 2011.
    @brief Comomon tools for unilateral contact and Coulomb friction bricks.
 */
#ifndef GETFEM_CONTACT_AND_FRICTION_COMMON_H__
#define GETFEM_CONTACT_AND_FRICTION_COMMON_H__

#include "getfem_models.h"
#include "getfem_assembling_tensors.h"
#include "getfem/bgeot_rtree.h"
#include <getfem/getfem_mesher.h>


#include <getfem/getfem_arch_config.h>
#if GETFEM_HAVE_MUPARSER_MUPARSER_H
#include <muParser/muParser.h>
#elif GETFEM_HAVE_MUPARSER_H
#include <muParser.h>
#endif


namespace getfem {

  //=========================================================================
  //
  //  Projection on a ball and gradient of the projection.
  //
  //=========================================================================

  template<typename VEC> void ball_projection(const VEC &x,
					      scalar_type radius) {
    if (radius <= scalar_type(0))
      gmm::clear(const_cast<VEC&>(x));
    else {
      scalar_type a = gmm::vect_norm2(x);
      if (a > radius) gmm::scale(const_cast<VEC&>(x), radius/a);
    }
  }
  
  template<typename VEC, typename VECR>
  void ball_projection_grad_r(const VEC &x, scalar_type radius,
                              VECR &g) {
    if (radius > scalar_type(0)) {
      scalar_type a = gmm::vect_norm2(x);
      if (a >= radius) {
        gmm::copy(x, g); gmm::scale(g, scalar_type(1)/a);
        return;
      }
    }
    gmm::clear(g);
  }

  template <typename VEC, typename MAT>
  void ball_projection_grad(const VEC &x, scalar_type radius, MAT &g) {
    if (radius <= scalar_type(0)) { gmm::clear(g); return; }
    gmm::copy(gmm::identity_matrix(), g);
    scalar_type a = gmm::vect_norm2(x);
    if (a >= radius) {
      gmm::scale(g, radius/a);
      // gmm::rank_one_update(g, gmm::scaled(x, -radius/(a*a*a)), x);
      for (size_type i = 0; i < x.size(); ++i)
        for (size_type j = 0; j < x.size(); ++j)
          g(i,j) -= radius*x[i]*x[j] / (a*a*a);
    }
  }



  template <typename VEC, typename VECR>
  void coupled_projection(const VEC &x, const VEC &n,
			  scalar_type f, VECR &g) {
    scalar_type xn = gmm::vect_sp(x, n);
    scalar_type xnm = gmm::neg(xn);
    scalar_type th = f * xnm;
    scalar_type xtn = gmm::sqrt(gmm::vect_norm2_sqr(x) - xn*xn);
    
    gmm::copy(gmm::scaled(n, -xnm), g);
    if (th > scalar_type(0)) {
      if (xtn <= th) {
	gmm::add(x, g);
	gmm::add(gmm::scaled(n, -xn), g);
      } else {
	gmm::add(gmm::scaled(x, f*xnm/xtn), g);
	gmm::add(gmm::scaled(n, -f*xnm*xn/xtn), g);
      }
    }
  }


  template <typename VEC, typename MAT>
  void coupled_projection_grad(const VEC &x, const VEC &n,
			       scalar_type f, MAT &g) {
    scalar_type xn = gmm::vect_sp(x, n);
    scalar_type xnm = gmm::neg(xn);
    scalar_type th = f * xnm;
    scalar_type xtn = gmm::sqrt(gmm::vect_norm2_sqr(x) - xn*xn);
    size_type N = gmm::vect_size(x);
    gmm::clear(g);
    
    if (th > scalar_type(0)) {
      if (xtn <= th) {
	gmm::copy(gmm::identity_matrix(), g);
	gmm::rank_one_update(g, gmm::scaled(n, -scalar_type(1)), n);
      } else if (xn < scalar_type(0)) {
	static base_small_vector t; gmm::resize(t, N);
	gmm::add(x, gmm::scaled(n, -xn), t);
        gmm::scale(t, scalar_type(1)/xtn);
	if (N > 2) {
	  gmm::copy(gmm::identity_matrix(), g);
	  gmm::rank_one_update(g, gmm::scaled(t, -scalar_type(1)), t);
	  gmm::rank_one_update(g, gmm::scaled(n, -scalar_type(1)), n);
	  gmm::scale(g, -xn*th/xtn);
        }
        gmm::rank_one_update(g, gmm::scaled(t, -f), n);
      }
    }

    if (xn < scalar_type(0)) gmm::rank_one_update(g, n, n);
  }

  //=========================================================================
  //
  //  De Saxce projection and its gradients.
  //
  //=========================================================================


  template<typename VEC>
  void De_Saxce_projection(const VEC &x, const VEC &n_, scalar_type f) {
    static base_small_vector n; // For more robustness, n_ is not supposed unitary
    size_type N = gmm::vect_size(x);
    gmm::resize(n, N);
    gmm::copy(gmm::scaled(n_, scalar_type(1)/gmm::vect_norm2(n_)), n);
    scalar_type xn = gmm::vect_sp(x, n);
    scalar_type nxt = sqrt(gmm::abs(gmm::vect_norm2_sqr(x) - xn*xn));
    if (xn >= scalar_type(0) && f * nxt <= xn) {
      gmm::clear(const_cast<VEC&>(x));
    } else if (xn > scalar_type(0) || nxt > -f*xn) {
      gmm::add(gmm::scaled(n, -xn), const_cast<VEC&>(x));
      gmm::scale(const_cast<VEC&>(x), -f / nxt);
      gmm::add(n, const_cast<VEC&>(x));
      gmm::scale(const_cast<VEC&>(x), (xn - f * nxt) / (f*f+scalar_type(1)));
    }
  }

  template<typename VEC, typename MAT>
  void De_Saxce_projection_grad(const VEC &x, const VEC &n_,
				scalar_type f, MAT &g) {
    static base_small_vector n;
    size_type N = gmm::vect_size(x);
    gmm::resize(n, N);
    gmm::copy(gmm::scaled(n_, scalar_type(1)/gmm::vect_norm2(n_)), n);
    scalar_type xn = gmm::vect_sp(x, n);
    scalar_type nxt = sqrt(gmm::abs(gmm::vect_norm2_sqr(x) - xn*xn));


    if (xn > scalar_type(0) && f * nxt <= xn) {
      gmm::clear(g);
    } else if (xn > scalar_type(0) || nxt > -f*xn) {
      static base_small_vector xt;
      gmm::resize(xt, N);
      gmm::add(x, gmm::scaled(n, -xn), xt);
      gmm::scale(xt, scalar_type(1)/nxt);

      if (N > 2) {
        gmm::copy(gmm::identity_matrix(), g);
        gmm::rank_one_update(g, gmm::scaled(n, -scalar_type(1)), n);
        gmm::rank_one_update(g, gmm::scaled(xt, -scalar_type(1)), xt);
        gmm::scale(g, f*(f - xn/nxt));
      } else {
        gmm::clear(g);
      }

      gmm::scale(xt, -f); gmm::add(n, xt);
      gmm::rank_one_update(g, xt, xt);
      gmm::scale(g, scalar_type(1) / (f*f+scalar_type(1)));
    } else {
      gmm::copy(gmm::identity_matrix(), g);
    }
  }


  template<typename VEC, typename MAT>
  static void De_Saxce_projection_gradn(const VEC &x, const VEC &n_,
					scalar_type f, MAT &g) {
    static base_small_vector n;
    size_type N = gmm::vect_size(x);
    scalar_type nn = gmm::vect_norm2(n_);
    gmm::resize(n, N);
    gmm::copy(gmm::scaled(n_, scalar_type(1)/nn), n);
    scalar_type xn = gmm::vect_sp(x, n);
    scalar_type nxt = sqrt(gmm::abs(gmm::vect_norm2_sqr(x) - xn*xn));
    gmm::clear(g);

    if (!(xn > scalar_type(0) && f * nxt <= xn)
	&& (xn > scalar_type(0) || nxt > -f*xn)) {
      static base_small_vector xt, aux;
      gmm::resize(xt, N); gmm::resize(aux, N);
      gmm::add(x, gmm::scaled(n, -xn), xt);
      gmm::scale(xt, scalar_type(1)/nxt);

      scalar_type c = (scalar_type(1) + f*xn/nxt)/nn;
      for (size_type i = 0; i < N; ++i) g(i,i) = c;
      gmm::rank_one_update(g, gmm::scaled(n, -c), n);
      gmm::rank_one_update(g, gmm::scaled(n, f/nn), xt);
      gmm::rank_one_update(g, gmm::scaled(xt, -f*xn/(nn*nxt)), xt);
      gmm::scale(g, xn - f*nxt);

      gmm::add(gmm::scaled(xt, -f), n, aux);
      gmm::rank_one_update(g, aux, gmm::scaled(xt, (nxt+f*xn)/nn));
     
      gmm::scale(g, scalar_type(1) / (f*f+scalar_type(1)));
    }
  }


#if 0

  //=========================================================================
  //
  //  Structure which store the contact boundaries, rigid obstacles and
  //  computes the contact pairs in large sliding/large deformation
  //
  //=========================================================================


  
  // TODO:
  // - algo de detection des paires de contact dans les deux cas (boites et
  //   delaunay)
  // - algo de présélection des zones d'influences/pts de gauss
  // - algo de projection
  // - Dans le cas Delaunay, un passage éventuel à l'élément voisin est-il
  //   envisageable ?
  // - algo de selection des paires de contact valides



  // - penser à eliminer le stockage non necessaire après la detection des
  //   paires de contact
  // - A la fin, revoir éventuellement si la stratégie de calcul des
  //   extensions des déplacement est la bonne (calcul en début de
  //   construction de la liste des boites ou points des bords de contact).
  // - Dans la structure finale, penser à séparer mieux les données (bords
  //   de contact, obstacles) et les résultats stockés. Ajouter les zones
  //   d'influence d'éléments, les distances de coupure, l'algo de ...
  // - Gerer le cas configuration de référence
  // - Dans le cas Delaunay, gérer les points coincidents ... ou voir si le
  //   delaunay les gère correctement, ou les perturber infinitesimalement ....




  class multi_contact_frame {

    // Structure describing a contact boundary
    struct contact_boundary {
      size_type region;                 // Boundary number
      const getfem::mesh_fem *mfu;      // F.e.m. for the displacement.
      const getfem::mesh_im *mim;       // Integration method for the boundary.
      size_type ind_U;                  // Index of displacement.
      contact_boundary(void) {}
      contact_boundary(size_type r, const mesh_fem &mf,
                       const mesh_im &mi, size_type i)
        : region(r), mfu(&mf), mim(&mi), ind_U(i) {}
    };


    size_type N;          // Meshes dimensions
    bool self_contact;    // self-contact is searched or not.
    bool ref_conf;        // contact in reference configuration
                          // for linear elasticity small sliding contact.
    bool use_delaunay;    // Use delaunay to detect the contact pairs instead
                          // of influence boxes.
    int fem_nodes_mode;   // 0 = Use Gauss points for both slave and master
                          // 1 = Use finite element nodes for slave and
                          //     Gauss points for master.
                          // 2 = Use finite element nodes for both slave
                          //     and master
                          
    scalar_type release_distance;  // Limit distance beyond which the contact
    // will not be considered. CAUTION: should be comparable to the element
    // size (if it is too large, a too large set of influence boxes will be
    // detected and the computation will be slow, except for delaunay option) 

    typedef model_real_plain_vector VECTOR;
    std::vector<const VECTOR *> Us;  // Displacement vectors
    std::vector<VECTOR> ext_Us;      // Unreduced displacement vectors
                                     // CAUTION : they have to be updated

    // Contact pairs are seached for a certain boundary (master or slave,
    // depending on the contact algorithm) on the master ones. If contact pairs
    // are searched for a master boundary, self-contact is taken into account
    // if the flag 'self_contact' is set to 'true'. Self-contact is never taken
    // into account for a slave boundary.
    dal::bit_vector slave_boundaries;
    std::vector<contact_boundary> contact_boundaries;

    std::vector<std::string> coordinates;
    base_node pt_eval;
#if GETFEM_HAVE_MUPARSER_MUPARSER_H || GETFEM_HAVE_MUPARSER_H
    std::vector<mu::Parser> obstacles_parsers;
#endif
    std::vector<std::string> obstacles;
    std::vector<std::string> obstacles_velocities;

      

    // 
    // Influence boxes part
    //
    bgeot::rtree element_boxes;                  // influence boxes

    // 
    // Delaunay part
    //
    std::vector<base_node> pts; // Il faudrait éliminer les doublons ...

    // 
    // Common part
    //
    struct normal_cone {
      std::vector<base_node> unit_normals;
      
      void add_normal(const base_node &n) { unit_normals.push_back(n); }
      // et faire des fonctions de tests de coincidence, de reduction si tout
      // le cone est conpris dans un certain angle ...
      void normal_cone(const base_node &n) : unit_normals(1, n) { }
    };
    std::vector<size_type> boundary_of_elements; // boundary number each box
    std::vector<size_type> ind_of_elements;  // element index in the boundary
    std::vector<size_type> face_of_elements; // face of the element
    std::vector<normal_cone> normal_cones;     // mean outward unit normal
    
    



    size_type add_U(const getfem::mesh_fem &mfu,
                    const model_real_plain_vector &U) {
      size_type i = 0;
      for (; i < Us.size(); ++i) if (Us[i] == &U) return i;
      Us.push_back(&U);
      model_real_plain_vector ext_U(mfu.nb_basic_dof());
      mfu.extend_vector(U, ext_U);
      ext_Us.push_back(ext_U);
      return i;
    }

    void extend_vectors(void) {
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

  public:

    const getfem::mesh_fem &mfu_of_boundary(size_type n) const
    { return *(contact_boundaries[n].mfu); }
    const getfem::mesh_im  &mim_of_boundary(size_type n) const
    { return *(contact_boundaries[n].mim); }
    const model_real_plain_vector &disp_of_boundary(size_type n) const
    { return ext_Us[contact_boundaries[n].ind_U]; }
    size_type region_of_boundary(size_type n) const
    { return contact_boundaries[n].region; }

    multi_contact_frame(size_type NN, scalar_type r_dist,
                        int fem_nodes = 0, bool dela = true,
                        bool refc = false, bool selfc = true)
      : N(NN), self_contact(selfc), ref_conf(refc), use_delaunay(dela), 
        fem_nodes_mode(fem_nodes), release_distance(r_dist),
        coordinates(N), pt_eval(N) {
      if (N > 0) coordinates[0] = "x";
      if (N > 1) coordinates[1] = "y";
      if (N > 2) coordinates[2] = "z";
      if (N > 3) coordinates[3] = "w";
      GMM_ASSERT1(N <= 4, "Complete the definition for contact in "
                  "dimension greater than 4");
    }

    size_type add_obstacle(const std::string &obs) {
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

    size_type add_slave_boundary(const getfem::mesh_im &mim,
                                  const getfem::mesh_fem &mfu,
                                  const model_real_plain_vector &U,
                                  size_type reg) {
      GMM_ASSERT1(mfu.linked_mesh().dim() == N,
                  "Mesh dimension is " << mfu.linked_mesh().dim()
                  << "should be " << N << ".");
      contact_boundary cb(reg, mfu, mim, add_U(mfu, U));
      size_type ind = contact_boundaries.size();
      contact_boundaries.push_back(cb);
      slave_boundaries.add(ind);
      return ind;
    }

    size_type add_master_boundary(const getfem::mesh_im &mim,
                                 const getfem::mesh_fem &mfu,
                                 const model_real_plain_vector &U,
                                 size_type reg) {
      GMM_ASSERT1(mfu.linked_mesh().dim() == N,
                  "Mesh dimension is " << mfu.linked_mesh().dim()
                  << "should be " << N << ".");
      contact_boundary cb(reg, mfu, mim, add_U(mfu, U));
      contact_boundaries.push_back(cb);
      return size_type(contact_boundaries.size() - 1);
    }

    // Compute the influence boxes of master boundary elements. To be run
    // before the detection of contact pairs. The influence box is the
    // bounding box extended by a distance equal to the release distance.
    void compute_influence_boxes(void) {
      fem_precomp_pool fppool;
      base_matrix G;
      model_real_plain_vector coeff;

      element_boxes.clear();
      normal_cones.resize(0);
      boundary_of_elements.resize(0);
      ind_of_elements.resize(0);
      face_of_elements.resize(0);
      
      for (size_type i = 0; i < contact_boundaries.size(); ++i)
        if (!(slave_boundaries[i])) {
          size_type bnum = region_of_boundary(i);
          const mesh_fem &mfu = mfu_of_boundary(i);
          const model_real_plain_vector &U = disp_of_boundary(i);
          const mesh &m = mfu.linked_mesh();
          
          base_node val(N), bmin(N), bmax(N), n0(N), n(N), n_mean(N);
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
            element_boxes.add_box(bmin, bmax, normal_cones.size());
            n_mean /= gmm::vect_norm2(n_mean);
            unit_normals.push_back(normal_cone(n_mean));
            boundary_of_elements.push_back(i);
            ind_of_elements.push_back(cv);
            face_of_elements.push_back(v.f());
          }
        }
    }
    
    // For delaunay triangulation. Advantages compared to influence boxes:
    // No degeneration of the algorithm complexity with refinement and
    // more easy to extend to fictitious domain with contact.
    void compute_boundary_points(void) {
      fem_precomp_pool fppool;
      base_matrix G;
      model_real_plain_vector coeff;

      pts.clear();
      normal_cones.resize(0);
      boundary_of_elements.resize(0);
      ind_of_elements.resize(0);
      face_of_elements.resize(0);
      
      // should store all the boundary deformed points relatively to
      // an integration method (for the master boundaries) and
      // relatively to an integration method or to finite element
      // nodes for a slave boundary.
      // Storing sufficient information to perform a Delaunay triangulation
      // and to be able to recover the boundary number, element number,
      // face number, mean normal ...
      
      for (size_type i = 0; i < contact_boundaries.size(); ++i) {

        size_type bnum = region_of_boundary(i);
        const mesh_fem &mfu = mfu_of_boundary(i);
        const mesh_fem &mim = mim_of_boundary(i);
        const model_real_plain_vector &U = disp_of_boundary(i);
        const mesh &m = mfu.linked_mesh();
        bool on_fem_nodes = (fem_nodes_mode == 2
                             || (fem_nodes_mode == 1 && slave_boundaries[i]));
        
        base_node val(N), bmin(N), bmax(N), n0(N), n(N), n_mean(N);
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
          size_type nbd_t = pgt->nb_points();
          slice_vector_on_basic_dof_of_element(mfu, U, cv, coeff);
          bgeot::vectors_to_base_matrix
            (G, mfu.linked_mesh().points_of_convex(cv));
          dim_type qqdim = mfu.get_qdim() / pf_s->target_dim();;
          
          pfem_precomp pfp(0); size_type nbptf(0);
          if (on_fem_nodes) {
            pfp = fppool(pf_s, pf_s->node_tab(cv));
            nbptf = pf_s->node_convex->structure()->nb_points_on_face(v.f());
          }
          else {
            pfp = fppool(pf_s,&(pim->approx_method()->integration_points()));
            nbptf = pim->approx_method()->nb_points_on_face(v.f());
          }
          fem_interpolation_context ctx(pgt,pfp,size_type(-1), G, cv, v.f());
          
          
          for (short_type ip = 0; ip < nbptf; ++ip) {

            size_type ind(0), inddof(0);
            if (on_fem_nodes) {
              inddof = mfu.ind_dof_of_element(cv)[ip*qqdim];
              ind =
                pf_s->node_convex->structure()->ind_points_of_face(v.f())[ip];
            }
            else
              ind = pim->approx_method()->ind_first_point_on_face(v.f())+ip;
            ctx.set_ii(ind);

            if (!(on_fem_nodes && dof_already_interpolated[inddof])) {
              pf_s->interpolation(ctx, coeff, val, dim_type(N));
              val += ctx.xreal();            
              if (on_fem_nodes)  dof_ind[inddof] = pts.size();
              pts.push_back(val);
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
              normal_cones[dof_ind[inddof]].add_normal(n);
            } else {
              normal_cones.push_back(normal_cone(n));
              boundary_of_elements.push_back(i);
              ind_of_elements.push_back(cv);
              face_of_elements.push_back(v.f());
            }

            if (on_fem_nodes)  dof_already_interpolated.add(inddof);
          }
        } 
      } 
    }



    // The whole process of the computation of contact pairs
    void compute_contact_pairs(void) {
      extend_vectors();

      if (use_delaunay)
        compute_boundary_points();
      else
        compute_influence_boxes();

      // Loop over contact points on slave surfaces
      for (size_type i = 0; i < contact_boundaries.size(); ++i)
        if (slave_boundaries[i] || self_contact) {
          
          size_type bnum = region_of_boundary(i);
          const mesh_fem &mfu = mfu_of_boundary(i);
          const mesh_fem &mim = mim_of_boundary(i);
          const model_real_plain_vector &U = disp_of_boundary(i);
          const mesh &m = mfu.linked_mesh();
          
          base_node val(N), bmin(N), bmax(N), n0(N), n(N), n_mean(N);
          base_matrix grad(N,N);
          mesh_region region = m.region(bnum);
          GMM_ASSERT1(mfu.get_qdim() == N, "Wrong mesh_fem qdim");
          
          dal::bit_vector points_already_interpolated;
          for (getfem::mr_visitor v(region,m); !v.finished(); ++v) {
            size_type cv = v.cv();
            bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
            pfem pf_s = mfu.fem_of_element(cv);
            pintegration_method pim = mim.int_method_of_element(cv);
            size_type nbd_t = pgt->nb_points();
            slice_vector_on_basic_dof_of_element(mfu, U, cv, coeff);
            bgeot::vectors_to_base_matrix
              (G, mfu.linked_mesh().points_of_convex(cv));
            
            pfem_precomp pfp(0); size_type nbptf(0);
            if (slave_boundaries[i] && use_fem_nodes) {
              pfp = fppool(pf_s, pf_s->node_tab(cv));
              nbptf = pf_s->node_convex->structure()->nb_points_on_face(v.f());
            }
            else {
              pfp = fppool(pf_s,&(pim->approx_method()->integration_points()));
              nbptf = pim->approx_method()->nb_points_on_face(v.f());
            }
            fem_interpolation_context ctx(pgt,pfp,size_type(-1), G, cv, v.f());
            
            for (short_type ip = 0; ip < nbptf; ++ip) {

              // ici la liste des points ...


            }
          }
        }



      // boucle sur les noeuds/points de Gauss surface esclaves et
      // maitres eventuellement

      //   - simplification éventuelle des normales si nodal 

      //   - calcul des paires potentielles par l'un ou l'autre moyen

      //   - calcul de la projection

      //   - éventuel passage au voisin

      //   - discrimination


    }





  };















#endif






}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONTACT_AND_FRICTION_COMMON_H__ */
