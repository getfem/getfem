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



  //=========================================================================
  //
  //  Structure which store the contact boundaries, rigid obstacles and
  //  computes the contact pairs in large sliding/large deformation
  //
  //=========================================================================

  class multi_contact_frame {

    // Structure describing a contact boundary
    struct contact_boundary {
      size_type region;            // Boundary number
      const getfem::mesh_fem *mfu; // F.e.m. for the displacement.
      const getfem::mesh_fem *mflambda; // F.e.m. for the displacement.
      const getfem::mesh_im *mim;  // Integration method for the boundary.
      std::string multname;        // Name of the optional contact stress
                                   // multiplier when linked to a model.
      size_type ind_U;             // Index of displacement.
      size_type ind_lambda;        // Index of multiplier (if any).
      contact_boundary(void) {}
      contact_boundary(size_type r, const mesh_fem *mf,
                       const mesh_im &mi, size_type i, const mesh_fem *mfl,
                       size_type j = size_type(-1))
        : region(r), mfu(mf), mflambda(mfl), mim(&mi), ind_U(i),
          ind_lambda(j) {}
    };


    size_type N;          // Meshes dimensions
    bool self_contact;    // Self-contact is searched or not.
    bool ref_conf;        // Contact in reference configuration
                          // for linear elasticity small sliding contact.
    bool use_delaunay;    // Use delaunay to detect the contact pairs instead
                          // of influence boxes.
    int fem_nodes_mode;   // 0 = Use Gauss points for both slave and master
                          // 1 = Use finite element nodes for slave and
                          //     Gauss points for master.
                          // 2 = Use finite element nodes for both slave
                          //     and master
    bool raytrace;        // Use raytrace instead of projection.
                          
    scalar_type release_distance;  // Limit distance beyond which the contact
    // will not be considered. CAUTION: should be comparable to the element
    // size (if it is too large, a too large set of influence boxes will be
    // detected and the computation will be slow, except for delaunay option) 

    scalar_type cut_angle; // Cut angle (in radian) for normal cones
    scalar_type EPS;       // Should be typically hmin/1000 (for computing
                           // gradients with finite differences 
    const model *md;       // The model if the structure is linked to a model.

    typedef model_real_plain_vector VECTOR;
    std::vector<const VECTOR *> Us;  // Displacement vectors
    std::vector<const VECTOR *> Ws;  // "Velocity" vectors
    std::vector<std::string> Unames; // Displacement vectors names. 
    std::vector<std::string> Wnames; // "Velocity" vectors names. 
    std::vector<VECTOR> ext_Us;      // Unreduced displacement vectors
    std::vector<VECTOR> ext_Ws;      // Unreduced "velocity" vectors
    std::vector<const VECTOR *> lambdas;  // Displacement vectors
    std::vector<std::string> lambdanames; // Displacement vectors names. 
    std::vector<VECTOR> ext_lambdas;      // Unreduced displacement vectors

    dal::bit_vector slave_boundaries;
    std::vector<contact_boundary> contact_boundaries;

    std::vector<std::string> coordinates;
    base_node pt_eval;
#if GETFEM_HAVE_MUPARSER_MUPARSER_H || GETFEM_HAVE_MUPARSER_H
    std::vector<mu::Parser> obstacles_parsers;
#endif
    std::vector<std::string> obstacles;
    std::vector<std::string> obstacles_velocities;

      
    struct normal_cone : public std::vector<base_small_vector> {
      
      void add_normal(const base_small_vector &n)
      { std::vector<base_small_vector>::push_back(n);}
      normal_cone(void) {}
      normal_cone(const base_small_vector &n)
        : std::vector<base_small_vector>(1, n) { }
    };

    // 
    // Influence boxes
    //
    struct influence_box {     // Additional information for an influence box
      size_type ind_boundary;  // Boundary number
      size_type ind_element;   // Element number
      short_type ind_face;     // Face number in element
      base_small_vector mean_normal;   // Mean outward normal unit vector
      influence_box(void) {}
      influence_box(size_type ib, size_type ie,
                    short_type i_f, const base_small_vector &n)
        : ind_boundary(ib), ind_element(ie), ind_face(i_f), mean_normal(n) {}
    };

    bgeot::rtree element_boxes;                  // influence boxes
    std::vector<influence_box> element_boxes_info;

    // 
    // Stored points (for Delaunay and slave nodal boundaries)
    //
    
    struct boundary_point {    // Additional information for a boundary point
      base_node ref_point;     // Point coordinate in reference configuration
      size_type ind_boundary;  // Boundary number
      size_type ind_element;   // Element number
      short_type ind_face;     // Face number in element
      size_type ind_pt;        // Dof number for fem_nodes_mode or point number
                               // of integration method 
      normal_cone normals;     // Set of outward unit normal vectors
      boundary_point(void) {}
      boundary_point(const base_node &rp, size_type ib, size_type ie,
                     short_type i_f, size_type id, const base_small_vector &n)
        : ref_point(rp), ind_boundary(ib), ind_element(ie), ind_face(i_f),
          ind_pt(id), normals(n) {}
    };
      
    std::vector<base_node> boundary_points;
    std::vector<boundary_point> boundary_points_info;   


    size_type add_U(const model_real_plain_vector *U, const std::string &name,
                    const model_real_plain_vector *w,const std::string &wname);
    size_type add_lambda(const model_real_plain_vector &lambda,
                         const std::string &name);

    void extend_vectors(void);

    void normal_cone_simplicication(void);

    bool test_normal_cones_compatibility(const normal_cone &nc1,
                                         const normal_cone &nc2);

    bool test_normal_cones_compatibility(const base_small_vector &n,
                                         const normal_cone &nc2);

    dal::bit_vector aux_dof_cv; // An auxilliary variable for are_dof_linked
    // function (in order to be of constant complexity).

    bool are_dof_linked(size_type ib1, size_type idof1,
                        size_type ib2, size_type idof2);

    bool is_dof_linked(size_type ib1, size_type idof1,
                       size_type ib2, size_type cv);
  public:

    struct face_info {
      size_type ind_boundary;
      size_type ind_element;
      short_type ind_face;
      face_info(void) {}
      face_info(size_type ib, size_type ie, short_type iff)
        : ind_boundary(ib), ind_element(ie), ind_face(iff) {}
    };
    
  protected:

    std::vector<std::vector<face_info> > potential_pairs;

    void add_potential_contact_face(size_type ip, size_type ib, size_type ie,
                                    short_type i_f);
  public:

    // stored information for contact pair
    struct contact_pair {
      
      base_node slave_point;         // The transformed slave point
      base_small_vector slave_n;     // Normal unit vector to slave surface
      size_type slave_ind_boundary;  // Boundary number
      size_type slave_ind_element;   // Element number
      short_type slave_ind_face;     // Face number in element
      size_type slave_ind_pt;        // Dof number for fem_nodes_mode
                                     // or integration point number otherwise

      base_node master_point_ref;    // The master point on ref element
      base_node master_point;        // The transformed master point
      base_small_vector master_n;    // Normal unit vector to master surface
      face_info master_face_info;

      scalar_type signed_dist;
      
      size_type irigid_obstacle;

      contact_pair(void) {}
      contact_pair(const base_node &spt, const base_small_vector &nx, 
                   const boundary_point &bp,
                   const base_node &mptr,  const base_node &mpt,
                   const base_small_vector &ny, 
                   const face_info &mfi, scalar_type sd)
        : slave_point(spt), slave_n(nx), slave_ind_boundary(bp.ind_boundary),
          slave_ind_element(bp.ind_element), slave_ind_face(bp.ind_face),
          slave_ind_pt(bp.ind_pt), master_point_ref(mptr),
          master_point(mpt), master_n(ny), master_face_info(mfi),
          signed_dist(sd), irigid_obstacle(-1) {}
      contact_pair(const base_node &spt, const base_small_vector &nx,
                   const boundary_point &bp,
                   const base_node &mpt, const base_small_vector &ny,
                   size_type ir, scalar_type sd)
        : slave_point(spt), slave_n(nx), slave_ind_boundary(bp.ind_boundary),
          slave_ind_element(bp.ind_element), slave_ind_face(bp.ind_face),
          slave_ind_pt(bp.ind_pt), master_point(mpt), master_n(ny),
          signed_dist(sd),
          irigid_obstacle(ir) {}
      
    };


    // Compute the influence boxes of master boundary elements. To be run
    // before the detection of contact pairs. The influence box is the
    // bounding box extended by a distance equal to the release distance.
    void compute_influence_boxes(void);
    
    // For delaunay triangulation. Advantages compared to influence boxes:
    // No degeneration of the algorithm complexity with refinement and
    // more easy to extend to fictitious domain with contact.
    // Stores all the boundary deformed points relatively to
    // an integration method or to finite element nodes (depending on
    // fem_nodes_mode flag). Storing sufficient information to perform
    // a Delaunay triangulation and to be able to recover the boundary
    // number, element number, face number, unit normal vector ...
    void compute_boundary_points(bool slave_only = false);
    void compute_potential_contact_pairs_delaunay(void);
    void compute_potential_contact_pairs_influence_boxes(void);

  protected:
      
    std::vector<contact_pair> contact_pairs;

    void clear_aux_info(void); // Delete auxillairy information

  public:

    size_type dim(void) const { return N; }
    const std::vector<contact_pair> &ct_pairs(void) const
    { return contact_pairs; }


    const getfem::mesh_fem &mfu_of_boundary(size_type n) const
    { return *(contact_boundaries[n].mfu); }
    const getfem::mesh_fem &mflambda_of_boundary(size_type n) const
    { return *(contact_boundaries[n].mflambda); }
    const getfem::mesh_im  &mim_of_boundary(size_type n) const
    { return *(contact_boundaries[n].mim); }
    size_type nb_variables(void) const { return Us.size(); }
    size_type nb_multipliers(void) const { return lambdas.size(); }
    const std::string &varname(size_type i) const { return Unames[i]; }
    const std::string &multname(size_type i) const { return lambdanames[i]; }
    const model_real_plain_vector &disp_of_boundary(size_type n) const
    { return ext_Us[contact_boundaries[n].ind_U]; }
    const model_real_plain_vector &w_of_boundary(size_type n) const
    { return ext_Ws[contact_boundaries[n].ind_U]; }
    const model_real_plain_vector &mult_of_boundary(size_type n) const
    { return ext_lambdas[contact_boundaries[n].ind_lambda]; }
    size_type region_of_boundary(size_type n) const
    { return contact_boundaries[n].region; }
    const std::string &varname_of_boundary(size_type n) const
    { return Unames[contact_boundaries[n].ind_U]; }
    size_type ind_varname_of_boundary(size_type n) const
    { return contact_boundaries[n].ind_U; }
    const std::string &multname_of_boundary(size_type n) const {
      static const std::string vname;
      size_type ind = contact_boundaries[n].ind_lambda;
      return (ind == size_type(-1)) ? vname : lambdanames[ind];
    }
    size_type ind_multname_of_boundary(size_type n) const
    { return contact_boundaries[n].ind_lambda; }
    size_type nb_boundaries(void) const { return contact_boundaries.size(); }
    bool is_self_contact(void) const { return self_contact; }
    bool is_slave_boundary(size_type n) const { return slave_boundaries[n]; }
    void set_raytrace(bool b) { raytrace = b; }
    void set_fem_nodes_mode(int m) { fem_nodes_mode = m; }
    size_type nb_contact_pairs(void) const { return contact_pairs.size(); }
    const contact_pair &get_contact_pair(size_type i)
    { return contact_pairs[i]; }

    multi_contact_frame(size_type NN, scalar_type r_dist,
                        bool dela = true, bool selfc = true,
                        scalar_type cut_a = 0.3, bool rayt = false,
                        int fem_nodes = 0, bool refc = false);
    multi_contact_frame(const model &md, size_type NN, scalar_type r_dist,
                        bool dela = true, bool selfc = true,
                        scalar_type cut_a = 0.3, bool rayt = false,
                        int fem_nodes = 0, bool refc = false);

    size_type add_obstacle(const std::string &obs);

    size_type add_slave_boundary(const getfem::mesh_im &mim,
                                 const getfem::mesh_fem *mfu,
                                 const model_real_plain_vector *U,
                                 size_type reg,
                                 const getfem::mesh_fem *mflambda = 0,
                                 const model_real_plain_vector *lambda = 0,
                                 const model_real_plain_vector *w = 0,
                                 const std::string &varname = std::string(),
                                 const std::string &multname = std::string(),
                                 const std::string &wname = std::string());

    size_type add_slave_boundary(const getfem::mesh_im &mim, size_type reg,
                                 const std::string &varname,
                                 const std::string &multname = std::string(),
                                 const std::string &wname = std::string());
    

    size_type add_master_boundary(const getfem::mesh_im &mim,
                                  const getfem::mesh_fem *mfu,
                                  const model_real_plain_vector *U,
                                  size_type reg,
                                  const getfem::mesh_fem *mflambda = 0,
                                  const model_real_plain_vector *lambda = 0,   
                                  const model_real_plain_vector *w = 0,
                                  const std::string &varname = std::string(),
                                  const std::string &multname = std::string(),
                                  const std::string &wname = std::string());

    size_type add_master_boundary(const getfem::mesh_im &mim, size_type reg,
                                  const std::string &varname,
                                  const std::string &multname = std::string(),
                                  const std::string &wname = std::string());
    


    // The whole process of the computation of contact pairs
    // Contact pairs are seached for a certain boundary (master or slave,
    // depending on the contact algorithm) on the master ones. If contact pairs
    // are searched for a master boundary, self-contact is taken into account
    // if the flag 'self_contact' is set to 'true'. Self-contact is never taken
    // into account for a slave boundary.
    void compute_contact_pairs(void);

  };










}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONTACT_AND_FRICTION_COMMON_H__ */
