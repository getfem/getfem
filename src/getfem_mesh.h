// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 1999-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

/**@file getfem_mesh.h
   @author  Yves Renard <Yves.Renard@insa-toulouse.fr>
   @date November 05, 1999.
   @brief Define a getfem::getfem_mesh object.
*/

#ifndef GETFEM_MESH_H__
#define GETFEM_MESH_H__

#include <bitset>
#include <dal_shared_ptr.h>
#include <ftool.h>
#include <bgeot_mesh.h>
#include <bgeot_convex_ref.h>
#include <bgeot_geotrans_inv.h>
#include <linkmsg.h>
#include <getfem_context.h>
#include <getfem_mesh_region.h>

namespace getfem {

  /* ********************************************************************* */
  /*								   	   */
  /*	I. Classes of message                                		   */
  /*									   */
  /* ********************************************************************* */

  class mesh;

  struct MESH_CLEAR { /* clear message for the structure.                  */
    enum { ID = 0 }; 
  };
  struct MESH_DELETE { /* clear message for the structure.                  */
    enum { ID = 1 }; 
    const mesh *msh;
    MESH_DELETE(const mesh &m) : msh(&m) {}
  };
  struct MESH_ADD_CONVEX { 
    enum { ID = 5 }; 
    size_t icv;
    MESH_ADD_CONVEX(size_t i) { icv = i; }
    MESH_ADD_CONVEX(void) {}
  };
  struct MESH_SUP_CONVEX { 
    enum { ID = 6 }; 
    size_t icv;
    MESH_SUP_CONVEX(size_t i) { icv = i; }
    MESH_SUP_CONVEX(void) {}
  };
  struct MESH_SWAP_CONVEX { 
    enum { ID = 7 }; 
    size_t icv1, icv2;
    MESH_SWAP_CONVEX(size_t i, size_t j) { icv1 = i; icv2 = j; }
    MESH_SWAP_CONVEX(void) {}
  };
  struct MESH_REFINE_CONVEX { 
    enum { ID = 8 }; 
    size_t icv;
    std::vector<size_type> &sub_cv_list;
    bool is_refine; // true : refine, false unrefine.
    MESH_REFINE_CONVEX(size_t i, std::vector<size_type> &s, bool r)
      : icv(i), sub_cv_list(s), is_refine(r) { }
    // MESH_SWAP_CONVEX(void) {}
  };

  /** base class for objects that receive notification messages from a
      mesh. */
  class mesh_receiver : public lmsg::virtual_linkmsg_receiver {
    public :

      virtual void receipt(const MESH_CLEAR           &)
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_DELETE          &)
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_ADD_CONVEX      &) 
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_SUP_CONVEX      &) 
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_SWAP_CONVEX     &)
      { DAL_THROW(internal_error, "internal error");}
      virtual void receipt(const MESH_REFINE_CONVEX   &)
      { DAL_THROW(internal_error, "internal error");}

      virtual ~mesh_receiver() {}
  };

  class approx_integration;

  /* ********************************************************************* */
  /*								   	   */
  /*	II. Class getfem::mesh                                 		   */
  /*									   */
  /* ********************************************************************* */
  
  /**@addtogroup mesh*/
  /**@{*/
  
  /** Describe a mesh (collection of convexes and points).  Note that
      mesh object have no copy constructor, use
      mesh::copy_from instead.  This class inherits from
      bgeot::mesh<base_node> and bgeot::mesh_structure. Compared to
      bgeot::mesh, It provides some additional methods for convenience
      (add_simplex_by_points etc...), and it is able to send
      "messages" to objects that depend on it (for example
      getfem::mesh_fem) when something in the mesh has been
      modified. For example, when a convex is removed from the mesh,
      the mesh_fem that uses this mesh automagically removes any finite
      element on the convex, and renumbers the dof.

      Points and convex numbering:

      The numbering of points and convexes is not dynamical, which
      means that when a point or a convex has been removed from the mesh,
      there might be holes in the numbering. To loop over the set of
      valid points in the mesh, one should use

      @code
      for (dal::bv_visitor ip(mesh.points_index()); !ip.finished(); ++ip) {
      ...
      }
      @endcode

      instead of 
      
      @code
      for (size_type ip = 0; ip < mesh.nb_points(); ++ip) {
      }
      @endcode
      (same thing for the convexes, always use convex_index()).
  */

  class mesh : virtual public dal::static_stored_object,
		      public bgeot::basic_mesh,
		      public context_dependencies {
  public :
    
    typedef lmsg::linkmsg_sender<mesh_receiver> msg_sender;
    typedef dal::approx_less<scalar_type> val_comp;
    typedef bgeot::basic_mesh::pt_comp pt_comp;
    typedef bgeot::basic_mesh::PT_TAB PT_TAB;
    typedef bgeot::mesh_structure::ind_cv_ct ind_cv_ct;
    typedef bgeot::mesh_structure::ind_set ind_set;
    typedef bgeot::mesh_structure::ind_pt_face_ct ind_pt_face_ct;
    typedef bgeot::mesh_structure::point_ct point_ct;
    typedef bgeot::basic_mesh::ref_mesh_pt_ct ref_mesh_pt_ct;
    typedef bgeot::basic_mesh::ref_mesh_face_pt_ct ref_mesh_face_pt_ct;
    typedef bgeot::convex<base_node, ref_mesh_pt_ct> ref_convex;
    
    
  protected :
    /*
     * When a new field is added, do NOT forget to add it in copy_from method!
     */
    
    mutable msg_sender lkmsg; /* gestionnaire de msg.                    */
    
    dal::dynamic_array<mesh_region> cvf_sets;
    dal::bit_vector valid_cvf_sets;
    void handle_region_refinement(size_type, const std::vector<size_type> &,
				  bool);
    
    mutable bool cuthill_mckee_uptodate;
    mutable std::vector<size_type> cmk_order; // cuthill-mckee
    void init(void);

#if GETFEM_PARA_LEVEL > 1
    mutable bool modified;
    mutable mesh_region mpi_region;
    mutable dal::dynamic_array<mesh_region> mpi_sub_region;
    mutable dal::bit_vector valid_sub_regions;

    void touch(void) { 
      modified = true; cuthill_mckee_uptodate = false;
      context_dependencies::touch();
    }    
    void compute_mpi_region(void) const ;
    void compute_mpi_sub_region(size_type) const;

  public :
    
    const mesh_region& get_mpi_region(void) const
    { if (modified) compute_mpi_region(); return mpi_region; }
    const mesh_region& get_mpi_sub_region(size_type n) const {
      if (modified) compute_mpi_region();
      if (n == size_type(-1)) return mpi_region;
      if (!(valid_sub_regions.is_in(n))) compute_mpi_sub_region(n);
      return mpi_sub_region[n];
    }
    void intersect_with_mpi_region(mesh_region &rg) const; 
#else
    void touch(void)
    { cuthill_mckee_uptodate = false; context_dependencies::touch(); }
  public :
    const mesh_region get_mpi_region(void) const 
    { return mesh_region::all_convexes(); }
    const mesh_region get_mpi_sub_region(size_type n) const {
      if (n == size_type(-1)) return get_mpi_region();
      return cvf_sets[n];
    }
    void intersect_with_mpi_region(mesh_region &) const {}
#endif
    
    /// Constructor.
    mesh(dim_type NN = dim_type(-1)); 
    mesh(const bgeot::basic_mesh &m);
    double eps(void) const { return eps_p; }
    msg_sender &lmsg_sender(void) const { return lkmsg; }
    void update_from_context(void) const {}
    /// Mesh dimension.
    dim_type dim(void) const { return dimension; }
    /// Return the array of PT.
    const PT_TAB &points(void) const { return pts; }
    /// Return the array of PT.
    PT_TAB &points(void) { return pts; }

    /// Return a (pseudo)container of the points of a given convex 
    ref_mesh_pt_ct points_of_convex(size_type ic) const {
      const ind_cv_ct &rct = ind_points_of_convex(ic);
      return ref_mesh_pt_ct(pts.begin(), rct.begin(), rct.end());
    } 
    
    /// Return a (pseudo)container of points of face of a given convex 
    ref_mesh_face_pt_ct points_of_face_of_convex(size_type ic,
						 size_type f) const {
      ind_pt_face_ct rct = ind_points_of_face_of_convex(ic,f);
      return ref_mesh_face_pt_ct(pts.begin(), rct.begin(), rct.end());
    }

    /// return a bgeot::convex object for the convex number ic. 
    ref_convex convex(size_type ic) const
    { return ref_convex(structure_of_convex(ic), points_of_convex(ic)); }

    /** Add the point pt to the mesh and return the index of the
	point. 

	If the point is too close to an existing point, the
	function does not create a new point, and returns the index of the
	already existing point.
	@param pt the point coordinates.
	@param norepeat should not be used.
    */
    size_type add_point(const base_node &pt, bool norepeat = true);
    /// Give the number of geometrical nodes in the mesh.
    size_type nb_points(void) const { return pts.card(); }
    /// Return the points index
    const dal::bit_vector &points_index(void) const { return pts.index(); }
    /// Delete the point of index i from the mesh.
    void sup_point(size_type i) { if (!is_point_valid(i)) pts.sup(i); }
    /// Swap the indexes of points of index i and j in the whole structure.
    void swap_points(size_type i, size_type j)
    { if (i != j) { pts.swap(i,j); mesh_structure::swap_points(i,j); } }
    /** Search a point given its coordinates.
	@param pt the point that is searched.
	@return the point index if pt was found in (or approximatively in)
	the mesh nodes, size_type(-1) if not found.
    */
    size_type search_point(const base_node &pt) const
    { return pts.search(pt); }
    /** Return the bgeot::geometric_trans attached to a convex.
	@param ic the convex number.
    */
    bgeot::pgeometric_trans trans_of_convex(size_type ic) const { 
      if (!(trans_exists[ic]))
	DAL_THROW(internal_error, "internal error");
      return gtab[ic]; 
    }
    
    /** Add a convex to the mesh.
	This methods assume that the convex nodes have already been 
	added to the mesh.
	@param pgt the geometric transformation of the convex.
	@param ipts an iterator to a set of point index.
	@return the number of the new convex.
     */
    template<class ITER>
    size_type add_convex(bgeot::pgeometric_trans pgt, ITER ipts) { 
      bool present;
      size_type i = bgeot::mesh_structure::add_convex(pgt->structure(),
						       ipts, &present);
      gtab[i] = pgt; trans_exists[i] = true;
      if (!present) { lmsg_sender().send(MESH_ADD_CONVEX(i)); touch(); }
      return i;
    }
    
    /** Add a convex to the mesh, given a geometric transformation and a
	list of point coordinates. 

	As a side-effect, the points are also added to the mesh (if
	they were not already in the mesh).

	@param pgt the geometric transformation of the convex.
	@param ipts an iterator on a set of getfem::base_node.
	@return the number of the new convex.
    */
    template<class ITER>
    size_type add_convex_by_points(bgeot::pgeometric_trans pgt, ITER ipts);
    
    /** Add a simplex to the mesh, given its dimension and point numbers.
     */
    template<class ITER>
    size_type add_simplex(dim_type di, ITER ipts)
    { return add_convex(bgeot::simplex_geotrans(di, 1), ipts); }
    /** Add a simplex to the mesh, given its dimension and point coordinates.
	@see add_convex_by_points.
     */
    template<class ITER>
    size_type add_simplex_by_points(dim_type dim, ITER ipts);
    /** Add a segment to the mesh, given the point id of its vertices. */
    size_type add_segment(size_type a, size_type b);
    /** Add a segment to the mesh, given the coordinates of its vertices. */
    size_type add_segment_by_points(const base_node &pt1,
				    const base_node &pt2) {
      size_type ipt1 = add_point(pt1);
      size_type ipt2 = add_point(pt2);
      return add_segment(ipt1, ipt2);
    }
    /** Add a triangle to the mesh, given the point id of its vertices. */
    size_type add_triangle(size_type a,size_type b, size_type c);
    /** Add a triangle to the mesh, given the coordinates of its vertices. */
    size_type add_triangle_by_points(const base_node &p1,
				     const base_node &p2,
				     const base_node &p3);
    /** Add a tetrahedron to the mesh, given the point id of its vertices. */
    size_type add_tetrahedron(size_type a,
			      size_type b, size_type c, size_type d);
    /** Add a tetrahedron to the mesh, given the coordinates of its vertices. */
    size_type add_tetrahedron_by_points(const base_node &p1,
					const base_node &p2,
					const base_node &p3,
					const base_node &p4);
    /** Add a parallelepiped to the mesh. 
	@param di dimension of the parallelepiped
	@param ipts iterator on the list of point id.
	@return the number of the new convex.
    */
    template<class ITER>
    size_type add_parallelepiped(dim_type di, const ITER &ipts);
    /** Add a parallelepiped to the mesh. 
	@param di dimension of the parallelepiped
	@param ps iterator on the list of point coordinates.
	@return the number of the new convex.
    */
    template<class ITER>
    size_type add_parallelepiped_by_points(dim_type di, const ITER &ps);
    /* Add a parallelepiped of dimension dim to the
     *          mesh. org is the point of type base_node representing
     *          the origine and "it" is an iterator on a list of
     *          vectors of type base_vector.
     *          Return the index of the convex in the mesh.

    template<class ITER>
    size_type add_parallelepiped_by_vectors(dim_type di,
				  const base_node &org, const ITER &vects);
     */    
    /** Add a prism to the mesh. 
	@param di dimension of the prism
	@param ipts iterator on the list of point id.
	@return the number of the new convex.
    */
    template<class ITER>
    size_type add_prism(dim_type di, const ITER &ipts);
    
    /** Add a prism to the mesh. 
	@param di dimension of the prism
	@param ps iterator on the list of point coordinates.
	@return the number of the new convex.
    */
    template<class ITER>
    size_type add_prism_by_points(dim_type di, const ITER &ps);
    
    /// Delete the convex of index ic from the mesh.
    void sup_convex(size_type ic, bool sup_points = false);
    /** Swap the indexes of the convex of indexes i and j 
     *          in the whole structure.
     */
    void swap_convex(size_type i, size_type j);
    
    /** Return the normal of the given convex face, evaluated at the point pt.
	@param ic the convex number.
	@param f the face number.
	@param pt the point at which the normal is taken, in the
	reference convex. This point should of course be on the
	correspounding face of the reference convex, except if the
	geometric transformation is linear: in that case, the normal
	constant.
	@return the face normal.
    */
    base_small_vector normal_of_face_of_convex(size_type ic, short_type f,
					       const base_node &pt) const;
    /* Return the normal of the given convex face.

       @param ic the convex number.
       @param f the face number.
       @param n local point index in the reference convex, at which the normal will be evaluated -- should be faster than the function above since it uses a bgeot::geotrans_precomp .
       
       @return the face normal.
    */
    base_small_vector normal_of_face_of_convex(size_type ic, short_type f,
					       size_type n=0) const;
    /** Return a local basis for the convex face.
       @param ic the convex number.
       @param f the face number.
       @param pt the point coordinates.
       @return a matrix with the normal vector in its first column, and the
       tangent vectors in the other columns.
    */ 
    base_matrix local_basis_of_face_of_convex(size_type ic, short_type f,
					      const base_node &pt) const;
    /** Return a local basis for the convex face.
       @param ic the convex number.
       @param f the face number.
       @param n local point index in the reference convex, at which the basis will be evaluated -- should be faster than the function above since it uses a bgeot::geotrans_precomp .
       @return a matrix with the normal vector in its first column, and the
       tangent vectors in the other columns.
    */ 
    base_matrix local_basis_of_face_of_convex(size_type ic, short_type f,
					      size_type n) const;
    /** Return an estimate of the convex quality (0 <= Q <= 1). */
    scalar_type convex_quality_estimate(size_type ic) const;
    /** Return an estimate of the convex area. 
	@param ic the convex number
	@param degree the degree of the approximation integration
	method used to compute the area.
    */
    scalar_type convex_area_estimate(size_type ic, size_type degree=2) const;
    /** Return an estimate of the convex largest dimension. @see getfem::convex_quality_estimate */
    scalar_type convex_radius_estimate(size_type ic) const;
    /** Return an estimate of the convex smallest dimension. @see getfem::convex_radius_estimate */
    scalar_type minimal_convex_radius_estimate() const;
    /** Apply the given translation to each mesh node. */
    void translation(base_small_vector);
    /** apply the given matrix transformation to each mesh node. */
    void transformation(base_matrix);
    /** Return the bounding box [Pmin - Pmax] of the mesh. */
    void bounding_box(base_node& Pmin, base_node &Pmax) const;
    /** Return the region of index 'id'. Regions stored in mesh are
	automagically created on their first use. Moreover, they are
	updated when the mesh is modified (i.e. when a convex is
	removed from the mesh, it is also removed from all regions of
	the mesh.
    */
    const mesh_region region(size_type id) const { 
      if (has_region(id)) return cvf_sets[id]; 
      else return mesh_region(const_cast<mesh&>(*this),id);
    }
    /* Return a reference such that operator= works as expected */
    mesh_region &region(size_type id) { 
      if (!has_region(id)) 
	/* will be added into valid_cvf_sets 
	   later via mesh_region::maybe_notify_parent_mesh */
	cvf_sets[id] = mesh_region(*this,id); 
      return cvf_sets[id];
    }
    /** Return true if the region of index 's' exists in the mesh */
    bool has_region(size_type s) const { return valid_cvf_sets[s]; }
    void add_face_to_set(size_type s, size_type c, short_type f) IS_DEPRECATED;
    const dal::bit_vector &regions_index() const { return valid_cvf_sets; }
    
    /*
    void sup_face_from_set(size_type b, size_type c, short_type f)
    { if (valid_cvf_sets[b]) { cvf_sets[b].sup(c,f); touch(); } }*/
    /** Remove the region of index b. */
    void sup_region(size_type b) {
      if (valid_cvf_sets[b])
	{ cvf_sets[b].clear(); valid_cvf_sets.sup(b); touch(); }
    }
    /** Remove all references to a convex from all regions stored in the mesh.
     @param cv the convex number.*/
    void sup_convex_from_regions(size_type cv);
    /** Pack the mesh : renumber convexes and nodes such that there
	is no holes in their numbering. Do NOT do the Cuthill-McKee. */
    void optimize_structure(void);
    /// Return the list of convex IDs for a Cuthill-McKee ordering
    const std::vector<size_type> &cuthill_mckee_ordering() const;
    /// Erase the mesh.
    void clear(void);
    /** Write the mesh to a file. The format is getfem-specific.
	@param name the file name.
    */
    void write_to_file(const std::string &name) const;
    /** Write the mesh to a text stream.
	@param ost the stream.
    */
    void write_to_file(std::ostream &ost) const;
    /** Load the mesh from a file. 
	@param name the file name.
	@see getfem::import_mesh.
    */
    void read_from_file(const std::string &name);
    /** Load the mesh from a text stream.
	@param ist the text stream.
	@see getfem::import_mesh.
    */
    void read_from_file(std::istream &ist);
    /** Clone a mesh */
    void copy_from(const mesh& m); /* might be the copy constructor */
    size_type memsize(void) const;
    ~mesh() {
      lmsg_sender().send(MESH_DELETE(*this));
      if (Bank_info) delete Bank_info;
    }

    friend class mesh_region;
  private:
    void swap_convex_in_regions(size_type c1, size_type c2);
    void touch_from_region(size_type id) { valid_cvf_sets.add(id); touch(); }
    void to_edges() {} /* to be done, the to_edges of mesh_structure does   */
                       /* not handle geotrans */

    //
    // Bank refinement
    //

    // TODO : - mise a jour des mesh_regions ;
    //        - dialogue avec mesh_im et mesh_fem

    struct green_simplex {
      bgeot::pgeometric_trans pgt;
      std::vector<size_type> sub_simplices;
      bgeot::convex<base_node> cv;
      std::vector<size_type> ipt_loc;
    };
    struct edge { 
      size_type i0, i1, i2;
      bool operator <(const edge &e) const;
      edge(size_type a, size_type b) : i0(0), i1(a), i2(b) {}
      edge(size_type a, size_type b, size_type c) : i0(a), i1(b), i2(c) {}
    };
    typedef std::set<edge> edge_set;

    struct Bank_info_struct {
      dal::bit_vector is_green_simplex; // indices of green simplices.
      std::map<size_type, size_type> num_green_simplex;
      dal::dynamic_tas<green_simplex> green_simplices;
      edge_set edges;
    };

    Bank_info_struct *Bank_info;

    void Bank_convex_with_edge(size_type, size_type,
			       std::vector<size_type> &);
    bool Bank_is_convex_having_points(size_type,
				      const std::vector<size_type> &);
    void Bank_sup_convex_from_green(size_type);
    void Bank_swap_convex(size_type, size_type);
    void Bank_build_first_mesh(mesh &, size_type);
    void Bank_basic_refine_convex(size_type);
    void Bank_refine_normal_convex(size_type);
    size_type Bank_test_and_refine_convex(size_type, dal::bit_vector &,
					  bool = true);
    void Bank_build_green_simplexes(size_type, std::vector<size_type> &);

  public :
    
    /** Use the Bank strategy to refine some convexes. */
    void Bank_refine(dal::bit_vector);

  };

  inline void mesh::add_face_to_set(size_type s, size_type c, short_type f) {
    region(s).add(c,f);
  }

  /**
   * build a N+1 dimensions mesh from a N-dimensions mesh by extrusion. 
   */
  void extrude(const mesh& in, mesh& out, unsigned nb_layers);

 template<class ITER>
    size_type mesh::add_convex_by_points(bgeot::pgeometric_trans pgt,
					                           ITER ipts)
  {
    short_type nb = pgt->nb_points();
    std::vector<size_type> ind(nb);
    for (short_type i = 0; i < nb; ++ipts, ++i) ind[i] = add_point(*ipts);
    return add_convex(pgt, ind.begin());
  }

  template<class ITER>
   size_type mesh::add_simplex_by_points(dim_type di, ITER ipts)
  { return add_convex_by_points(bgeot::simplex_geotrans(di, 1), ipts); }

  template<class ITER>
    size_type mesh::add_parallelepiped(dim_type di, const ITER &ipts)
  { return add_convex(bgeot::parallelepiped_geotrans(di, 1), ipts); }

  template<class ITER>
    size_type mesh::add_parallelepiped_by_points
    (dim_type di, const ITER &ps)
  { return add_convex_by_points(bgeot::parallelepiped_geotrans(di, 1), ps); }

  template<class ITER>
    size_type mesh::add_prism(dim_type di, const ITER &ipts)
  { return add_convex(bgeot::prism_geotrans(di, 1), ipts); }

  template<class ITER>
    size_type mesh::add_prism_by_points
    (dim_type di, const ITER &ps)
  { return add_convex_by_points(bgeot::prism_geotrans(di, 1), ps); }

  typedef mesh *pmesh;

  /** rough estimate of the convex area.
      @param pgt the geometric transformation.
      @param pts the convex nodes.
      @param pai the approximate integration used for the computation
      of the convex area.
   */
  scalar_type convex_area_estimate(bgeot::pgeometric_trans pgt,
				   const base_matrix& pts, 
				   const approx_integration* pai);

  /** rough estimate of the maximum value of the condition 
   * number of the jacobian of the geometric transformation */
  scalar_type convex_quality_estimate(bgeot::pgeometric_trans pgt,
				      const base_matrix& pts);

  /** rough estimate of the radius of the convex using the largest eigenvalue
   * of the jacobian of the geometric transformation */
  scalar_type convex_radius_estimate(bgeot::pgeometric_trans pgt,
				     const base_matrix& pts);

  /* stores a convex face. if f == -1, it is the whole convex.             */
  struct convex_face  {
    size_type cv;
    size_type f;    
    inline bool operator < (const convex_face &e) const
    {
      if (cv < e.cv) return true; if (cv > e.cv) return false; 
      if (f < e.f) return true; else if (f > e.f) return false;
      return false;
    }
    bool is_face() const { return f != size_type(-1); }
    convex_face(size_type cv_, size_type f_=size_type(-1)) : cv(cv_), f(f_) {}
    convex_face() : cv(size_type(-1)), f(size_type(-1)) {}
  } IS_DEPRECATED;
  typedef std::vector<convex_face> convex_face_ct;

  /** returns a list of "exterior" faces of a mesh
   * (i.e. faces which are not shared by two convexes) 
   * + convexes whose dimension is smaller that m.dim()
   */
  void  outer_faces_of_mesh(const mesh &m, 
			const dal::bit_vector& cvlst, convex_face_ct& flist);

  inline void outer_faces_of_mesh(const mesh &m, convex_face_ct& flist)
  { outer_faces_of_mesh(m,m.convex_index(),flist); }

  void  outer_faces_of_mesh(const mesh &m, const mesh_region &cvlst,
			    mesh_region &flist);

  inline void  outer_faces_of_mesh(const mesh &m, mesh_region &flist)
  { outer_faces_of_mesh(m,m.convex_index(),flist); }
  ///@}
}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MESH_H__  */
