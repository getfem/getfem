/* -*- c++ -*- (enables emacs c++ mode)                                    */
#ifndef GETFEM_MESH_SLICES_H
#define GETFEM_MESH_SLICES_H

#include <bitset>
#include <getfem_mesh.h>
#include <getfem_fem.h>
#include <getfem_poly_composite.h>

namespace getfem {

  /* this should not stay here .. */
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
    convex_face(size_type cv_, size_type f_ = size_type(-1)) : cv(cv_), f(f_) {}
    convex_face() : cv(size_type(-1)), f(size_type(-1)) {}
  };
  typedef std::vector<convex_face> convex_face_ct;


  /**
     node data in a slice: contains both real position, and position in the reference convex
   */
  
  struct slice_node {
    typedef std::bitset<32> faces_ct; /// broken for convexes with more than 32 faces
    base_node pt, pt_ref;
    faces_ct faces; 
    slice_node() {}
    slice_node(const base_node& pt_, const base_node& pt_ref_) : pt(pt_), pt_ref(pt_ref_) {}
  };

  /**
     simplex data in a slice: just a list of slice_node ids.
   */
  struct slice_simplex {
    std::vector<size_type> inodes;
    size_type dim() const { return inodes.size()-1; }
    slice_simplex(size_type n) : inodes(n) {} 
    slice_simplex() : inodes(4) {}
    bool operator==(const slice_simplex& o) const { return inodes == o.inodes; }
    bool operator!=(const slice_simplex& o) const { return inodes != o.inodes; }
  };

  /**
     mesh slice: a list of nodes/simplexes, which can be seen as a P1 discontinuous
     mesh_fem on which the interpolation is very fast
  */
  class mesh_slice_cv_dof_data_base;
  class slicer;

  class mesh_slice {
  public:
    typedef std::vector<slice_node> cs_nodes_ct;
    typedef std::vector<slice_simplex> cs_simplexes_ct;
  private:
    /* nodes lists and simplexes lists for each convex of the original mesh */
    struct convex_slice {
      size_type cv_num; 
      dim_type cv_dim;
      dim_type fcnt, cv_nbfaces; // number of faces of the convex (fcnt also counts the faces created by the slicing of the convex)
      cs_nodes_ct nodes;
      cs_simplexes_ct simplexes;
    };
    typedef std::deque<convex_slice> cvlst_ct;

    /* keep track of the original mesh (hence it should not be destroyed before the slice) */
    const getfem_mesh& m;
    std::vector<size_type> simplex_cnt; // count simplexes of dimension 0,1,...,dim
    size_type points_cnt;
    cvlst_ct cvlst;
    size_type _dim;
    void do_slicing(size_type cv, bgeot::pconvex_ref cvr, slicer *ms, cs_nodes_ct cv_nodes, 
		    cs_simplexes_ct cv_simplexes, dal::bit_vector& splx_in);
  protected:
    size_type set_convex(size_type pos, size_type cv, bgeot::pconvex_ref cvr, 
                         cs_nodes_ct cv_nodes, cs_simplexes_ct cv_simplexes, 
                         dim_type fcnt, dal::bit_vector& splx_in);

  public:
    mesh_slice(const getfem_mesh& m);
    virtual ~mesh_slice() {}
    /**
       build a new mesh_slice given:
       @param m the mesh that is to be sliced
       @param ms a slicer
       @param nrefine number of refinments for each convex of the original mesh
       @param cvlst the list of convex numbers (or convex faces) of m that will 
       be taken into account for the slice
       @param def_mf_data an optional pointer to a mesh_slice_cv_dof_data object 
       describing a deformation to apply to the mesh before slicing it.
    */
    void build(slicer* ms, size_type nrefine, convex_face_ct& cvlst, 
	       mesh_slice_cv_dof_data_base *def_mf_data=0);

    /**
       build a new mesh slice given:
       @param sl an initial mesh_slice
       @param ms a slicer
    */
    void build_from_slice(const mesh_slice& sl, slicer* ms);

    /**
       build a new mesh_slice than can be used to interpolate a field
       on a fixed set of points.

       @param m the mesh that is to be sliced
       @param ms a slicer
       @param cvlst the list of convex numbers (or convex faces) of m that will 
       be taken into account for the slice       
       @param pts the list of points
     */
    void build_from_points(const std::vector<base_node>& pts, slicer* ms);
		      
    size_type nb_convex() const { return cvlst.size(); }
    size_type convex_num(size_type ic) const { return cvlst[ic].cv_num; }
    size_type dim() const { return _dim; }
    const getfem_mesh& linked_mesh() const { return m; }
    void nb_simplexes(std::vector<size_type>& c) const { c = simplex_cnt; }
    size_type nb_simplexes(size_type sdim) const { return simplex_cnt[sdim]; }
    size_type nb_points() const { return points_cnt; }
    const cs_nodes_ct& nodes(size_type ic) const { return cvlst[ic].nodes; }
    cs_nodes_ct& nodes(size_type ic) { return cvlst[ic].nodes; }
    void edges_mesh(getfem_mesh &edges_m, dal::bit_vector& slice_edges) const;
    const cs_simplexes_ct& simplexes(size_type ic) const { return cvlst[ic].simplexes; }
    size_type memsize() const;

    /** merges with another mesh slice */
    void merge(const mesh_slice& sl);

    /** interpolation of a mesh_fem on a slice (the mesh_fem
	and the slice must share the same mesh, of course)
    */
    template<typename V1, typename V2> void 
    interpolate(const getfem::mesh_fem &mf, const V1& U, V2& V) {
      _fem_precomp fprecomp;
      bgeot::stored_point_tab refpts;
      base_vector coeff;
      base_matrix G;
      base_node val(mf.get_qdim());
      typename V2::iterator out = V.begin();
      for (size_type i=0; i < nb_convex(); ++i) {
        size_type cv = convex_num(i);
        refpts.resize(nodes(i).size());
        for (size_type j=0; j < refpts.size(); ++j) refpts[j] = nodes(i)[j].pt_ref;
        pfem pf = mf.fem_of_element(cv);
	if (pf->need_G())transfert_to_G(G, mf.linked_mesh().points_of_convex(cv));
        fem_precomp_not_stored(pf, &refpts, fprecomp);
        
        ref_mesh_dof_ind_ct dof = mf.ind_dof_of_element(cv);
        coeff.resize(mf.nb_dof_of_element(cv));
        base_vector::iterator cit = coeff.begin();
        for (ref_mesh_dof_ind_ct::iterator it=dof.begin(); it != dof.end(); ++it, ++cit)
          *cit = U[*it];

        for (size_type j=0; j < refpts.size(); ++j) {
          pf->interpolation(&fprecomp, j,
                            G, mf.linked_mesh().trans_of_convex(cv),
                            coeff, val, mf.get_qdim());
          out = std::copy(val.begin(), val.end(), out);
        }
      }
    }
  };


  /**
     returns a list of "exterior" faces of a mesh (i.e. faces which are not shared by two convexes)
  */
  void  outer_faces_of_mesh(const getfem::getfem_mesh &m, const dal::bit_vector& cvlst, convex_face_ct& flist);



  /* stupid class in order to use any vector type for field data associated to mesh_fems
     in slices (used for slice deformation and isovalues) */
  class mesh_slice_cv_dof_data_base {
  public:
    const mesh_fem *pmf;    
    virtual void copy(size_type cv, base_vector& coeff) const = 0;
    virtual scalar_type maxval() const = 0;
    virtual ~mesh_slice_cv_dof_data_base() {}
  };

  /**
     use this structure to specify that the mesh must be deformed 
     (with a mesh_fem and an associated field)
     before the slicing 
  */
  template<typename VEC> class mesh_slice_cv_dof_data : public mesh_slice_cv_dof_data_base {
    const VEC &u;
  public:
    mesh_slice_cv_dof_data(const mesh_fem &mf_, VEC &u_) : u(u_) { pmf = &mf_; }
    virtual void copy(size_type cv, base_vector& coeff) const {
      coeff.resize(pmf->nb_dof_of_element(cv));
      ref_mesh_dof_ind_ct dof = pmf->ind_dof_of_element(cv);
      base_vector::iterator out = coeff.begin();
      for (ref_mesh_dof_ind_ct::iterator it=dof.begin(); it != dof.end(); ++it, ++out)
        *out = u[*it];
    }
    scalar_type maxval() const { return bgeot::vect_norminf(u); }
  };
  
  /**
     generic slicer class: given a list of slice_simplex/slice_node, it 
     build a news list of slice_simplex/slice_node, indicating which ones 
     are in the slice with the bit_vector splx_in
  */
  class slicer {
  public:
    static const float EPS;
    virtual void slice(size_type cv, dim_type& fcnt,
                       mesh_slice::cs_nodes_ct& nodes, mesh_slice::cs_simplexes_ct& splxs, 
                       dal::bit_vector& splx_in) = 0;    
    virtual ~slicer() {}
  };

  /**
     extraction of the boundary of a slice
  */
  class slicer_boundary : public slicer {
    slicer *A;
    std::vector<slice_node::faces_ct> convex_faces;
    bool test_bound(const slice_simplex& s, slice_node::faces_ct& fmask, 
                    const mesh_slice::cs_nodes_ct& nodes) const;
  public:
    slicer_boundary(const getfem_mesh& m, slicer *sA, const convex_face_ct& fbound);
    void slice(size_type cv, dim_type& fcnt,
	       mesh_slice::cs_nodes_ct& nodes, mesh_slice::cs_simplexes_ct& splxs, 
               dal::bit_vector& splx_in);
  };

  /**
     base class for general slices of a mesh (planar, sphere, cylinder,isosurface)
  */
  class slicer_volume : public slicer {
  protected:
    enum {VOLIN=-1, VOLBOUND=0, VOLOUT=+1}; 
    int orient;
    dal::bit_vector pt_in, pt_bin;    
    
    /* overload either 'prepare' or 'test_point' */
    virtual void prepare(size_type /*cv*/, const mesh_slice::cs_nodes_ct& nodes) {
      pt_in.clear(); pt_bin.clear();
      for (size_type i=0; i < nodes.size(); ++i) {
	bool in, bin; test_point(nodes[i].pt, in, bin);        
	pt_in[i] = bin || ((orient > 0) ? !in : in);
        pt_bin[i] = bin;
      }
    }
    virtual void test_point(const base_node&, bool& in, bool& bound) const { in=true; bound=true; }

    /* edge_intersect should always be overloaded */
    virtual scalar_type edge_intersect(size_type /*i*/, size_type /*j*/, 
				       const mesh_slice::cs_nodes_ct& /*nodes*/) const = 0;

    slicer_volume(int orient_) : orient(orient_) {}

    /* utility function */
    static scalar_type trinom(scalar_type a, scalar_type b, scalar_type c) {
      scalar_type delta = b*b - 4*a*c;
      if (delta < 0.) return 1./EPS;
      delta = sqrt(delta);
      scalar_type s1 = (-b - delta) / (2*a);
      scalar_type s2 = (-b + delta) / (2*a);
      if (dal::abs(s1-.5) < dal::abs(s2-.5)) return s1; else return s2;
    }
    void split_simplex(mesh_slice::cs_nodes_ct& nodes, 
                       mesh_slice::cs_simplexes_ct& splxs, dal::bit_vector& splx_in,
                       const slice_simplex s, /* s is NOT a reference, it is on purpose (push_back in the function)*/
                       size_type sstart);
  public:
    void slice(size_type cv, dim_type& fcnt,
	       mesh_slice::cs_nodes_ct& nodes, mesh_slice::cs_simplexes_ct& splxs, 
	       dal::bit_vector& splx_in);
  };

  /**
     slices a mesh with a half-space (or its boundary)
  */
  class slicer_half_space : public slicer_volume {
    base_node x0, n; /* normal directed from inside toward outside */
    void test_point(const base_node& P, bool& in, bool& bound) const {
      scalar_type s = 0.;
      for (unsigned i=0; i < P.size(); ++i) s += (P[i] - x0[i])*n[i];
      in = (s <= 0); bound = (s*s <= dal::sqr(EPS*bgeot::vect_norm2_sqr(P)));
    }
    scalar_type edge_intersect(size_type iA, size_type iB, const mesh_slice::cs_nodes_ct& nodes) const {
      const base_node& A=nodes[iA].pt;
      const base_node& B=nodes[iB].pt;
      scalar_type s1 = 0., s2 = 0.;
      for (unsigned i=0; i < A.size(); ++i) { s1 += (A[i] - B[i])*n[i]; s2 += (A[i]-x0[i])*n[i]; }
      if (dal::abs(s1) < EPS) return 1./EPS;
      else return s2/s1;
    }
  public:
    slicer_half_space(base_node x0_, base_node n_, int orient_) : 
      slicer_volume(orient_), x0(x0_), n(n_) {
      n *= (1./bgeot::vect_norm2(n));
    }
  };

  /**
     slices a mesh with a sphere (or its boundary)
  */
  class slicer_sphere : public slicer_volume {
    base_node x0;
    scalar_type R;
    void test_point(const base_node& P, bool& in, bool& bound) const {
      scalar_type R2 = bgeot::vect_dist2_sqr(P,x0);
      bound = (R2 >= (1-EPS)*R*R && R2 <= (1+EPS)*R*R);
      in = R2 <= R*R;
    }
    scalar_type edge_intersect(size_type iA, size_type iB, const mesh_slice::cs_nodes_ct& nodes) const {
      const base_node& A=nodes[iA].pt;
      const base_node& B=nodes[iB].pt;
      scalar_type a,b,c; // a*x^2 + b*x + c = 0
      a = bgeot::vect_norm2_sqr(B-A); if (a < EPS) return pt_bin[iA] ? 0. : 1./EPS;
      b = 2*bgeot::vect_sp(A-x0,B-A);
      c = bgeot::vect_norm2_sqr(A-x0)-R*R;
      return slicer_volume::trinom(a,b,c);
    }
  public:
    slicer_sphere(base_node x0_, scalar_type R_, int orient_) : 
      slicer_volume(orient_), x0(x0_), R(R_) {} //cerr << "slicer_volume, x0=" << x0 << ", R=" << R << endl; }
  };
  
  /**
     slices a mesh with a cylinder (or its boundary)
  */
  class slicer_cylinder : public slicer_volume {
    base_node x0, d;
    scalar_type R;
    void test_point(const base_node& P, bool& in, bool& bound) const {
      base_node N = P-x0;
      scalar_type axpos = bgeot::vect_sp(d, N);
      scalar_type dist2 = bgeot::vect_norm2_sqr(N) - dal::sqr(axpos);
      bound = dal::abs(dist2-R*R) < EPS;
      in = dist2 < R*R;
    }
    scalar_type edge_intersect(size_type iA, size_type iB, const mesh_slice::cs_nodes_ct& nodes) const {
      base_node F=nodes[iA].pt-x0; scalar_type Fd = bgeot::vect_sp(F,d);
      base_node D=nodes[iB].pt-nodes[iA].pt; scalar_type Dd = bgeot::vect_sp(D,d);
      scalar_type a = bgeot::vect_norm2_sqr(D) - dal::sqr(Dd); if (a < EPS) return pt_bin[iA] ? 0. : 1./EPS; assert(a> -EPS);
      scalar_type b = 2*(bgeot::vect_sp(F,D) - Fd*Dd);
      scalar_type c = bgeot::vect_norm2_sqr(F) - dal::sqr(Fd) - dal::sqr(R);
      return slicer_volume::trinom(a,b,c);
    }
  public:
    slicer_cylinder(base_node x0_, base_node x1_, scalar_type R_, int orient_) : 
      slicer_volume(orient_), x0(x0_), d(x1_-x0_), R(R_) {
      d /= bgeot::vect_norm2(d);
    }
  };


  /**
     extract an isosurface
  */
  class slicer_isovalues : public slicer_volume {
    const mesh_slice_cv_dof_data_base& mfU;
    scalar_type val;
    scalar_type val_scaling; /* = max(abs(U)) */
    std::vector<scalar_type> Uval;
    void prepare(size_type cv, const mesh_slice::cs_nodes_ct& nodes);
    scalar_type edge_intersect(size_type iA, size_type iB, const mesh_slice::cs_nodes_ct& /*nodes*/) const {
      assert(iA < Uval.size() && iB < Uval.size()); /* this should never happen! */
      if (((Uval[iA] < val) && (Uval[iB] > val)) ||
	  ((Uval[iA] > val) && (Uval[iB] < val)))
	return (val-Uval[iA])/(Uval[iB]-Uval[iA]);
      else return 1./EPS;
    }
  public:
    /* orient = -1: u(x) <= val, 0: u(x) == val, +1: u(x) >= val */
    slicer_isovalues(const mesh_slice_cv_dof_data_base& mfU_, scalar_type val_, int orient_) : 
      slicer_volume(orient_), mfU(mfU_), val(val_) {
      if (mfU.pmf->get_qdim() != 1) DAL_THROW(dal::failure_error, "can't compute isovalues of a vector field !");
      val_scaling = mfU.maxval();
    }
  };

  /**
     union of two slices
  */
  class slicer_union : public slicer {
    slicer *A, *B;
  public:
    slicer_union(slicer *sA, slicer *sB) : A(sA), B(sB) {}
    /*    void test_point(const base_node& P, bool& in, bool& bound) const {
      bool inA, boundA, inB, boundB;
      A->test_point(P,inA,boundA); B->test_point(P,inB,boundB);
      bound = (boundA && !inB) || (boundB && !inA) || (boundA && boundB);
      in = inA || inB;
      }*/
    void slice(size_type cv, dim_type& fcnt,
	       mesh_slice::cs_nodes_ct& nodes, mesh_slice::cs_simplexes_ct& splxs, 
	       dal::bit_vector& splx_in);
  };

  /**
     intersection of two slices
  */
  class slicer_intersect : public slicer {
    slicer *A, *B;
  public:
    slicer_intersect(slicer *sA, slicer *sB) : A(sA), B(sB) {}
    /*    void test_point(const base_node& P, bool& in, bool& bound) const {
      bool inA, boundA, inB, boundB;
      A->test_point(P,inA,boundA); B->test_point(P,inB,boundB);
      bound = (boundA && inB) || (boundB && inA) || (boundA && boundB);
      in = inA && inB;
      }*/
    void slice(size_type cv, dim_type& fcnt,
	       mesh_slice::cs_nodes_ct& nodes, mesh_slice::cs_simplexes_ct& splxs, 
	       dal::bit_vector& splx_in);
  };

  /**
     complementary of a slice
  */
  class slicer_complementary : public slicer {
    slicer *A;
  public:
    slicer_complementary(slicer *sA) : A(sA) {}
    /*    void test_point(const base_node& P, bool& in, bool& bound) const {
      A->test_point(P,in,bound); in = !in;
      }*/
    void slice(size_type cv, dim_type& fcnt,
	       mesh_slice::cs_nodes_ct& nodes, mesh_slice::cs_simplexes_ct& splxs, 
	       dal::bit_vector& splx_in);
  };
  

  /** 
      this slicer does nothing! 
  */
  class slicer_none : public slicer {
  public:
    slicer_none() {}
    void slice(size_type /*cv*/, dim_type& /*fcnt*/, mesh_slice::cs_nodes_ct& /*nodes*/, 
	       mesh_slice::cs_simplexes_ct& /*splxs*/, dal::bit_vector& /*splx_in*/) {}
  };
  

}

#endif /*GETFEM_MESH_SLICES_H*/
