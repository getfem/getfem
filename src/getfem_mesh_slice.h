// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_mesh_slice.h : 
//           
// Date    : May 2003.
// Author  : Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2003-2005 Julien Pommier
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================

#ifndef GETFEM_MESH_SLICE_H
#define GETFEM_MESH_SLICE_H

#include <getfem_mesh_slicers.h>

namespace getfem {
  class slicer_build_stored_mesh_slice;

  class stored_mesh_slice {
  protected:
    /* nodes lists and simplexes lists for each convex of the original mesh */
    struct convex_slice {
      size_type cv_num; 
      dim_type cv_dim;
      dim_type fcnt, cv_nbfaces; // number of faces of the convex (fcnt also counts the faces created by the slicing of the convex)
      mesh_slicer::cs_nodes_ct nodes;
      mesh_slicer::cs_simplexes_ct simplexes;
      size_type global_points_count;
    };
    typedef std::deque<convex_slice> cvlst_ct;

    struct merged_node_t {
      const slice_node *P;
      unsigned pos; /* [0..nb_points()-1] */
    };
    /* these three arrays allow the association of a global point
       number to a merged point number, and of a merged point to a
       list of global points */
    mutable std::vector<merged_node_t> merged_nodes;     /* size = points_cnt */
    mutable std::vector<size_type> merged_nodes_idx;     /* size = nb of merged points + 1 */
    mutable std::vector<size_type> to_merged_index;      /* size = points_cnt */
    mutable bool merged_nodes_available;

    /* keep track of the original mesh (hence it should not be destroyed before the slice) */
    const getfem_mesh *poriginal_mesh;
    std::vector<size_type> simplex_cnt; // count simplexes of dimension 0,1,...,dim
    size_type points_cnt;
    cvlst_ct cvlst;
    size_type dim_;
    std::vector<size_type> cv2pos; // convex id -> pos in cvlst
    friend class slicer_build_stored_mesh_slice;
    friend class mesh_slicer;
  public:
    stored_mesh_slice() : poriginal_mesh(0), points_cnt(0), dim_(size_type(-1)) { }
    stored_mesh_slice(const getfem::getfem_mesh& m, size_type nrefine = 1) 
      : poriginal_mesh(0), points_cnt(0), dim_(size_type(-1)) { 
      this->build(m, slicer_none(), nrefine);
    };
    virtual ~stored_mesh_slice() {}
    size_type nb_convex() const { return cvlst.size(); }
    size_type convex_num(size_type ic) const { return cvlst[ic].cv_num; }
    void set_dim(size_type newdim);
    size_type dim() const { return dim_; }
    const getfem_mesh& linked_mesh() const { return *poriginal_mesh; }
    void nb_simplexes(std::vector<size_type>& c) const { c = simplex_cnt; }
    size_type nb_simplexes(size_type sdim) const { return simplex_cnt[sdim]; }
    size_type nb_points() const { return points_cnt; }
    const mesh_slicer::cs_nodes_ct& nodes(size_type ic) const { return cvlst[ic].nodes; }
    mesh_slicer::cs_nodes_ct& nodes(size_type ic) { return cvlst[ic].nodes; }
    const mesh_slicer::cs_simplexes_ct& simplexes(size_type ic) const { return cvlst[ic].simplexes; }
    size_type memsize() const;
    void clear() { poriginal_mesh = 0; cvlst.clear(); points_cnt = 0; 
      dim_ = size_type(-1); cv2pos.clear(); simplex_cnt.clear(); clear_merged_nodes(); }
    /** merges with another mesh slice */
    void merge(const stored_mesh_slice& sl);

    /** build a list of merged nodes, i.e. nodes which have the same
	geometrical location but were extracted from two different
	convexes will be considered as one same node. Use for
	exportation purposes, as VTK and OpenDX do not like
	'discontinuous' meshes */
    void merge_nodes() const;
    size_type nb_merged_nodes() const 
      { return merged_nodes_idx.size() - 1; }
    const base_node merged_point(size_type i_merged) const 
      { return merged_nodes[merged_nodes_idx[i_merged]].P->pt; }
    size_type merged_index(size_type ic, size_type ipt) const
      { return to_merged_index[global_index(ic,ipt)]; }
    size_type global_index(size_type ic, size_type ipt) const 
      { return cvlst[ic].global_points_count+ipt; }
    size_type merged_point_cnt(size_type i_merged) const 
      { return merged_nodes_idx[i_merged+1] - merged_nodes_idx[i_merged]; }
    std::vector<merged_node_t>::const_iterator 
    merged_point_nodes(size_type i_merged) const { 
      return merged_nodes.begin() + merged_nodes_idx[i_merged]; 
    }
    void clear_merged_nodes() const;

    /* extract the list of mesh edges into 'edges' (size = 2* number
       of edges). 'slice_edges' indicates which one were created after
       slicing.  The from_merged_nodes flag may be used if you want to
       use (and merge common edges according to) the merged points */
    void get_edges(std::vector<size_type> &edges,
		   dal::bit_vector &slice_edges,
		   bool from_merged_nodes) const;

    void set_convex(size_type cv, bgeot::pconvex_ref cvr, 
	       mesh_slicer::cs_nodes_ct cv_nodes, 
	       mesh_slicer::cs_simplexes_ct cv_simplexes, 
		    dim_type fcnt, const dal::bit_vector& splx_in);

    void build(const getfem::getfem_mesh& m, const slicer_action &a, 
	       size_type nrefine = 1) { build(m,&a,0,0,nrefine); }
    void build(const getfem::getfem_mesh& m, const slicer_action &a, const slicer_action &b, 
	       size_type nrefine = 1) { build(m,&a,&b,0,nrefine); }
    void build(const getfem::getfem_mesh& m, const slicer_action &a, const slicer_action &b, const slicer_action &c, 
	       size_type nrefine = 1) { build(m,&a,&b,&c,nrefine); }
    void build(const getfem::getfem_mesh& m, const slicer_action *a, const slicer_action *b, const slicer_action *c, 
	       size_type nrefine);
    
      
    /* apply the listed slicer_action(s) to the slice object -- 
       the stored_mesh_slice is not modified */
    void replay(slicer_action &a) const { replay(&a,0,0); }
    void replay(slicer_action &a, slicer_action &b) const { replay(&a, &b, 0); }
    void replay(slicer_action &a, slicer_action &b, slicer_action &c) const { replay(&a, &b, &c); }
    void replay(slicer_action *a, slicer_action *b, slicer_action *c) const;

    /** interpolation of a mesh_fem on a slice (the mesh_fem
	and the slice must share the same mesh, of course)
    */
    template<typename V1, typename V2> void 
    interpolate(const getfem::mesh_fem &mf, const V1& U, V2& V) const {
      typedef typename gmm::linalg_traits<V2>::value_type T;
      bgeot::stored_point_tab refpts;
      std::vector<std::vector<T> > coeff;
      base_matrix G;
      size_type qdim = mf.get_qdim();
      size_type qqdim = gmm::vect_size(U)/mf.nb_dof();
      size_type pos = 0;
      coeff.resize(qqdim);
      for (size_type i=0; i < nb_convex(); ++i) {
        size_type cv = convex_num(i);
        refpts.resize(nodes(i).size());
        for (size_type j=0; j < refpts.size(); ++j) refpts[j] = nodes(i)[j].pt_ref;
        pfem pf = mf.fem_of_element(cv);
	if (pf->need_G()) 
	  bgeot::vectors_to_base_matrix(G, mf.linked_mesh().points_of_convex(cv));
	fem_precomp_pool fppool;
        pfem_precomp pfp = fppool(pf, &refpts);
        
        ref_mesh_dof_ind_ct dof = mf.ind_dof_of_element(cv);
        for (size_type qq=0; qq < qqdim; ++qq) {
          coeff[qq].resize(mf.nb_dof_of_element(cv));
          typename std::vector<T>::iterator cit = coeff[qq].begin();
          for (ref_mesh_dof_ind_ct::iterator it=dof.begin(); it != dof.end(); ++it, ++cit)
            *cit = U[(*it)*qqdim+qq];
        }
	fem_interpolation_context ctx(mf.linked_mesh().trans_of_convex(cv),pfp,0,G,cv);
        for (size_type j=0; j < refpts.size(); ++j) {
	  ctx.set_ii(j);
          for (size_type qq = 0; qq < qqdim; ++qq) {
            typename gmm::sub_vector_type<V2*, gmm::sub_interval>::vector_type dest = 
              gmm::sub_vector(V,gmm::sub_interval(pos,qdim));
            pf->interpolation(ctx,coeff[qq],dest,qdim);
            pos += qdim;
          }
        }
      }
      if (pos != V.size()) DAL_THROW(dal::failure_error, "bad dimensions");
    }
  };

  /**
     slicer whose side effect is to build a stored_mesh_slice object
   */
  class slicer_build_stored_mesh_slice : public slicer_action {
    stored_mesh_slice &sl;
  public:
    slicer_build_stored_mesh_slice(stored_mesh_slice& sl_) : sl(sl_) {
      if (sl.cvlst.size()) 
	DAL_THROW(dal::failure_error, "the stored_mesh_slice already contains data");
    }
    void exec(mesh_slicer& ms);
  };

  std::ostream& operator<<(std::ostream& o, const stored_mesh_slice& m);
}

#endif /*GETFEM_MESH_SLICE_H*/
