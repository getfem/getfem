#include <getfem_mesh_fem.h>
#include <getfem_mesh_slice.h>
#include <bgeot_geotrans_inv.h>

namespace getfem {

  std::ostream& operator<<(std::ostream& o, const stored_mesh_slice& m) {
    o << "stored_mesh_slice, containing " << m.nb_convex() << " convexes\n";
    for (size_type ic = 0; ic < m.nb_convex(); ++ic) {
      o << "slice convex #" << ic << " (original = " << m.convex_num(ic) << ")\n";
      for (size_type i = 0; i < m.nodes(ic).size(); ++i) {
        o << "node " << i << ": " << m.nodes(ic)[i].pt << ", ref=" << m.nodes(ic)[i].pt_ref << " flist=" << m.nodes(ic)[i].faces << endl;
      }
      for (size_type i = 0; i < m.simplexes(ic).size(); ++i) {
        o << "simplex " << i << ", inodes=";
        for (size_type j=0;j< m.simplexes(ic)[i].dim()+1;++j)
          o << m.simplexes(ic)[i].inodes[j] << " ";
        o << endl;
      }
      o << endl;
    }
    return o;
  }

  void slicer_build_stored_mesh_slice::exec(mesh_slicer &ms) {
    if (!sl.poriginal_mesh) {
      sl.poriginal_mesh = &ms.m; sl.dim_ = sl.linked_mesh().dim();
      sl.cv2pos.clear(); sl.cv2pos.resize(sl.linked_mesh().convex_index().last_true() + 1, size_type(-1));
    } else if (sl.poriginal_mesh != &ms.m) DAL_THROW(dal::failure_error, "wrong mesh..");
    sl.set_convex(ms.cv, ms.cvr, ms.nodes, ms.simplexes, ms.fcnt, ms.splx_in);
  }

  void stored_mesh_slice::set_convex(size_type cv, bgeot::pconvex_ref cvr, 
				     mesh_slicer::cs_nodes_ct cv_nodes, 
				     mesh_slicer::cs_simplexes_ct cv_simplexes, 
				     dim_type fcnt, dal::bit_vector& splx_in) {
    /* push the used nodes and simplexes in the final list */
    if (splx_in.card() == 0) return;
    std::vector<size_type> nused(cv_nodes.size(), size_type(-1));
    convex_slice *sc = 0;
    if (cv >= cv2pos.size()) DAL_THROW(dal::internal_error, "");
    if (cv2pos[cv] == size_type(-1)) {
      cv2pos[cv] = cvlst.size();
      cvlst.push_back(convex_slice());
      sc = &cvlst.back();
      sc->cv_num = cv;
      sc->cv_dim = cvr->structure()->dim();
      sc->cv_nbfaces = cvr->structure()->nb_faces();
      sc->fcnt = fcnt;
    } else {
      sc = &cvlst[cv2pos[cv]];
      assert(sc->cv_num == cv);
    }
    for (size_type snum = splx_in.take_first(); snum != size_type(-1); snum << splx_in) {
      for (size_type i=0; i < cv_simplexes[snum].dim()+1; ++i) {
        size_type lnum = cv_simplexes[snum].inodes[i];
        if (nused[lnum] == size_type(-1)) {
          nused[lnum] = sc->nodes.size(); sc->nodes.push_back(cv_nodes[lnum]);
	  dim_ = std::max(int(dim_), int(cv_nodes[lnum].pt.size())); // dim_ may be equal to size_type(-1)
          points_cnt++;
        }
        cv_simplexes[snum].inodes[i] = nused[lnum];
      }
      simplex_cnt.resize(dim_+1, 0);
      simplex_cnt[cv_simplexes[snum].dim()]++;
      sc->simplexes.push_back(cv_simplexes[snum]);
    }
  }

  void stored_mesh_slice::build(const getfem::getfem_mesh& m, 
				const slicer_action *a, const slicer_action *b, const slicer_action *c, 
				size_type nrefine) {
    clear();
    mesh_slicer slicer(m);
    slicer.push_back_action(*const_cast<slicer_action*>(a));
    if (b) slicer.push_back_action(*const_cast<slicer_action*>(b));
    if (c) slicer.push_back_action(*const_cast<slicer_action*>(c));
    slicer_build_stored_mesh_slice sbuild(*this);
    slicer.push_back_action(sbuild);
    slicer.exec(nrefine);
  }

  void stored_mesh_slice::replay(slicer_action *a, slicer_action *b, slicer_action *c) const {
    mesh_slicer slicer(linked_mesh());
    slicer.push_back_action(*a); 
    if (b) slicer.push_back_action(*b);
    if (c) slicer.push_back_action(*c);
    slicer.exec(*this);
  }

  void stored_mesh_slice::set_dim(size_type newdim) {
    dim_ = newdim;
    for (size_type ic=0; ic < nb_convex(); ++ic) {
      for (mesh_slicer::cs_nodes_ct::iterator it=nodes(ic).begin(); it != nodes(ic).end(); ++it) {
	it->pt.resize(newdim);
      }
    }
  }

  void stored_mesh_slice::merge(const stored_mesh_slice& sl) {
    if (dim() != sl.dim()) DAL_THROW(dal::dimension_error, "inconsistent dimensions for slice merging");
    cv2pos.resize(std::max(cv2pos.size(), sl.cv2pos.size()), size_type(-1));
    for (size_type i=0; i < sl.nb_convex(); ++i) 
      if (cv2pos[sl.convex_num(i)] != size_type(-1) &&
	  cvlst[cv2pos[sl.convex_num(i)]].cv_dim != sl.cvlst[i].cv_num)
	DAL_THROW(dal::dimension_error, "inconsistent dimensions for convex " << sl.cvlst[i].cv_num << " on the slices");

    for (size_type i=0; i < sl.nb_convex(); ++i) {
      size_type cv = sl.convex_num(i);
      if (cv2pos[cv] == size_type(-1)) {
	cv2pos[cv] = cvlst.size(); cvlst.push_back(convex_slice());
      }
      const stored_mesh_slice::convex_slice *src = &sl.cvlst[i];
      stored_mesh_slice::convex_slice *dst = &cvlst[cv2pos[cv]];
      size_type n = dst->nodes.size();
      dst->nodes.insert(dst->nodes.end(), src->nodes.begin(), src->nodes.end());
      for (mesh_slicer::cs_simplexes_ct::const_iterator it = src->simplexes.begin(); it != src->simplexes.end(); ++it) {
	dst->simplexes.push_back(*it);
	for (size_type j = 0; j < (*it).dim()+1; ++j) dst->simplexes.back().inodes[j] += n;	
	simplex_cnt[dst->simplexes.back().dim()]++;
      }
      points_cnt += src->nodes.size();
    }
  }

  size_type stored_mesh_slice::memsize() const {
    size_type sz = sizeof(stored_mesh_slice);
    for (cvlst_ct::const_iterator it = cvlst.begin(); it != cvlst.end(); ++it) {
      sz += sizeof(size_type);
      /*cerr << "memsize: convex " << it->cv_num << " nodes:" 
	<< it->nodes.size() << ", splxs:" << it->simplexes.size() << ", sz=" << sz << "\n";*/
      for (size_type i=0; i < it->nodes.size(); ++i) {
	/*cerr << "  point " << i << ": size+= " << sizeof(slice_node) << "+" <<  
          it->nodes[i].pt.memsize() << "+" << it->nodes[i].pt_ref.memsize() << "-" << sizeof(it->nodes[i].pt)*2 << "\n";*/
        sz += sizeof(slice_node) + 
          (it->nodes[i].pt.memsize()+it->nodes[i].pt_ref.memsize()) - sizeof(it->nodes[i].pt)*2;
      }
      for (size_type i=0; i < it->simplexes.size(); ++i) {
	/*cerr << "  simplex " << i << ": size+= " << sizeof(slice_simplex) << "+" << 
          it->simplexes[i].inodes.size()*sizeof(size_type) << "\n";*/
        sz += sizeof(slice_simplex) + 
          it->simplexes[i].inodes.size()*sizeof(size_type);
      }
    }
    sz += cv2pos.size() * sizeof(size_type);
    return sz;
  }

}
