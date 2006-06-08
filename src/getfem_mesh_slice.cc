// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2003-2006 Julien Pommier
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

  void stored_mesh_slice::write_to_file(const std::string &name, 
			       bool with_mesh) const {
    std::ofstream o(name.c_str());
    if (!o)
      DAL_THROW(failure_error, "impossible to open file '" << name << "'");
    o << "% GETFEM SLICE FILE " << '\n';
    o << "% GETFEM VERSION " << GETFEM_VERSION << '\n' << '\n' << '\n';
    if (with_mesh) linked_mesh().write_to_file(o);
    write_to_file(o);
  }


  void stored_mesh_slice::write_to_file(std::ostream &os) const {
    os << "\nBEGIN MESH_SLICE\n";
    os << " DIM " << int(dim()) << "\n"; 
    for (unsigned i=0; i < cvlst.size(); ++i) {
      const convex_slice &cs = cvlst[i];
      os << " CONVEX " << cs.cv_num 
	 << " " << int(cs.fcnt)
	 << " " << int(cs.discont) << "\n"
	 << " " << cs.nodes.size() << " " << cs.simplexes.size() << "\n";
      for (unsigned j=0; j < cs.nodes.size(); ++j) {
	os << "\t";
	for (unsigned k=0; k < cs.nodes[j].pt.size(); ++k) {
	  if (k) os << " ";
	  os << cs.nodes[j].pt[k];
	}
	os << ";";
	for (unsigned k=0; k < cs.nodes[j].pt_ref.size(); ++k)
	  os << " " << cs.nodes[j].pt_ref[k];
	os << "; "; os << cs.nodes[j].faces.to_ulong();;
	os << "\n";
      }
      for (unsigned j=0; j < cs.simplexes.size(); ++j) {
	os << "\t" << cs.simplexes[j].inodes.size() << ":";
	for (unsigned k=0; k < cs.simplexes[j].inodes.size(); ++k) {
	  os << " " << cs.simplexes[j].inodes[k];
	}
	os << "\n";
      }
    }
    os << "END MESH_SLICE\n";
  }

  void stored_mesh_slice::read_from_file(const std::string &name, const getfem::mesh &m) {
    std::ifstream o(name.c_str());
    if (!o) DAL_THROW(file_not_found_error,
		      "slice file '" << name << "' does not exist");
    read_from_file(o,m);
  }

  void stored_mesh_slice::read_from_file(std::istream &ist, const getfem::mesh &m) {
    if (!poriginal_mesh) {
      poriginal_mesh = &m; 
    } else if (poriginal_mesh != &m) 
      DAL_THROW(dal::failure_error, "wrong mesh..");

    dim_ = m.dim();
    cv2pos.clear(); 
    cv2pos.resize(m.convex_index().last_true() + 1, 
		  size_type(-1));

    std::string tmp;
    ist.precision(16);
    ist.seekg(0);ist.clear();
    ftool::read_until(ist, "BEGIN MESH_SLICE");

    mesh_slicer::cs_nodes_ct nodes;
    mesh_slicer::cs_simplexes_ct simplexes;
    

    while (true) {
      ist >> std::ws; ftool::get_token(ist, tmp);
      if (ftool::casecmp(tmp, "END")==0) {
	break;
      } else if (ftool::casecmp(tmp, "DIM")==0) {
	int d; ist >> d;
	dim_ = d;
      } else if (ftool::casecmp(tmp, "CONVEX")==0) {
	ftool::get_token(ist,tmp);
	size_type ic = atoi(tmp.c_str());
	if (!m.convex_index().is_in(ic)) {
	  DAL_THROW(failure_error, "Convex " << ic <<
		    " does not exist, are you sure "
		    "that the mesh attached to this object is right one ?");
	}
	bgeot::pconvex_ref cvr = m.trans_of_convex(ic)->convex_ref();
	unsigned fcnt, discont, nbn, nbs;
	ist >> fcnt >> discont >> nbn >> nbs;
	nodes.resize(nbn); 
	simplexes.resize(nbs);
	for (unsigned i=0; i < nbn; ++i) {
	  nodes[i].pt.resize(dim()); 
	  nodes[i].pt_ref.resize(cvr->structure()->dim());
	  for (unsigned j=0; j < dim(); ++j) 
	    ist >> nodes[i].pt[j];
	  ist >> ftool::skip(";");
	  for (unsigned j=0; j < cvr->structure()->dim(); ++j) 
	    ist >> nodes[i].pt_ref[j];
	  ist >> ftool::skip(";");
	  unsigned long ul; ist >> ul;
	  nodes[i].faces = slice_node::faces_ct(ul);
	}
	for (unsigned i=0; i < nbs; ++i) {
	  unsigned np(0);
	  ist >> np >> ftool::skip(":");
	  if (np > dim()+1) 
	    DAL_THROW(dal::failure_error, "invalid simplex..");
	  simplexes[i].inodes.resize(np);
	  for (unsigned j=0; j < np; ++j) 
	    ist >> simplexes[i].inodes[j];
	}
	dal::bit_vector bv; bv.add(0, nbs);
	set_convex(ic, cvr, nodes, simplexes, dim_type(fcnt), bv, discont);
      } else if (tmp.size()) {
	DAL_THROW(failure_error, "Unexpected token '" << tmp <<
		  "' [pos=" << std::streamoff(ist.tellg()) << "]");
      } else if (ist.eof()) {
	DAL_THROW(failure_error, "Unexpected end of stream "
		  << "(missing BEGIN MESH_SLICE/END MESH_SLICE ?)");	
      }
    }
  }

  void slicer_build_stored_mesh_slice::exec(mesh_slicer &ms) {
    if (!sl.poriginal_mesh) {
      sl.poriginal_mesh = &ms.m; sl.dim_ = sl.linked_mesh().dim();
      sl.cv2pos.clear(); sl.cv2pos.resize(sl.linked_mesh().convex_index().last_true() + 1, size_type(-1));
    } else if (sl.poriginal_mesh != &ms.m) DAL_THROW(dal::failure_error, "wrong mesh..");
    sl.set_convex(ms.cv, ms.cvr, ms.nodes, ms.simplexes, ms.fcnt, ms.splx_in, ms.discont);
  }

  void stored_mesh_slice::set_convex(size_type cv, bgeot::pconvex_ref cvr, 
				     mesh_slicer::cs_nodes_ct cv_nodes, 
				     mesh_slicer::cs_simplexes_ct cv_simplexes, 
				     dim_type fcnt, const dal::bit_vector& splx_in, bool discont) {
    /* push the used nodes and simplexes in the final list */
    if (splx_in.card() == 0) return;
    merged_nodes_available = false;
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
      sc->global_points_count = points_cnt;
      sc->discont = discont;
    } else {
      sc = &cvlst[cv2pos[cv]];
      assert(sc->cv_num == cv);
    }
    for (dal::bv_visitor snum(splx_in); !snum.finished(); ++snum) {
      slice_simplex& s = cv_simplexes[snum];
      for (size_type i=0; i < s.dim()+1; ++i) {
        size_type lnum = s.inodes[i];
        if (nused[lnum] == size_type(-1)) {
          nused[lnum] = sc->nodes.size(); sc->nodes.push_back(cv_nodes[lnum]);
	  dim_ = std::max(int(dim_), int(cv_nodes[lnum].pt.size())); // dim_ may be equal to size_type(-1)
          points_cnt++;
        }
        s.inodes[i] = nused[lnum];
      }
      simplex_cnt.resize(dim_+1, 0);
      simplex_cnt[cv_simplexes[snum].dim()]++;
      sc->simplexes.push_back(cv_simplexes[snum]);
    }
  }

  struct get_edges_aux {
    size_type iA, iB;
    mutable bool slice_edge;
    get_edges_aux(size_type a, size_type b, bool slice_edge_) :
      iA(std::min(a,b)), iB(std::max(a,b)), slice_edge(slice_edge_) {}
    bool operator<(const get_edges_aux& other) const {
      /* ignore the slice_edge on purpose */
      return (iA < other.iA || (iA == other.iA && iB < other.iB));
    }
  };

  void stored_mesh_slice::get_edges(std::vector<size_type> &edges,
				    dal::bit_vector &slice_edges,
				    bool from_merged_nodes) const {
    if (from_merged_nodes && !merged_nodes_available) merge_nodes();
    std::set<get_edges_aux> e;
    for (cvlst_ct::const_iterator it = cvlst.begin(); it != cvlst.end(); ++it) {
      for (size_type is=0; is < it->simplexes.size(); ++is) {
	const slice_simplex &s = it->simplexes[is];
	for (size_type i=0; i < s.dim(); ++i) {
	  for (size_type j=i+1; j <= s.dim(); ++j) {
	    const slice_node& A = it->nodes[s.inodes[i]];
	    const slice_node& B = it->nodes[s.inodes[j]];
	    /* duplicate with slicer_build_edges_mesh which also 
	       builds a list of edges */
	    if ((A.faces & B.faces).count() >= unsigned(it->cv_dim-1)) {
	      slice_node::faces_ct fmask((1 << it->cv_nbfaces)-1); fmask.flip();
	      size_type iA, iB;
	      iA = it->global_points_count + s.inodes[i];
	      iB = it->global_points_count + s.inodes[j];
	      if (from_merged_nodes) {
		iA = to_merged_index[iA]; iB = to_merged_index[iB];
	      }
	      get_edges_aux a(iA,iB,((A.faces & B.faces) & fmask).any());
	      std::set<get_edges_aux>::iterator p=e.find(a);
	      if (p != e.end()) {
		if (p->slice_edge && !a.slice_edge) p->slice_edge = false;
	      } else e.insert(a);
	    }
	  }
	}
      }
    }
    slice_edges.clear(); slice_edges.sup(0, e.size());
    edges.clear(); edges.reserve(2*e.size());
    for (std::set<get_edges_aux>::const_iterator p=e.begin(); p != e.end(); ++p) {
      if (p->slice_edge) slice_edges.add(edges.size()/2);
      edges.push_back(p->iA);edges.push_back(p->iB);
    }
  }

  void stored_mesh_slice::build(const getfem::mesh& m, 
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
    clear_merged_nodes();
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
    size_type count = 0;
    for (size_type ic=0; ic < nb_convex(); ++ic) {
      cvlst[ic].global_points_count = count; count += nodes(ic).size();
    }
    assert(count == points_cnt);
  }

  void stored_mesh_slice::clear_merged_nodes() const { 
    merged_nodes_idx.clear(); merged_nodes.clear(); 
    to_merged_index.clear();
    merged_nodes_available = false; 
  }

  void stored_mesh_slice::merge_nodes() const {
    size_type count = 0;
    mesh mp;
    clear_merged_nodes();
    std::vector<size_type> iv;
    std::vector<const slice_node*> nv(nb_points());
    to_merged_index.resize(nb_points());
    for (cvlst_ct::const_iterator it = cvlst.begin(); it != cvlst.end(); ++it) {
      for (size_type i=0; i < it->nodes.size(); ++i) {
	nv[count] = &it->nodes[i];
	to_merged_index[count++] = mp.add_point(it->nodes[i].pt);
	//cout << "orig[" << count-1 << "] = " << nv[count-1]->pt << ", idx=" << to_merged_index[count-1] << "\n";
      }
    }
    dal::sorted_indexes(to_merged_index,iv);
    //cout << "to_merged_index = " << iv << "\n";
    merged_nodes.resize(nb_points());
    merged_nodes_idx.reserve(nb_points()/8);
    
    merged_nodes_idx.push_back(0);
    for (size_type i=0; i < nb_points(); ++i) {
      merged_nodes[i].P = nv[iv[i]];
      merged_nodes[i].pos = iv[i];
      //cout << "i=" << i << " -> {" << merged_nodes[i].P->pt << "," << merged_nodes[i].pos << "}\n";
      if (i == nb_points()-1 || to_merged_index[iv[i+1]] != to_merged_index[iv[i]]) 
	merged_nodes_idx.push_back(i+1);
    }
    //cout << "merged_nodes_idx = " << merged_nodes_idx << "\n";
    merged_nodes_available = true;
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
