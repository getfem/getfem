#include <map>
#include <getfem_mesh_fem.h>
#include <getfem_mesh_slice.h>

namespace getfem {
  const float slicer::EPS = 1e-13;


  /* ---------------------------- extraction of  outer faces of a mesh --------------------- */

  /* identify convex faces by their mesh point ids */
  struct mesh_faces_by_pts_list_elt  {
    std::vector<size_type> ptid; // point numbers of faces
    int cnt; // number of convexes sharing that face
    int cv, f;
    bool operator<(const mesh_faces_by_pts_list_elt &e) const {
      return ptid < e.ptid;
    }
    template<typename CONT> 
    mesh_faces_by_pts_list_elt(size_type _cv, size_type _f, const CONT& p) 
      : ptid(p.size()), cnt(0), cv(_cv), f(_f)  {
      if (p.size() == 0) DAL_THROW(dal::internal_error, "internal error");
      std::partial_sort_copy(p.begin(), p.end(), ptid.begin(), ptid.end());
    }
    mesh_faces_by_pts_list_elt() {}
  };
  typedef dal::dynamic_tree_sorted<mesh_faces_by_pts_list_elt> mesh_faces_by_pts_list;

  /**
     returns a list of "exterior" faces of a mesh (i.e. faces which are not shared by two convexes)
  */
  void
  outer_faces_of_mesh(const getfem::getfem_mesh &m, dal::bit_vector cvlst, convex_face_ct& flist) {
    mesh_faces_by_pts_list lst;
    dal::bit_vector convex_tested;
  
    for (dal::bit_vector::const_iterator it_cv = cvlst.begin(); it_cv != cvlst.end(); ++it_cv) {
      if (!(*it_cv)) continue;

      size_type ic = it_cv.index();
      if (m.structure_of_convex(ic)->dim() == m.dim()) {
	for (size_type f = 0; f < m.structure_of_convex(ic)->nb_faces(); f++) {
	  size_type idx = lst.add_norepeat(mesh_faces_by_pts_list_elt(ic,f,m.ind_points_of_face_of_convex(ic, f)));
	  lst[idx].cnt++;
	}
      } else { /* les objets de dim inferieure sont considérés comme "exterieurs" 
		(c'ets plus pratique pour faire des dessins)
	       */
	size_type idx = lst.add_norepeat(mesh_faces_by_pts_list_elt(ic,size_type(-1),m.ind_points_of_convex(ic)));
	lst[idx].cnt++;
      }
    }

    size_type fcnt = 0;
    for (size_type i = 0; i < lst.size(); i++)
      if (lst[i].cnt == 1) ++fcnt;
    flist.resize(fcnt); fcnt = 0;
    for (size_type i = 0; i < lst.size(); i++) if (lst[i].cnt == 1) { 
      flist[fcnt].cv = lst[i].cv; flist[fcnt].f = lst[i].f; ++fcnt; 
    }
  }



  /* -------------------------------------- slicers --------------------------------------*/

  /* nodes : list of nodes (new nodes may be added)
     splxs : list of simplexes (new simplexes may be added)
     splx_in : input: simplexes to take into account, output: list of simplexes inside the slice

     note that the simplexes in the list may have different dimensions
  */
  void slicer_volume::slice(std::deque<slice_node>& nodes, 
                            std::deque<slice_simplex>& splxs, dal::bit_vector& splx_in) const {    
    /*size_type cnt=0;
    for (std::deque<slice_simplex>::iterator it = splxs.begin();
         it != splxs.end(); ++it, ++cnt) {    
      if (!splx_in[cnt]) continue;
    */
    //cerr << "\n\n------------------------------------------------\nslicer_volume::slice : entree, splx_in=" << splx_in << endl;
    if (splx_in.card() == 0) return;
    dal::bit_vector pt_in; pt_in.sup(0,nodes.size());
    for (size_type i=0; i < nodes.size(); ++i) if (is_in(nodes[i].pt)) pt_in.add(i);
    dal::bit_vector pt_bin; pt_bin.sup(0,nodes.size());
    for (size_type i=0; i < nodes.size(); ++i) if (is_in(nodes[i].pt,true)) pt_bin.add(i);

    dal::bit_vector bv = splx_in;
    for (size_type cnt=bv.take_first(); cnt != size_type(-1); cnt << bv) {
      slice_simplex& s = splxs[cnt];
      /*cerr << "\n--------slicer_volume::slice : slicing convex " << cnt << endl;
      for (size_type i=0; i < s.dim()+1; ++i)
        cerr << "   * pt[" << i << "]=" << nodes[s.inodes[i]].pt << ", is_in=" << 
          is_in(nodes[s.inodes[i]].pt) << ", is_bin=" << is_in(nodes[s.inodes[i]].pt,true) << endl;
      */
      size_type in_cnt = 0, in_bcnt = 0;
      for (size_type i=0; i < s.dim()+1; ++i) {
	if (pt_in[s.inodes[i]]) ++in_cnt;
        if (pt_bin[s.inodes[i]]) ++in_bcnt;
      }

      if (in_cnt == 0) {
        splx_in.sup(cnt);
      } else if (in_cnt != s.dim()+1 || is_boundary_slice) {           /* the simplex crosses the slice boundary */
        splx_in.sup(cnt);
        size_type l = splxs.size();//, n = nodes.size();
        //cerr << "slicer_volume::slice : convex " << cnt << " will be splited" << endl;
	split_simplex(nodes, pt_in, pt_bin, splxs, s, splxs.size(), is_boundary_slice);
        splx_in.add(l,splxs.size()-l);
      }
    }
  }

  void slicer_union::slice(std::deque<slice_node>& nodes, std::deque<slice_simplex>& splxs, dal::bit_vector& splx_in) const {
    dal::bit_vector idx = splx_in;
    for (size_type cnt = idx.take_first(); cnt != size_type(-1); cnt << idx) {
      size_type in_cntA = A->is_in(nodes, splxs[cnt]);

      if (in_cntA == splxs[cnt].dim()+1) { /* the simplex is inside A */
        continue;
      }
      size_type in_cntB = B->is_in(nodes, splxs[cnt]);
      if (in_cntB == splxs[cnt].dim()+1) { /* the simplex is inside B */
        continue;
      }      
      dal::bit_vector bv; bv.add(cnt);
      if (in_cntA == 0) {
        B->slice(nodes, splxs, bv);
      } else if (in_cntB) {
        A->slice(nodes, splxs, bv);
      } else splx_in.sup(cnt);
      splx_in |= bv;
    }
    
  }

  void slicer_intersect::slice(std::deque<slice_node>& nodes, std::deque<slice_simplex>& splxs, dal::bit_vector& splx_in) const {
    dal::bit_vector idx = splx_in;
    for (size_type cnt = idx.take_first(); cnt != size_type(-1); cnt << idx) {
      size_type in_cntA = A->is_in(nodes, splxs[cnt]);
      size_type in_cntB = B->is_in(nodes, splxs[cnt]);
      if (in_cntA == splxs[cnt].dim()+1 && in_cntB == splxs[cnt].dim()+1) {
        continue;
      } else {
        dal::bit_vector bv; bv.add(cnt);
        splx_in.sup(cnt);
        A->slice(nodes, splxs, bv);
        B->slice(nodes, splxs, bv);
        splx_in |= bv;
      }
    }
  }



  struct sorted_order_aux {
    const std::vector<size_type>& w;
    bool operator()(size_type i, size_type j) { return w[i] > w[j]; }
    sorted_order_aux(const std::vector<size_type>& W) : w(W) {}
  };
  static void sorted_order(const std::vector<size_type>& w, std::vector<size_type>&order) {
    order.resize(w.size());
    for (size_type i=0; i < w.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), sorted_order_aux(w));
  }

  /* intersects the simplex with the slice, and (recusively) decomposes it
     into sub-simplices, which are added to the list 'splxs' 
     if 'reduce_dimension' is true, then it is the faces of 
     sub-simplices which are added
  */
  void slicer::split_simplex(std::deque<slice_node>& nodes, dal::bit_vector& pt_in, dal::bit_vector& pt_bin,
                             std::deque<slice_simplex>& splxs, 
                             const slice_simplex& s, size_type sstart, bool reduce_dimension) const {
    scalar_type alpha = 0; size_type iA=0, iB = 0;
    bool intersection = false;
    static int level = 0;

    level++;    
    /*cerr << "split_simplex: level=" << level << " simplex: ";
    for (iA=0; iA < s.dim()+1 && !intersection; ++iA) 
      cerr << "node#" << s.inodes[iA] << "=" << nodes[s.inodes[iA]].pt 
           << ", in=" << pt_in[s.inodes[iA]] << ", bin=" << pt_bin[s.inodes[iA]] << "; "; cerr << endl;
    */
    assert(level < 100);
    for (iA=0; iA < s.dim()+1; ++iA) {
      for (iB=iA+1; iB < s.dim()+1; ++iB) {
        if (pt_in[s.inodes[iA]] != pt_in[s.inodes[iB]] && !pt_bin[s.inodes[iA]] && !pt_bin[s.inodes[iB]]) {
          alpha=edge_intersect(nodes[s.inodes[iA]].pt,nodes[s.inodes[iB]].pt);
          //cerr << " : intersection #" << iA << ":"<< nodes[s.inodes[iA]].pt << "-#"<<iB<<":"<<nodes[s.inodes[iB]].pt<<", f=" << alpha << endl;
          if (alpha >= 1e-12 && alpha <= 1-1e-12) { intersection = true; break; }
        }
      }
      if (intersection) break;
    }
    if (intersection) {
      const slice_node& A = nodes[s.inodes[iA]]; 
      const slice_node& B = nodes[s.inodes[iB]]; 
      slice_node n; 
      n.pt = A.pt + alpha*(B.pt-A.pt);
      n.pt_ref = A.pt_ref + alpha*(B.pt_ref-A.pt_ref);
      n.faces = A.faces & B.faces;
      size_type nn = nodes.size(); 
      nodes.push_back(n); 
      pt_bin.add(nn); pt_in.add(nn);
      //cerr << " -> intersection at alpha=" << alpha << ", EPS=" << EPS << ", n=" << n.pt << endl;

      slice_simplex s1(s.dim()+1); 
      for (size_type k=0; k < s.dim()+1; k++)
        s1.inodes[k] = (k != iA) ? s.inodes[k] : nn;
      split_simplex(nodes,pt_in,pt_bin,splxs,s1,sstart, reduce_dimension);
      for (size_type k=0; k < s.dim()+1; k++)
        s1.inodes[k] = (k != iB) ? s.inodes[k] : nn;
      split_simplex(nodes,pt_in,pt_bin,splxs,s1,sstart, reduce_dimension);
    } else {
      if (!reduce_dimension) {
        bool all_in = true;
        for (size_type i=0; i < s.dim()+1; ++i) if (!pt_in[s.inodes[i]]) all_in = false;
        //cerr << " -> no intersection , all_in=" << all_in << endl;
        if (all_in)
          splxs.push_back(s);
      } else {
        slice_simplex face(s.dim());
        for (size_type f=0; f < s.dim()+1; ++f) {
          bool all_in = true;
          for (size_type i=0; i < s.dim(); ++i) {
            size_type p = s.inodes[i + (i<f?0:1)];
            if (!pt_bin[p]) { all_in = false; break; }
            else face.inodes[i] = s.inodes[i + (i<f?0:1)];
          }
          if (all_in) {
            /* prevent addition of a face twice */
            std::sort(face.inodes.begin(), face.inodes.end());
            if (std::find(splxs.begin()+sstart, splxs.end(), face) == splxs.end()) {
              splxs.push_back(face);
            }
          }
        }
      }
    }
    level--;
  }

  struct slice_node_compare_pt_ref : public std::binary_function<slice_node,slice_node,bool> {
    bool operator()(const slice_node& a, const slice_node& b) const { return a.pt_ref < b.pt_ref; }
  };

  /** 
      build the list of edges and store them in a mesh object (hence all common nodes/edges are eliminated)
  */
  void mesh_slice::edges_mesh(getfem::getfem_mesh &m) const {
    cs_nodes_ct n;
    std::vector<size_type> fpts;
    std::vector<dim_type> bits;
    std::vector<dim_type> cnt;
    for (size_type ic = 0; ic < cvlst.size(); ++ic) {
      //cerr << "examen du convexe " << ic << " (" << nodes(ic).size() << " nodes, " << simplexes(ic).size() << " simplexes" << endl;
      n.resize(nodes(ic).size());
      /* the trick is here: the nodes are sorted according to their
	 reference point. Since the edges of the reference convex are
	 assumed to be straight, the edges points will be also
	 sorted */
      std::partial_sort_copy(cvlst[ic].nodes.begin(), cvlst[ic].nodes.end(), 
			     n.begin(), n.end(), slice_node_compare_pt_ref());
      typedef std::map<unsigned long,dal::bit_vector> fmap_t;
      fmap_t fmap;
      dim_type nface = 0;//cvlst[ic].cv_dim-1;
      for (size_type j=0; j < cvlst[ic].simplexes.size(); ++j) 
        if (cvlst[ic].simplexes[j].dim()) 
          nface = std::max(nface, dim_type(cvlst[ic].simplexes[j].dim()-1));
      fpts.resize(n.size());
      slice_node::faces_ct fmask((unsigned long)(-1));
      for (size_type ip = 0; ip < n.size(); ++ip) fmask &= n[ip].faces;
      for (size_type ip = 0; ip < n.size(); ++ip) {
	slice_node::faces_ct f = n[ip].faces & (~fmask);
        dim_type nbits = f.count();
	//cerr << "  nface = " << int(nface) << ", nbits=" << int(nbits) << endl;
	//cerr << "  noeud " << ip << ": " << n[ip].pt << ", " << n[ip].pt_ref << ", f=" << n[ip].faces << endl;
        if (nface <= nbits) { /* not sure that is the right test for
				 dimension > 3 .. */
          /* add the point to the mesh */
          fpts[ip] = m.add_point(n[ip].pt);

          /* now generate all combinations of (cvlst[ic].cv_dim-1)
	     bits in the f.count() that are set in f .. */

          /* build list of faces */
          bits.clear(); 
	  for (dim_type bcnt = 0; f.any(); ++bcnt) { 
	    if (f[0]) bits.push_back(bcnt); 
	    f >>= 1; 
	  }
	  //cerr << "   bits="; std::copy(bits.begin(), bits.end(), std::ostream_iterator<size_type>(cerr," ")); cerr << endl;
	  
          /* init counter */
          cnt.resize(nface+1); 
	  for (size_type i=0; i < nface; ++i) cnt[i] = i; 
	  cnt[nface] = nbits+1;
	  //cerr << "    cnt="; std::copy(cnt.begin(), cnt.end(), std::ostream_iterator<size_type>(cerr," ")); cerr << endl;
	  
          while (cnt[nface-1] != nbits) {
            f.reset(); for (size_type i=0; i < nface; ++i) f[bits[cnt[i]]] = 1;
	    //cerr << "        ---> f <= " << f << endl;
            if (int(f.count()) == nface) {
	      //cerr << "    -> ajout du noeud sur la face " << f << endl;
              fmap[f.to_ulong()].add(ip);
            }
            /* next combination */
            for (dim_type i=0, pcnt=0; i < nface; ++i) {
              if (++cnt[i] < cnt[i+1]) break; else cnt[i] = pcnt;
              pcnt = cnt[i]+1;
            }
	    //cerr << "    cnt="; std::copy(cnt.begin(), cnt.end(), std::ostream_iterator<size_type>(cerr," ")); cerr << endl;
          }
        }
      }
      for (fmap_t::iterator it = fmap.begin(); it != fmap.end(); ++it) {
        dal::bit_vector &bv = (*it).second;
        size_type pip = bv.take_first();
        for (size_type ip = bv.take_first(); ip != size_type(-1); ip << bv) {
	  if (fpts[pip] != fpts[ip])
	    m.add_segment(fpts[pip], fpts[ip]); 
	  pip = ip;
	}
      }
    }
  }

  static void flag_points_on_faces(const bgeot::pconvex_ref& cvr, 
                                   const bgeot::stored_point_tab& pts, 
                                   std::vector<slice_node::faces_ct>& faces) {
    if (cvr->structure()->nb_faces() > 32) DAL_THROW(std::out_of_range, "won't work for convexes with more than 32 faces (hardcoded limit)");
    faces.resize(pts.size());
    for (size_type i=0; i < pts.size(); ++i) {
      faces[i].reset();      
      for (size_type f=0; f < cvr->structure()->nb_faces(); ++f)
        faces[i][f] = (dal::abs(cvr->is_in_face(f, pts[i])) < 1e-10);
    }
  }

  std::ostream& operator<<(std::ostream& o, const mesh_slice& m) {
    o << "mesh_slice, containing " << m.nb_convex() << " convexes\n";
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

  class mesh_slice_pre_deform {
    mesh_slice_cv_dof_data_base *defdata;
    pfem pf;
    bgeot::pstored_point_tab pspt;
    _fem_precomp fprecomp;
    base_matrix G;
    base_vector coeff;
    size_type cv;
  public:
    mesh_slice_pre_deform(mesh_slice_cv_dof_data_base *_defdata) 
      : defdata(_defdata), pf(0), pspt(0), cv(size_type(-1)) {
      if (defdata && defdata->pmf->get_qdim() != defdata->pmf->linked_mesh().dim()) 
        DAL_THROW(dal::dimension_error, "wrong Q(=" << int(defdata->pmf->get_qdim()) 
                  << ") dimension for slice deformation: should be equal to "
                  "the mesh dimension which is " << int(defdata->pmf->linked_mesh().dim()));
    }
    void prepare(size_type _cv, bgeot::stored_point_tab& refpts, bool force_update) {
      cv = _cv;
      if (defdata) {
        if (force_update || defdata->pmf->fem_of_element(cv) != pf) {
          pf = defdata->pmf->fem_of_element(cv);
          fem_precomp_not_stored(pf, &refpts, fprecomp);
        }
        defdata->copy(cv, coeff);
      }
    }

    /* apply the geometric transformation and an optional mesh_fem deformation */
    void apply(const getfem_mesh& m, const bgeot::mesh_structure *cvms, 
               _geotrans_precomp& gp, 
               const bgeot::stored_point_tab &cvm_pts, std::vector<slice_node::faces_ct> &points_on_faces,
               mesh_slice::cs_nodes_ct& cv_nodes, mesh_slice::cs_simplexes_ct& cv_simplexes) const {
      std::vector<size_type> ptsid(cvm_pts.size(), size_type(-1));
      
      cv_simplexes.resize(cvms->nb_convex());
      cv_nodes.resize(0);
      size_type pcnt = 0;
      for (size_type snum = 0; snum < cvms->nb_convex(); ++snum) { /* cvms should not contain holes in its convex index.. */
        cv_simplexes[snum].inodes.resize(cvms->nb_points_of_convex(snum));
        std::copy(cvms->ind_points_of_convex(snum).begin(),
                  cvms->ind_points_of_convex(snum).end(), cv_simplexes[snum].inodes.begin());
        /* store indices of points which are really used , and renumbers them */
        for (std::vector<size_type>::iterator itp = cv_simplexes[snum].inodes.begin();
             itp != cv_simplexes[snum].inodes.end(); ++itp) {
          if (ptsid[*itp] == size_type(-1)) {
            cv_nodes.push_back(slice_node());
            ptsid[*itp] = pcnt;
            cv_nodes[pcnt].pt_ref = cvm_pts[*itp];
            cv_nodes[pcnt].faces = points_on_faces[*itp];
            cv_nodes[pcnt].pt.resize(m.dim()); cv_nodes[pcnt].pt.fill(0.);
            if (defdata) {
              cv_nodes[pcnt].pt.resize(defdata->pmf->get_qdim());
              pf->interpolation(&fprecomp, *itp,
                                G, m.trans_of_convex(cv),
                                coeff, cv_nodes[pcnt].pt, defdata->pmf->get_qdim());
            }
            gp.transform(m.points_of_convex(cv), *itp, cv_nodes[pcnt].pt);
            pcnt++;
          }
          *itp = ptsid[*itp];
        }
      }
    }
  };

  void mesh_slice::do_slicing(size_type cv, bgeot::pconvex_ref cvr, const slicer *ms, cs_nodes_ct cv_nodes, 
			      cs_simplexes_ct cv_simplexes, dal::bit_vector& splx_in) {
    /* do the slices */
    if (ms) ms->slice(cv_nodes, cv_simplexes, splx_in);
    
    /* push the used nodes and simplexes in the final list */         
    if (splx_in.card()) {
      std::vector<size_type> nused(cv_nodes.size(), size_type(-1));      
      cvlst.push_back(convex_slice());
      cvlst.back().cv_num = cv;
      cvlst.back().cv_dim = cvr->structure()->dim();
      for (size_type snum = splx_in.take_first(); snum != size_type(-1); snum << splx_in) {
	for (size_type i=0; i < cv_simplexes[snum].dim()+1; ++i) {
	  size_type lnum = cv_simplexes[snum].inodes[i];
	  if (nused[lnum] == size_type(-1)) {
	    nused[lnum] = cvlst.back().nodes.size(); cvlst.back().nodes.push_back(cv_nodes[lnum]);
	    points_cnt++;
	  }
	  cv_simplexes[snum].inodes[i] = nused[lnum];
	}
	simplex_cnt[cv_simplexes[snum].dim()]++;
	cvlst.back().simplexes.push_back(cv_simplexes[snum]);
      }
    }
  }

  /* of course, nodes created from edge/slice intersection are almost always duplicated */
  mesh_slice::mesh_slice(const getfem_mesh& m, const slicer* ms, size_type nrefine, 
                         convex_face_ct& in_cvlst, mesh_slice_cv_dof_data_base *def_mf_data) 
    : simplex_cnt(m.dim()+1, size_type(0)), points_cnt(0), _dim(m.dim()) {
    _geotrans_precomp gp;
    bgeot::stored_point_tab cvm_pts;
    bgeot::pconvex_ref prev_cvr = 0;
    const getfem_mesh *cvm = 0;
    const bgeot::mesh_structure *cvms = 0;

    cs_nodes_ct cv_nodes;
    cs_simplexes_ct cv_simplexes;

    std::auto_ptr<mesh_slice_pre_deform> def( new mesh_slice_pre_deform(def_mf_data) );
    
    /* list of nodes and simplexes for the current convex 
       (including those who are outside of the slice)
     */
    std::vector<slice_node::faces_ct> points_on_faces;
    for (convex_face_ct::const_iterator it = in_cvlst.begin(); it != in_cvlst.end(); ++it) {
      size_type cv = (*it).cv;
      size_type face = (*it).f;
      bgeot::pconvex_ref cvr = m.trans_of_convex(cv)->convex_ref();

      //cerr << "mesh_slice::mesh_slice: doing convex " << cv << endl;
      /* update structure-dependent data */
      if (prev_cvr != cvr) {
	prev_cvr = cvr;
	cvm = getfem::refined_simplex_mesh_for_convex(cvr, nrefine);
	cvm_pts.resize(cvm->nb_points());
	std::copy(cvm->points().begin(), cvm->points().end(), cvm_pts.begin());
	geotrans_precomp_not_stored(m.trans_of_convex(cv), 
				    &cvm_pts, gp);
        flag_points_on_faces(cvr, cvm_pts, points_on_faces);
        def->prepare(cv, cvm_pts, true);
      } else def->prepare(cv, cvm_pts, false);

      if (face < dim_type(-1))
        cvms = getfem::refined_simplex_mesh_for_convex_faces(cvr, nrefine)[face];
      else
        cvms = cvm; 

      /*cerr << "doing convex " << cv << ", face=" << face << ", cvm->nb_convex=" << cvm->nb_convex() 
	<< ", cvms->nb_convex=" << cvms->nb_convex() << endl;*/


      def->apply(m, cvms, gp, cvm_pts, points_on_faces, cv_nodes, cv_simplexes);

      dal::bit_vector splx_in; splx_in.add(0, cv_simplexes.size());
      do_slicing(cv, cvr, ms, cv_nodes, cv_simplexes, splx_in);
    }
    //cerr << *this << endl;
  }

  void mesh_slice::merge(const mesh_slice& sl) {
    if (dim() != sl.dim()) DAL_THROW(dal::dimension_error, "inconsistent dimensions for slice merging");
    size_type maxcvid=0;
    for (size_type i=0; i < nb_convex(); ++i) maxcvid=std::max(maxcvid,convex_num(i));
    for (size_type i=0; i < sl.nb_convex(); ++i) maxcvid=std::max(maxcvid,sl.convex_num(i));
    std::vector<size_type> cvs(maxcvid,size_type(-1));
    for (size_type i=0; i < nb_convex(); ++i) cvs[convex_num(i)] = i;
    for (size_type i=0; i < sl.nb_convex(); ++i) 
      if (cvs[sl.convex_num(i)] != size_type(-1) &&
	  cvlst[cvs[sl.convex_num(i)]].cv_dim != sl.cvlst[i].cv_num)
	DAL_THROW(dal::dimension_error, "inconsistent dimensions for convex " << sl.cvlst[i].cv_num << " on the slices");

    for (size_type i=0; i < sl.nb_convex(); ++i) {
      size_type cv = sl.convex_num(i);
      if (cvs[cv] == size_type(-1)) {
	cvs[cv] = cvlst.size(); cvlst.push_back(convex_slice());
      }
      const mesh_slice::convex_slice *src = &sl.cvlst[i];
      mesh_slice::convex_slice *dst = &cvlst[cvs[cv]];
      size_type n = dst->nodes.size();
      dst->nodes.insert(dst->nodes.end(), src->nodes.begin(), src->nodes.end());
      for (cs_simplexes_ct::const_iterator it = src->simplexes.begin(); it != src->simplexes.end(); ++it) {
	dst->simplexes.push_back(*it);
	for (size_type j = 0; j < (*it).dim()+1; ++j) dst->simplexes.back().inodes[j] += n;	
	simplex_cnt[dst->simplexes.back().dim()]++;
      }
      points_cnt += src->nodes.size();
    }
  }

  size_type mesh_slice::memsize() const {
    size_type sz = sizeof(mesh_slice);
    for (cvlst_ct::const_iterator it = cvlst.begin(); it != cvlst.end(); ++it) {
      sz += sizeof(size_type);
      for (size_type i=0; i < it->nodes.size(); ++i) 
        sz += sizeof(slice_node) + 
          (it->nodes[i].pt.memsize()+it->nodes[i].pt_ref.memsize()) - sizeof(it->nodes[i].pt)*2;
      for (size_type i=0; i < it->simplexes.size(); ++i) 
        sz += sizeof(slice_simplex) + 
          it->simplexes[i].inodes.size()*sizeof(size_type);
    }
    return sz;
  }

}
