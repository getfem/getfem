#include <getfem_mesh_fem.h>
#include <getfem_mesh_slice.h>
#include <bgeot_geotrans_inv.h>

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
    mesh_faces_by_pts_list_elt(size_type cv_, size_type f_, const CONT& p) 
      : ptid(p.size()), cnt(0), cv(cv_), f(f_)  {
      if (p.size() == 0) DAL_THROW(dal::internal_error, "internal error");
      std::partial_sort_copy(p.begin(), p.end(), ptid.begin(), ptid.end());
    }
    mesh_faces_by_pts_list_elt() {}
  };
  typedef dal::dynamic_tree_sorted<mesh_faces_by_pts_list_elt> mesh_faces_by_pts_list;

  void
  outer_faces_of_mesh(const getfem::getfem_mesh &m, const dal::bit_vector& cvlst, convex_face_ct& flist) {
    mesh_faces_by_pts_list lst;
    dal::bit_vector convex_tested;
  
    for (dal::bv_visitor ic(cvlst); !ic.finished(); ++ic) {
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
  
  slicer_boundary::slicer_boundary(const getfem_mesh& m, slicer *sA, 
				   const convex_face_ct& cvflst) : A(sA) {
    if (m.convex_index().card()==0) return;
    convex_faces.resize(m.convex_index().last()+1, slice_node::faces_ct(0L));
    for (size_type i=0; i < cvflst.size(); ++i) 
      if (cvflst[i].is_face() && cvflst[i].f<32) convex_faces[cvflst[i].cv][cvflst[i].f]=1;
      else convex_faces[cvflst[i].cv].set();
    /* set the mask to 1 for all other possible faces of the convexes, which may 
       appear after slicing the convex, hence they will be part of the "boundary" */
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv) {
      for (size_type f=m.structure_of_convex(cv)->nb_faces(); f < convex_faces[cv].size(); ++f)
	convex_faces[cv][f]=1;
    }
  }

  bool slicer_boundary::test_bound(const slice_simplex& s, slice_node::faces_ct& fmask, const mesh_slice::cs_nodes_ct& nodes) const {
    slice_node::faces_ct f; f.set();
    for (size_type i=0; i < s.dim()+1; ++i) {
      f &= nodes[s.inodes[i]].faces;
    }
    f &= fmask;
    return (f.any());
  }

  void slicer_boundary::slice(size_type cv, dim_type& fcnt,
                              mesh_slice::cs_nodes_ct& nodes, 
                              mesh_slice::cs_simplexes_ct& splxs, 
                              dal::bit_vector& splx_in) {
    A->slice(cv, fcnt, nodes, splxs, splx_in);
    if (splx_in.card() == 0) return;
    slice_node::faces_ct fmask(cv < convex_faces.size() ? convex_faces[cv] : 0);
    //cerr << "slicer_boundary::slice(cv=" << cv << ", fmask=" << fmask << ")\n";
    /* quickly check if the convex have any chance to be part of the boundary */
    if (!convex_faces[cv].any()) { splx_in.clear(); return; }

    for (dal::bv_visitor_c cnt(splx_in); !cnt.finished(); ++cnt) {
      const slice_simplex& s = splxs[cnt]; 
      /*cerr << "slicer_boundary::slice, cnt= " << cnt << ", fmask=" << fmask << endl;
      for (size_type iA=0; iA < s.dim()+1; ++iA) 
	cerr << " node#" << s.inodes[iA] << "=" << nodes[s.inodes[iA]].pt << ", f=" << nodes[s.inodes[iA]].faces << endl;
      */
      if (s.dim() < nodes[0].pt.size()) {
	//cerr << " -> splx_in[cnt]=" << test_bound(s, fmask, nodes) << endl;
        if (!test_bound(s, fmask, nodes)) splx_in.sup(cnt);
      } else if (s.dim() == 2) {
	splx_in.sup(cnt);
        slice_simplex s2(2);
        for (size_type j=0; j < 3; ++j) {
          /* usage of s forbidden in this loop since push_back happens .. */
	  static unsigned ord[][2] = {{0,1},{1,2},{2,0}}; /* keep orientation of faces */
          for (size_type k=0; k < 2; ++k) { s2.inodes[k] = splxs[cnt].inodes[ord[j][k]]; }
          if (test_bound(s2, fmask, nodes)) {
            splx_in.add(splxs.size()); splxs.push_back(s2); 
          }
        }
      } else if (s.dim() == 3) {
        splx_in.sup(cnt);
        slice_simplex s2(3);
        for (size_type j=0; j < 4; ++j) {
          /* usage of s forbidden in this loop since push_back happens .. */
	  static unsigned ord[][3] = {{0,2,1},{1,2,3},{1,3,0},{0,3,2}}; /* keep orientation of faces */
          for (size_type k=0; k < 3; ++k) { s2.inodes[k] = splxs[cnt].inodes[ord[j][k]]; } //k + (k<j ? 0 : 1)]; }
	  /*cerr << " -> testing "; for (size_type iA=0; iA < s2.dim()+1; ++iA) cerr << s2.inodes[iA] << " "; 
	    cerr << " : " << test_bound(s2, fmask, nodes) << endl;*/
          if (test_bound(s2, fmask, nodes)) {
            splx_in.add(splxs.size()); splxs.push_back(s2); 
          }
        }
      } /* simplexes of higher dimension are ignored... */
    }
  }

  void slicer_isovalues::prepare(size_type cv, const mesh_slice::cs_nodes_ct& nodes) 
  {
    pt_in.clear(); pt_bin.clear();
    bgeot::stored_point_tab refpts(nodes.size());
    Uval.resize(nodes.size());
    base_vector coeff;
    base_matrix G;
    pfem pf = mfU->pmf->fem_of_element(cv);
    fem_precomp_ fprecomp;
    if (pf->need_G()) bgeot::vectors_to_base_matrix(G, mfU->pmf->linked_mesh().points_of_convex(cv));
    for (size_type i=0; i < nodes.size(); ++i) refpts[i] = nodes[i].pt_ref;
    fem_precomp_not_stored(pf, &refpts, fprecomp);
    mfU->copy(cv, coeff);
    //cerr << "cv=" << cv << ", val=" << val << ", coeff=" << coeff << endl;
    base_vector v(1); 
    for (size_type i=0; i < nodes.size(); ++i) {
      v[0] = 0;
      pf->interpolation(&fprecomp, i,
			G, mfU->pmf->linked_mesh().trans_of_convex(cv),
			coeff, v, mfU->pmf->get_qdim());
      Uval[i] = v[0];
      pt_bin[i] = (dal::abs(Uval[i] - val) < EPS * val_scaling);
      pt_in[i] = (Uval[i] - val < 0); if (orient>0) pt_in[i] = !pt_in[i]; 
      pt_in[i] = pt_in[i] || pt_bin[i];
      //cerr << "cv=" << cv << ", node["<< i << "]=" << nodes[i].pt << ", Uval[i]=" << Uval[i] << ", pt_in[i]=" << pt_in[i] << ", pt_bin[i]=" << pt_bin[i] << endl;
    }
  }



  /* intersects the simplex with the slice, and (recursively) decomposes it
     into sub-simplices, which are added to the list 'splxs' 
     if 'reduce_dimension' is true, then it is the faces of 
     sub-simplices which are added

     assertion: when called, it will always push *at least* one new simplex on the stack
  */
  void slicer_volume::split_simplex(mesh_slice::cs_nodes_ct& nodes, 
				    mesh_slice::cs_simplexes_ct& splxs, dal::bit_vector& splx_in, 
				    const slice_simplex s, size_type sstart) {
    scalar_type alpha = 0; size_type iA=0, iB = 0;
    bool intersection = false;
    static int level = 0;

    level++;    
    /*
      cerr << "split_simplex: level=" << level << " simplex: ";
      for (iA=0; iA < s.dim()+1 && !intersection; ++iA) 
      cerr << "node#" << s.inodes[iA] << "=" << nodes[s.inodes[iA]].pt 
      << ", in=" << pt_in[s.inodes[iA]] << ", bin=" << pt_bin[s.inodes[iA]] << "; "; cerr << endl;
    */
    assert(level < 100);
    for (iA=0; iA < s.dim()+1; ++iA) {
      for (iB=iA+1; iB < s.dim()+1; ++iB) {
        if (pt_in[s.inodes[iA]] != pt_in[s.inodes[iB]] && !pt_bin[s.inodes[iA]] && !pt_bin[s.inodes[iB]]) {
          alpha=edge_intersect(s.inodes[iA],s.inodes[iB],nodes);
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
      split_simplex(nodes,splxs,splx_in,s1,sstart);
      for (size_type k=0; k < s.dim()+1; k++)
        s1.inodes[k] = (k != iB) ? s.inodes[k] : nn;
      split_simplex(nodes,splxs,splx_in,s1,sstart);
    } else {
      bool all_in = true;
      for (size_type i=0; i < s.dim()+1; ++i) if (!pt_in[s.inodes[i]]) all_in = false;
      //cerr << " -> no intersection , all_in=" << all_in << endl;
      splxs.push_back(s); // even simplexes "outside" are pushed, in case of a slicer_complementary op
      if (all_in && orient != 0) splx_in.add(splxs.size()-1);
      if (orient==0) { /* reduce dimension */
        slice_simplex face(s.dim());
        for (size_type f=0; f < s.dim()+1; ++f) {
          all_in = true;
          for (size_type i=0; i < s.dim(); ++i) {
            size_type p = s.inodes[i + (i<f?0:1)];
            if (!pt_bin[p]) { all_in = false; break; }
            else face.inodes[i] = s.inodes[i + (i<f?0:1)];
          }
          if (all_in) {
            /* prevent addition of a face twice */
            std::sort(face.inodes.begin(), face.inodes.end());
            if (std::find(splxs.begin()+sstart, splxs.end(), face) == splxs.end()) {
              splxs.push_back(face); splx_in.add(splxs.size()-1);
            }
          }
        }
      }
    }
    level--;
  }


  /* nodes : list of nodes (new nodes may be added)
     splxs : list of simplexes (new simplexes may be added)
     splx_in : input: simplexes to take into account, output: list of simplexes inside the slice

     note that the simplexes in the list may have different dimensions
  */
  void slicer_volume::slice(size_type cv, dim_type& fcnt,
                            mesh_slice::cs_nodes_ct& nodes, 
                            mesh_slice::cs_simplexes_ct& splxs, dal::bit_vector& splx_in) {
    /*size_type cnt=0;
    for (mesh_slice::cs_simplexes_ct::iterator it = splxs.begin();
         it != splxs.end(); ++it, ++cnt) {    
      if (!splx_in[cnt]) continue;
    */
    //cerr << "\n\n------------------------------------------------\nslicer_volume::slice : entree, splx_in=" << splx_in << endl;
    if (splx_in.card() == 0) return;
    prepare(cv,nodes);
    for (dal::bv_visitor_c cnt(splx_in); !cnt.finished(); ++cnt) {
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
      } else if (in_cnt != s.dim()+1 || orient==0) {           /* the simplex crosses the slice boundary */
        splx_in.sup(cnt);
        //size_type l = splxs.size();//, n = nodes.size();
        //cerr << "slicer_volume::slice : convex " << cnt << " will be splited" << endl;
	split_simplex(nodes, splxs, splx_in, slice_simplex(s), splxs.size());
        splxs[cnt] = splxs.back(); splxs.pop_back(); splx_in.swap(cnt,splxs.size()); // replace the sliced simplex by one of its slices
        //splx_in.add(l,splxs.size()-l);
      }
    }

    /* signalement des points qui se trouvent pile-poil sur la bordure */
    if (pt_bin.card()) {
      if (fcnt == dim_type(-1)) DAL_THROW(dal::internal_error, 
					  "too much {faces}/{slices faces} in the convex " << cv 
					  << " (nbfaces=" << fcnt << ")");
      for (size_type cnt=pt_bin.take_first(); cnt != size_type(-1); cnt << pt_bin) {
	nodes[cnt].faces[fcnt] = 1;
      }
      fcnt++;
    }
  }

  void slicer_union::slice(size_type cv, dim_type& fcnt,
                           mesh_slice::cs_nodes_ct& nodes, mesh_slice::cs_simplexes_ct& splxs, 
                           dal::bit_vector& splx_in) {
    dal::bit_vector splx_inA = splx_in;
    A->slice(cv,fcnt,nodes,splxs,splx_inA);
    B->slice(cv,fcnt,nodes,splxs,splx_in);
    splx_in |= splx_inA;
  }

  void slicer_intersect::slice(size_type cv, dim_type& fcnt,
			       mesh_slice::cs_nodes_ct& nodes, mesh_slice::cs_simplexes_ct& splxs, 
                               dal::bit_vector& splx_in) {
    A->slice(cv,fcnt,nodes,splxs,splx_in);
    B->slice(cv,fcnt,nodes,splxs,splx_in);
  }

  void slicer_complementary::slice(size_type cv, dim_type& fcnt,
                                   mesh_slice::cs_nodes_ct& nodes, mesh_slice::cs_simplexes_ct& splxs, 
                                   dal::bit_vector& splx_in) {
    dal::bit_vector splx_inA = splx_in;
    size_type sz = splxs.size();
    A->slice(cv, fcnt, nodes, splxs, splx_inA);
    dal::bit_vector bv = splx_in; bv.add(sz, splxs.size()-sz);
    for (dal::bv_visitor_c i(bv); !i.finished(); ++i) {
      /*cerr << "convex " << cv << ",examining simplex #" << i << ": {";
      for (size_type j=0; j < splxs[i].inodes.size(); ++j) cerr << nodes[splxs[i].inodes[j]].pt << " ";
      cerr << "}, splx_in=" << splx_in[i] << "splx_inA=" << splx_inA[i] << endl;*/
      splx_in[i] = !splx_inA[i];
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

  struct slice_node_compare_pt_ref : public std::binary_function<slice_node,slice_node,bool> {
    bool operator()(const slice_node& a, const slice_node& b) const { return a.pt_ref < b.pt_ref; }
  };

  /** 
      build the list of edges (i.e. segments) and store them in a mesh
      object (hence all common nodes/edges are eliminated) slice_edges
      contains the list of edges with where not part of the original
      mesh, but come from the slice of faces.
  */
  void mesh_slice::edges_mesh(getfem_mesh &edges_m, dal::bit_vector& slice_edges) const {
    
    for (size_type ic = 0; ic < cvlst.size(); ++ic) {
      const cs_nodes_ct& n = cvlst[ic].nodes;
      for (cs_simplexes_ct::const_iterator is = cvlst[ic].simplexes.begin();
	   is != cvlst[ic].simplexes.end(); ++is) {
	for (size_type i=0; i < is->dim(); ++i) {
	  for (size_type j=i+1; j <= is->dim(); ++j) {
	    const slice_node& A = n[is->inodes[i]];
	    const slice_node& B = n[is->inodes[j]];
	    if ((A.faces & B.faces).count() >= unsigned(cvlst[ic].cv_dim-1)) {
	      slice_node::faces_ct fmask((1 << cvlst[ic].cv_nbfaces)-1); fmask.flip();
	      size_type e = edges_m.add_segment_by_points(A.pt,B.pt);
	      if (((A.faces & B.faces) & fmask).any()) slice_edges.add(e);
	    }
	  }
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

  /* apply a precompted deformation before slicing */
  class mesh_slice_pre_deform {
    mesh_slice_cv_dof_data_base *defdata;
    pfem pf;
    bgeot::pstored_point_tab pspt;
    fem_precomp_ fprecomp;
    base_matrix G;
    base_vector coeff;
    size_type cv;
  public:
    mesh_slice_pre_deform(mesh_slice_cv_dof_data_base *defdata_) 
      : defdata(defdata_), pf(0), pspt(0), cv(size_type(-1)) {
      if (defdata && defdata->pmf->get_qdim() != defdata->pmf->linked_mesh().dim()) 
        DAL_THROW(dal::dimension_error, "wrong Q(=" << int(defdata->pmf->get_qdim()) 
                  << ") dimension for slice deformation: should be equal to "
                  "the mesh dimension which is " << int(defdata->pmf->linked_mesh().dim()));
    }
    void prepare(size_type cv_, bgeot::stored_point_tab& refpts, bool force_update) {
      cv = cv_;
      if (defdata) {
        if (force_update || defdata->pmf->fem_of_element(cv) != pf) {
          pf = defdata->pmf->fem_of_element(cv);
	  if (pf->need_G()) bgeot::vectors_to_base_matrix(G, defdata->pmf->linked_mesh().points_of_convex(cv));
          fem_precomp_not_stored(pf, &refpts, fprecomp);
        }
        defdata->copy(cv, coeff);
      }
    }

    /* apply the geometric transformation and an optional mesh_fem deformation */
    void apply(const getfem_mesh& m, const bgeot::mesh_structure *cvms, 
               geotrans_precomp_& gp, 
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
	base_vector val(m.dim());
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
                                coeff, val, defdata->pmf->get_qdim());
	      std::copy(val.begin(), val.end(), cv_nodes[pcnt].pt.begin());
            }
            gp.transform(m.points_of_convex(cv), *itp, cv_nodes[pcnt].pt);
            pcnt++;
          }
          *itp = ptsid[*itp];
        }
      }
    }
  }; /* mesh_slice_pre_deform  */

  void mesh_slice::do_slicing(size_type cv, bgeot::pconvex_ref cvr, slicer *ms, cs_nodes_ct cv_nodes, 
			      cs_simplexes_ct cv_simplexes, dal::bit_vector& splx_in) {
    dim_type fcnt = cvr->structure()->nb_faces();
    /* do the slices */
    if (ms) ms->slice(cv, fcnt, cv_nodes, cv_simplexes, splx_in);
    
    set_convex(size_type(-1), cv, cvr, cv_nodes, cv_simplexes, fcnt, splx_in);
  }

  size_type mesh_slice::set_convex(size_type pos, size_type cv, bgeot::pconvex_ref cvr, 
                                   cs_nodes_ct cv_nodes, cs_simplexes_ct cv_simplexes, 
                                   dim_type fcnt, dal::bit_vector& splx_in) {
    /* push the used nodes and simplexes in the final list */
    if (splx_in.card() == 0) return size_type(-1);
    std::vector<size_type> nused(cv_nodes.size(), size_type(-1));
    convex_slice *sc = 0;
    if (pos == size_type(-1)) {
      pos = cvlst.size();
      cvlst.push_back(convex_slice());
      sc = &cvlst.back();
      sc->cv_num = cv;
      sc->cv_dim = cvr->structure()->dim();
      sc->cv_nbfaces = cvr->structure()->nb_faces();
      sc->fcnt = fcnt;
    } else {
      sc = &cvlst[pos];
      assert(sc->cv_num == cv);
    }
    for (size_type snum = splx_in.take_first(); snum != size_type(-1); snum << splx_in) {
      for (size_type i=0; i < cv_simplexes[snum].dim()+1; ++i) {
        size_type lnum = cv_simplexes[snum].inodes[i];
        if (nused[lnum] == size_type(-1)) {
          nused[lnum] = sc->nodes.size(); sc->nodes.push_back(cv_nodes[lnum]);
	  for (size_type k=0; k < sc->nodes.back().pt_ref.size(); ++k) {
	    assert(!isnan(sc->nodes.back().pt[k]));
	    assert(!isnan(sc->nodes.back().pt_ref[k]));
	  }
          points_cnt++;
        }
        cv_simplexes[snum].inodes[i] = nused[lnum];
      }
      simplex_cnt[cv_simplexes[snum].dim()]++;
      sc->simplexes.push_back(cv_simplexes[snum]);
    }
    return pos;
  }

  mesh_slice::mesh_slice(const getfem_mesh& m_) :
    m(m_), simplex_cnt(m.dim()+1, size_type(0)), points_cnt(0), dim_(m.dim()) {
  }

  /* of course, nodes created from edge/slice intersection are almost always duplicated */
  void mesh_slice::build(slicer* ms, size_type nrefine, 
			 convex_face_ct& in_cvlst, mesh_slice_cv_dof_data_base *def_mf_data) 
  {
    if (cvlst.size()) DAL_THROW(dal::failure_error, "non empty slice: should use mesh_slice::merge");
    geotrans_precomp_ gp;
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
      /* apply the geometric transformation and an optional mesh_fem deformation */
      def->apply(m, cvms, gp, cvm_pts, points_on_faces, cv_nodes, cv_simplexes);
      /* and now do the slicing */
      dal::bit_vector splx_in; splx_in.add(0, cv_simplexes.size());
      do_slicing(cv, cvr, ms, cv_nodes, cv_simplexes, splx_in);
    }
  }

  void mesh_slice::build_from_slice(const mesh_slice& sl, slicer* ms) 
  {
    if (cvlst.size()) DAL_THROW(dal::failure_error, "non empty slice: should use mesh_slice::merge");
    if (&sl.linked_mesh() != &m) DAL_THROW(dal::failure_error, "wrong mesh");
    cs_nodes_ct cv_nodes;
    cs_simplexes_ct cv_simplexes;
    for (cvlst_ct::const_iterator it = sl.cvlst.begin(); it != sl.cvlst.end(); ++it) {
      cv_nodes = it->nodes;
      cv_simplexes = it->simplexes;
      dal::bit_vector splx_in; splx_in.add(0, cv_simplexes.size());
      do_slicing(it->cv_num, m.trans_of_convex(it->cv_num)->convex_ref(), ms, cv_nodes, cv_simplexes, splx_in);
    }
  }
  
  void mesh_slice::build_from_points(const std::vector<base_node>& pts, slicer* ms) {
    if (cvlst.size()) DAL_THROW(dal::failure_error,"non empty slice: should use mesh_slice::merge");
    bgeot::geotrans_inv gti;
    gti.add_points(pts);
    dal::dynamic_array<base_node> ptab;
    dal::dynamic_array<size_type> itab;
    cs_nodes_ct cv_nodes;
    cs_simplexes_ct cv_simplexes;
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv) {
      bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
      size_type nb = gti.points_in_convex(m.convex(cv), pgt, ptab, itab);
      if (nb) {
	for (size_type i=0; i < nb; ++i) {
	  cerr << "point " << itab[i] << "(" << pts[itab[i]] << ") trouve dans le convex " << cv << " [pt_ref=" << ptab[i] << "]\n";
	  cv_nodes.push_back(slice_node(pts[itab[i]],ptab[i])); cv_nodes.back().faces=0;
	  cv_simplexes.push_back(slice_simplex(1)); cv_simplexes.back().inodes[0] = cv_nodes.size()-1;
	}
	dal::bit_vector splx_in; splx_in.add(0, cv_simplexes.size());
	do_slicing(cv, pgt->convex_ref(), ms, cv_nodes, cv_simplexes, splx_in);
      }
    }
  }

  void mesh_slice::set_dim(size_type newdim) {
    dim_ = newdim;
    for (size_type ic=0; ic < nb_convex(); ++ic) {
      for (cs_nodes_ct::iterator it=nodes(ic).begin(); it != nodes(ic).end(); ++it) {
	it->pt.resize(newdim);
      }
    }
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

  void mesh_slice::to_mesh(getfem_mesh& mm) const {
    std::vector<size_type> pnums;
    for (cvlst_ct::const_iterator it=cvlst.begin(); it != cvlst.end(); ++it) {
      pnums.clear();
      for (cs_nodes_ct::const_iterator itn=(*it).nodes.begin(); 
	   itn != (*it).nodes.end(); ++itn)
	pnums.push_back(mm.add_point((*itn).pt));
      for (cs_simplexes_ct::const_iterator its=(*it).simplexes.begin();
	   its != (*it).simplexes.end(); ++its)
	mm.add_convex(bgeot::simplex_geotrans((*its).dim(),1), 
		      dal::index_ref_iterator(pnums.begin(),
					      (*its).inodes.begin()));
    }
  }

  size_type mesh_slice::memsize() const {
    size_type sz = sizeof(mesh_slice);
    for (cvlst_ct::const_iterator it = cvlst.begin(); it != cvlst.end(); ++it) {
      sz += sizeof(size_type);
      cerr << "memsize: convex " << it->cv_num << "\n";
      for (size_type i=0; i < it->nodes.size(); ++i) {
	cerr << "  point " << i << ": size+= " << sizeof(slice_node) << "+" <<  
          it->nodes[i].pt.memsize() << "+" << it->nodes[i].pt_ref.memsize() << "-" << sizeof(it->nodes[i].pt)*2 << "\n";
        sz += sizeof(slice_node) + 
          (it->nodes[i].pt.memsize()+it->nodes[i].pt_ref.memsize()) - sizeof(it->nodes[i].pt)*2;
      }
      for (size_type i=0; i < it->simplexes.size(); ++i) {
	cerr << "  simplex " << i << ": size+= " << sizeof(slice_simplex) << "+" << 
          it->simplexes[i].inodes.size()*sizeof(size_type) << "\n";
        sz += sizeof(slice_simplex) + 
          it->simplexes[i].inodes.size()*sizeof(size_type);
      }
    }
    return sz;
  }

}
