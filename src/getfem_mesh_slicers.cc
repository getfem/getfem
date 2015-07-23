/*===========================================================================

 Copyright (C) 2004-2015 Julien Pommier

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

#include "getfem/dal_singleton.h"
#include "getfem/getfem_mesh_slicers.h"
#include "getfem/getfem_mesh_slice.h"
#include "getfem/bgeot_geotrans_inv.h"
#include "getfem/getfem_mesh_level_set.h"

namespace getfem {
  const float slicer_action::EPS = float(1e-13);

  /* ------------------------------ slicers ------------------------------- */

  slicer_none& slicer_none::static_instance() {
    return dal::singleton<slicer_none>::instance();
  }

  /* boundary extraction */
  slicer_boundary::slicer_boundary(const mesh& m, slicer_action &sA, 
                                   const mesh_region& cvflst) : A(&sA) {
    build_from(m,cvflst);
  }

  slicer_boundary::slicer_boundary(const mesh& m, slicer_action &sA) : A(&sA) {
    mesh_region cvflist;
    outer_faces_of_mesh(m, cvflist);
    build_from(m,cvflist);
  }

  void slicer_boundary::build_from(const mesh& m, const mesh_region& cvflst) {
    if (m.convex_index().card()==0) return;
    convex_faces.resize(m.convex_index().last()+1, slice_node::faces_ct(0L));
    for (mr_visitor i(cvflst); !i.finished(); ++i) 
      if (i.is_face()) convex_faces[i.cv()][i.f()]=1;
      else convex_faces[i.cv()].set();
    /* set the mask to 1 for all other possible faces of the convexes, which may 
       appear after slicing the convex, hence they will be part of the "boundary" */
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv) {
      for (short_type f=m.structure_of_convex(cv)->nb_faces(); f < convex_faces[cv].size(); ++f)
        convex_faces[cv][f]=1;
    }
  }

  bool slicer_boundary::test_bound(const slice_simplex& s, slice_node::faces_ct& fmask, const mesh_slicer::cs_nodes_ct& nodes) const {
    slice_node::faces_ct f; f.set();
    for (size_type i=0; i < s.dim()+1; ++i) {
      f &= nodes[s.inodes[i]].faces;
    }
    f &= fmask;
    return (f.any());
  }

  void slicer_boundary::exec(mesh_slicer& ms) {
    if (A) A->exec(ms);
    if (ms.splx_in.card() == 0) return;
    slice_node::faces_ct fmask(ms.cv < convex_faces.size() ? convex_faces[ms.cv] : 0);
    /* quickly check if the convex have any chance to be part of the boundary */
    if (!convex_faces[ms.cv].any()) { ms.splx_in.clear(); return; }

    for (dal::bv_visitor_c cnt(ms.splx_in); !cnt.finished(); ++cnt) {
      const slice_simplex& s = ms.simplexes[cnt];
      if (s.dim() < ms.nodes[0].pt.size()) {
        if (!test_bound(s, fmask, ms.nodes)) ms.splx_in.sup(cnt);
      } else if (s.dim() == 2) {
        ms.sup_simplex(cnt);
        slice_simplex s2(2);
        for (size_type j=0; j < 3; ++j) {
          /* usage of s forbidden in this loop since push_back happens .. */
          static unsigned ord[][2] = {{0,1},{1,2},{2,0}}; /* keep orientation of faces */
          for (size_type k=0; k < 2; ++k) { s2.inodes[k] = ms.simplexes[cnt].inodes[ord[j][k]]; }
          if (test_bound(s2, fmask, ms.nodes)) {
            ms.add_simplex(s2, true);
          }
        }
      } else if (s.dim() == 3) {
        //ms.simplex_orientation(ms.simplexes[cnt]);
        ms.sup_simplex(cnt);
        slice_simplex s2(3);
        for (size_type j=0; j < 4; ++j) {
          /* usage of s forbidden in this loop since push_back happens .. */
          static unsigned ord[][3] = {{0,2,1},{1,2,3},{1,3,0},{0,3,2}}; /* keep orientation of faces */
          for (size_type k=0; k < 3; ++k) { s2.inodes[k] = ms.simplexes[cnt].inodes[ord[j][k]]; } //k + (k<j ? 0 : 1)]; }
          /*cerr << " -> testing "; for (size_type iA=0; iA < s2.dim()+1; ++iA) cerr << s2.inodes[iA] << " "; 
            cerr << " : " << test_bound(s2, fmask, nodes) << endl;*/
          if (test_bound(s2, fmask, ms.nodes)) {
            ms.add_simplex(s2, true);
          }
        }
      } /* simplexes of higher dimension are ignored... */
    }
    ms.update_nodes_index();
  }

  /* apply deformation from a mesh_fem to the nodes */
  void slicer_apply_deformation::exec(mesh_slicer& ms) {
    base_vector coeff;
    base_matrix G;
    bool ref_pts_changed = false;
    if (ms.cvr != ms.prev_cvr
        || defdata->pmf->fem_of_element(ms.cv) != pf) {
      pf = defdata->pmf->fem_of_element(ms.cv);
      if (pf->need_G()) 
        bgeot::vectors_to_base_matrix
          (G, defdata->pmf->linked_mesh().points_of_convex(ms.cv));
    }
    /* check that the points are still the same 
     * -- or recompute the fem_precomp */
    std::vector<base_node> ref_pts2; ref_pts2.reserve(ms.nodes_index.card());
    for (dal::bv_visitor i(ms.nodes_index); !i.finished(); ++i) {
      ref_pts2.push_back(ms.nodes[i].pt_ref);
      if (ref_pts2.size() > ref_pts.size()
          || gmm::vect_dist2_sqr(ref_pts2[i],ref_pts[i])>1e-20) 
        ref_pts_changed = true;
    }
    if (ref_pts2.size() != ref_pts.size()) ref_pts_changed = true;
    if (ref_pts_changed) {
      ref_pts.swap(ref_pts2);
      fprecomp.clear();
    }
    bgeot::pstored_point_tab pspt = store_point_tab(ref_pts);
    pfem_precomp pfp = fprecomp(pf, pspt);
    defdata->copy(ms.cv, coeff);
    
    base_vector val(ms.m.dim());
    size_type cnt = 0;
    fem_interpolation_context ctx(ms.pgt, pfp, 0, G, ms.cv, short_type(-1));
    for (dal::bv_visitor i(ms.nodes_index); !i.finished(); ++i, ++cnt) {
      ms.nodes[i].pt.resize(defdata->pmf->get_qdim());
      ctx.set_ii(cnt);
      pf->interpolation(ctx, coeff, val, defdata->pmf->get_qdim());
      gmm::add(val, ms.nodes[cnt].pt);
    }
  }


  //static bool check_flat_simplex(mesh_slicer::cs_nodes_ct& /*nodes*/, const slice_simplex /*s*/) {
    /*base_matrix M(s.dim(),s.dim());
    for (size_type i=0; i < s.dim(); ++i) {
      for (size_type j=0; j < s.dim(); ++j) {
        M(i,j) = nodes[s.inodes[i+1]].pt_ref[j] - nodes[s.inodes[0]].pt_ref[j];
      }
    }
    scalar_type d = gmm::lu_det(M);
    if (gmm::abs(d) < pow(1e-6,s.dim())) {
      cout.precision(10);
      cout << "!!Flat (" << d << ") :";
      for (size_type i=0; i < s.dim()+1; ++i) cout << " " << nodes[s.inodes[i]].pt;
      cout << "\n";
      return false;
      }*/
  //    return true;
  //}

  /* 
     intersects the simplex with the slice, and (recursively)
     decomposes it into sub-simplices, which are added to the list
     'splxs'. If orient == 0, then it is the faces of sub-simplices
     which are added

     assertion: when called, it will always push *at least* one new
     simplex on the stack.

     remark: s is not reference: on purpose.
  */
  void slicer_volume::split_simplex(mesh_slicer& ms,
                                    slice_simplex s, size_type sstart, 
                                    std::bitset<32> spin, std::bitset<32> spbin) {
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
    for (iA=0; iA < s.dim(); ++iA) {
      if (spbin[iA]) continue;
      for (iB=iA+1; iB < s.dim()+1; ++iB) {
        if (!spbin[iB] && spin[iA] != spin[iB]) {
          alpha=edge_intersect(s.inodes[iA],s.inodes[iB],ms.nodes);
          if (alpha >= 1e-8 && alpha <= 1-1e-8) { intersection = true; break; }
        }
      }
      if (intersection) break;
    }
    if (intersection) {
      /* will call split_simplex recursively */
      const slice_node& A = ms.nodes[s.inodes[iA]]; 
      const slice_node& B = ms.nodes[s.inodes[iB]]; 
      slice_node n(A.pt + alpha*(B.pt-A.pt), A.pt_ref + alpha*(B.pt_ref-A.pt_ref));
      n.faces = A.faces & B.faces;
      size_type nn = ms.nodes.size();
      ms.nodes.push_back(n); /* invalidate A and B.. */
      pt_bin.add(nn); pt_in.add(nn);
      
      std::bitset<32> spin2(spin), spbin2(spbin); 
      std::swap(s.inodes[iA],nn);
      spin2.set(iA); spbin2.set(iA);
      split_simplex(ms, s, sstart, spin2, spbin2);

      std::swap(s.inodes[iA],nn); std::swap(s.inodes[iB],nn);
      spin2 = spin; spbin2 = spbin; spin2.set(iB); spbin2.set(iB);
      split_simplex(ms, s, sstart, spin2, spbin2);

    } else {
      /* end of recursion .. */
      bool all_in = true;
      for (size_type i=0; i < s.dim()+1; ++i) if (!spin[i]) { all_in = false; break; }
      //check_flat_simplex(ms.nodes,s);    

      // even simplexes "outside" are pushed, in case of a slicer_complementary op
      ms.add_simplex(s,(all_in && orient != VOLBOUND) || orient == VOLSPLIT); 
      if (orient==0) { /* reduce dimension */
        slice_simplex face(s.dim());
        for (size_type f=0; f < s.dim()+1; ++f) {
          all_in = true;
          for (size_type i=0; i < s.dim(); ++i) {
            size_type p = i + (i<f?0:1);
            if (!spbin[p]) { all_in = false; break; }
            else face.inodes[i] = s.inodes[p];
          }
          if (all_in) {
            /* prevent addition of a face twice */
            std::sort(face.inodes.begin(), face.inodes.end());
            if (std::find(ms.simplexes.begin()+sstart, ms.simplexes.end(), face) == ms.simplexes.end()) {
              ms.add_simplex(face,true);
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
  void slicer_volume::exec(mesh_slicer& ms) {
    //cerr << "\n----\nslicer_volume::slice : entree, splx_in=" << splx_in << endl;
    if (ms.splx_in.card() == 0) return;
    prepare(ms.cv,ms.nodes,ms.nodes_index);
    for (dal::bv_visitor_c cnt(ms.splx_in); !cnt.finished(); ++cnt) {
      slice_simplex& s = ms.simplexes[cnt];
      /*cerr << "\n--------slicer_volume::slice : slicing convex " << cnt << endl;
      for (size_type i=0; i < s.dim()+1; ++i)
        cerr << "   * pt[" << i << "]=" << nodes[s.inodes[i]].pt << ", is_in=" << 
          is_in(nodes[s.inodes[i]].pt) << ", is_bin=" << is_in(nodes[s.inodes[i]].pt,true) << endl;
      */
      size_type in_cnt = 0, in_bcnt = 0;
      std::bitset<32> spin, spbin;
      for (size_type i=0; i < s.dim()+1; ++i) {
        if (pt_in.is_in(s.inodes[i])) { ++in_cnt; spin.set(i); }
        if (pt_bin.is_in(s.inodes[i])) { ++in_bcnt; spbin.set(i); }
      }

      if (in_cnt == 0) {
        if (orient != VOLSPLIT) ms.splx_in.sup(cnt);
      } else if (in_cnt != s.dim()+1 || orient==VOLBOUND) {           /* the simplex crosses the slice boundary */
        ms.sup_simplex(cnt);
        split_simplex(ms, s, ms.simplexes.size(), spin, spbin);
      }
    }

    /* signalement des points qui se trouvent pile-poil sur la bordure */
    if (pt_bin.card()) {
      GMM_ASSERT1(ms.fcnt != dim_type(-1), 
                  "too much {faces}/{slices faces} in the convex " << ms.cv 
                  << " (nbfaces=" << ms.fcnt << ")");
      for (dal::bv_visitor cnt(pt_bin); !cnt.finished(); ++cnt) {
        ms.nodes[cnt].faces.set(ms.fcnt);
      }
      ms.fcnt++;
    }
    ms.update_nodes_index();
  }

  slicer_mesh_with_mesh::slicer_mesh_with_mesh(const mesh& slm_) :  slm(slm_) { 
    base_node min,max;
    for (dal::bv_visitor cv(slm.convex_index()); !cv.finished(); ++cv) {
      bgeot::bounding_box(min,max,slm.points_of_convex(cv),slm.trans_of_convex(cv));
      tree.add_box(min, max, cv);
    }
  }

  void slicer_mesh_with_mesh::exec(mesh_slicer &ms) {
    /* identify potientially intersecting convexes of slm */
    base_node min(ms.nodes[0].pt),max(ms.nodes[0].pt);
    for (size_type i=1; i < ms.nodes.size(); ++i) {
      for (size_type k=0; k < min.size(); ++k) {
        min[k] = std::min(min[k], ms.nodes[i].pt[k]);
        max[k] = std::max(max[k], ms.nodes[i].pt[k]);
      }
    }
    std::vector<size_type> slmcvs;
    tree.find_intersecting_boxes(min, max, slmcvs);
    /* save context */
    mesh_slicer::cs_simplexes_ct simplexes_final(ms.simplexes); 
    dal::bit_vector splx_in_save(ms.splx_in), 
      simplex_index_save(ms.simplex_index), nodes_index_save(ms.nodes_index); 
    size_type scnt0 = ms.simplexes.size();
    /* loop over candidate convexes of slm */
    //cout << "slicer_mesh_with_mesh: convex " << ms.cv << ", " << ms.splx_in.card() << " candidates\n";
    for (size_type i=0; i < slmcvs.size(); ++i) {
      size_type slmcv = slmcvs[i];
      dim_type fcnt_save = dim_type(ms.fcnt);
      ms.simplexes.resize(scnt0); 
      ms.splx_in = splx_in_save; ms.simplex_index = simplex_index_save; ms.nodes_index = nodes_index_save;
      //cout << "test intersection of " << ms.cv << " and " << slmcv << "\n";
      /* loop over the faces and apply slicer_half_space */
      for (short_type f=0; f<slm.structure_of_convex(slmcv)->nb_faces(); ++f) {
        base_node x0,n;
        if (slm.structure_of_convex(slmcv)->dim() == 3 && slm.dim() == 3) {
          x0 = slm.points_of_face_of_convex(slmcv,f)[0];
          base_node A = slm.points_of_face_of_convex(slmcv,f)[1] - x0;
          base_node B = slm.points_of_face_of_convex(slmcv,f)[2] - x0;
          base_node G = gmm::mean_value(slm.points_of_convex(slmcv).begin(),slm.points_of_convex(slmcv).end());
          n.resize(3);
          n[0] = A[1]*B[2] - A[2]*B[1];
          n[1] = A[2]*B[0] - A[0]*B[2];
          n[2] = A[0]*B[1] - A[1]*B[0];
          if (gmm::vect_sp(n,G-x0) > 0) n *= -1.;
        } else {
          size_type ip = slm.structure_of_convex(slmcv)->nb_points_of_face(f) / 2;
          x0 = slm.points_of_face_of_convex(slmcv,f)[ip];
          n = slm.normal_of_face_of_convex(slmcv,f, x0);
        }
        slicer_half_space slf(x0,n,slicer_volume::VOLIN);
        slf.exec(ms);
        if (ms.splx_in.card() == 0) break;
      }
      if (ms.splx_in.card()) intersection_callback(ms, slmcv);
      for (dal::bv_visitor is(ms.splx_in); !is.finished(); ++is) {
        simplexes_final.push_back(ms.simplexes[is]);
      }
      ms.fcnt=fcnt_save;
    }
    ms.splx_in.clear(); ms.splx_in.add(scnt0, simplexes_final.size()-scnt0); 
    ms.simplexes.swap(simplexes_final);
    ms.simplex_index = ms.splx_in;
    ms.update_nodes_index();
    /*cout << "convex " << ms.cv << "was sliced into " << ms.splx_in.card() 
         << " simplexes, nodes.size=" << ms.nodes.size() 
         << ", used nodes=" << ms.nodes_index.card() << "\n";*/
  }

  /* isosurface computations */
  void slicer_isovalues::prepare(size_type cv,
                                 const mesh_slicer::cs_nodes_ct& nodes, 
                                 const dal::bit_vector& nodes_index) {
    pt_in.clear(); pt_bin.clear();
    std::vector<base_node> refpts(nodes.size());
    Uval.resize(nodes.size());
    base_vector coeff;
    base_matrix G;
    pfem pf = mfU->pmf->fem_of_element(cv);
    if (pf == 0) return;
    fem_precomp_pool fprecomp;
    if (pf->need_G()) 
      bgeot::vectors_to_base_matrix
        (G,mfU->pmf->linked_mesh().points_of_convex(cv));
    for (size_type i=0; i < nodes.size(); ++i) refpts[i] = nodes[i].pt_ref;
    pfem_precomp pfp = fprecomp(pf, store_point_tab(refpts));
    mfU->copy(cv, coeff);
    //cerr << "cv=" << cv << ", val=" << val << ", coeff=" << coeff << endl;
    base_vector v(1); 
    fem_interpolation_context ctx(mfU->pmf->linked_mesh().trans_of_convex(cv),
                                  pfp, 0, G, cv, short_type(-1));
    for (dal::bv_visitor i(nodes_index); !i.finished(); ++i) {
      v[0] = 0;
      ctx.set_ii(i);
      pf->interpolation(ctx, coeff, v, mfU->pmf->get_qdim());
      Uval[i] = v[0];
      // optimisable -- les bit_vectors sont lents..
      pt_bin[i] = (gmm::abs(Uval[i] - val) < EPS * val_scaling);
      pt_in[i] = (Uval[i] - val < 0); if (orient>0) pt_in[i] = !pt_in[i]; 
      pt_in[i] = pt_in[i] || pt_bin[i];
      // cerr << "cv=" << cv << ", node["<< i << "]=" << nodes[i].pt
      //      << ", Uval[i]=" << Uval[i] << ", pt_in[i]=" << pt_in[i]
      //      << ", pt_bin[i]=" << pt_bin[i] << endl;
    }
  }


  void slicer_union::exec(mesh_slicer &ms) {
    dal::bit_vector splx_in_base = ms.splx_in;
    size_type c = ms.simplexes.size();
    dim_type fcnt_0 = dim_type(ms.fcnt);
    A->exec(ms); 
    dal::bit_vector splx_inA(ms.splx_in);
    dim_type fcnt_1 = dim_type(ms.fcnt);

    dal::bit_vector splx_inB = splx_in_base;
    splx_inB.add(c, ms.simplexes.size()-c);
    splx_inB.setminus(splx_inA);
    for (dal::bv_visitor_c i(splx_inB); !i.finished(); ++i) {
      if (!ms.simplex_index[i]) splx_inB.sup(i);
    }
    //splx_inB &= ms.simplex_index;
    ms.splx_in = splx_inB;
    B->exec(ms); splx_inB = ms.splx_in;
    ms.splx_in |= splx_inA;

    /* 
       the boring part : making sure that the "slice face" node markers
       are correctly set
    */
    for (unsigned f=fcnt_0; f < fcnt_1; ++f) {
      for (dal::bv_visitor i(splx_inB); !i.finished(); ++i) {
        for (unsigned j=0; j < ms.simplexes[i].dim()+1; ++j) {
          bool face_boundA = true;
          for (unsigned k=0; k < ms.simplexes[i].dim()+1; ++k) {
            if (j != k && !ms.nodes[ms.simplexes[i].inodes[k]].faces[f]) {
              face_boundA = false; break;
            }
          }
          if (face_boundA) {
            /* now we know that the face was on slice A boundary, and
               that the convex is inside slice B. The conclusion: the
               face is not on a face of A union B.
            */
            for (unsigned k=0; k < ms.simplexes[i].dim()+1; ++k)
              if (j != k) ms.nodes[ms.simplexes[i].inodes[k]].faces[f] = false;            
          }
        }
      }
    }
    ms.update_nodes_index();
  }

  void slicer_intersect::exec(mesh_slicer& ms) {
    A->exec(ms);
    B->exec(ms);
  }

  void slicer_complementary::exec(mesh_slicer& ms) {
    dal::bit_vector splx_inA = ms.splx_in;
    size_type sz = ms.simplexes.size();
    A->exec(ms); splx_inA.swap(ms.splx_in);
    ms.splx_in &= ms.simplex_index;
    dal::bit_vector bv = ms.splx_in; bv.add(sz, ms.simplexes.size()-sz); bv &= ms.simplex_index;
    for (dal::bv_visitor_c i(bv); !i.finished(); ++i) {
      /*cerr << "convex " << cv << ",examining simplex #" << i << ": {";
      for (size_type j=0; j < simplexes[i].inodes.size(); ++j) cerr << nodes[simplexes[i].inodes[j]].pt << " ";
      cerr << "}, splx_in=" << splx_in[i] << "splx_inA=" << splx_inA[i] << endl;*/
      ms.splx_in[i] = !splx_inA.is_in(i);
    }
  }

  void slicer_compute_area::exec(mesh_slicer &ms) {
    for (dal::bv_visitor is(ms.splx_in); !is.finished(); ++is) {
      const slice_simplex& s = ms.simplexes[is];
        base_matrix M(s.dim(),s.dim());
      for (size_type i=0; i < s.dim(); ++i) 
        for (size_type j=0; j < s.dim(); ++j)
          M(i,j) = ms.nodes[s.inodes[i+1]].pt[j] - ms.nodes[s.inodes[0]].pt[j];
      scalar_type v = gmm::abs(gmm::lu_det(M));
      for (size_type d=2; d <= s.dim(); ++d) v /= scalar_type(d);
      a += v;
    }
  }

  void slicer_build_edges_mesh::exec(mesh_slicer &ms) {
    for (dal::bv_visitor is(ms.splx_in); !is.finished(); ++is) {
      const slice_simplex& s = ms.simplexes[is];
      for (size_type i=0; i < s.dim(); ++i) {
        for (size_type j=i+1; j <= s.dim(); ++j) {
          const slice_node& A = ms.nodes[s.inodes[i]];
          const slice_node& B = ms.nodes[s.inodes[j]];
          /* duplicate with stored_mesh_slice which also 
             builds a list of edges */
          if ((A.faces & B.faces).count() >= unsigned(ms.cv_dim-1)) {
            slice_node::faces_ct fmask((1 << ms.cv_nbfaces)-1); fmask.flip();
            size_type e = edges_m.add_segment_by_points(A.pt,B.pt);
            if (pslice_edges && (((A.faces & B.faces) & fmask).any())) pslice_edges->add(e);
          }
        }
      }
    }
  }

  void slicer_build_mesh::exec(mesh_slicer &ms) {
    std::vector<size_type> pid(ms.nodes_index.last_true()+1);
    for (dal::bv_visitor i(ms.nodes_index); !i.finished(); ++i)
      pid[i] = m.add_point(ms.nodes[i].pt);
    for (dal::bv_visitor i(ms.splx_in); !i.finished(); ++i) {
      for (unsigned j=0; j < ms.simplexes.at(i).inodes.size(); ++j) {
        assert(m.points_index().is_in(pid.at(ms.simplexes.at(i).inodes[j])));
      }
      m.add_convex(bgeot::simplex_geotrans(ms.simplexes[i].dim(),1),
                   gmm::index_ref_iterator(pid.begin(),
                                           ms.simplexes[i].inodes.begin()));
    }
  }

  void slicer_explode::exec(mesh_slicer &ms) {
    if (ms.nodes_index.card() == 0) return;    

    base_node G;
    if (ms.face < dim_type(-1))
      G = gmm::mean_value(ms.m.points_of_face_of_convex(ms.cv, ms.face).begin(), 
                          ms.m.points_of_face_of_convex(ms.cv, ms.face).end());
    else
      G = gmm::mean_value(ms.m.points_of_convex(ms.cv).begin(), 
                          ms.m.points_of_convex(ms.cv).end());    
    for (dal::bv_visitor i(ms.nodes_index); !i.finished(); ++i)
      ms.nodes[i].pt = G + coef*(ms.nodes[i].pt - G);

    for (dal::bv_visitor cnt(ms.splx_in); !cnt.finished(); ++cnt) {
      const slice_simplex& s = ms.simplexes[cnt];
      if (s.dim() == 3) { // keep only faces
              ms.sup_simplex(cnt);
        slice_simplex s2(3);
        for (size_type j=0; j < 4; ++j) {
          /* usage of s forbidden in this loop since push_back happens .. */
          static unsigned ord[][3] = {{0,2,1},{1,2,3},{1,3,0},{0,3,2}}; /* keep orientation of faces */
          for (size_type k=0; k < 3; ++k) { s2.inodes[k] = ms.simplexes[cnt].inodes[ord[j][k]]; } //k + (k<j ? 0 : 1)]; }

          slice_node::faces_ct f; f.set();
          for (size_type i=0; i < s2.dim()+1; ++i) {
            f &= ms.nodes[s2.inodes[i]].faces;
          }
          if (f.any()) {
            ms.add_simplex(s2, true);
          }
        }
      }
    }
  }

  /* -------------------- member functions of mesh_slicer -------------- */

  mesh_slicer::mesh_slicer(const mesh_level_set &mls_) :
    m(mls_.linked_mesh()), mls(&mls_), pgt(0), cvr(0) {}
  mesh_slicer::mesh_slicer(const mesh& m_) : 
    m(m_), mls(0), pgt(0), cvr(0) {}

  void mesh_slicer::using_mesh_level_set(const mesh_level_set &mls_) { 
    mls = &mls_;
    GMM_ASSERT1(&m == &mls->linked_mesh(), "different meshes");
  }

  void mesh_slicer::pack() {
    std::vector<size_type> new_id(nodes.size());
    size_type ncnt = 0;
    for (dal::bv_visitor i(nodes_index); !i.finished(); ++i) {
      if (i != ncnt) {
        nodes[i].swap(nodes[ncnt]);
      }
      new_id[i] = ncnt++;
    }
    nodes.resize(ncnt);
    size_type scnt = 0;
    for (dal::bv_visitor j(simplex_index); !j.finished(); ++j) {
      if (j != scnt) { simplexes[scnt] = simplexes[j]; }
      for (std::vector<size_type>::iterator it = simplexes[scnt].inodes.begin(); 
           it != simplexes[scnt].inodes.end(); ++it) {
        *it = new_id[*it];
      }
    }
    simplexes.resize(scnt);
  }

  void mesh_slicer::update_nodes_index() {
    nodes_index.clear();
    for (dal::bv_visitor j(simplex_index); !j.finished(); ++j) {
      assert(j < simplexes.size());
      for (std::vector<size_type>::iterator it = simplexes[j].inodes.begin(); 
           it != simplexes[j].inodes.end(); ++it) {
        assert(*it < nodes.size());
        nodes_index.add(*it);
      }
    }
  }

  static void flag_points_on_faces(const bgeot::pconvex_ref& cvr, 
                                   const std::vector<base_node>& pts, 
                                   std::vector<slice_node::faces_ct>& faces) {
    GMM_ASSERT1(cvr->structure()->nb_faces() <= 32,
                "won't work for convexes with more than 32 faces "
                "(hardcoded limit)");
    faces.resize(pts.size());
    for (size_type i=0; i < pts.size(); ++i) {
      faces[i].reset();      
      for (short_type f=0; f < cvr->structure()->nb_faces(); ++f)
        faces[i][f] = (gmm::abs(cvr->is_in_face(f, pts[i])) < 1e-10);
    }
  }

  void mesh_slicer::update_cv_data(size_type cv_, short_type f_) {
    cv = cv_;
    face = f_;
    pgt = m.trans_of_convex(cv);
    prev_cvr = cvr; cvr = pgt->convex_ref();      
    cv_dim = cvr->structure()->dim();
    cv_nbfaces = cvr->structure()->nb_faces();
    fcnt = cvr->structure()->nb_faces();
    discont = (mls && mls->is_convex_cut(cv));
  }

  void mesh_slicer::apply_slicers() {
    simplex_index.clear(); simplex_index.add(0, simplexes.size());
    splx_in = simplex_index;
    nodes_index.clear(); nodes_index.add(0, nodes.size());      
    for (size_type i=0; i < action.size(); ++i) {
      action[i]->exec(*this);
      //cout << "simplex_index=" << simplex_index << "\n   splx_in=" << splx_in << "\n";
      assert(simplex_index.contains(splx_in));
    }
  }

  void mesh_slicer::simplex_orientation(slice_simplex& s) {
    size_type N = m.dim();
    if (s.dim() == N) {
      base_matrix M(N,N);
      for (size_type i=1; i < N+1; ++i) {
        base_small_vector d = nodes[s.inodes[i]].pt - nodes[s.inodes[0]].pt;
        gmm::copy_n(d.const_begin(), N, M.begin() + (i-1)*N);
      }
      scalar_type J = lu_det(M);
      //cout << " lu_det = " << J << "\n";        
      if (J < 0) {
        std::swap(s.inodes[1],s.inodes[0]);
      }
    }
  }

  void mesh_slicer::exec(size_type nrefine, const mesh_region& cvlst) {
    short_type n = short_type(nrefine);
    exec_(&n, 0, cvlst);
  }

  void mesh_slicer::exec(const std::vector<short_type> &nrefine, const mesh_region& cvlst) {
    exec_(&nrefine[0], 1, cvlst);
  }

  static bool check_orient(size_type cv, bgeot::pgeometric_trans pgt, const mesh& m) {
    if (pgt->dim() == m.dim() && m.dim()>=2) { /* no orient check for 
                                                  convexes of lower dim */
      base_matrix G; bgeot::vectors_to_base_matrix(G,m.points_of_convex(cv));
      base_node g(pgt->dim()); g.fill(.5); 
      base_matrix pc; pgt->poly_vector_grad(g,pc);
      base_matrix K(pgt->dim(),pgt->dim());
      gmm::mult(G,pc,K);
      scalar_type J = gmm::lu_det(K);
      // bgeot::geotrans_interpolation_context ctx(pgp,0,G);
      // scalar_type J = gmm::lu_det(ctx.B()); // pb car inverse K même
      if (J < 0) return true;
      //cout << "cv = " << cv << ", J = " << J << "\n";
    }
    return false;
  }

#if OLD_MESH_SLICER
  void mesh_slicer::exec_(const short_type *pnrefine, int nref_stride, const mesh_region& cvlst) {
    std::vector<base_node> cvm_pts;
    const bgeot::basic_mesh *cvm = 0;
    const bgeot::mesh_structure *cvms = 0;
    bgeot::geotrans_precomp_pool gppool;
    bgeot::pgeotrans_precomp pgp = 0;
    std::vector<slice_node::faces_ct> points_on_faces;

    cvlst.from_mesh(m);
    size_type prev_nrefine = 0;
    for (mr_visitor it(cvlst); !it.finished(); ++it) {
      size_type nrefine = pnrefine[it.cv()*nref_stride];
      update_cv_data(it.cv(),it.f());      
      bool revert_orientation = check_orient(cv, pgt,m);

      /* update structure-dependent data */
      if (prev_cvr != cvr || nrefine != prev_nrefine) {
        cvm = bgeot::refined_simplex_mesh_for_convex(cvr, nrefine);
        cvm_pts.resize(cvm->points().card());
        std::copy(cvm->points().begin(), cvm->points().end(), cvm_pts.begin());
        pgp = gppool(pgt, store_point_tab(cvm_pts));
        flag_points_on_faces(cvr, cvm_pts, points_on_faces);
        prev_nrefine = nrefine;
      }
      if (face < dim_type(-1))
        cvms = bgeot::refined_simplex_mesh_for_convex_faces(cvr, nrefine)[face];
      else
        cvms = cvm; 

      /* apply the initial geometric transformation */
      std::vector<size_type> ptsid(cvm_pts.size()); std::fill(ptsid.begin(), ptsid.end(), size_type(-1));
      simplexes.resize(cvms->nb_convex());
      nodes.resize(0);
      for (size_type snum = 0; snum < cvms->nb_convex(); ++snum) { 
        /* cvms should not contain holes in its convex index.. */
        simplexes[snum].inodes.resize(cvms->nb_points_of_convex(snum));
        std::copy(cvms->ind_points_of_convex(snum).begin(),
                  cvms->ind_points_of_convex(snum).end(), simplexes[snum].inodes.begin());
        if (revert_orientation) std::swap(simplexes[snum].inodes[0],simplexes[snum].inodes[1]);
        /* store indices of points which are really used , and renumbers them */
        for (std::vector<size_type>::iterator itp = simplexes[snum].inodes.begin();
             itp != simplexes[snum].inodes.end(); ++itp) {
          if (ptsid[*itp] == size_type(-1)) {
            ptsid[*itp] = nodes.size();
            nodes.push_back(slice_node());
            nodes.back().pt_ref = cvm_pts[*itp];
            nodes.back().faces = points_on_faces[*itp];
            nodes.back().pt.resize(m.dim()); nodes.back().pt.fill(0.);
            pgp->transform(m.points_of_convex(cv), *itp, nodes.back().pt);
          }
          *itp = ptsid[*itp];
        }
        if (0) { 
          static int once = 0;
          if (once++ < 3) {
            cout << "check orient cv " << cv << ", snum=" << snum << "/" << cvms->nb_convex();
          }
          simplex_orientation(simplexes[snum]);
        }
      }
      apply_slicers();
    }
  }
#endif // OLD_MESH_SLICER

  template <typename CONT> 
  static void add_degree1_convex(bgeot::pgeometric_trans pgt, const CONT &pts, 
                                 mesh &m) {
    unsigned N = pgt->dim();
    std::vector<base_node> v; v.reserve(N+1);
    for (unsigned i=0; i < pgt->nb_points(); ++i) {
      const base_node P = pgt->convex_ref()->points()[i];
      scalar_type s = 0; 
      for (unsigned j=0; j < N; ++j) { 
        s += P[j]; if (P[j] == 1) { v.push_back(pts[i]); break; }
      }
      if (s == 0) v.push_back(pts[i]);
    }
    assert(v.size() == N+1);
    base_node G = gmm::mean_value(v);
    /*for (unsigned i=0; i < v.size();++i) 
      v[i] = v[i] + 0.1 * (G - v[i]);*/
    m.add_convex_by_points(bgeot::simplex_geotrans(N,1), v.begin());
  }

  const mesh& mesh_slicer::refined_simplex_mesh_for_convex_cut_by_level_set
  (const mesh &cvm, unsigned nrefine) {
    mesh mm; mm.copy_from(cvm);
    while (nrefine > 1) { 
      mm.Bank_refine(mm.convex_index());
      nrefine /= 2;
    }

    std::vector<size_type> idx;
    tmp_mesh.clear();
    //cerr << "nb cv = " << tmp_mesh.convex_index().card() << "\n";
    for (dal::bv_visitor_c ic(mm.convex_index()); !ic.finished(); ++ic) {
      add_degree1_convex(mm.trans_of_convex(ic), mm.points_of_convex(ic), tmp_mesh);
    }
    /*tmp_mesh.write_to_file(std::cerr);
      cerr << "\n";*/
    return tmp_mesh;
  }
  
  const bgeot::mesh_structure &
  mesh_slicer::refined_simplex_mesh_for_convex_faces_cut_by_level_set
  (short_type f) {
    mesh &cvm = tmp_mesh;
    tmp_mesh_struct.clear();
    unsigned N = cvm.dim();
    
    dal::bit_vector pt_in_face; pt_in_face.sup(0, cvm.points().index().last_true()+1);
    for (dal::bv_visitor ip(cvm.points().index()); !ip.finished(); ++ip)
      if (gmm::abs(cvr->is_in_face(short_type(f), cvm.points()[ip]))) pt_in_face.add(ip);

    for (dal::bv_visitor_c ic(cvm.convex_index()); !ic.finished(); ++ic) {
      for (short_type ff=0; ff < cvm.nb_faces_of_convex(ic); ++ff) {
        bool face_ok = true;
        for (unsigned i=0; i < N; ++i) {
          if (!pt_in_face.is_in(cvm.ind_points_of_face_of_convex(ic,ff)[i])) {
            face_ok = false; break;
          }
        }
        if (face_ok) {
          tmp_mesh_struct.add_convex(bgeot::simplex_structure(dim_type(N-1)), 
                                     cvm.ind_points_of_face_of_convex(ic, ff).begin());
        }
      }
    }
    return tmp_mesh_struct;
  }

  void mesh_slicer::exec_(const short_type *pnrefine, 
                          int nref_stride, 
                          const mesh_region& cvlst) {
    std::vector<base_node> cvm_pts;
    const bgeot::basic_mesh *cvm = 0;
    const bgeot::mesh_structure *cvms = 0;
    bgeot::geotrans_precomp_pool gppool;
    bgeot::pgeotrans_precomp pgp = 0;
    std::vector<slice_node::faces_ct> points_on_faces;
    bool prev_discont = true;

    cvlst.from_mesh(m);
    size_type prev_nrefine = 0;
    // size_type prev_cv = size_type(-1);
    for (mr_visitor it(cvlst); !it.finished(); ++it) {
      size_type nrefine = pnrefine[it.cv()*nref_stride];
      update_cv_data(it.cv(),it.f());
      bool revert_orientation = check_orient(cv, pgt,m);

      /* update structure-dependent data */
      
      /* TODO : fix levelset handling when slicing faces .. */
      if (prev_cvr != cvr || nrefine != prev_nrefine || discont || prev_discont) {
        if (discont) {
          //if (prev_cv != it.cv())
          cvm = &refined_simplex_mesh_for_convex_cut_by_level_set
            (mls->mesh_of_convex(cv), unsigned(nrefine));
        } else {
          cvm = bgeot::refined_simplex_mesh_for_convex(cvr,
                                                       short_type(nrefine));
        }
        cvm_pts.resize(cvm->points().card());
        std::copy(cvm->points().begin(), cvm->points().end(), cvm_pts.begin());
        pgp = gppool(pgt, store_point_tab(cvm_pts));
        flag_points_on_faces(cvr, cvm_pts, points_on_faces);
        prev_nrefine = nrefine;
      }
      if (face < dim_type(-1)) {
        if (!discont) {
          cvms = bgeot::refined_simplex_mesh_for_convex_faces(cvr, short_type(nrefine))[face];
        } else {
          cvms = &refined_simplex_mesh_for_convex_faces_cut_by_level_set(face);
        }
      } else {
        cvms = cvm; 
      }



      /* apply the initial geometric transformation */
      std::vector<size_type> ptsid(cvm_pts.size()); std::fill(ptsid.begin(), ptsid.end(), size_type(-1));
      simplexes.resize(cvms->nb_convex());
      nodes.resize(0);

      base_node G;
      for (size_type snum = 0; snum < cvms->nb_convex(); ++snum) { 
        /* cvms should not contain holes in its convex index.. */
        simplexes[snum].inodes.resize(cvms->nb_points_of_convex(snum));
        std::copy(cvms->ind_points_of_convex(snum).begin(),
                  cvms->ind_points_of_convex(snum).end(), simplexes[snum].inodes.begin());
        if (revert_orientation) std::swap(simplexes[snum].inodes[0],simplexes[snum].inodes[1]);
        /* store indices of points which are really used , and renumbers them */
        if (discont) {
          G.resize(m.dim()); G.fill(0.);
          for (std::vector<size_type>::iterator itp = 
                 simplexes[snum].inodes.begin();
               itp != simplexes[snum].inodes.end(); ++itp) {
            G += cvm_pts[*itp];
          }
          G /= scalar_type(simplexes[snum].inodes.size());
        }
          

        for (std::vector<size_type>::iterator itp = 
               simplexes[snum].inodes.begin();
             itp != simplexes[snum].inodes.end(); ++itp) {
          if (discont || ptsid[*itp] == size_type(-1)) {
            ptsid[*itp] = nodes.size();
            nodes.push_back(slice_node());
            if (!discont) {
              nodes.back().pt_ref = cvm_pts[*itp];
            } else {
              /* displace the ref point such that one will not interpolate
                 on the discontinuity (yes this is quite ugly and not 
                 robust) 
              */
              nodes.back().pt_ref = cvm_pts[*itp] + 0.01*(G - cvm_pts[*itp]);
            }
            nodes.back().faces = points_on_faces[*itp];
            nodes.back().pt.resize(m.dim()); nodes.back().pt.fill(0.);
            pgp->transform(m.points_of_convex(cv), *itp, nodes.back().pt);
            //nodes.back().pt = pgt->transform(G, m.points_of_convex(cv));
            //cerr << "G = " << G << " -> pt = " << nodes.back().pt << "\n";
          }
          *itp = ptsid[*itp];
        }
      }
      //cerr << "cv = " << cv << ", cvm.nb_points_ = "<< cvm->points().size() << ", nbnodes = " << nodes.size() << ", nb_simpl=" << simplexes.size() << "\n";

      apply_slicers();
      // prev_cv = it.cv();
      prev_discont = discont;
    }
  }


  void mesh_slicer::exec(size_type nrefine) {
    exec(nrefine,mesh_region(m.convex_index()));
  }
  
  /* apply slice ops to an already stored slice object */
  void mesh_slicer::exec(const stored_mesh_slice& sl) {
    GMM_ASSERT1(&sl.linked_mesh() == &m, "wrong mesh");
    for (stored_mesh_slice::cvlst_ct::const_iterator it = sl.cvlst.begin(); it != sl.cvlst.end(); ++it) {
      update_cv_data((*it).cv_num);
      nodes = (*it).nodes;
      simplexes = (*it).simplexes;
      apply_slicers();
    }
  }
  
  /* apply slice ops to a set of nodes */
  void mesh_slicer::exec(const std::vector<base_node>& pts) {
    bgeot::geotrans_inv gti;
    gti.add_points(pts);
    dal::dynamic_array<base_node> ptab;
    dal::dynamic_array<size_type> itab;
    for (dal::bv_visitor ic(m.convex_index()); !ic.finished(); ++ic) {
      size_type nb = gti.points_in_convex(m.convex(ic), m.trans_of_convex(ic), ptab, itab);
      if (!nb) continue;
      update_cv_data(ic);
      nodes.resize(0); simplexes.resize(0);
      for (size_type i=0; i < nb; ++i) {
        //cerr << "point " << itab[i] << "(" << pts[itab[i]] 
        //<< ") trouve dans le convex " << ic << " [pt_ref=" << ptab[i] << "]\n";
        nodes.push_back(slice_node(pts[itab[i]],ptab[i])); nodes.back().faces=0;
        slice_simplex s(1); s.inodes[0] = nodes.size()-1;
        simplexes.push_back(s);
      }
      apply_slicers();
    }
  }
}
