/*===========================================================================

 Copyright (C) 1999-2015 Yves Renard

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



#include "getfem/bgeot_mesh_structure.h"

namespace bgeot {


  dal::bit_vector mesh_structure::convex_index(dim_type n) const {
    dal::bit_vector res = convex_tab.index();
    for (dal::bv_visitor cv(convex_tab.index()); !cv.finished(); ++cv)
      if (structure_of_convex(cv)->dim() != n) res.sup(cv);
    return res;
  }

  size_type mesh_structure::local_ind_of_convex_point(size_type ic,
                                                      size_type ip) const {
    mesh_convex_structure::ind_pt_ct::const_iterator it;
    size_type ind = 0;
    for (it=convex_tab[ic].pts.begin();
         it != convex_tab[ic].pts.end() && (*it) != ip; ++it) ++ind;
    GMM_ASSERT1(it != convex_tab[ic].pts.end(),
                "This point does not exist on this convex.");
    return ind;
  }

  void mesh_structure::swap_points(size_type i, size_type j) {
    if (i == j) return;
    std::vector<size_type> doubles;

    for (size_type k = 0; k < points_tab[i].size(); ++k) {
      size_type cv = points_tab[i][k];
      for (size_type l = 0; l < convex_tab[cv].pts.size(); ++l) {
        size_type &ind = convex_tab[cv].pts[l];
        if (ind == i) ind = j;
        else if (ind == j) { ind = i; doubles.push_back(cv); }
      }
    }
    for (size_type k = 0; k < points_tab[j].size(); ++k) {
      size_type cv = points_tab[j][k];
      if (std::find(doubles.begin(), doubles.end(), cv) == doubles.end()) {
        for (size_type l = 0; l < convex_tab[cv].pts.size(); ++l)
          if (convex_tab[cv].pts[l] == j) convex_tab[cv].pts[l] = i;
      }
    }
    points_tab.swap(i,j);
  }

  void mesh_structure::swap_convex(size_type i, size_type j) {
    if (i == j) return;
    std::vector<size_type> doubles;

    if (is_convex_valid(i))
      for (size_type k = 0; k < convex_tab[i].pts.size(); ++k) {
        size_type ip = convex_tab[i].pts[k];
        for (size_type l = 0; l < points_tab[ip].size(); ++l) {
          size_type &ind = points_tab[ip][l];
          if (ind == i) ind = j;
          else if (ind == j) { ind = i; doubles.push_back(ip); }
        }
      }
    if (is_convex_valid(j))
      for (size_type k = 0; k < convex_tab[j].pts.size(); ++k) {
        size_type ip = convex_tab[j].pts[k];
        if (std::find(doubles.begin(), doubles.end(), ip) == doubles.end()) {
          for (size_type l = 0; l < points_tab[ip].size(); ++l)
            if (points_tab[ip][l] == j) points_tab[ip][l] = i;
        }
      }
    convex_tab.swap(i,j);
  }

  size_type mesh_structure::add_segment(size_type a, size_type b) {
    static pconvex_structure cs = NULL;
    if (!cs) cs = simplex_structure(1);
    size_type t[2]; t[0] = a; t[1] = b;
    return add_convex(cs, &t[0]);
  }

  void mesh_structure::sup_convex(size_type ic) {
    if (!(is_convex_valid(ic))) return;
    for (size_type l = 0; l < convex_tab[ic].pts.size(); ++l) {
      size_type &ind = convex_tab[ic].pts[l];
      std::vector<size_type>::iterator it1= points_tab[ind].begin(), it2 = it1;
      std::vector<size_type>::iterator ite= points_tab[ind].end();
      for (; it2 != ite; ++it2) { *it1 = *it2; if (*it1 != ic) ++it1; }
      points_tab[ind].pop_back();
    }
    convex_tab.sup(ic);
  }

  size_type mesh_structure::add_face_of_convex(size_type ic, short_type f) {
    return add_convex( (structure_of_convex(ic)->faces_structure())[f],
                  ind_points_of_face_of_convex(ic, f).begin());
  }

  void mesh_structure::add_faces_of_convex(size_type ic) {
    pconvex_structure ps = structure_of_convex(ic);
    for (short_type iff = 0; iff < ps->nb_faces(); ++iff)
      add_convex( (ps->faces_structure())[iff],
                  ind_points_of_face_of_convex(ic, iff).begin());
  }

  void mesh_structure::to_faces(dim_type n) { // not efficient
    dal::bit_vector nn = convex_index();
    for (dal::bv_visitor i(nn); !i.finished(); ++i)
      if (convex_tab[i].cstruct->dim() == n)
        { add_faces_of_convex(i); sup_convex(i); }
  }

  void mesh_structure::to_edges(void) { // not efficient at all !
    dim_type dmax = 0;
    dal::bit_vector nn = convex_index();
    for (dal::bv_visitor i(nn); !i.finished(); ++i)
      dmax = std::max(dmax, convex_tab[i].cstruct->dim());
    for ( ; dmax > 1; --dmax) to_faces(dmax);
  }

  size_type mesh_structure::nb_convex_with_edge(size_type i1, size_type i2) {
    size_type nb = 0;
    for (size_type k = 0; k < points_tab[i1].size(); ++k) {
      size_type cv = points_tab[i1][k];
      for (size_type l = 0; l < convex_tab[cv].pts.size(); ++l)
        if (convex_tab[cv].pts[l] == i2) { ++nb; break; }
    }
    return nb;
  }

  void mesh_structure::convex_with_edge(size_type i1, size_type i2,
                                        std::vector<size_type> &ipt) {
    ipt.resize(0);
    for (size_type k = 0; k < points_tab[i1].size(); ++k) {
      size_type cv = points_tab[i1][k];
      for (size_type l = 0; l < convex_tab[cv].pts.size(); ++l)
        if (convex_tab[cv].pts[l] == i2) { ipt.push_back(cv); break; }
    }
  }

  void mesh_structure::ind_points_to_point(size_type ip, ind_set &s) const {
    // not efficient
    s.resize(0);
    for (size_type k = 0; k < points_tab[ip].size(); ++k) {
      size_type cv = points_tab[ip][k];
      for (size_type l = 0; l < convex_tab[cv].pts.size(); ++l)
        if (ip != convex_tab[cv].pts[l]) {
          size_type ind = convex_tab[cv].pts[l];
          if (std::find(s.begin(), s.end(), ind) != s.end()) s.push_back(ind);
        }
    }
  }

  mesh_structure::ind_pt_face_ct
     mesh_structure::ind_points_of_face_of_convex(size_type ic,
                                                  short_type iff) const {
    const mesh_convex_structure *q = &(convex_tab[ic]);
    std::vector<size_type>::const_iterator r = q->pts.begin();
    const convex_ind_ct *p = &(q->cstruct->ind_points_of_face(iff));
    return ind_pt_face_ct(r, p->begin(), p->end());
  }

  size_type mesh_structure::memsize(void) const {
    size_type mems =  sizeof(mesh_structure) + points_tab.memsize()
      + convex_tab.memsize();
    for (size_type i = 0; i < convex_tab.size(); ++i)
      mems += convex_tab[i].pts.size() * sizeof(size_type);
    for (size_type i = 0; i < points_tab.size(); ++i)
      mems += points_tab[i].size() * sizeof(size_type);
    return mems;
  }

  void mesh_structure::optimize_structure() {
    size_type i, j = nb_convex();
    for (i = 0; i < j; i++)
      if (!convex_tab.index_valid(i))
        swap_convex(i, convex_tab.ind_last());

    if (points_tab.size())
      for (i=0, j = (points_tab.end()-points_tab.begin())-1; i < j; ++i, --j){
        while (i < j && !(points_tab[i].empty())) ++i;
        while (i < j && points_tab[j].empty()) --j;
        if (i < j) swap_points(i, j);
      }
  }

  void mesh_structure::clear(void) {
    points_tab = dal::dynamic_tas<ind_cv_ct, 8>();
    convex_tab = dal::dynamic_tas<mesh_convex_structure, 8>();

  }

  void mesh_structure::stat(void) {
    cout << "mesh structure with " << nb_convex() << " valid convex, "
         << "for a total memory size of " << memsize() << " bytes.\n";
  }


  void mesh_structure::neighbours_of_convex(size_type ic, short_type iff,
                                            ind_set &s) const {
    s.resize(0);
    ind_pt_face_ct pt = ind_points_of_face_of_convex(ic, iff);

    for (size_type i = 0; i < points_tab[pt[0]].size(); ++i) {
      size_type icv = points_tab[pt[0]][i];
      if (icv != ic && is_convex_having_points(icv, short_type(pt.size()),
                                               pt.begin())
          && (convex_tab[ic].cstruct->dim()==convex_tab[icv].cstruct->dim()))
        s.push_back(icv);
    }
  }


  void mesh_structure::neighbours_of_convex(size_type ic,
					    const std::vector<short_type> &ftab,
					    ind_set &s) const {
    s.resize(0);
    size_type nb = nb_points_of_convex(ic);
    const mesh_convex_structure &q = convex_tab[ic];
    std::vector<size_type> cpt(nb, size_type(0)), ipt(nb);
    for (short_type iff : ftab)
      for (short_type i : q.cstruct->ind_points_of_face(iff))
	cpt[i]++;
    ipt.resize(0);
    for (size_type i = 0; i < nb; ++i)
      if (cpt[i] == ftab.size()) ipt.push_back(q.pts[i]);

    for (size_type i = 0; i < points_tab[ipt[0]].size(); ++i) {
      size_type icv = points_tab[ipt[0]][i];
      if (icv != ic && is_convex_having_points(icv, short_type(ipt.size()),
                                               ipt.begin())
          && (convex_tab[ic].cstruct->dim()==convex_tab[icv].cstruct->dim()))
        s.push_back(icv);
    }
  }

  void mesh_structure::neighbours_of_convex(size_type ic, ind_set &s) const {
    s.resize(0);
    unsigned nbf = nb_faces_of_convex(ic);
    for (short_type iff = 0; iff < nbf; ++iff) {
      ind_pt_face_ct pt = ind_points_of_face_of_convex(ic, iff);

      for (size_type i = 0; i < points_tab[pt[0]].size(); ++i) {
        size_type icv = points_tab[pt[0]][i];
        if (icv != ic && is_convex_having_points(icv, short_type(pt.size()),
                                                 pt.begin())
            && (convex_tab[ic].cstruct->dim()==convex_tab[icv].cstruct->dim()))
          if (std::find(s.begin(), s.end(), icv) == s.end())
            s.push_back(icv);
      }
    }
  }

  size_type mesh_structure::neighbour_of_convex(size_type ic,
                                                short_type iff) const {
    ind_pt_face_ct pt = ind_points_of_face_of_convex(ic, iff);

    for (size_type i = 0; i < points_tab[pt[0]].size(); ++i) {
      size_type icv = points_tab[pt[0]][i];
      if (icv != ic && is_convex_having_points(icv, short_type(pt.size()),
                                               pt.begin())
          && (convex_tab[ic].cstruct->dim()==convex_tab[icv].cstruct->dim()))
        return icv;
    }
    return size_type(-1);
  }

  convex_face mesh_structure::adjacent_face(size_type cv, short_type f) const {
    size_type neighbour_element = neighbour_of_convex(cv, f);
    if (neighbour_element == size_type(-1)) return convex_face::invalid_face();
    auto pcs = structure_of_convex(neighbour_element);
    auto face_points = ind_points_of_face_of_convex(cv, f);
    auto nNeighbourElementFaces = pcs->nb_faces();
    for (short_type iff = 0; iff < nNeighbourElementFaces; ++iff) {
      auto nPointsOnFace = pcs->nb_points_of_face(iff);
      if (is_convex_face_having_points(neighbour_element, iff,
				       nPointsOnFace, face_points.begin()))
        return {neighbour_element, iff};
    }
    GMM_ASSERT2(false, "failed to determine neighbouring face");
    return convex_face::invalid_face();
  }

  size_type mesh_structure::ind_in_convex_of_point(size_type ic,
                                                   size_type ip) const {
    const ind_cv_ct &ct = ind_points_of_convex(ic);
    ind_cv_ct::const_iterator it = std::find(ct.begin(), ct.end(), ip);
    return (it != ct.end()) ? (it - ct.begin()): size_type(-1);
  }

  /* ********************************************************************* */
  /*                                                                       */
  /*  Gives the list of edges of a convex face.                            */
  /*                                                                       */
  /* ********************************************************************* */

  void mesh_edge_list_convex(pconvex_structure cvs, std::vector<size_type> points_of_convex,
                             size_type cv_id, edge_list &el, bool merge_convex)
  { // a tester ... optimisable.
    dim_type n = cvs->dim();
    short_type nbp = cvs->nb_points();
    size_type ncv = merge_convex ? 0 : cv_id;

    if (nbp == n+1 && cvs == simplex_structure(n)) {
      for (dim_type k = 0; k < n; ++k)
        for (dim_type l = dim_type(k+1); l <= n; ++l)
          el.add(edge_list_elt(points_of_convex[k],
                               points_of_convex[l], ncv));
    }
    else if (nbp == (size_type(1) << n)
             && cvs == parallelepiped_structure(n)) {
      for (size_type k = 0; k < (size_type(1) << n); ++k)
        for (dim_type j = 0; j < n; ++j)
          if ((k & (1 << j)) == 0)
            el.add(edge_list_elt(points_of_convex[k],
                                 points_of_convex[k | (1 << j)], ncv));
    }
    else if (nbp == 2 * n && cvs == prism_structure(n)) {
      for (dim_type k = 0; k < n - 1; ++k)
        for (dim_type l = dim_type(k+1); l < n; ++l) {
          el.add(edge_list_elt(points_of_convex[k],
                               points_of_convex[l], ncv));
          el.add(edge_list_elt(points_of_convex[k+n],
                               points_of_convex[l+n], ncv));
        }
      for (dim_type k = 0; k < n; ++k)
        el.add(edge_list_elt(points_of_convex[k],
                             points_of_convex[k+n], ncv));
    }
    else {
      dal::dynamic_array<pconvex_structure> cvstab;
      dal::dynamic_array< std::vector<size_type> > indpttab;
      size_type ncs = 1;
      cvstab[0] = cvs;
      indpttab[0].resize(cvstab[0]->nb_points());
      std::copy(points_of_convex.begin(),
                points_of_convex.end(), indpttab[0].begin());

      /* pseudo recursive decomposition of the initial convex into
         face of face of ... of face of convex , cvstab and indpttab
         being the "stack" of convexes to treat
      */
      while (ncs != 0) {
        ncs--;
        cvs = cvstab[ncs];
        std::vector< size_type > pts = indpttab[ncs];
        if (cvs->dim() == 1) { // il faudrait étendre aux autres cas classiques.

          for (size_type j = 1; j < cvs->nb_points(); ++j) {
            //cerr << "ncs=" << ncs << "j=" << j << ", ajout de " << (indpttab[ncs])[j-1] << "," << (indpttab[ncs])[j] << endl;
            el.add(edge_list_elt((indpttab[ncs])[j-1],(indpttab[ncs])[j],ncv));
          }
        }
        else {
          short_type nf = cvs->nb_faces();
          //cerr << "ncs = " << ncs << ",cvs->dim=" << int(cvs->dim()) << ", nb_faces=" << nf << endl;
          for (short_type f = 0; f < nf; ++f) {
            cvstab[ncs+f] = (cvs->faces_structure())[f];
            indpttab[ncs+f].resize(cvs->nb_points_of_face(f));
            /*
            cerr << "   -> f=" << f << ", cvs->nb_points_of_face(f)=" <<
              cvs->nb_points_of_face(f) << "=[";
            for (size_type k = 0; k < cvs->nb_points_of_face(f); ++k)
              cerr << (std::string(k==0?"":","))
                   << int(cvs->ind_points_of_face(f)[k]) << "{" << pts[cvs->ind_points_of_face(f)[k]] << "}";
            cerr << "]" << endl;
            */
            for (size_type k = 0; k < cvs->nb_points_of_face(f); ++k)
              (indpttab[ncs+f])[k] = pts[cvs->ind_points_of_face(f)[k]];
          }
          //          cvstab[ncs] = cvstab[ncs + nf - 1];
          //indpttab[ncs] = indpttab[ncs + nf - 1];
          ncs += nf;
          //cerr << "on empile les " << nf << " faces, -> ncs=" << ncs << endl;
        }
      }
    }
  }


  void mesh_edge_list(const mesh_structure &m, edge_list &el,
                      bool merge_convex) {
    std::vector<size_type> p;
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv) {
      p.resize(m.nb_points_of_convex(cv));
      std::copy(m.ind_points_of_convex(cv).begin(),
                m.ind_points_of_convex(cv).end(), p.begin());
      mesh_edge_list_convex(m.structure_of_convex(cv), p, cv,
                            el, merge_convex);
    }
  }


  // static double enumerate_dof_time = 0;
  void cuthill_mckee_on_convexes(const bgeot::mesh_structure &ms,
                                 std::vector<size_type> &cmk) {
    const dal::bit_vector &cvlst = ms.convex_index();
    cmk.reserve(cvlst.card()); cmk.resize(0);
    if (ms.convex_index().card() == 0) return;
    std::deque<int> pile;
    std::vector<size_type> tab, connectivity(cvlst.last_true()+1),
      temp(cvlst.last_true()+1);

    size_type cv = cvlst.first_true();

    std::fill(connectivity.begin(), connectivity.end(), size_type(-1));
    // double t = dal::uclock_sec();

    /* count neighbours for each convex */
    for (dal::bv_visitor j(cvlst); !j.finished(); ++j) {
      const mesh_structure::ind_cv_ct &ct = ms.ind_points_of_convex(j);
      mesh_structure::ind_cv_ct::const_iterator itp = ct.begin(),
        itpe = ct.end();
      size_type nei = 0;

      for (; itp != itpe; ++itp) {
        size_type ip = *itp;
        mesh_structure::ind_cv_ct::const_iterator
          it = ms.convex_to_point(ip).begin(),
          ite =  ms.convex_to_point(ip).end();
        for ( ; it != ite; ++it)
          if (temp[*it] != j+1) { temp[*it] = j+1; nei++; }
      }
      connectivity[j] = nei-1;
      if (nei < connectivity[cv]) cv = j;
    }

    /* do the cuthill mckee */
    while (cv != size_type(-1)) {
      connectivity[cv] = size_type(-1);
      cmk.push_back(cv);
      size_type nbp = ms.nb_points_of_convex(cv);

      for (size_type i = 0; i < nbp; i++) {
        size_type ip = ms.ind_points_of_convex(cv)[i];
        mesh_structure::ind_cv_ct::const_iterator
          it = ms.convex_to_point(ip).begin(),
          ite =  ms.convex_to_point(ip).end();
        for ( ; it != ite; ++it)
          if (connectivity[*it] != size_type(-1)) {
            connectivity[*it] = size_type(-1);
            pile.push_front(int(*it));
          }
      }

      if (pile.empty()) {
        cv = std::min_element(connectivity.begin(), connectivity.end()) - connectivity.begin();
        if (connectivity[cv] == size_type(-1)) break;
      }
      else { cv = pile.back(); pile.pop_back(); }
    }

    // enumerate_dof_time += dal::uclock_sec() - t;
    // cerr << "cuthill_mckee_on_convexes: " << enumerate_dof_time << " sec\n";
  }

}  /* end of namespace bgeot.                                              */
