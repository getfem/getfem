/*===========================================================================

 Copyright (C) 2001-2020 Yves Renard

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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
#include "getfem/bgeot_convex_ref.h"
#include "getfem/bgeot_mesh_structure.h"
#include "getfem/bgeot_comma_init.h"

namespace bgeot {
  
  // ******************************************************************
  //    Interface with qhull
  // ******************************************************************
  
# if !defined(GETFEM_HAVE_LIBQHULL_R_QHULL_RA_H)
  void qhull_delaunay(const std::vector<base_node> &,
                gmm::dense_matrix<size_type>&) {
    GMM_ASSERT1(false, "Qhull header files not installed. "
                "Install qhull library and reinstall GetFEM library.");
  }
# else

  extern "C" { // old versions of this header lack a __cplusplus guard
# include <libqhull_r/qhull_ra.h>
  }
  
  void qhull_delaunay(const std::vector<base_node> &pts,
		      gmm::dense_matrix<size_type>& simplexes) {
    // cout << "running delaunay with " << pts.size() << " points\n";
    size_type dim = pts[0].size();   /* points dimension.           */
    if (pts.size() <= dim) { gmm::resize(simplexes, dim+1, 0); return; }
    if (pts.size() == dim+1) {
      gmm::resize(simplexes, dim+1, 1);
      for (size_type i=0; i <= dim; ++i) simplexes(i, 0) = i;
      return;
    }
    std::vector<coordT> Pts(dim * pts.size());
    for (size_type i=0; i < pts.size(); ++i)
      gmm::copy(pts[i], gmm::sub_vector(Pts, gmm::sub_interval(i*dim, dim)));
    boolT ismalloc=0;  /* True if qhull should free points in
                        * qh_freeqhull() or reallocation      */
    /* Be Aware: option QJ could destabilizate all, it can break everything. */
    /* option Qbb -> QbB (????) */
    /* option flags for qhull, see qh_opt.htm */
    char flags[]= "qhull QJ d Qbb Pp T0"; //QJ s i TO";//"qhull Tv";
    FILE *outfile= 0;    /* output from qh_produce_output()
                          *  use NULL to skip qh_produce_output() */
    FILE *errfile= stderr;    /* error messages from qhull code */
    int exitcode;             /* 0 if no error from qhull */
    facetT *facet;                  /* set by FORALLfacets */
    int curlong, totlong;          /* memory remaining after qh_memfreeshort */
    vertexT *vertex, **vertexp;
    qhT context = {};
    qhT* qh = &context;
    exitcode = qh_new_qhull (qh, int(dim), int(pts.size()), &Pts[0], ismalloc,
                             flags, outfile, errfile);
    if (!exitcode) { /* if no error */
      size_type nbf=0;
      FORALLfacets { if (!facet->upperdelaunay) nbf++; }
      gmm::resize(simplexes, dim+1, nbf);
        /* 'qh facet_list' contains the convex hull */
      nbf=0;
      FORALLfacets {
        if (!facet->upperdelaunay) {
          size_type s=0;
          FOREACHvertex_(facet->vertices) {
            assert(s < (unsigned)(dim+1));
            simplexes(s++,nbf) = qh_pointid(qh, vertex->point);
          }
          nbf++;
        }
      }
    }
    qh_freeqhull(qh, !qh_ALL);
    qh_memfreeshort(qh, &curlong, &totlong);
    if (curlong || totlong)
      cerr << "qhull internal warning (main): did not free " << totlong <<
        " bytes of long memory (" << curlong << " pieces)\n";

  }

#endif

  

  size_type simplexified_tab(pconvex_structure cvs, size_type **tab);

  static void simplexify_convex(bgeot::convex_of_reference *cvr,
                                mesh_structure &m) {
    pconvex_structure cvs = cvr->structure();
    m.clear();
    auto basic_cvs = basic_structure(cvs);
    dim_type n = basic_cvs->dim();
    std::vector<size_type> ipts(n+1);
    if (basic_cvs->nb_points() == n + 1) {
      for (size_type i = 0; i <= n; ++i) ipts[i] = i;
      m.add_simplex(n, ipts.begin());
    }
    else {
      size_type *tab;
      size_type nb = simplexified_tab(basic_cvs, &tab);
      if (nb) {
	for (size_type nc = 0; nc < nb; ++nc) {
	  for (size_type i = 0; i <= n; ++i) ipts[i] = *tab++;
	  m.add_simplex(n, ipts.begin());
	}
      }	else {
#       if defined(GETFEM_HAVE_LIBQHULL_R_QHULL_RA_H)
	gmm::dense_matrix<size_type> t;
	qhull_delaunay(cvr->points(), t);
	for (size_type nc = 0; nc < gmm::mat_ncols(t); ++nc) {
	  for (size_type i = 0; i <= n; ++i) ipts[i] = t(i,nc);
	  m.add_simplex(n, ipts.begin());
	}
#       endif
      }
    }
  }

  /* ********************************************************************* */
  /*       Point tab storage.                                              */
  /* ********************************************************************* */

  struct stored_point_tab_key : virtual public dal::static_stored_object_key  {
    const stored_point_tab *pspt;
    bool compare(const static_stored_object_key &oo) const override {
      const stored_point_tab_key &o
        = dynamic_cast<const stored_point_tab_key &>(oo);
      const stored_point_tab &x = *pspt;
      const stored_point_tab &y = *(o.pspt);
      if (&x == &y) return false;
      std::vector<base_node>::const_iterator it1 = x.begin(), it2 = y.begin();
      base_node::const_iterator itn1, itn2, itne;
      for ( ; it1 != x.end() && it2 != y.end() ; ++it1, ++it2) {
        if ((*it1).size() < (*it2).size()) return true;
        if ((*it1).size() > (*it2).size()) return false;
        itn1 = (*it1).begin(); itne = (*it1).end(); itn2 = (*it2).begin();
        for ( ; itn1 != itne ; ++itn1, ++itn2)
          if (*itn1 < *itn2) return true;
          else if (*itn1 > *itn2) return false;
      }
      if (it2 != y.end()) return true;
      return false;
    }
    bool equal(const static_stored_object_key &oo) const override {
      auto &o = dynamic_cast<const stored_point_tab_key &>(oo);
      auto &x = *pspt;
      auto &y = *(o.pspt);
      if (&x == &y) return true;
      if (x.size() != y.size()) return false;
      auto it1 = x.begin();
      auto it2 = y.begin();
      base_node::const_iterator itn1, itn2, itne;
      for ( ; it1 != x.end() && it2 != y.end() ; ++it1, ++it2) {
        if ((*it1).size() != (*it2).size()) return false;
        itn1 = (*it1).begin(); itne = (*it1).end(); itn2 = (*it2).begin();
        for ( ; itn1 != itne ; ++itn1, ++itn2)
          if (*itn1 != *itn2) return false;
      }
      return true;
    }
    stored_point_tab_key(const stored_point_tab *p) : pspt(p) {}
  };


  pstored_point_tab store_point_tab(const stored_point_tab &spt) {
    dal::pstatic_stored_object_key
      pk = std::make_shared<stored_point_tab_key>(&spt);
    dal::pstatic_stored_object o = dal::search_stored_object(pk);
    if (o) return std::dynamic_pointer_cast<const stored_point_tab>(o);
    pstored_point_tab p = std::make_shared<stored_point_tab>(spt);
    DAL_STORED_OBJECT_DEBUG_CREATED(p.get(), "Stored point tab");
    dal::pstatic_stored_object_key
      psp = std::make_shared<stored_point_tab_key>(p.get());
    dal::add_stored_object(psp, p, dal::AUTODELETE_STATIC_OBJECT);
    return p;
  }

  convex_of_reference::convex_of_reference
  (pconvex_structure cvs_, bool auto_basic_) :
    convex<base_node>(move(cvs_)), basic_convex_ref_(0),
    auto_basic(auto_basic_) {
    DAL_STORED_OBJECT_DEBUG_CREATED(this, "convex of refrence");
    psimplexified_convex = std::make_shared<mesh_structure>();
    // dal::singleton<cleanup_simplexified_convexes>::instance()
    //        .push_back(psimplexified_convex);
  }

  /* should be called on the basic_convex_ref */
  const mesh_structure* convex_of_reference::simplexified_convex() const {
    GMM_ASSERT1(auto_basic,
                "always use simplexified_convex on the basic_convex_ref() "
                "[this=" << nb_points() << ", basic="
                << basic_convex_ref_->nb_points());
    return psimplexified_convex.get();
  }

  // Key type for static storing
  class convex_of_reference_key : virtual public dal::static_stored_object_key{
    int type; // 0 = simplex structure degree K
              // 1 = equilateral simplex of ref.
              // 2 = dummy
    dim_type N; short_type K; short_type nf;
  public :
    bool compare(const static_stored_object_key &oo) const override{
      const convex_of_reference_key &o
        = dynamic_cast<const convex_of_reference_key &>(oo);
      if (type < o.type) return true;
      if (type > o.type) return false;
      if (N < o.N) return true;
      if (N > o.N) return false;
      if (K < o.K) return true;
      if (K > o.K) return false;
      if (nf < o.nf) return true;
      return false;
    }
    bool equal(const static_stored_object_key &oo) const override{
      auto &o = dynamic_cast<const convex_of_reference_key &>(oo);
      if (type != o.type) return false;
      if (N != o.N) return false;
      if (K != o.K) return false;
      if (nf != o.nf) return false;
      return true;
    }
    convex_of_reference_key(int t, dim_type NN, short_type KK = 0,
                            short_type nnf = 0)
      : type(t), N(NN), K(KK), nf(nnf) {}
  };


  /* simplexes.                                                            */

  class K_simplex_of_ref_ : public convex_of_reference {
  public :
    scalar_type is_in(const base_node &pt) const {
      // return a negative number if pt is in the convex
      GMM_ASSERT1(pt.size() == cvs->dim(),
                  "K_simplex_of_ref_::is_in: Dimensions mismatch");
      scalar_type e = -1.0, r = (pt.size() > 0) ? -pt[0] : 0.0;
      base_node::const_iterator it = pt.begin(), ite = pt.end();
      for (; it != ite; e += *it, ++it) r = std::max(r, -(*it));
      return std::max(r, e);
    }

    scalar_type is_in_face(short_type f, const base_node &pt) const {
      // return zero if pt is in the face of the convex
      // negative if the point is on the side of the face where the element is
      GMM_ASSERT1(pt.size() == cvs->dim(),
                  "K_simplex_of_ref_::is_in_face: Dimensions mismatch");
      if (f > 0) return -pt[f-1];
      scalar_type e = -1.0;
      base_node::const_iterator it = pt.begin(), ite = pt.end();
      for (; it != ite; e += *it, ++it) {};
      return e / sqrt(scalar_type(pt.size()));
    }

    void project_into(base_node &pt) const {
      if (auto_basic) {
        GMM_ASSERT1(pt.size() == cvs->dim(),
                    "K_simplex_of_ref_::project_into: Dimensions mismatch");
        scalar_type sum_coordinates = 0.0;
        for (const auto &coord : pt) sum_coordinates += coord;
        if (sum_coordinates > 1.0) gmm::scale(pt, 1.0 / sum_coordinates);
        for (auto &coord : pt) {
          if (coord < 0.0) coord = 0.0;
          if (coord > 1.0) coord = 1.0;
        }
      } else
        basic_convex_ref_->project_into(pt);
    }

    K_simplex_of_ref_(dim_type NN, short_type KK) :
      convex_of_reference(simplex_structure(NN, KK), (KK == 1) || (NN == 0))
    {
      if ((KK != 1) && (NN != 0)) basic_convex_ref_ = simplex_of_reference(NN, 1);
      size_type R = cvs->nb_points();
      convex<base_node>::points().resize(R);
      normals_.resize(NN+1);
      base_node null(NN); null.fill(0.0);
      std::fill(normals_.begin(), normals_.end(), null);
      std::fill(convex<base_node>::points().begin(),
                convex<base_node>::points().end(), null);
      for (size_type i = 1; i <= NN; ++i)
        normals_[i][i-1] = -1.0;
      if (NN > 0)
        std::fill(normals_[0].begin(), normals_[0].end(),
                  scalar_type(1.0)/sqrt(scalar_type(NN)));
      base_node c(NN);  c.fill(0.0);

      if (KK == 0) {
        c.fill(1.0/(NN+1));
        convex<base_node>::points()[0] = c;
      }
      else {
        size_type sum = 0, l;
        for (size_type r = 0; r < R; ++r) {
          convex<base_node>::points()[r] = c;
          if (KK != 0 && NN > 0) {
            l = 0; c[l] += 1.0 / scalar_type(KK); sum++;
            while (sum > KK) {
              sum -= int(floor(0.5+(c[l] * KK)));
              c[l] = 0.0; l++; if (l == NN) break;
              c[l] += 1.0 / scalar_type(KK); sum++;
            }
          }
        }
      }
      ppoints = store_point_tab(convex<base_node>::points());
      if (auto_basic) simplexify_convex(this, *psimplexified_convex);
    }
  };

  pconvex_ref simplex_of_reference(dim_type nc, short_type K) {
    dal::pstatic_stored_object_key
      pk = std::make_shared<convex_of_reference_key>(0, nc, K);
    dal::pstatic_stored_object o = dal::search_stored_object(pk);
    if (o) return std::dynamic_pointer_cast<const convex_of_reference>(o);
    pconvex_ref p = std::make_shared<K_simplex_of_ref_>(nc, K);
    dal::add_stored_object(pk, p, p->structure(), p->pspt(),
                           dal::PERMANENT_STATIC_OBJECT);
    pconvex_ref p1 = basic_convex_ref(p);
    if (p != p1) add_dependency(p, p1);
    return p;
  }

  /* ******************************************************************** */
  /*    Incomplete Q2 quadrilateral or hexahedral of reference.           */
  /* ******************************************************************** */
  /* By Yao Koutsawa  <yao.koutsawa@tudor.lu> 2012-12-10                  */

  class Q2_incomplete_of_ref_ : public convex_of_reference {
  public :
    scalar_type is_in(const base_node& pt) const
    { return basic_convex_ref_->is_in(pt); }
    scalar_type is_in_face(short_type f, const base_node& pt) const
    { return basic_convex_ref_->is_in_face(f, pt); }

    Q2_incomplete_of_ref_(dim_type nc) :
      convex_of_reference(Q2_incomplete_structure(nc), false)
    {
      GMM_ASSERT1(nc == 2 || nc == 3, "Sorry exist only in dimension 2 or 3");
      convex<base_node>::points().resize(cvs->nb_points());
      normals_.resize(nc == 2 ? 4: 6);
      basic_convex_ref_ = parallelepiped_of_reference(nc);

      if(nc==2) {
        normals_[0] = { 1, 0};
        normals_[1] = {-1, 0};
        normals_[2] = { 0, 1};
        normals_[3] = { 0,-1};

        convex<base_node>::points()[0] = base_node(0.0, 0.0);
        convex<base_node>::points()[1] = base_node(0.5, 0.0);
        convex<base_node>::points()[2] = base_node(1.0, 0.0);
        convex<base_node>::points()[3] = base_node(0.0, 0.5);
        convex<base_node>::points()[4] = base_node(1.0, 0.5);
        convex<base_node>::points()[5] = base_node(0.0, 1.0);
        convex<base_node>::points()[6] = base_node(0.5, 1.0);
        convex<base_node>::points()[7] = base_node(1.0, 1.0);

      } else {
        normals_[0] = { 1, 0, 0};
        normals_[1] = {-1, 0, 0};
        normals_[2] = { 0, 1, 0};
        normals_[3] = { 0,-1, 0};
        normals_[4] = { 0, 0, 1};
        normals_[5] = { 0, 0,-1};

        convex<base_node>::points()[0] = base_node(0.0, 0.0, 0.0);
        convex<base_node>::points()[1] = base_node(0.5, 0.0, 0.0);
        convex<base_node>::points()[2] = base_node(1.0, 0.0, 0.0);
        convex<base_node>::points()[3] = base_node(0.0, 0.5, 0.0);
        convex<base_node>::points()[4] = base_node(1.0, 0.5, 0.0);
        convex<base_node>::points()[5] = base_node(0.0, 1.0, 0.0);
        convex<base_node>::points()[6] = base_node(0.5, 1.0, 0.0);
        convex<base_node>::points()[7] = base_node(1.0, 1.0, 0.0);

        convex<base_node>::points()[8] = base_node(0.0, 0.0, 0.5);
        convex<base_node>::points()[9] = base_node(1.0, 0.0, 0.5);
        convex<base_node>::points()[10] = base_node(0.0, 1.0, 0.5);
        convex<base_node>::points()[11] = base_node(1.0, 1.0, 0.5);

        convex<base_node>::points()[12] = base_node(0.0, 0.0, 1.0);
        convex<base_node>::points()[13] = base_node(0.5, 0.0, 1.0);
        convex<base_node>::points()[14] = base_node(1.0, 0.0, 1.0);
        convex<base_node>::points()[15] = base_node(0.0, 0.5, 1.0);
        convex<base_node>::points()[16] = base_node(1.0, 0.5, 1.0);
        convex<base_node>::points()[17] = base_node(0.0, 1.0, 1.0);
        convex<base_node>::points()[18] = base_node(0.5, 1.0, 1.0);
        convex<base_node>::points()[19] = base_node(1.0, 1.0, 1.0);
      }
      ppoints = store_point_tab(convex<base_node>::points());
    }
  };


  DAL_SIMPLE_KEY(Q2_incomplete_of_reference_key_, dim_type);

  pconvex_ref Q2_incomplete_of_reference(dim_type nc) {
     dal::pstatic_stored_object_key
      pk = std::make_shared<Q2_incomplete_of_reference_key_>(nc);
    dal::pstatic_stored_object o = dal::search_stored_object(pk);
    if (o) return std::dynamic_pointer_cast<const convex_of_reference>(o);
    pconvex_ref p = std::make_shared<Q2_incomplete_of_ref_>(nc);
    dal::add_stored_object(pk, p, p->structure(), p->pspt(),
                           dal::PERMANENT_STATIC_OBJECT);
     pconvex_ref p1 = basic_convex_ref(p);
    if (p != p1) add_dependency(p, p1);
    return p;
  }

  /* ******************************************************************** */
  /*    Pyramidal element of reference.                                   */
  /* ******************************************************************** */

  class pyramid_QK_of_ref_ : public convex_of_reference {
  public :
    scalar_type is_in_face(short_type f, const base_node& pt) const {
      // return zero if pt is in the face of the convex
      // negative if the point is on the side of the face where the element is
      GMM_ASSERT1(pt.size() == 3, "Dimensions mismatch");
      if (f == 0)
        return -pt[2];
      else
        return gmm::vect_sp(normals_[f], pt) - sqrt(2.)/2.;
    }

    scalar_type is_in(const base_node& pt) const {
      // return a negative number if pt is in the convex
      scalar_type r = is_in_face(0, pt);
      for (short_type f = 1; f < 5; ++f) r = std::max(r, is_in_face(f, pt));
      return r;
    }

    void project_into(base_node &pt) const {
      if (auto_basic) {
        GMM_ASSERT1(pt.size() == 3, "Dimensions mismatch");
        if (pt[2] < .0) pt[2] = 0.;
        for (short_type f = 1; f < 5; ++f) {
          scalar_type reldist = gmm::vect_sp(normals_[f], pt)*sqrt(2.);
          if (reldist > 1.)
            gmm::scale(pt, 1./reldist);
        }
      } else
        basic_convex_ref_->project_into(pt);
    }

    pyramid_QK_of_ref_(dim_type k) : convex_of_reference(pyramid_QK_structure(k), k == 1) {
      GMM_ASSERT1(k == 1 || k == 2,
                  "Sorry exist only in degree 1 or 2, not " << k);

      convex<base_node>::points().resize(cvs->nb_points());
      normals_.resize(cvs->nb_faces());
      if (k != 1) basic_convex_ref_ = pyramid_QK_of_reference(1);

      normals_[0] = { 0., 0., -1.};
      normals_[1] = { 0.,-1.,  1.};
      normals_[2] = { 1., 0.,  1.};
      normals_[3] = { 0., 1.,  1.};
      normals_[4] = {-1., 0.,  1.};

      for (size_type i = 0; i < normals_.size(); ++i)
        gmm::scale(normals_[i], 1. / gmm::vect_norm2(normals_[i]));

      if (k==1) {
        convex<base_node>::points()[0]  = base_node(-1.0, -1.0, 0.0);
        convex<base_node>::points()[1]  = base_node( 1.0, -1.0, 0.0);
        convex<base_node>::points()[2]  = base_node(-1.0,  1.0, 0.0);
        convex<base_node>::points()[3]  = base_node( 1.0,  1.0, 0.0);
        convex<base_node>::points()[4]  = base_node( 0.0,  0.0, 1.0);
      } else if (k==2) {
        convex<base_node>::points()[0]  = base_node(-1.0, -1.0, 0.0);
        convex<base_node>::points()[1]  = base_node( 0.0, -1.0, 0.0);
        convex<base_node>::points()[2]  = base_node( 1.0, -1.0, 0.0);
        convex<base_node>::points()[3]  = base_node(-1.0,  0.0, 0.0);
        convex<base_node>::points()[4]  = base_node( 0.0,  0.0, 0.0);
        convex<base_node>::points()[5]  = base_node( 1.0,  0.0, 0.0);
        convex<base_node>::points()[6]  = base_node(-1.0,  1.0, 0.0);
        convex<base_node>::points()[7]  = base_node( 0.0,  1.0, 0.0);
        convex<base_node>::points()[8]  = base_node( 1.0,  1.0, 0.0);
        convex<base_node>::points()[9]  = base_node(-0.5, -0.5, 0.5);
        convex<base_node>::points()[10] = base_node( 0.5, -0.5, 0.5);
        convex<base_node>::points()[11] = base_node(-0.5,  0.5, 0.5);
        convex<base_node>::points()[12] = base_node( 0.5,  0.5, 0.5);
        convex<base_node>::points()[13] = base_node( 0.0,  0.0, 1.0);
      }
      ppoints = store_point_tab(convex<base_node>::points());
      if (auto_basic) simplexify_convex(this, *psimplexified_convex);
    }
  };


  DAL_SIMPLE_KEY(pyramid_QK_reference_key_, dim_type);

  pconvex_ref pyramid_QK_of_reference(dim_type k) {
     dal::pstatic_stored_object_key
      pk = std::make_shared<pyramid_QK_reference_key_>(k);
    dal::pstatic_stored_object o = dal::search_stored_object(pk);
    if (o) return std::dynamic_pointer_cast<const convex_of_reference>(o);
    pconvex_ref p = std::make_shared<pyramid_QK_of_ref_>(k);
    dal::add_stored_object(pk, p, p->structure(), p->pspt(),
                           dal::PERMANENT_STATIC_OBJECT);
    pconvex_ref p1 = basic_convex_ref(p);
    if (p != p1) add_dependency(p, p1);
    return p;
  }


  /* ******************************************************************** */
  /*    Incomplete quadratic pyramidal element of reference.              */
  /* ******************************************************************** */

  class pyramid_Q2_incomplete_of_ref_ : public convex_of_reference {
  public :
    scalar_type is_in(const base_node& pt) const
    { return basic_convex_ref_->is_in(pt); }
    scalar_type is_in_face(short_type f, const base_node& pt) const
    { return basic_convex_ref_->is_in_face(f, pt); }

    pyramid_Q2_incomplete_of_ref_() : convex_of_reference(pyramid_Q2_incomplete_structure(), false) {
      convex<base_node>::points().resize(cvs->nb_points());
      normals_.resize(cvs->nb_faces());
      basic_convex_ref_ = pyramid_QK_of_reference(1);

      normals_ = basic_convex_ref_->normals();

      convex<base_node>::points()[0]  = base_node(-1.0, -1.0, 0.0);
      convex<base_node>::points()[1]  = base_node( 0.0, -1.0, 0.0);
      convex<base_node>::points()[2]  = base_node( 1.0, -1.0, 0.0);
      convex<base_node>::points()[3]  = base_node(-1.0,  0.0, 0.0);
      convex<base_node>::points()[4]  = base_node( 1.0,  0.0, 0.0);
      convex<base_node>::points()[5]  = base_node(-1.0,  1.0, 0.0);
      convex<base_node>::points()[6]  = base_node( 0.0,  1.0, 0.0);
      convex<base_node>::points()[7]  = base_node( 1.0,  1.0, 0.0);
      convex<base_node>::points()[8]  = base_node(-0.5, -0.5, 0.5);
      convex<base_node>::points()[9]  = base_node( 0.5, -0.5, 0.5);
      convex<base_node>::points()[10] = base_node(-0.5,  0.5, 0.5);
      convex<base_node>::points()[11] = base_node( 0.5,  0.5, 0.5);
      convex<base_node>::points()[12] = base_node( 0.0,  0.0, 1.0);

      ppoints = store_point_tab(convex<base_node>::points());
      if (auto_basic) simplexify_convex(this, *psimplexified_convex);
    }
  };


  DAL_SIMPLE_KEY(pyramid_Q2_incomplete_reference_key_, dim_type);

  pconvex_ref pyramid_Q2_incomplete_of_reference() {
    dal::pstatic_stored_object_key
      pk = std::make_shared<pyramid_Q2_incomplete_reference_key_>(0);
    dal::pstatic_stored_object o = dal::search_stored_object(pk);
    if (o)
      return std::dynamic_pointer_cast<const convex_of_reference>(o);
    else {
      pconvex_ref p = std::make_shared<pyramid_Q2_incomplete_of_ref_>();
      dal::add_stored_object(pk, p, p->structure(), p->pspt(),
                             dal::PERMANENT_STATIC_OBJECT);
      pconvex_ref p1 = basic_convex_ref(p);
      if (p != p1) add_dependency(p, p1);
      return p;
    }
  }


  /* ******************************************************************** */
  /*    Incomplete quadratic triangular prism element of reference.       */
  /* ******************************************************************** */

  class prism_incomplete_P2_of_ref_ : public convex_of_reference {
  public :
    scalar_type is_in(const base_node& pt) const
    { return basic_convex_ref_->is_in(pt); }
    scalar_type is_in_face(short_type f, const base_node& pt) const
    { return basic_convex_ref_->is_in_face(f, pt); }

    prism_incomplete_P2_of_ref_() : convex_of_reference(prism_incomplete_P2_structure(), false) {
      convex<base_node>::points().resize(cvs->nb_points());
      normals_.resize(cvs->nb_faces());
      basic_convex_ref_ = prism_of_reference(3);

      normals_ = basic_convex_ref_->normals();

      convex<base_node>::points()[0]  = base_node(0.0, 0.0, 0.0);
      convex<base_node>::points()[1]  = base_node(0.5, 0.0, 0.0);
      convex<base_node>::points()[2]  = base_node(1.0, 0.0, 0.0);
      convex<base_node>::points()[3]  = base_node(0.0, 0.5, 0.0);
      convex<base_node>::points()[4]  = base_node(0.5, 0.5, 0.0);
      convex<base_node>::points()[5]  = base_node(0.0, 1.0, 0.0);
      convex<base_node>::points()[6]  = base_node(0.0, 0.0, 0.5);
      convex<base_node>::points()[7]  = base_node(1.0, 0.0, 0.5);
      convex<base_node>::points()[8]  = base_node(0.0, 1.0, 0.5);
      convex<base_node>::points()[9]  = base_node(0.0, 0.0, 1.0);
      convex<base_node>::points()[10] = base_node(0.5, 0.0, 1.0);
      convex<base_node>::points()[11] = base_node(1.0, 0.0, 1.0);
      convex<base_node>::points()[12] = base_node(0.0, 0.5, 1.0);
      convex<base_node>::points()[13] = base_node(0.5, 0.5, 1.0);
      convex<base_node>::points()[14] = base_node(0.0, 1.0, 1.0);

      ppoints = store_point_tab(convex<base_node>::points());
      if (auto_basic) simplexify_convex(this, *psimplexified_convex);
    }
  };


  DAL_SIMPLE_KEY(prism_incomplete_P2_reference_key_, dim_type);

  pconvex_ref prism_incomplete_P2_of_reference() {
    dal::pstatic_stored_object_key
      pk = std::make_shared<prism_incomplete_P2_reference_key_>(0);
    dal::pstatic_stored_object o = dal::search_stored_object(pk);
    if (o)
      return std::dynamic_pointer_cast<const convex_of_reference>(o);
    else {
      pconvex_ref p = std::make_shared<prism_incomplete_P2_of_ref_>();
      dal::add_stored_object(pk, p, p->structure(), p->pspt(),
                             dal::PERMANENT_STATIC_OBJECT);
      pconvex_ref p1 = basic_convex_ref(p);
      if (p != p1) add_dependency(p, p1);
      return p;
    }
  }


  /* ******************************************************************** */
  /*    Products.                                                         */
  /* ******************************************************************** */

  DAL_DOUBLE_KEY(product_ref_key_, pconvex_ref, pconvex_ref);

  struct product_ref_ : public convex_of_reference {
    pconvex_ref cvr1, cvr2;

    scalar_type is_in(const base_node &pt) const {
      dim_type n1 = cvr1->structure()->dim(), n2 = cvr2->structure()->dim();
      base_node pt1(n1), pt2(n2);
      GMM_ASSERT1(pt.size() == cvs->dim(),
                  "product_ref_::is_in: Dimension does not match");
      std::copy(pt.begin(), pt.begin()+n1, pt1.begin());
      std::copy(pt.begin()+n1,   pt.end(), pt2.begin());
      return std::max(cvr1->is_in(pt1), cvr2->is_in(pt2));
    }

    scalar_type is_in_face(short_type f, const base_node &pt) const {
      // ne controle pas si le point est dans le convexe mais si un point
      // suppos\E9 appartenir au convexe est dans une face donn\E9e
      dim_type n1 = cvr1->structure()->dim(), n2 = cvr2->structure()->dim();
      base_node pt1(n1), pt2(n2);
      GMM_ASSERT1(pt.size() == cvs->dim(), "Dimensions mismatch");
      std::copy(pt.begin(), pt.begin()+n1, pt1.begin());
      std::copy(pt.begin()+n1,   pt.end(), pt2.begin());
      if (f < cvr1->structure()->nb_faces()) return cvr1->is_in_face(f, pt1);
      else return cvr2->is_in_face(short_type(f - cvr1->structure()->nb_faces()), pt2);
    }

    void project_into(base_node &pt) const {
      if (auto_basic) {
        GMM_ASSERT1(pt.size() == cvs->dim(), "Dimensions mismatch");
        dim_type n1 = cvr1->structure()->dim(), n2 = cvr2->structure()->dim();
        base_node pt1(n1), pt2(n2);
        std::copy(pt.begin(), pt.begin()+n1, pt1.begin());
        std::copy(pt.begin()+n1,   pt.end(), pt2.begin());
        cvr1->project_into(pt1);
        cvr2->project_into(pt2);
        std::copy(pt1.begin(), pt1.end(), pt.begin());
        std::copy(pt2.begin(), pt2.end(), pt.begin()+n1);
      } else
        basic_convex_ref_->project_into(pt);
    }

    product_ref_(pconvex_ref a, pconvex_ref b) :
      convex_of_reference(
        convex_direct_product(*a, *b).structure(),
        basic_convex_ref(a) == a && basic_convex_ref(b) == b)
    {
      if (a->structure()->dim() < b->structure()->dim())
        GMM_WARNING1("Illegal convex: swap your operands: dim(cv1)=" <<
                    int(a->structure()->dim()) << " < dim(cv2)=" <<
                    int(b->structure()->dim()));
      cvr1 = a; cvr2 = b;
      *((convex<base_node> *)(this)) = convex_direct_product(*a, *b);
      normals_.resize(cvs->nb_faces());
      base_small_vector null(cvs->dim()); null.fill(0.0);
      std::fill(normals_.begin(), normals_.end(), null);
      for (size_type r = 0; r < cvr1->structure()->nb_faces(); r++)
        std::copy(cvr1->normals()[r].begin(), cvr1->normals()[r].end(),
                  normals_[r].begin());
      for (size_type r = 0; r < cvr2->structure()->nb_faces(); r++)
        std::copy(cvr2->normals()[r].begin(), cvr2->normals()[r].end(),
                  normals_[r+cvr1->structure()->nb_faces()].begin()
                  + cvr1->structure()->dim());
      ppoints = store_point_tab(convex<base_node>::points());

      if (basic_convex_ref(a) != a || basic_convex_ref(b) != b)
        basic_convex_ref_ = convex_ref_product(basic_convex_ref(a),
                                               basic_convex_ref(b));
      if (auto_basic) simplexify_convex(this, *psimplexified_convex);
    }
  };


  pconvex_ref convex_ref_product(pconvex_ref a, pconvex_ref b) {
    dal::pstatic_stored_object_key
      pk = std::make_shared<product_ref_key_>(a, b);
    dal::pstatic_stored_object o = dal::search_stored_object(pk);
    if (o)
      return std::dynamic_pointer_cast<const convex_of_reference>(o);
    else {
      pconvex_ref p = std::make_shared<product_ref_>(a, b);
      dal::add_stored_object(pk, p, a, b,
                             convex_product_structure(a->structure(),
                                                      b->structure()),
                             p->pspt(), dal::PERMANENT_STATIC_OBJECT);
      pconvex_ref p1 = basic_convex_ref(p);
      if (p != p1) add_dependency(p, p1);
      return p;
    }
  }

  pconvex_ref parallelepiped_of_reference(dim_type nc, dim_type k) {
    if (nc <= 1) return simplex_of_reference(nc,k);
    return convex_ref_product(parallelepiped_of_reference(dim_type(nc-1),k),
                              simplex_of_reference(k));
  }

  pconvex_ref prism_of_reference(dim_type nc) {
    if (nc <= 2) return parallelepiped_of_reference(nc);
    else return convex_ref_product(simplex_of_reference(dim_type(nc-1)),
                                   simplex_of_reference(1));
  }

  /* equilateral ref convexes are used for estimation of convex quality */
  class equilateral_simplex_of_ref_ : public convex_of_reference {
    scalar_type r_inscr;
  public:
    scalar_type is_in(const base_node &pt) const {
      GMM_ASSERT1(pt.size() == cvs->dim(), "Dimension does not match");
      scalar_type d(0);
      for (size_type f = 0; f < normals().size(); ++f) {
        const base_node &x0 = (f ? convex<base_node>::points()[f-1]
                               : convex<base_node>::points().back());
        scalar_type v = gmm::vect_sp(pt-x0, normals()[f]);
        if (f == 0) d = v; else d = std::max(d,v);
      }
      return d;
    }

    scalar_type is_in_face(short_type f, const base_node &pt) const {
      GMM_ASSERT1(pt.size() == cvs->dim(), "Dimension does not match");
      const base_node &x0 = (f ? convex<base_node>::points()[f-1]
                             : convex<base_node>::points().back());
      return gmm::vect_sp(pt-x0, normals()[f]);
    }

    void project_into(base_node &pt) const {
      dim_type N = cvs->dim();
      GMM_ASSERT1(pt.size() == N, "Dimension does not match");
      base_node G(N); G.fill(0.);
      for (const base_node &x : convex<base_node>::points())
        G += x;
      gmm::scale(G, scalar_type(1)/scalar_type(N+1));
      for (size_type f = 0; f < normals().size(); ++f) {
        scalar_type r = gmm::vect_sp(pt-G, normals()[f]);
        if (r > r_inscr)
          pt = G + r_inscr/r*(pt-G);
      }

    }

    equilateral_simplex_of_ref_(size_type N) :
      convex_of_reference(simplex_structure(dim_type(N), 1), true)
    {
      //https://math.stackexchange.com/questions/2739915/radius-of-inscribed-sphere-of-n-simplex
      r_inscr = scalar_type(1)/sqrt(scalar_type(2*N)*scalar_type(N+1));

      pconvex_ref prev = equilateral_simplex_of_reference(dim_type(N-1));
      convex<base_node>::points().resize(N+1);
      normals_.resize(N+1);
      base_node G(N); G.fill(0.);
      for (short_type i=0; i < N+1; ++i) {
        convex<base_node>::points()[i].resize(N);
        if (i != N) {
          std::copy(prev->convex<base_node>::points()[i].begin(),
                    prev->convex<base_node>::points()[i].end(),
                    convex<base_node>::points()[i].begin());
          convex<base_node>::points()[i][N-1] = 0.;
        } else {
          convex<base_node>::points()[i] = scalar_type(1)/scalar_type(N) * G;
          convex<base_node>::points()[i][N-1]
            = sqrt(1. - gmm::vect_norm2_sqr(convex<base_node>::points()[i]));
        }
        G += convex<base_node>::points()[i];
      }
      gmm::scale(G, scalar_type(1)/scalar_type(N+1));
      for (size_type f=0; f < N+1; ++f) {
        normals_[f] = G - convex<base_node>::points()[f];
        gmm::scale(normals_[f], 1/gmm::vect_norm2(normals_[f]));
      }
      ppoints = store_point_tab(convex<base_node>::points());
      if (auto_basic) simplexify_convex(this, *psimplexified_convex);
    }
  };

  pconvex_ref equilateral_simplex_of_reference(dim_type nc) {
    if (nc <= 1) return simplex_of_reference(nc);
     dal::pstatic_stored_object_key
      pk = std::make_shared<convex_of_reference_key>(1, nc);
    dal::pstatic_stored_object o  = dal::search_stored_object(pk);
    if (o) return std::dynamic_pointer_cast<const convex_of_reference>(o);
    pconvex_ref p = std::make_shared<equilateral_simplex_of_ref_>(nc);
    dal::add_stored_object(pk, p, p->structure(), p->pspt(),
                           dal::PERMANENT_STATIC_OBJECT);
    return p;
  }


  /* generic convex with n global nodes      */

  class generic_dummy_ : public convex_of_reference {
  public:
    scalar_type is_in(const base_node &) const
    { GMM_ASSERT1(false, "Information not available here"); }
    scalar_type is_in_face(short_type, const base_node &) const
    { GMM_ASSERT1(false, "Information not available here"); }
    void project_into(base_node &) const
    { GMM_ASSERT1(false, "Operation not available here"); }

    generic_dummy_(dim_type d, size_type n, short_type nf) :
      convex_of_reference(generic_dummy_structure(d, n, nf), true)
    {
      convex<base_node>::points().resize(n);
      normals_.resize(0);
      base_node P(d);
      std::fill(P.begin(), P.end(), scalar_type(1)/scalar_type(20));
      std::fill(convex<base_node>::points().begin(), convex<base_node>::points().end(), P);
      ppoints = store_point_tab(convex<base_node>::points());
    }
  };

  pconvex_ref generic_dummy_convex_ref(dim_type nc, size_type n,
                                       short_type nf) {
    dal::pstatic_stored_object_key
      pk = std::make_shared<convex_of_reference_key>(2, nc, short_type(n), nf);
    dal::pstatic_stored_object o = dal::search_stored_object(pk);
    if (o) return std::dynamic_pointer_cast<const convex_of_reference>(o);
    pconvex_ref p = std::make_shared<generic_dummy_>(nc, n, nf);
    dal::add_stored_object(pk, p, p->structure(), p->pspt(),
                           dal::PERMANENT_STATIC_OBJECT);
    return p;
  }


}  /* end of namespace bgeot.                                              */
