// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_convex_structure.cc : convex structures.
//           
// Date    : December 20, 1999.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1999-2005 Yves Renard
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


#include <dal_singleton.h>
#include <dal_static_stored_objects.h>
#include <bgeot_convex_structure.h>


namespace bgeot {

  /* ******************************************************************** */
  /*									  */
  /* class   convex_structure                                             */
  /*									  */
  /* ******************************************************************** */

  void convex_structure::add_point_adaptative(short_type i, short_type f) {
    if (nbpt < i) throw dal::internal_error(
		   "convex_structure::add_point_adaptative : internal error");
    if (i == nbpt) nbpt++;
    if (f != short_type(-1)) {
      faces[f].resize(faces[f].size() + 1);
      (faces[f])[faces[f].size() - 1] = i;
    }
  }

  void convex_structure::init_for_adaptative(pconvex_structure cvs) {
    *this = *(cvs->basic_structure());
    std::fill(faces_struct.begin(),faces_struct.end(),
	      (const convex_structure *)(0));
    std::fill(faces.begin(),faces.end(), convex_ind_ct());
     dir_points_ = convex_ind_ct();
    nbpt = 0;
  }

  std::ostream &operator <<(std::ostream &o, const convex_structure &cv) {
    o << "convex structure of dimension " << int(cv.dim()) << " with "
      << cv.nb_points() << " points and " << cv.nb_faces() << " faces "
      << endl;
    // a completer au besoin
    return o;
  }

  // Key type for static storing
  class convex_structure_key : public dal::static_stored_object_key {
    int type; // 0 = simplex structure degree K
              // 1 = polygon (N = nb of points, K = 0)
              // 2 = dummy (N = dimension, K = nbpt)
    dim_type N; short_type K; short_type nf;
  public :
    virtual bool compare(const static_stored_object_key &oo) const {
      const convex_structure_key &o
	= dynamic_cast<const convex_structure_key &>(oo);
      if (type < o.type) return true;
      if (type > o.type) return false;
      if (N < o.N) return true;
      if (N > o.N) return false;
      if (K < o.K) return true;
      if (K > o.K) return false;
      if (nf < o.nf) return true;
      return false;
    }
    convex_structure_key(int t, dim_type NN, short_type KK = 0,
			 short_type nnf = 0)
      : type(t), N(NN), K(KK), nf(nnf)  {}
  };
    
  /* ******************************************************************** */
  /* simplex structures                                                   */
  /* ******************************************************************** */

  class simplex_structure_ : public convex_structure
  { friend pconvex_structure simplex_structure(dim_type nc); };

#ifdef GETFEM_HAVE_QDLIB
#  include <qd/fpu.h>
#endif

  pconvex_structure simplex_structure(dim_type nc) {
#ifdef GETFEM_HAVE_QDLIB
    /* initialisation for QD on intel CPUs */
    static bool fpu_init = false;
    if (!fpu_init) {
      unsigned int old_cw;
      fpu_fix_start(&old_cw);
      fpu_init = true;
    }
#endif
    dal::pstatic_stored_object o
      = dal::search_stored_object(convex_structure_key(0, nc, 1));
    if (o) return dal::stored_cast<convex_structure>(o);
    
    simplex_structure_ *p = new simplex_structure_;
    p->Nc = nc; p->nbpt = nc+1; p->nbf = nc+1;
    p->faces_struct.resize(p->nbf);
    p->faces.resize(p->nbf);
    p->dir_points_.resize(p->Nc + 1);
    p->basic_pcvs = p;
    if (nc == 0) {	
      p->faces_struct[0] = p;
      p->faces[0].resize(1);
      (p->faces[0])[0] = 0;
    }
    else
      for (int i = 0; i < p->nbf; i++) { 
	p->dir_points_[i] = i;
	p->faces_struct[i] = simplex_structure(nc-1).get();
	p->faces[i].resize(nc);
	for (int j = 0; j < nc; j++)
	  (p->faces[i])[j] = (j >= i) ? (j + 1) : j;
      }
    if (nc == 0)
      dal::add_stored_object(new convex_structure_key(0, nc, 1), p,
			     dal::PERMANENT_STATIC_OBJECT);
    else
      dal::add_stored_object(new convex_structure_key(0, nc, 1), p, 
			     simplex_structure(nc-1),
			     dal::PERMANENT_STATIC_OBJECT);
    return p;
  }

  /* ******************************************************************** */
  /* K-simplex structures                                                 */
  /* ******************************************************************** */

  struct K_simplex_structure_ : public convex_structure {
    
    K_simplex_structure_(dim_type NN, short_type KK) {
      Nc = NN; nbpt = alpha(Nc, KK); nbf = Nc+1;
      basic_pcvs = simplex_structure(NN).get();
      faces_struct.resize(nbf);
      faces.resize(nbf);
      dir_points_.resize(Nc+1);
      
      for (int i = 0; i < nbf; i++) { 
	if (KK > 0) {
	  faces_struct[i] = simplex_structure(Nc-1, KK).get(); 
	  faces[i].resize(faces_struct[i]->nb_points());
	}
	else {
	  faces_struct[i] = NULL; 
	  faces[i].resize(0);
	}
      }
      
      base_node c(Nc); c.fill(0.0);
      std::vector<int> pf(Nc+1); std::fill(pf.begin(), pf.end(), 0);
      size_type l, sum = 0, pd = 0;
      if (KK == 0) c.fill(scalar_type(1.0) / scalar_type(Nc+1));
      else {
	for (l = 1; l <= Nc; ++l) (faces[l])[(pf[l])++] = 0;
	dir_points_[pd++] = 0;
      }
      
      for (size_type r = 1; r < nbpt; ++r) {
	l = 0;
	c[l] += scalar_type(1.0) / scalar_type(KK); ++sum;
	while (sum > KK) {
	  sum -= size_type(floor(0.5+(c[l] * KK)));
	  c[l] = 0.0; ++l; c[l] += scalar_type(1.0) / scalar_type(KK);
	  ++sum;
	}
	for (l = 1; l <= Nc; ++l)
	  if (c[l-1] == scalar_type(0.0)) (faces[l])[(pf[l])++] = r;
	if (sum == KK) {
	  (faces[0])[(pf[0])++] = r;
	  if (*(std::max_element(c.begin(), c.end())) == scalar_type(1.0))
	    dir_points_[pd++] = r;
	}
      }
    }
  };
  
  pconvex_structure simplex_structure(dim_type nc, short_type K) {
    if (nc == 0) K = 1;
    if (K == 1) return simplex_structure(nc);
    dal::pstatic_stored_object o
      = dal::search_stored_object(convex_structure_key(0, nc, K));
    if (o) return dal::stored_cast<convex_structure>(o);
    pconvex_structure p = new K_simplex_structure_(nc, K);
    dal::add_stored_object(new convex_structure_key(0, nc, K), p,
			   simplex_structure(nc-1, K),
			   dal::PERMANENT_STATIC_OBJECT);
    return p;
  }
  
  /* ******************************************************************** */
  /* polygon structures                                                   */
  /* ******************************************************************** */

  struct polygon_structure_ : public convex_structure {
    friend pconvex_structure polygon_structure(short_type nc);
  };
  
  pconvex_structure polygon_structure(short_type nbt) {
    if (nbt <= 1) return simplex_structure(0);
    if (nbt <= 3) return simplex_structure(nbt-1);
    dal::pstatic_stored_object o
      = dal::search_stored_object(convex_structure_key(1, nbt));
    if (o) return dal::stored_cast<convex_structure>(o);

    polygon_structure_ *p = new polygon_structure_;
    p->Nc = 2; p->nbpt = nbt; p->nbf = nbt;
    p->basic_pcvs = p;
    p->faces_struct = std::vector<const convex_structure *>(p->nbf);
    p->faces = std::vector< std::vector<short_type> >(p->nbf);
    p->dir_points_ = std::vector<short_type>(p->Nc + 1);
    
    for (int i = 0; i < p->nbf; i++) { 
      p->faces_struct[i] = simplex_structure(1).get();
      p->faces[i] = std::vector<short_type>(2);
      for (int j = 0; j < 2; j++)
	(p->faces[i])[j] = ((i+j) % nbt);
    }
    
    p->dir_points_[0] = 0;
    p->dir_points_[1] = 1;
    p->dir_points_[2] = nbt - 1;
    
    dal::add_stored_object(new convex_structure_key(1, nbt), p,
			   simplex_structure(1),
			   dal::PERMANENT_STATIC_OBJECT);
    return p;
  }

  /* ******************************************************************** */
  /* direct product of convex structures                                  */
  /* ******************************************************************** */

  struct cv_pr_key_ : public dal::static_stored_object_key {
    pconvex_structure cv1, cv2;
    virtual bool compare(const static_stored_object_key &oo) const {
      const cv_pr_key_ &o = dynamic_cast<const cv_pr_key_ &>(oo);
      if (cv1 < o.cv1) return true;
      if (o.cv1 < cv1) return false; 
      if (cv2 < o.cv2) return true;
      return false;
    }
    cv_pr_key_(pconvex_structure a, pconvex_structure b) { cv1=a; cv2=b; }
  };

  struct cv_pr_structure_ : public convex_structure {
    cv_pr_structure_(pconvex_structure cv1, pconvex_structure cv2) {
      Nc = cv1->dim() + cv2->dim();
      prod_a = cv1; prod_b = cv2;
      nbpt = cv1->nb_points() * cv2->nb_points();
      nbf = cv1->nb_faces() + cv2->nb_faces();
      if (cv1->basic_structure() != cv1 || cv2->basic_structure() != cv2)
	basic_pcvs = convex_product_structure(cv1->basic_structure(),
					      cv2->basic_structure()).get();
      else
	basic_pcvs = this;
      faces_struct = std::vector<const convex_structure *>(nbf);
      faces = std::vector< std::vector<short_type> >(nbf);

      if (cv1->ind_dir_points().size() && cv2->ind_dir_points().size()) {
	dir_points_ = std::vector<short_type>(Nc + 1);

	for (int i = 0; i <= cv1->dim(); i++)
	  dir_points_[i] = cv1->ind_dir_points()[i]
	    + cv2->ind_dir_points()[0] * cv1->nb_points();
	for (int i = 1; i <= cv2->dim(); i++)
	  dir_points_[cv1->dim()+i] = cv1->ind_dir_points()[0]
	    + cv2->ind_dir_points()[i] * cv1->nb_points();
      }

      for (int i = 0; i < cv1->nb_faces(); i++) { 
	if (cv1->nb_points_of_face(i) == 1)
	  faces_struct[i] = cv2.get();
	else
	  faces_struct[i]
	    = (cv1->faces_structure()[i] == NULL) ? NULL
	    : convex_product_structure(cv1->faces_structure()[i], cv2).get();

	faces[i] = std::vector<short_type>(cv1->nb_points_of_face(i)
					      * cv2->nb_points());

	for (int j = 0; j < cv1->nb_points_of_face(i); j++)
	  for (int l = 0; l < cv2->nb_points(); l++) {
	    (faces[i])[l*cv1->nb_points_of_face(i)+j]
	      = (cv1->ind_points_of_face(i))[j] + l * cv1->nb_points();
	  }
      }
      for (int i = 0; i < cv2->nb_faces(); i++) { 
	int k = cv1->nb_faces();
	if (cv2->nb_points_of_face(i) == 1)
	  faces_struct[i+k] = cv1.get();
	else
	  faces_struct[i+k]
	    = (cv2->faces_structure()[i] == NULL) ? NULL
	    : convex_product_structure(cv1, cv2->faces_structure()[i]).get();

	faces[i+k] = std::vector<short_type>(cv2->nb_points_of_face(i)
					      * cv1->nb_points());

	for (int j = 0; j < cv2->nb_points_of_face(i); j++)
	  for (int l = 0; l < cv1->nb_points(); l++) {
	    (faces[i+k])[j*cv1->nb_points()+l]
	      = l + (cv2->ind_points_of_face(i))[j] * cv1->nb_points(); 
	  }
      }
    }
  };

  pconvex_structure convex_product_structure(pconvex_structure a,
					     pconvex_structure b) {
    dal::pstatic_stored_object o
      = dal::search_stored_object(cv_pr_key_(a, b));
    if (o) return dal::stored_cast<convex_structure>(o);
    pconvex_structure p = new cv_pr_structure_(a, b);
    dal::add_stored_object(new cv_pr_key_(a, b), p, a, b,
			   dal::PERMANENT_STATIC_OBJECT);
    for (size_type k = 0; k < p->nb_faces(); ++k)
      dal::add_dependency(p, p->faces_structure()[k]);
    return p;
  }

  /* ******************************************************************** */
  /* parallelepiped structures.                                           */
  /* ******************************************************************** */

  struct parallelepiped_ : public dal::static_stored_object {
    pconvex_structure p;
  };

  struct parallelepiped_key_ : public dal::static_stored_object_key {
    dim_type n;
    virtual bool compare(const static_stored_object_key &oo) const {
      const parallelepiped_key_ &o
	= dynamic_cast<const parallelepiped_key_ &>(oo);
      if (n < o.n) return true;
      return false;
    }
    parallelepiped_key_(dim_type nn) : n(nn) { }
  };

  pconvex_structure parallelepiped_structure(dim_type nc) {
    if (nc <= 1) return simplex_structure(nc);
    dal::pstatic_stored_object o
      = dal::search_stored_object(parallelepiped_key_(nc));
    if (o) return (dal::stored_cast<parallelepiped_>(o))->p;
    parallelepiped_ *p = new parallelepiped_;
    p->p = convex_product_structure(parallelepiped_structure(nc-1),
				    simplex_structure(1));
    dal::add_stored_object(new parallelepiped_key_(nc), p, p->p,
			   dal::PERMANENT_STATIC_OBJECT);
    return p->p;
  }

  // generic convex with n global nodes

  struct dummy_structure_ : public convex_structure {
    friend pconvex_structure generic_dummy_structure(dim_type, size_type,
						     size_type);
  };
  
  pconvex_structure generic_dummy_structure(dim_type nc, size_type n,
					    size_type nf) {
    dal::pstatic_stored_object o
      = dal::search_stored_object(convex_structure_key(2, nc, n, nf));
    if (o) return dal::stored_cast<convex_structure>(o);
    dummy_structure_ *p = new dummy_structure_;
    p->Nc = nc; p->nbpt = n; p->nbf = 0;
    p->faces_struct.resize(nf);
    p->faces.resize(nf);
    for (size_type j = 0; j < nf; ++j) {
      if (nc == 0) p->faces_struct[j] = p;
      else p->faces_struct[j] = generic_dummy_structure(nc-1, n, nc).get();
      p->faces[j].resize(n);
      for (size_type k = 0; k < n; ++k) p->faces[j][k] = k;
    }
    p->dir_points_.resize(0);
    p->basic_pcvs = p;
    if (nc == 0)
      dal::add_stored_object(new convex_structure_key(2, nc, n, nf), p,
			     dal::PERMANENT_STATIC_OBJECT);
    else
      dal::add_stored_object(new convex_structure_key(2, nc, n, nf), p,
			     generic_dummy_structure(nc-1, n, nc),
			     dal::PERMANENT_STATIC_OBJECT);
    return p;
  }

}  /* end of namespace bgeot.                                            */
