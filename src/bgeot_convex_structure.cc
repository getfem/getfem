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


#include "getfem/dal_singleton.h"
#include "getfem/dal_static_stored_objects.h"
#include "getfem/bgeot_convex_structure.h"
#include "getfem/bgeot_comma_init.h"

namespace bgeot {

  /* ******************************************************************** */
  /*									  */
  /* class   convex_structure                                             */
  /*									  */
  /* ******************************************************************** */

  void convex_structure::add_point_adaptative(short_type i, short_type f) {
    GMM_ASSERT1(i <= nbpt,  "convex_structure::add_point_adaptative: "
		"internal error");
    if (i == nbpt) nbpt++;
    if (f != short_type(-1)) {
      faces[f].resize(faces[f].size() + 1);
      (faces[f])[faces[f].size() - 1] = i;
    }
  }

  void convex_structure::init_for_adaptative(pconvex_structure cvs) {
    *this = *(basic_structure(cvs));
    std::fill(faces_struct.begin(),faces_struct.end(),
	      pconvex_structure());
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
  class convex_structure_key : virtual public dal::static_stored_object_key {
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
    dal::pstatic_stored_object_key
      pcsk = std::make_shared<convex_structure_key>(0, nc, 1);
    dal::pstatic_stored_object o = dal::search_stored_object(pcsk);
    if (o) return std::dynamic_pointer_cast<const convex_structure>(o);
    
    auto p = std::make_shared<simplex_structure_>();
    pconvex_structure pcvs = pconvex_structure(p);
    p->Nc = dim_type(nc); p->nbpt = short_type(nc+1);
    p->nbf = short_type(nc ? nc+1 : 0);
    p->faces_struct.resize(p->nbf);
    p->faces.resize(p->nbf);
    p->dir_points_.resize(p->Nc + 1);
    p->auto_basic = true;
    for (short_type i = 0; i < p->nbf; i++) { 
      p->dir_points_[i] = i;
      p->faces_struct[i] = simplex_structure(dim_type(nc-1));
      p->faces[i].resize(nc);
      for (short_type j = 0; j < nc; j++)
	(p->faces[i])[j] = (j >= i) ? short_type(j + 1) : j;
    }
    if (nc == 0)
      dal::add_stored_object(pcsk, pcvs, dal::PERMANENT_STATIC_OBJECT);
    else
      dal::add_stored_object(pcsk, pcvs, simplex_structure(dim_type(nc-1)),
			     dal::PERMANENT_STATIC_OBJECT);
    return pcvs;
  }

  /* ******************************************************************** */
  /* K-simplex structures                                                 */
  /* ******************************************************************** */

  struct K_simplex_structure_ : public convex_structure {
    
    K_simplex_structure_(dim_type NN, short_type KK) {
      Nc = NN; nbpt = short_type(alpha(Nc, KK)); nbf = short_type(Nc+1);
      basic_pcvs = simplex_structure(NN);
      faces_struct.resize(nbf);
      faces.resize(nbf);
      dir_points_.resize(Nc+1);
      
      for (int i = 0; i < nbf; i++) { 
	if (KK > 0) {
	  faces_struct[i] = simplex_structure(dim_type(Nc-1), KK); 
	  faces[i].resize(faces_struct[i]->nb_points());
	}
	else {
	  faces_struct[i] = pconvex_structure(); 
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
      
      for (short_type r = 1; r < nbpt; ++r) {
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
    dal::pstatic_stored_object_key
      pcsk = std::make_shared<convex_structure_key>(0, nc, K);
    dal::pstatic_stored_object o = dal::search_stored_object(pcsk);
    if (o) return std::dynamic_pointer_cast<const convex_structure>(o);

    pconvex_structure p = std::make_shared<K_simplex_structure_>(nc, K);
    dal::add_stored_object(pcsk, p, simplex_structure(dim_type(nc-1), K),
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
    if (nbt <= 3) return simplex_structure(dim_type(nbt-1));
    
    dal::pstatic_stored_object_key
      pcsk = std::make_shared<convex_structure_key>(1, dim_type(nbt));
    dal::pstatic_stored_object o = dal::search_stored_object(pcsk);
    if (o) return std::dynamic_pointer_cast<const convex_structure>(o);

    auto p = std::make_shared<polygon_structure_>();
    pconvex_structure pcvs(p);
    p->Nc = 2; p->nbpt = nbt; p->nbf = nbt;
    p->auto_basic = true;
    p->faces_struct.resize(p->nbf);
    p->faces = std::vector< std::vector<short_type> >(p->nbf);
    p->dir_points_ = std::vector<short_type>(p->Nc + 1);
    
    for (int i = 0; i < p->nbf; i++) { 
      p->faces_struct[i] = simplex_structure(1);
      p->faces[i] = std::vector<short_type>(2);
      for (int j = 0; j < 2; j++)
	(p->faces[i])[j] = short_type((i+j) % nbt);
    }
    
    p->dir_points_[0] = 0;
    p->dir_points_[1] = 1;
    p->dir_points_[2] = short_type(nbt - 1);
    
    dal::add_stored_object(pcsk, pcvs, simplex_structure(1),
			   dal::PERMANENT_STATIC_OBJECT);
    return pcvs;
  }

  /* ******************************************************************** */
  /* Direct product of convex structures                                  */
  /* ******************************************************************** */

  DAL_DOUBLE_KEY(cv_pr_key_, pconvex_structure, pconvex_structure);

  struct cv_pr_structure_ : public convex_structure {
    cv_pr_structure_(pconvex_structure cv1, pconvex_structure cv2) {
      Nc = dim_type(cv1->dim() + cv2->dim());
      prod_a = cv1; prod_b = cv2;
      nbpt = short_type(cv1->nb_points() * cv2->nb_points());
      nbf = short_type(cv1->nb_faces() + cv2->nb_faces());
      if (basic_structure(cv1) != cv1 || basic_structure(cv2) != cv2)
	basic_pcvs = convex_product_structure(basic_structure(cv1),
					      basic_structure(cv2));
      else
	auto_basic = true;
      faces_struct.resize(nbf);
      faces = std::vector< std::vector<short_type> >(nbf);

      if (cv1->ind_dir_points().size() && cv2->ind_dir_points().size()) {
	dir_points_ = std::vector<short_type>(Nc + 1);

	for (int i = 0; i <= cv1->dim(); i++)
	  dir_points_[i]
	    = short_type(cv1->ind_dir_points()[i]
			 + cv2->ind_dir_points()[0] * cv1->nb_points());
	for (int i = 1; i <= cv2->dim(); i++)
	  dir_points_[cv1->dim()+i]
	    = short_type(cv1->ind_dir_points()[0]
			 + cv2->ind_dir_points()[i] * cv1->nb_points());
      }

      for (short_type i = 0; i < cv1->nb_faces(); i++) { 
	if (cv1->nb_points_of_face(i) == 1)
	  faces_struct[i] = cv2;
	else
	  faces_struct[i]
	    = (cv1->faces_structure()[i] == pconvex_structure()) ?
	    pconvex_structure()
	    : convex_product_structure(cv1->faces_structure()[i], cv2);

	faces[i] = std::vector<short_type>(cv1->nb_points_of_face(i)
					      * cv2->nb_points());

	for (short_type j = 0; j < cv1->nb_points_of_face(i); j++)
	  for (short_type l = 0; l < cv2->nb_points(); l++) {
	    (faces[i])[l*cv1->nb_points_of_face(i)+j]
	      = short_type((cv1->ind_points_of_face(i))[j]
			   + l * cv1->nb_points());
	  }
      }
      for (short_type i = 0; i < cv2->nb_faces(); i++) { 
	short_type k = cv1->nb_faces();
	if (cv2->nb_points_of_face(i) == 1)
	  faces_struct[i+k] = cv1;
	else
	  faces_struct[i+k]
	    = (cv2->faces_structure()[i] == pconvex_structure()) ? 
	    pconvex_structure()
	    : convex_product_structure(cv1, cv2->faces_structure()[i]);

	faces[i+k] = std::vector<short_type>(cv2->nb_points_of_face(i)
					      * cv1->nb_points());

	for (short_type j = 0; j < cv2->nb_points_of_face(i); j++)
	  for (short_type l = 0; l < cv1->nb_points(); l++) {
	    (faces[i+k])[j*cv1->nb_points()+l]
	      = short_type(l + (cv2->ind_points_of_face(i))[j]
			   * cv1->nb_points()); 
	  }
      }
    }
  };

  pconvex_structure convex_product_structure(pconvex_structure a,
					     pconvex_structure b) {
    
    dal::pstatic_stored_object_key pcsk = std::make_shared<cv_pr_key_>(a, b);
    dal::pstatic_stored_object o = dal::search_stored_object(pcsk);
    if (o) return std::dynamic_pointer_cast<const convex_structure>(o);
    pconvex_structure p = std::make_shared<cv_pr_structure_>(a, b);
    dal::add_stored_object(pcsk, p, a, b, dal::PERMANENT_STATIC_OBJECT);
    for (size_type k = 0; k < p->nb_faces(); ++k) {
      if (exists_stored_object(p->faces_structure()[k]))
	dal::add_dependency(p, p->faces_structure()[k]);
    }
    return p;
  }

  /* ******************************************************************** */
  /* parallelepiped structures.                                           */
  /* ******************************************************************** */

  struct parallelepiped_ : virtual public dal::static_stored_object {
    pconvex_structure p;
    parallelepiped_()
    { DAL_STORED_OBJECT_DEBUG_CREATED(this, "parallelepiped structure"); }
    ~parallelepiped_()
    { DAL_STORED_OBJECT_DEBUG_DESTROYED(this, "parallelepiped structure"); }
  };

  DAL_SIMPLE_KEY(parallelepiped_key_, dim_type);

  pconvex_structure parallelepiped_structure(dim_type nc) {
    if (nc <= 1) return simplex_structure(nc);
    dal::pstatic_stored_object_key
      pcsk = std::make_shared<parallelepiped_key_>(nc);

    dal::pstatic_stored_object o = dal::search_stored_object(pcsk);
    if (o) return ((std::dynamic_pointer_cast<const parallelepiped_>(o))->p);
    auto p = std::make_shared<parallelepiped_>();
    p->p = convex_product_structure(parallelepiped_structure(dim_type(nc-1)),
				    simplex_structure(1));
    dal::add_stored_object(pcsk, dal::pstatic_stored_object(p), p->p,
			   dal::PERMANENT_STATIC_OBJECT);
    return p->p;
  }


  /* ******************************************************************** */
  /*	Incomplete Q2 structure for n=2 or 3.                             */
  /* ******************************************************************** */
  /* By Yao Koutsawa  <yao.koutsawa@tudor.lu> 2012-12-10                  */

  struct Q2_incomplete_structure_ : public convex_structure {
    friend pconvex_structure Q2_incomplete_structure(dim_type nc);
  };
  
  DAL_SIMPLE_KEY(Q2_incomplete_structure_key_, dim_type);
  
  pconvex_structure Q2_incomplete_structure(dim_type nc) {
    GMM_ASSERT1(nc == 2 || nc == 3, "Bad parameter, expected value 2 or 3");
    dal::pstatic_stored_object_key
      pcsk = std::make_shared<Q2_incomplete_structure_key_>(nc);
    dal::pstatic_stored_object o = dal::search_stored_object(pcsk);
    if (o) return std::dynamic_pointer_cast<const convex_structure>(o);
    
    auto p = std::make_shared<Q2_incomplete_structure_>();
    pconvex_structure pcvs(p);
    p->Nc = nc;
    p->nbpt = (nc == 2) ? 8 : 20;
    p->nbf =  (nc == 2) ? 4 : 6;
    p->basic_pcvs =  parallelepiped_structure(nc);
    p->faces_struct.resize(p->nbf);
    p->faces = std::vector< std::vector<short_type> >(p->nbf);
    p->dir_points_ = std::vector<short_type>(p->Nc + 1);
    
    if (nc == 2) {
      // 5--6--7
      // |     |
      // 3     4
      // |     |
      // 0--1--2
      sc(p->faces[0]) = 2,4,7;
      sc(p->faces[1]) = 0,3,5;
      sc(p->faces[2]) = 5,6,7;
      sc(p->faces[3]) = 0,1,2;
      
      p->dir_points_[0] = 0;
      p->dir_points_[1] = 2;
      p->dir_points_[2] = 5;
    } else {
      //       17---18---19
      //      /|        /|
      //     / 10      / 11
      //   15  |      16 |
      //   /   5----6/---7
      //  /   /     /   /
      // 12---13---14  /
      // |  3      |  4
      // 8 /       9 /
      // |/        |/
      // 0----1----2
      sc(p->faces[0]) = 2,4,7,9,11,14,16,19;
      sc(p->faces[1]) = 0,3,5,8,10,12,15,17;
      
      sc(p->faces[2]) = 5,6,7,10,11,17,18,19;
      sc(p->faces[3]) = 0,1,2,8,9,12,13,14;
      
      sc(p->faces[4]) = 12,13,14,15,16,17,18,19;
      sc(p->faces[5]) = 0,1,2,3,4,5,6,7;
      
      p->dir_points_[0] = 0;
      p->dir_points_[1] = 2;
      p->dir_points_[2] = 5;
      p->dir_points_[3] = 12;
    }
    
    for (int i = 0; i < p->nbf; i++) {
      p->faces_struct[i] = (nc == 2) ? simplex_structure(1, 2)
        : Q2_incomplete_structure(2);
    }
    
    dal::add_stored_object(pcsk, pcvs, parallelepiped_structure(dim_type(nc-1)),
                           dal::PERMANENT_STATIC_OBJECT);
    return pcvs;
  }


  /* ******************************************************************** */
  /*	Generic dummy convex with n global nodes.                         */
  /* ******************************************************************** */

  struct dummy_structure_ : public convex_structure {
    friend pconvex_structure generic_dummy_structure(dim_type, size_type,
						     size_type);
  };
  
  pconvex_structure generic_dummy_structure(dim_type nc, size_type n,
					    size_type nf) {
    dal::pstatic_stored_object_key
      pcsk = std::make_shared<convex_structure_key>(2, nc, short_type(n),
						    short_type(nf));
    dal::pstatic_stored_object o = dal::search_stored_object(pcsk);
    if (o) return std::dynamic_pointer_cast<const convex_structure>(o);
    auto p = std::make_shared<dummy_structure_>();
    pconvex_structure pcvs(p);
    p->Nc = nc; p->nbpt = short_type(n); p->nbf = 0;
    p->faces_struct.resize(nf);
    p->faces.resize(nf);
    for (size_type j = 0; j < nf; ++j) {
      if (nc == 0)
	p->faces_struct[j] = simplex_structure(0);
      else p->faces_struct[j] = generic_dummy_structure(dim_type(nc-1), n, nc);
      p->faces[j].resize(n);
      for (short_type k = 0; k < n; ++k) p->faces[j][k] = k;
    }
    p->dir_points_.resize(0);
    p->auto_basic = true;
    if (nc == 0)
      dal::add_stored_object(pcsk, pcvs, dal::PERMANENT_STATIC_OBJECT);
    else
      dal::add_stored_object(pcsk, pcvs,
			     generic_dummy_structure(dim_type(nc-1), n, nc),
			     dal::PERMANENT_STATIC_OBJECT);
    return pcvs;
  }

}  /* end of namespace bgeot.                                            */
