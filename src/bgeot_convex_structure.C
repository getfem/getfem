/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_convex_structure.C : convex structures.                */
/*     									   */
/*                                                                         */
/* Date : December 20, 1999.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 1999-2002  Yves Renard.                                   */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */


#include <dal_singleton.h>
#include <bgeot_convex_structure.h>

namespace bgeot
{

  /* ******************************************************************** */
  /*									  */
  /* class   convex_structure                                             */
  /*									  */
  /* ******************************************************************** */

  void convex_structure::add_point_adaptative(short_type i, short_type f)
  {
    if (nbpt < i) throw dal::internal_error(
		   "convex_structure::add_point_adaptative : internal error");
    if (i == nbpt) nbpt++;
    if (f != short_type(-1))
    {
      faces[f].resize(faces[f].size() + 1);
      (faces[f])[faces[f].size() - 1] = i;
    }
  }

  void convex_structure::init_for_adaptative(pconvex_structure cvs)
  {
    *this = *(cvs->basic_structure());
    std::fill(faces_struct.begin(),faces_struct.end(),pconvex_structure(NULL));
    std::fill(faces.begin(),faces.end(), convex_ind_ct());
     dir_points_ = convex_ind_ct();
    nbpt = 0;
  }

  std::ostream &operator <<(std::ostream &o, const convex_structure &cv)
  {
    o << "convex structure of dimension " << int(cv.dim()) << " with "
      << cv.nb_points() << " points and " << cv.nb_faces() << " faces " << endl;
    // a completer au besoin
    return o;
  }
    
  /* ******************************************************************** */
  /* simplex structures                                                   */
  /* ******************************************************************** */

  class simplex_structure_ : public convex_structure
  { friend pconvex_structure simplex_structure(dim_type nc); };

#ifdef GETFEM_HAVE_QDLIB
#  include <qd/fpu.h>
#endif

  struct simplex_structure_data {
    dal::dynamic_array<simplex_structure_> simplex;
    int nb_simplex;
    simplex_structure_data() : nb_simplex(-1) {
#ifdef GETFEM_HAVE_QDLIB
      /* initialisation for QD on intel CPUs */
      unsigned int old_cw;
      fpu_fix_start(&old_cw);
#endif
    }
  };

  pconvex_structure simplex_structure(dim_type nc)
  {
    simplex_structure_data& ssd = dal::singleton<simplex_structure_data>::instance();
    simplex_structure_ *p = &(ssd.simplex[nc]);
    if (ssd.nb_simplex < nc)
    {
      p->Nc = nc; p->nbpt = nc+1; p->nbf = nc+1;
      p->faces_struct.resize(p->nbf);
      p->faces.resize(p->nbf);
      p->dir_points_.resize(p->Nc + 1);
      p->basic_pcvs = p;
      if (nc == 0)
      {	
	p->faces_struct[0] = p;
	p->faces[0].resize(1);
	(p->faces[0])[0] = 0;
      }
      else
	for (int i = 0; i < p->nbf; i++)
	{ 
	  p->dir_points_[i] = i;
	  p->faces_struct[i] = simplex_structure(nc-1);
	  p->faces[i].resize(nc);
	  for (int j = 0; j < nc; j++)
	    (p->faces[i])[j] = (j >= i) ? (j + 1) : j;
	}
      
      ssd.nb_simplex = nc;
    }
    return p;
  }

  /* ******************************************************************** */
  /* K-simplex structures                                                 */
  /* ******************************************************************** */

  struct K_simplex_light_
  {
    dim_type nc; short_type K;
    bool operator < (const K_simplex_light_ &l) const
    {
      if (nc < l.nc) return true; if (nc > l.nc) return false; 
      if (K < l.K) return true; return false;
    }
    K_simplex_light_(dim_type n, short_type k) { nc = n; K = k; }
    K_simplex_light_(void) { }
  };

  class K_simplex_structure_ : public convex_structure
  {
    public :

      K_simplex_structure_(const K_simplex_light_ &ls)
      {
	Nc = ls.nc; nbpt = alpha(Nc, ls.K); nbf = Nc+1;
	basic_pcvs = simplex_structure(ls.nc);
	faces_struct.resize(nbf);
	faces.resize(nbf);
	dir_points_.resize(Nc+1);
	
	for (int i = 0; i < nbf; i++)
	{ 
	  if (ls.K > 0) {
	    faces_struct[i] = simplex_structure(Nc-1, ls.K); 
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
	if (ls.K == 0) c.fill(scalar_type(1.0) / scalar_type(Nc+1));
	else {
	  for (l = 1; l <= Nc; ++l) (faces[l])[(pf[l])++] = 0;
	  dir_points_[pd++] = 0;
	}
	
	for (size_type r = 1; r < nbpt; ++r)
	{
	  l = 0;
	  c[l] += scalar_type(1.0) / scalar_type(ls.K); ++sum;
	  while (sum > ls.K)
	  {
	    sum -= size_type(floor(0.5+(c[l] * ls.K)));
	    c[l] = 0.0; ++l; c[l] += scalar_type(1.0) / scalar_type(ls.K);
	    ++sum;
	  }
	  for (l = 1; l <= Nc; ++l)
	    if (c[l-1] == scalar_type(0.0)) (faces[l])[(pf[l])++] = r;
	  if (sum == ls.K)
	  {
	    (faces[0])[(pf[0])++] = r;
	    if (*(std::max_element(c.begin(), c.end())) == scalar_type(1.0))
	      dir_points_[pd++] = r;
	  }
	}
      }

  };

  typedef dal::FONC_TABLE<K_simplex_light_, K_simplex_structure_> simplexes_FONC_TABLE;

  pconvex_structure simplex_structure(dim_type nc, short_type K) {
    if (nc == 0) return simplex_structure(0);
    if (K == 1) return simplex_structure(nc);
    return dal::singleton<simplexes_FONC_TABLE>::instance().add(K_simplex_light_(nc, K));
  }

  /* ******************************************************************** */
  /* polygon structures                                                   */
  /* ******************************************************************** */

  class polygon_structure_ : public convex_structure {
    friend pconvex_structure polygon_structure(short_type nc);
  };

  struct polygon_structure_data {
    dal::bit_vector ind_polygon;
    dal::dynamic_array<polygon_structure_> polygon;
  };
  
  pconvex_structure polygon_structure(short_type nbt)
  {
    polygon_structure_data &psd = dal::singleton<polygon_structure_data>::instance();

    polygon_structure_ *p = &(psd.polygon[nbt]);
    if (nbt < 4) return simplex_structure(nbt-1);

    if (!psd.ind_polygon.is_in(nbt))
    {
      p->Nc = 2; p->nbpt = nbt; p->nbf = nbt;
      p->basic_pcvs = p;
      p->faces_struct = std::vector<pconvex_structure>(p->nbf);
      p->faces = std::vector< std::vector<short_type> >(p->nbf);
      p->dir_points_ = std::vector<short_type>(p->Nc + 1);

      for (int i = 0; i < p->nbf; i++)
      { 
	p->faces_struct[i] = simplex_structure(1);
	p->faces[i] = std::vector<short_type>(2);
	for (int j = 0; j < 2; j++)
	  (p->faces[i])[j] = ((i+j) % nbt);
      }

      p->dir_points_[0] = 0;
      p->dir_points_[1] = 1;
      p->dir_points_[2] = nbt - 1;

      psd.ind_polygon.add(nbt);
    }
    return p;
  }

  /* ******************************************************************** */
  /* direct product of convex structures                                  */
  /* ******************************************************************** */

  struct cv_pr_light_ {
    pconvex_structure cv1, cv2;
    bool operator < (const cv_pr_light_ &ls) const {
      if (cv1 < ls.cv1) return true; if (cv1 > ls.cv1) return false; 
      if (cv2 < ls.cv2) return true; return false;
    }
    cv_pr_light_(pconvex_structure a, pconvex_structure b) { cv1=a; cv2=b; }
    cv_pr_light_(void) { }
  };


  struct cv_pr_structure_ : public convex_structure
  {
    cv_pr_structure_(const cv_pr_light_ &ls)
    {
      Nc = ls.cv1->dim() + ls.cv2->dim();
      prod_a = ls.cv1; prod_b = ls.cv2;
      nbpt = ls.cv1->nb_points() * ls.cv2->nb_points();
      nbf = ls.cv1->nb_faces() + ls.cv2->nb_faces();
      if (ls.cv1->basic_structure() != ls.cv1
	  || ls.cv2->basic_structure() != ls.cv2)
	basic_pcvs = convex_product_structure(ls.cv1->basic_structure(),
					      ls.cv2->basic_structure());
      else
	basic_pcvs = this;
      faces_struct = std::vector<pconvex_structure>(nbf);
      faces = std::vector< std::vector<short_type> >(nbf);

      if (ls.cv1->ind_dir_points().size() && ls.cv2->ind_dir_points().size()) {
	dir_points_ = std::vector<short_type>(Nc + 1);

	for (int i = 0; i <= ls.cv1->dim(); i++)
	  dir_points_[i] = ls.cv1->ind_dir_points()[i]
	    + ls.cv2->ind_dir_points()[0] * ls.cv1->nb_points();
	for (int i = 1; i <= ls.cv2->dim(); i++)
	  dir_points_[ls.cv1->dim()+i] = ls.cv1->ind_dir_points()[0]
	    + ls.cv2->ind_dir_points()[i] * ls.cv1->nb_points();
      }

      for (int i = 0; i < ls.cv1->nb_faces(); i++)
      { 
	if (ls.cv1->nb_points_of_face(i) == 1)
	  faces_struct[i] = ls.cv2;
	else
	  faces_struct[i]
	    = (ls.cv1->faces_structure()[i] == NULL) ? NULL
	    : convex_product_structure(ls.cv1->faces_structure()[i], ls.cv2);

	faces[i] = std::vector<short_type>(ls.cv1->nb_points_of_face(i)
					      * ls.cv2->nb_points());

	for (int j = 0; j < ls.cv1->nb_points_of_face(i); j++)
	  for (int l = 0; l < ls.cv2->nb_points(); l++)
	  {
	    (faces[i])[l*ls.cv1->nb_points_of_face(i)+j]
	      = (ls.cv1->ind_points_of_face(i))[j] + l * ls.cv1->nb_points();
	  }
      }
      for (int i = 0; i < ls.cv2->nb_faces(); i++)
      { 
	int k = ls.cv1->nb_faces();
	if (ls.cv2->nb_points_of_face(i) == 1)
	  faces_struct[i+k] = ls.cv1;
	else
	  faces_struct[i+k]
	    = (ls.cv2->faces_structure()[i] == NULL) ? NULL
	    : convex_product_structure(ls.cv1, ls.cv2->faces_structure()[i]);

	faces[i+k] = std::vector<short_type>(ls.cv2->nb_points_of_face(i)
					      * ls.cv1->nb_points());

	for (int j = 0; j < ls.cv2->nb_points_of_face(i); j++)
	  for (int l = 0; l < ls.cv1->nb_points(); l++)
	  {
	    (faces[i+k])[j*ls.cv1->nb_points()+l]
	      = l + (ls.cv2->ind_points_of_face(i))[j] * ls.cv1->nb_points(); 
	  }
      }
    }
  };

  pconvex_structure convex_product_structure(pconvex_structure a,
					     pconvex_structure b) {
    return dal::singleton< dal::FONC_TABLE<cv_pr_light_, cv_pr_structure_> >
      ::instance().add(cv_pr_light_(a, b));
  }

  /* ******************************************************************** */
  /* parallelepiped structures.                                           */
  /* ******************************************************************** */

  struct parallelepiped_structure_data {
    dal::dynamic_array<pconvex_structure> tab;
    int nb_parallelepiped;
    parallelepiped_structure_data() : nb_parallelepiped(-1) {}
  };

  pconvex_structure parallelepiped_structure(dim_type nc)
  {
    parallelepiped_structure_data &psd = dal::singleton<parallelepiped_structure_data>::instance();
    if (nc > psd.nb_parallelepiped) {
      if (psd.nb_parallelepiped < 0) {
	psd.tab[0] = simplex_structure(0);
	psd.tab[1] = simplex_structure(1);
      }
      for (int i = 1; i < nc; i++)
	psd.tab[i+1]
	  = convex_product_structure(psd.tab[i], simplex_structure(1));
      psd.nb_parallelepiped = nc;
    }
    return psd.tab[nc];
  }

  // generic convex with n global nodes

  class dummy_structure_ : public convex_structure
  { friend pconvex_structure generic_dummy_structure(dim_type, size_type); };


  pconvex_structure generic_dummy_structure(dim_type nc, size_type n) {
    static std::vector< std::vector<pconvex_structure> > tab;
    if (size_type(nc)+1 > tab.size()) tab.resize(nc+1);
    if (n+1 > tab[nc].size()) {
      size_type d = tab[nc].size();
      tab[nc].resize(n+1);
      for (size_type i = d; i <= n; ++i) {
	dummy_structure_ *p = new dummy_structure_;
	p->Nc = nc; p->nbpt = n; p->nbf = 0;
	p->faces_struct.resize(0);
	p->faces.resize(0);
	p->dir_points_.resize(0);
	p->basic_pcvs = p;
	tab[nc][i] = p;
      }
    }
    return tab[nc][n];
  }

}  /* end of namespace bgeot.                                            */

