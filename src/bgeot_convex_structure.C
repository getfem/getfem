/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot) version 1.0                    */
/* File    :  bgeot_convex_structure.C : convex structures.                */
/*     									   */
/*                                                                         */
/* Date : December 20, 1999.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; version 2 of the License.                 */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software Foundation, */
/* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.         */
/*                                                                         */
/* *********************************************************************** */


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
    assert(nbpt >= i);
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
     _dir_points = convex_ind_ct();
    nbpt = 0;
  }


  /* ******************************************************************** */
  /* Numeration of convex structures                                      */
  /* ******************************************************************** */

  static dal::dynamic_tree_sorted<pconvex_structure> *_numeration;

  static inline void init_tab(void) // because of problem with initialization
  {                                 // in dynamic libraries.
    static bool initialized = false;
    if (!initialized)
    { 
      initialized = true;
      _numeration = new dal::dynamic_tree_sorted<pconvex_structure>();
    }
  }
  
  void numerate_convex_structure(pconvex_structure cvs)
  { init_tab(); _numeration->add(cvs); }

  pconvex_structure convex_structure_with_number(size_t i)
  { init_tab(); return (*_numeration)[i]; }

  size_t number_of_convex_structure(pconvex_structure cvs)
  { init_tab(); return _numeration->search(cvs); }

  STD_NEEDED ostream &operator <<(STD_NEEDED ostream &o, const convex_structure &cv)
  {
    o << "convex structure of dimension " << int(cv.dim()) << " with "
      << cv.nb_points() << " points and " << cv.nb_faces() << " faces " << endl;
    // a completer au besoin
    return o;
  }
    
  /* ******************************************************************** */
  /* simplex structures                                                   */
  /* ******************************************************************** */

  class _simplex_structure : public convex_structure
  {
    friend pconvex_structure simplex_structure(dim_type nc);
  };

  static int __nb_simplex = -1;

  pconvex_structure simplex_structure(dim_type nc)
  {
    static dal::dynamic_array<_simplex_structure> *_simplex;
    static int _nb_simplex = -1;
    static bool isinit = false;
    if (!isinit) {
      _simplex = new dal::dynamic_array<_simplex_structure>();
      isinit = true;
    }

    _simplex_structure *p = &(*_simplex)[nc];
    if (_nb_simplex < nc)
    {
      p->Nc = nc; p->nbpt = nc+1; p->nbf = nc+1;
      p->faces_struct.resize(p->nbf);
      p->faces.resize(p->nbf);
      p->_dir_points.resize(p->Nc + 1);
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
	  p->_dir_points[i] = i;
	  p->faces_struct[i] = simplex_structure(nc-1);
	  p->faces[i].resize(nc);
	  for (int j = 0; j < nc; j++)
	    (p->faces[i])[j] = (j >= i) ? (j + 1) : j;
	}
      
      numerate_convex_structure(p);
      __nb_simplex=_nb_simplex = nc;
    }
    return p;
  }

  /* ******************************************************************** */
  /* K-simplex structures                                                 */
  /* ******************************************************************** */

  struct _K_simplex_light
  {
    dim_type nc; short_type K;
    bool operator < (const _K_simplex_light &l) const
    {
      if (nc < l.nc) return true; if (nc > l.nc) return false; 
      if (K < l.K) return true; return false;
    }
    _K_simplex_light(dim_type n, short_type k) { nc = n; K = k; }
    _K_simplex_light(void) { }
  };

  class _K_simplex_structure : public convex_structure
  {
    public :

      _K_simplex_structure(const _K_simplex_light &ls)
      {
	Nc = ls.nc; nbpt = alpha(Nc, ls.K); nbf = Nc+1;
	basic_pcvs = simplex_structure(ls.nc);
	faces_struct.resize(nbf);
	faces.resize(nbf);
	_dir_points.resize(Nc+1);
	
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
	  _dir_points[pd++] = 0;
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
	      _dir_points[pd++] = r;
	  }
	}
	
	numerate_convex_structure(this);
      }

  };
  
  pconvex_structure simplex_structure(dim_type nc, short_type K)
  {
    static dal::FONC_TABLE<_K_simplex_light, _K_simplex_structure> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_K_simplex_light, _K_simplex_structure>();
      isinit = true;
    }
    if (nc == 0) return simplex_structure(0);
    if (K == 1) return simplex_structure(nc);
    return tab->add(_K_simplex_light(nc, K));
  }


  /* ******************************************************************** */
  /* polygon structures                                                   */
  /* ******************************************************************** */

  class _polygon_structure : public convex_structure
  {
    friend pconvex_structure polygon_structure(short_type nc);
  };

  static dal::bit_vector *_ind_polygon = NULL;

  pconvex_structure polygon_structure(short_type nbt)
  {
    static dal::dynamic_array<_polygon_structure> *_polygon;
    static bool initialized = false;

    if (!initialized) {
      initialized = true; _ind_polygon = new dal::bit_vector();
      _polygon = new dal::dynamic_array<_polygon_structure>();
    }

    _polygon_structure *p = &(*_polygon)[nbt];
    if (nbt < 4) return simplex_structure(nbt-1);

    if (!_ind_polygon->is_in(nbt))
    {
      p->Nc = 2; p->nbpt = nbt; p->nbf = nbt;
      p->basic_pcvs = p;
      p->faces_struct = std::vector<pconvex_structure>(p->nbf);
      p->faces = std::vector< std::vector<short_type> >(p->nbf);
      p->_dir_points = std::vector<short_type>(p->Nc + 1);

      for (int i = 0; i < p->nbf; i++)
      { 
	p->faces_struct[i] = simplex_structure(1);
	p->faces[i] = std::vector<short_type>(2);
	for (int j = 0; j < 2; j++)
	  (p->faces[i])[j] = ((i+j) % nbt);
      }

      p->_dir_points[0] = 0;
      p->_dir_points[1] = 1;
      p->_dir_points[2] = nbt - 1;

      numerate_convex_structure(p);
      _ind_polygon->add(nbt);
    }
    return p;
  }

  /* ******************************************************************** */
  /* direct product of convex structures                                  */
  /* ******************************************************************** */

  struct _cv_pr_light
  {
    pconvex_structure cv1, cv2;
    bool operator < (const _cv_pr_light &ls) const
    {
      if (cv1 < ls.cv1) return true; if (cv1 > ls.cv1) return false; 
      if (cv2 < ls.cv2) return true; return false;
    }
    _cv_pr_light(pconvex_structure a, pconvex_structure b) { cv1=a; cv2=b; }
    _cv_pr_light(void) { }
   
  };


  struct _cv_pr_structure : public convex_structure
  {
    _cv_pr_structure(const _cv_pr_light &ls)
    {
      Nc = ls.cv1->dim() + ls.cv2->dim();
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
      _dir_points = std::vector<short_type>(Nc + 1);

      for (int i = 0; i <= ls.cv1->dim(); i++)
	_dir_points[i] = ls.cv1->ind_dir_points()[i]
	                   + ls.cv2->ind_dir_points()[0] * ls.cv1->nb_points();
      for (int i = 1; i <= ls.cv2->dim(); i++)
	_dir_points[ls.cv1->dim()+i] = ls.cv1->ind_dir_points()[0]
	  + ls.cv2->ind_dir_points()[i] * ls.cv1->nb_points();

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
      numerate_convex_structure(this);
    }
  };

  static dal::FONC_TABLE<_cv_pr_light, _cv_pr_structure> *_cv_pr_tab;
  

  pconvex_structure convex_product_structure(pconvex_structure a,
					     pconvex_structure b)
  {
    static bool initialized = false;
    if (!initialized) 
    {
      initialized = true;
      _cv_pr_tab = new dal::FONC_TABLE<_cv_pr_light, _cv_pr_structure>();
    }
    return _cv_pr_tab->add(_cv_pr_light(a, b));
  }

  /* ******************************************************************** */
  /* parallelepiped structures.                                           */
  /* ******************************************************************** */

  pconvex_structure parallelepiped_structure(dim_type nc)
  {
     static dal::dynamic_array<pconvex_structure> *tab;
     static int _nb_parallelepiped = -1;
     static bool isinit = false;
     if (!isinit) {
       tab = new dal::dynamic_array<pconvex_structure>();
       isinit = true;
     }

    if (nc > _nb_parallelepiped) {
      if (_nb_parallelepiped < 0) {
	(*tab)[0] = simplex_structure(0);
	(*tab)[1] = simplex_structure(1);
      }
      for (int i = 1; i < nc; i++)
	(*tab)[i+1]
	  = convex_product_structure((*tab)[i], simplex_structure(1));
      _nb_parallelepiped = nc;
    }
    return (*tab)[nc];
  }


  /* ******************************************************************** */
  /* read-write convex structures to files.                               */
  /* ******************************************************************** */

  void write_convex_structures_to_file(STD_NEEDED ostream &ost)
  {
    int i;

    ost << endl << "BEGIN CONVEX STRUCTURE DESCRIPTION" << endl << endl;
    
    simplex_structure(1); // to be sure that initialization is done.
    for (i = 0; i <= __nb_simplex; i++ )
      ost << "CONVEX "
	  << number_of_convex_structure(simplex_structure(i))
	  << "\t= SIMPLEX " << i << endl;

    polygon_structure(2); // to be sure that initialization is done.
    dal::bit_vector nn = (*_ind_polygon);
    for(i << nn; i >= 0; i << nn)
      ost << "CONVEX "
	  << number_of_convex_structure(polygon_structure(i))
	  << "\t= POLYGON " << i << endl;
    
    convex_product_structure(simplex_structure(1), simplex_structure(1));
    nn = _cv_pr_tab->index();
    for(i << nn; i >= 0; i << nn)
      ost << "CONVEX "
	  << number_of_convex_structure(_cv_pr_tab->table()[i])
	  << "\t= PRODUCT "
	  << number_of_convex_structure(_cv_pr_tab->light_table()[i].cv1)
	  << " \t"
	  << number_of_convex_structure(_cv_pr_tab->light_table()[i].cv2)
	  << endl;

    ost << endl << "END CONVEX STRUCTURE DESCRIPTION" << endl;

  }

  int read_convex_structures_from_file(STD_NEEDED istream &ist,
			dal::dynamic_array<pconvex_structure> &tab)
  {
    bool te = false, please_get = true; bool error = false;
    char tmp[100];
    dal::bit_vector nn;

    ist.seekg(0);
    ftool::read_untill(ist, "BEGIN CONVEX STRUCTURE DESCRIPTION");

    while (!te)
    {
      if (please_get) ftool::get_token(ist, tmp, 99); else please_get = true;
      
      if (!strcmp(tmp, "END"))
      { te = true; }
      else if (!strcmp(tmp, "CONVEX"))
      {
	ftool::get_token(ist, tmp, 99);
	int ics = dal::abs(atoi(tmp));
	ftool::get_token(ist, tmp, 99);
	ftool::get_token(ist, tmp, 99);
	if (!strcmp(tmp, "SIMPLEX"))
	{
	  ftool::get_token(ist, tmp, 99);
	  int dim = atoi(tmp);
	  tab[ics] = simplex_structure(dim);
	  nn.add(ics);
	}
	else if (!strcmp(tmp, "POLYGON"))
	{
	  ftool::get_token(ist, tmp, 99);
	  int nbt = dal::abs(atoi(tmp));
	  tab[ics] = polygon_structure(nbt);
	  nn.add(ics);
	}
	else if (!strcmp(tmp, "PRODUCT"))
	{
	  ftool::get_token(ist, tmp, 99);
	  int cv1 = dal::abs(atoi(tmp));
	  ftool::get_token(ist, tmp, 99);
	  int cv2 = dal::abs(atoi(tmp));
	  if (nn.is_in(cv1) && nn.is_in(cv2))
	    tab[ics] = convex_product_structure(tab[cv1], tab[cv2]);
	  else
	    { error = true; te = true; }
	  nn.add(ics);
	}
	else
	{
	  error = true; te = true;
	}

      }
      else
      {
	error = true; te = true;
      }
    }

    if (error)
    {
      STD_NEEDED cerr << "BGEOT_CONVEX : Warning, file not valid." << endl;
      return -1;
    }
    return 0;
  }

}  /* end of namespace bgeot.                                            */

