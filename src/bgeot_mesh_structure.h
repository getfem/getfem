/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    : bgeot_mesh_structure.h : mesh structures.                     */
/*     									   */
/*                                                                         */
/* Date :  November 5, 1999.                                               */
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


#ifndef __BGEOT_MESH_STRUCTURE_H
#define __BGEOT_MESH_STRUCTURE_H

#include <bgeot_convex_structure.h>

namespace bgeot
{

  /* ******************************************************************** */
  /* transfert containers and iterators.                                  */
  /* ******************************************************************** */

  class mesh_structure;

  struct mesh_point          /* structure for a point.                    */
  {
    size_type first;         /* first convex attached to point or ST_NIL. */
    short_type ind_in_first; /* index of the point in first convex.       */
    mesh_point(void) { first = ST_NIL; }
    bool is_valid(void) const { return (first != ST_NIL); }
  };

  struct mesh_convex_structure
  {
    pconvex_structure cstruct;       /* type of convexe.                  */
    size_type pts;                   /* index in  points list.            */
    pconvex_structure structure(void) const { return cstruct; }
    pconvex_structure &structure(void) { return cstruct; }
  };

  struct mesh_point_link      /* link between convex and points.          */
  {
    size_type next;           /* index of next convex linked to the point.*/
    short_type ind_in_next;   /* index of this point in next convex.      */
    mesh_point_link(void) { next = ST_NIL; }
  };

  typedef dal::dynamic_array<mesh_point, 8>          mesh_point_st_ct;
  typedef dal::dynamic_tas<mesh_convex_structure, 8> mesh_convex_ct;
  typedef dal::dynamic_alloc<mesh_point_link, 8>     mesh_link_ct;
  typedef dal::dynamic_array<size_type, 8>           mesh_point_ind_ct;
  typedef dal::tab_ref<mesh_point_ind_ct::const_iterator>
                                                     ref_mesh_point_ind_ct;
  typedef std::vector<size_type>                     mesh_point_search_ind_ct;
  typedef dal::tab_ref_index_ref<mesh_point_ind_ct::const_iterator,
                                 convex_ind_ct::const_iterator>
                                                    ind_ref_mesh_point_ind_ct;

  /* ********************************************************************* */
  /* Container : convex attached to a point.                               */
  /* ********************************************************************* */

  struct mc_const_iterator
  {
    typedef size_t              value_type;
    typedef const value_type*   pointer;
    typedef const value_type&   reference;
    typedef size_t              size_type;
    typedef ptrdiff_t           difference_type;
    typedef std::forward_iterator_tag iterator_category;
    
    const mesh_structure *p;
    const mesh_convex_structure *pc;
    size_type ind_cv;
    short_type ind_in_cv;
    
    mc_const_iterator(void) { ind_cv = ST_NIL; }
    mc_const_iterator(const mesh_structure &ms, size_type ip);
	
    mc_const_iterator &operator ++();
    
    mc_const_iterator operator ++(int)
    { mc_const_iterator tmp = *this; ++(*this); return tmp; }
    
    template<class ITER> bool is_convex_has_points(short_type nb,
						   const ITER &ipt) const;
    
    bool is_end(void) const { return (ind_cv == ST_NIL); }
    value_type operator *() const { return ind_cv; }
    pconvex_structure structure(void) const
    { return pc->structure(); }
    
    bool operator ==(const mc_const_iterator &i) const
    { return (i.ind_cv == ind_cv);}
    bool operator !=(const mc_const_iterator &i) const
    { return (i.ind_cv != ind_cv);}
  };

  class mesh_convex_ind_ct
  {
    public :
      
      typedef size_t               value_type;
      typedef value_type*          pointer;
      typedef const value_type*    const_pointer;
      typedef value_type&          reference;
      typedef const value_type&    const_reference;
      typedef size_t               size_type;
      typedef mc_const_iterator    const_iterator;
      typedef mc_const_iterator    iterator;
      
    protected :

      iterator _begin, _end;

    public :

      bool empty(void) const { return _begin == _end; }
    
      mesh_convex_ind_ct(const mesh_structure &ms, size_type ip)
      { _begin = iterator(ms, ip); }

      const iterator &begin() const { return _begin; }
      const iterator &end() const { return _end; }
      iterator &begin() { return _begin; }
      iterator &end() { return _end; }
      
      value_type front(void) const { return *_begin; }
      void pop_front(void) { ++_begin; }
  };

  /* ********************************************************************* */
  /* mesh_structure.                                                       */
  /* ********************************************************************* */


  class mesh_structure
  {
    protected :

      mesh_point_st_ct    points_tab;
      mesh_convex_ct      convex_tab;
      mesh_point_ind_ct   point_lists;
      mesh_link_ct        point_links;

    public :

      mesh_link_ct &links(void) { return point_links; }
      mesh_convex_ct &convex(void) { return convex_tab; }
      mesh_point_st_ct &points(void) { return points_tab; }

      const dal::bit_vector &convex_index(void) const
      { return convex_tab.index(); }
      dal::bit_vector convex_index(dim_type) const;
      const mesh_point_st_ct &point_structures(void) const
      { return points_tab; }
      const mesh_link_ct  &links(void) const { return point_links; }
      const mesh_convex_ct &convex(void) const { return convex_tab; }

      void swap_points(size_type i, size_type j);
      size_type nb_convex(void) const { return convex_tab.card(); }
      bool convex_is_valid(size_type i) { return (convex_tab.index())[i]; }
      bool point_is_valid(size_type i)
      { return (points_tab[i].first != size_type(-1)); }
      ref_mesh_point_ind_ct ind_points_of_convex(size_type ic) const;
      size_type ind_dir_point_of_convex(size_type ic, size_type j) const
      {
	return ind_points_of_convex(ic)
	  [convex_tab[ic].cstruct->ind_dir_points()[j] ];
      }

      pconvex_structure structure_of_convex(size_type i) const
      { return convex_tab[i].cstruct; }
      short_type nb_points_of_convex(size_type i) const
      { return convex_tab[i].cstruct->nb_points(); }
      void swap_convex(size_type i, size_type j);

      template<class ITER>
	size_type add_convex(pconvex_structure cs,
			     ITER ipts, bool *present = NULL);
      template<class ITER> size_type add_simplex(dim_type dim, ITER ipts)
      { return add_convex(simplex_structure(dim), ipts); }
      size_type add_segment(size_type a, size_type b);

      void sup_convex(size_type ic);
      template<class ITER> 
	void sup_convex_with_points(ITER ipts, short_type nb);
      void sup_segment(size_type a, size_type b)
      { size_type t[2]; t[0] = a; t[1] = b; sup_convex_with_points(&t[0], 2); }
      size_type add_face_of_convex(size_type ic, short_type f);
      void add_faces_of_convex(size_type ic);
      void to_faces(dim_type n);
      void to_edges(void);

      mesh_convex_ind_ct convex_to_point(size_type ip) const
      { return mesh_convex_ind_ct(*this, ip); }
      
      mesh_point_search_ind_ct ind_points_to_point(size_type ip) const;
      
      template<class ITER>
	bool is_convex_has_points(size_type ic,short_type nb, ITER pit) const;

      template<class ITER> 
        bool is_convex_face_has_points(size_type ic, size_type face_num,
                                  short_type nb, ITER pit) const;
      ind_ref_mesh_point_ind_ct ind_points_of_face_of_convex(size_type ic,
							 short_type f) const;

      size_type memsize(void) const;
      void optimize_structure(void);
      void clear(void);
      void stat(void);
      int write_to_file(std::ostream &ost) const;
      int read_from_file(std::istream &ist);
     
  };

  /* ******************************************************************** */
  /* template member function of search iterator.                         */
  /* ******************************************************************** */
    
  template<class ITER> bool mc_const_iterator::is_convex_has_points
                                       (short_type nb, const ITER &ipt) const
  { return p->is_convex_has_points(this->ind_cv, nb, ipt); }


  /* ******************************************************************** */
  /* search functions.                                                    */
  /* ******************************************************************** */

  template<class ITER> struct mesh_convex_has_points
  {
    ITER ipt;
    short_type nb;
    dim_type dim;
    bool tnb;
    size_type ict;

    bool operator() (const mesh_convex_ind_ct::iterator &it) const
    { 
      if (it.is_end()) return true;
      return (nb < it.structure()->nb_points()
	      && *it != ict
	      && (!tnb || nb == it.structure()->nb_points())
	      && (dim == dim_type(-1) || dim == it.structure()->dim())
	      && it.is_convex_has_points(nb, ipt));
    }
    mesh_convex_has_points(short_type n, ITER ip,
	     dim_type d = dim_type(-1), bool t = false, size_type ic = ST_NIL)
    { nb = n; ipt = ip; dim = d; tnb = t; ict = ic; }
  };

  template<class ITER> class mesh_convex_with_points_ind_ct
    : public dal::tab_ref_with_selection<mesh_convex_ind_ct::iterator,
                                         mesh_convex_has_points<ITER> >
  {
    public :

      mesh_convex_with_points_ind_ct(const mesh_convex_ind_ct::iterator &b,
				     const mesh_convex_ind_ct::iterator &e,
				     const mesh_convex_has_points<ITER> &c)
      : dal::tab_ref_with_selection<mesh_convex_ind_ct::iterator,
                                    mesh_convex_has_points<ITER> >(b,e,c)
      {}
  }; 

  template<class ITER> mesh_convex_with_points_ind_ct<ITER>
    convex_with_points(const mesh_structure &ms, short_type nb, ITER ipts)
  {
    mesh_convex_ind_ct ct = ms.convex_to_point(*ipts);
    return mesh_convex_with_points_ind_ct<ITER>(ct.begin(), ct.end(),
				  mesh_convex_has_points<ITER>(--nb, ++ipts));
  }

  typedef mesh_convex_with_points_ind_ct<ref_mesh_point_ind_ct::const_iterator>
    mesh_over_convex_ind_ct;

  mesh_over_convex_ind_ct over_convex(const mesh_structure &ms, size_type ic);

  typedef mesh_convex_with_points_ind_ct<
    ind_ref_mesh_point_ind_ct::const_iterator> mesh_face_convex_ind_ct;

  mesh_face_convex_ind_ct face_of_convex(const mesh_structure &ms, 
					 size_type ic, short_type iff);

  mesh_face_convex_ind_ct neighbour_of_convex(const mesh_structure &ms, 
					      size_type ic, short_type iff);

  /* ******************************************************************** */
  /* template member functions of mesh_structure.                         */
  /* ******************************************************************** */

  template<class ITER>
    bool mesh_structure::is_convex_has_points(size_type ic,
					      short_type nb, ITER pit) const
  {
    for (short_type i = 0; i < nb; ++i, ++pit)
    {
      ref_mesh_point_ind_ct pt = ind_points_of_convex(ic);
      if (std::find(pt.begin(), pt.end(), *pit) == pt.end()) return false;
    }
    return true;
  }
  

  template<class ITER>
  bool mesh_structure::is_convex_face_has_points(size_type ic, size_type face_num,
						 short_type nb, ITER pit) const
  {
    for (short_type i = 0; i < nb; ++i, ++pit)
    {
      ind_ref_mesh_point_ind_ct pt = ind_points_of_face_of_convex(ic, face_num);
      if (std::find(pt.begin(), pt.end(), *pit) == pt.end()) return false;
    }
    return true;
  }


  template<class ITER>
    size_type mesh_structure::add_convex(pconvex_structure cs,
					 ITER ipts, bool *present)
  {
    short_type nb = cs->nb_points();
    mesh_convex_with_points_ind_ct<ITER>
      ct = convex_with_points(*this, nb, ipts);
    size_type is = *(ct.begin());
    if (present != NULL) *present = false;

    if (!ct.empty() && structure_of_convex(is) == cs)
    { if (present != NULL) *present = true; return is; }
    
    mesh_convex_structure s; s.cstruct = cs; s.pts = point_links.alloc(nb);
    if (nb > 0)
      { point_lists[s.pts+nb-1] = 0; point_links[s.pts+nb-1].next = 0; }
    dal::copy_n(ipts, nb, point_lists.begin()+s.pts);
    is = convex_tab.add(s);

    mesh_link_ct::iterator ipl = point_links.begin(); ipl += s.pts;
    for (short_type i = 0; i < nb; ++i, ++ipts, ++ipl)
    {
      mesh_point *os = &(points_tab[*ipts]);
      (*ipl).next = os->first; os->first = is;
      (*ipl).ind_in_next = os->ind_in_first; os->ind_in_first = i;
    }
    
    return is;
  }

  template<class ITER>
    void mesh_structure::sup_convex_with_points(ITER ipts, short_type nb)
  {
    mesh_convex_with_points_ind_ct<ITER>
      ct = convex_with_points(*this, nb, ipts);
    typename mesh_convex_with_points_ind_ct<ITER>::iterator b = ct.begin();
    typename mesh_convex_with_points_ind_ct<ITER>::iterator e = ct.end();
    for ( ; b != e; ++b) sup_convex(*b);
  }



}  /* end of namespace bgeot.                                              */


#endif /* __BGEOT_MESH_STRUCTURE_H                                         */
