/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot) version 1.0                    */
/* File    :  bgeot_mesh_structure.C : mesh structures.                    */
/*     									   */
/*                                                                         */
/* Date : November 5, 1999 .                                               */
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


#include <bgeot_mesh_structure.h>

namespace bgeot
{
  mc_const_iterator::mc_const_iterator(const mesh_structure &ms,
				       size_type ip)
  { 
    p = &ms; const mesh_point *q = &(p->point_structures()[ip]);
    ind_cv = q->first; ind_in_cv = q->ind_in_first;
    pc = &(p->convex()[ind_cv]); 
  }
	
  mc_const_iterator &mc_const_iterator::operator ++()
  { 
    const mesh_point_link *q = &(p->links()[pc->pts+ind_in_cv]);
    ind_cv = q->next; ind_in_cv = q->ind_in_next;
    pc = &(p->convex()[ind_cv]); return *this;
  }
  
  dal::bit_vector mesh_structure::convex_index(dim_type n) const
  {
    dal::bit_vector res = convex_tab.index();
    dal::bit_vector::iterator it = res.begin(), ite = res.end();
    for ( ; it != ite; ++it)
      if (*it && structure_of_convex(it.index())->dim() != n) *it = false;
    return res;
  }

  ind_ref_mesh_point_ind_ct
     mesh_structure::ind_points_of_face_of_convex(size_type ic,
					      short_type iff) const
  {
    mesh_point_ind_ct::const_iterator r = point_lists.begin();
    const mesh_convex_structure *q = &(convex_tab[ic]); r += q->pts;
    const convex_ind_ct *p = &(q->cstruct->ind_points_of_face(iff));
    return ind_ref_mesh_point_ind_ct(r, p->begin(), p->end());  
  }

  ref_mesh_point_ind_ct mesh_structure::ind_points_of_convex(size_type ic) const
  {
    const mesh_convex_structure *q = &(convex_tab[ic]);
    mesh_point_ind_ct::const_iterator r = point_lists.begin();
    r += q->pts;
    return ref_mesh_point_ind_ct(r, r + q->cstruct->nb_points());
  }

  size_type mesh_structure::memsize(void) const
  {
    return sizeof(mesh_structure) + points_tab.memsize()
      + convex_tab.memsize() + point_lists.memsize()
      + point_links.memsize();
  }

  void mesh_structure::swap_points(size_type i, size_type j)
  {
    if (i != j)
    {
      mesh_convex_ind_ct ct = convex_to_point(i);
      for ( ; !ct.empty(); ct.pop_front())
	point_lists[ct.begin().pc->pts + ct.begin().ind_in_cv] = i;
      ct = convex_to_point(j);
      for ( ; !ct.empty(); ct.pop_front())
	point_lists[ct.begin().pc->pts + ct.begin().ind_in_cv] = j;
      points_tab.swap(i,j);
    }			 
  }

  void mesh_structure::sup_convex(size_type ic)
  {
    // STD_NEEDED cout << "sup convex" << ic << " !! " << endl;
    mesh_convex_structure *p = &(convex_tab[ic]);
    short_type nb = p->cstruct->nb_points();
    mesh_link_ct::iterator ipl = point_links.begin(); ipl += p->pts;
    mesh_point_ind_ct::const_iterator ipt = point_lists.begin(); ipt += p->pts;
    for (short_type i = 0; i < nb; ++i, ++ipl, ++ipt)
    {
      mesh_point *os = &(points_tab[*ipt]);
      mesh_convex_ind_ct::iterator b = convex_to_point(*ipt).begin(), c = b;
      if (*b == ic)
      { os->first = (*ipl).next; os->ind_in_first = (*ipl).ind_in_next; }
      else
      {
	++b; while (*b != ic) { c = b; ++b; }
	mesh_point_link *ssp = &(point_links[c.pc->pts + c.ind_in_cv]);
	ssp->next = (*ipl).next; ssp->ind_in_next = (*ipl).ind_in_next;
      }
    }
    point_links.free(p->pts, nb); convex_tab.sup(ic);
  }

  mesh_point_search_ind_ct
    mesh_structure::ind_points_to_point(size_type ip) const
  { // pas tres efficace.
    mesh_point_search_ind_ct r;
    dal::bit_vector nn;
    mesh_convex_ind_ct ct = convex_to_point(ip);
    for ( ; !ct.empty(); ct.pop_front())
    {
      ref_mesh_point_ind_ct pt = ind_points_of_convex(ct.front());
      for ( ; !pt.empty(); pt.pop_front())
      { if (ip != pt.front()) nn.add(pt.front()); }
    }
    r.resize(nn.card());
    size_type i, j;
    for (i = 0, j << nn; j != size_type(-1); j << nn, ++i) r[i] = j;
    return r;
  }

  size_type mesh_structure::add_face_of_convex(size_type ic, short_type f)
  {
    return add_convex( (structure_of_convex(ic)->faces_structure())[f],
		  ind_points_of_face_of_convex(ic, f).begin());
  }

  void mesh_structure::add_faces_of_convex(size_type ic)
  {
    // STD_NEEDED cout << "debut de add_faces(" << ic << ")" << endl;
    // write_to_file(STD_NEEDED cout); 
    pconvex_structure ps = structure_of_convex(ic);
    for (short_type iff = 0; iff < ps->nb_faces(); ++iff)
      add_convex( (ps->faces_structure())[iff],
		  ind_points_of_face_of_convex(ic, iff).begin());
    // STD_NEEDED cout << "fin de add_faces " << endl;
    // write_to_file(STD_NEEDED cout); getchar();
  }
    
  void mesh_structure::to_faces(dim_type n)
  {
    // STD_NEEDED cout << "debut de to_faces(" << int(n) << ")" << endl;
    mesh_convex_ct::tas_iterator b = convex_tab.tas_begin(),
                                 e = convex_tab.tas_end();
    for ( ; b != e; ++b)
      if ((*b).cstruct->dim() == n)
	{ add_faces_of_convex(b.index()); sup_convex(b.index()); }
    // STD_NEEDED cout << "fin de to_faces " << endl;
    // write_to_file(STD_NEEDED cout); getchar();
  }

  void mesh_structure::to_edges(void)
  { // inefficace, a refaire
    dim_type dmax = 0;
    mesh_convex_ct::tas_iterator b = convex_tab.tas_begin(),
                                 e = convex_tab.tas_end();
    for ( ; b != e; ++b) dmax = std::max(dmax, (*b).cstruct->dim());
    for ( ; dmax > 1; --dmax) to_faces(dmax);
  }

  void mesh_structure::swap_convex(size_type i, size_type j)
  {
    size_type ip, ic, *p;

    if (i != j)
    {
      if (convex_is_valid(i))
      {
	ref_mesh_point_ind_ct pt = ind_points_of_convex(i);
	for( ; !pt.empty(); pt.pop_front())
	{
	  ip = pt.front(); p = &(points_tab[ip].first);
	  mesh_convex_ind_ct ct = convex_to_point(ip);

	  for( ; !ct.empty(); ct.pop_front())
	  {
	    ic = ct.front();
	    if (*p == i) *p = j; else if (*p == j) *p = i;
	    p = &(point_links[convex_tab[ic].pts+ct.begin().ind_in_cv].next);
	  }
	}
      }
      if (convex_is_valid(j))
      {
	ref_mesh_point_ind_ct pt = ind_points_of_convex(j);
	for( ; !pt.empty(); pt.pop_front())
	{
	  ip = pt.front();
	  if (!convex_is_valid(i) || !is_convex_has_points( i, 1, &ip))
	  {
	    p = &(points_tab[ip].first);
	    mesh_convex_ind_ct ct = convex_to_point(ip);

	    for( ; !ct.empty(); ct.pop_front())
	    {
	      ic = ct.front();
	      if (*p == i) *p = j; else if (*p == j) *p = i;
	      p = &(point_links[convex_tab[ic].pts+ct.begin().ind_in_cv].next);
	    }
	  }
	}
      }
      convex_tab.swap(i,j);
    }
  }

  void mesh_structure::clear(void)
  { 
    points_tab.clear(); convex_tab.clear();
    point_links.clear(); point_lists.clear();
  }

  void mesh_structure::optimize_structure()
  {
    size_type i, j;
    j = nb_convex();
    for (i = 0; i < j; i++)
      if (!convex_tab.index_valid(i))
	swap_convex(i, convex_tab.ind_last());

    for (i = 0, j = (points_tab.end()-points_tab.begin())-1; i < j; ++i, --j)
    {
      while (i < j && points_tab[i].first != ST_NIL) ++i;
      while (i < j && points_tab[j].first == ST_NIL) --j;
      if (i < j) swap_points(i, j);
    }
  }

  void mesh_structure::stat(void)
  {
    STD_NEEDED cout << "mesh structure with " << nb_convex() << " valid convex, "
	 << "for a total memory size of " << memsize() << " bytes.\n";
  }

  int mesh_structure::read_from_file(STD_NEEDED istream &ist)
  {
    char tmp[100]; 
    bool tend = false, error = false;
    dal::dynamic_alloc<size_type> cv_pt;
    dal::dynamic_array<mesh_convex_structure> cv;
    dal::dynamic_array<pconvex_structure> cvs;
    dal::bit_vector ncv;
    int ic;

    if (read_convex_structures_from_file(ist, cvs) != 0) return -1;
    
    ist.seekg(0);
    if (!ftool::read_untill(ist, "BEGIN MESH STRUCTURE DESCRIPTION"))
      return -1;

    while (!tend)
    {
      tend = !ftool::get_token(ist, tmp, 99);
      if (!strcmp(tmp, "END"))
      { tend = true; }
      else if (!strcmp(tmp, "CONVEX"))
      {
	ftool::get_token(ist, tmp, 99);
        ic = atoi(tmp);
	if (ncv.is_in(ic))
	{
	  STD_NEEDED cerr << "CVMH : Fatal Error, negativ index or two convex with"
	       << "the same index. loading aborted" << endl;
	  return -1;
	}
	ncv.add(ic);
	ftool::get_token(ist, tmp, 99);
        size_type id = dal::abs(atoi(tmp));
	size_type nb = cvs[id]->nb_points();

	cv[ic].cstruct = cvs[id];
	cv[ic].pts = cv_pt.alloc(nb);
	for (size_type i = 0; i < nb; i++)
	{
	  ftool::get_token(ist, tmp, 99);
	  cv_pt[cv[ic].pts+i] = dal::abs(atoi(tmp));
	}
      }
      else
      { STD_NEEDED cerr << "CVMH : Syntax error reading mesh file" << endl; return -1; }
    }

    for (ic << ncv; ic >= 0; ic << ncv)
    {
      size_type i = add_convex(cv[ic].cstruct, cv_pt.begin() + cv[ic].pts);
      if (i != ic) swap_convex(i, ic);
    }

    return 0;
  }

  template<class ITER>
    static void _write_tab_to_file(STD_NEEDED ostream &ost, ITER b, ITER e)
  { for ( ; b != e; ++b) ost << "  " << *b; }

  template<class ITER>
    static void _write_convex_to_file(const mesh_structure &ms, STD_NEEDED ostream &ost,
				      ITER b, ITER e)
  {
    for ( ; b != e; ++b)
    {
      size_type i = b.index();
      ost << "CONVEX " << i << "    "
	  << number_of_convex_structure(ms.structure_of_convex(i)) << "    ";
      _write_tab_to_file(ost, ms.ind_points_of_convex(i).begin(),
			      ms.ind_points_of_convex(i).end()  );
      ost << endl;
    }
  }

  int mesh_structure::write_to_file(STD_NEEDED ostream &ost) const
  {
    write_convex_structures_to_file(ost);
    ost << endl << "BEGIN MESH STRUCTURE DESCRIPTION" << endl << endl;
    _write_convex_to_file(*this, ost, convex_tab.tas_begin(),
			              convex_tab.tas_end());
    ost << endl << "END MESH STRUCTURE DESCRIPTION" << endl;
    return 0; /* errors not suported */
  }

  mesh_over_convex_ind_ct over_convex(const mesh_structure &ms, size_type ic)
  {
    ref_mesh_point_ind_ct pt = ms.ind_points_of_convex(ic);
    mesh_convex_ind_ct ct = ms.convex_to_point(pt[0]);
    return mesh_over_convex_ind_ct(ct.begin(), ct.end(),
       mesh_convex_has_points<ref_mesh_point_ind_ct::const_iterator>
       (ms.structure_of_convex(ic)->nb_points() - 1, pt.begin() + 1,
	ms.structure_of_convex(ic)->dim()+1));
  }

  mesh_face_convex_ind_ct face_of_convex(const mesh_structure &ms, 
					 size_type ic, short_type iff)
  {
    pconvex_structure q = ms.structure_of_convex(ic);
    short_type n = q->nb_points_of_face(iff);
    ind_ref_mesh_point_ind_ct pt = ms.ind_points_of_face_of_convex(ic, iff);
    mesh_convex_ind_ct ct = ms.convex_to_point(pt[0]);
    return mesh_face_convex_ind_ct(ct.begin(), ct.end(), mesh_convex_has_points<ind_ref_mesh_point_ind_ct::const_iterator>(n-1, ++(pt.begin()), q->dim()-1, true));
  }

  mesh_face_convex_ind_ct neighbour_of_convex(const mesh_structure &ms, 
					  size_type ic, short_type iff)
  {
    pconvex_structure q = ms.structure_of_convex(ic);
    short_type n = q->nb_points_of_face(iff);
    ind_ref_mesh_point_ind_ct pt = ms.ind_points_of_face_of_convex(ic, iff);
    mesh_convex_ind_ct ct = ms.convex_to_point(pt[0]);
    return mesh_face_convex_ind_ct(ct.begin(), ct.end(), mesh_convex_has_points<ind_ref_mesh_point_ind_ct::const_iterator>(n-1, ++(pt.begin()), q->dim(), false, ic));
  }

}  /* end of namespace bgeot.                                              */
