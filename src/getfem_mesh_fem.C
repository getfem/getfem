/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Elements Methods (getfem) version 1.0*/
/* File    :  getfem_mesh_fem.C : finite element methods on convex meshes. */
/*     									   */
/*                                                                         */
/* Date : December 21, 1999.                                               */
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


#include <queue>
#include <getfem_mesh_fem.h>

namespace getfem
{

  void boundary_description::add_elt(size_type c, short_type f)
  { 
    faces[  (cvindex[c]) ? cv_in.search(c) : cv_in.add(c)  ].add(f); 
    cvindex.add(c);
  }
      
  void boundary_description::sup_elt(size_type c, short_type f)
  {
    if (cvindex[c]) 
    { 
      size_type i = cv_in.search(c);
      faces[i].sup(f);
      if (faces[i].card() == 0) cvindex.sup(c);
    }
  }

  void boundary_description::sup_convex(size_type c)
  {
    if (cvindex[c]) 
    { 
      size_type i = cv_in.search(c); faces[i].clear();
      cvindex.sup(c); cv_in.sup(i);
    }
  }

  const dal::bit_vector &
    boundary_description::faces_of_convex(size_type c) const
  {
    static dal::bit_vector *without;
    static bool isinit = false;
    if (!isinit) {
      without = new dal::bit_vector(); isinit = true;
    }
    return (cvindex[c]) ? faces[cv_in.search(c)] : *without;
  }

  void boundary_description::swap_convex(size_type c1, size_type c2)
  {
    size_type i1, i2;
    dal::bit_vector b1, b2;
    if (cvindex[c1])
    { i1=cv_in.search(c1); b1=faces[i1]; faces[i1].clear(); cv_in.sup(i1); }
    if (cvindex[c2])
    { i2=cv_in.search(c2); b2=faces[i2]; faces[i2].clear(); cv_in.sup(i2); }
    if (cvindex[c1])
    { i1 = cv_in.add(c2);  faces[i1] = b1; }
    if (cvindex[c2])
    { i2 = cv_in.add(c1);  faces[i2] = b2; }
    cvindex.swap(c1, c2);
  }

  bool intfem::operator < (const intfem &l) const {
    if (pf < l.pf) return true; if (pf > l.pf) return false; 
    if (pi.method.ppi < l.pi.method.ppi) return true;
    return false;
  }

  pintfem give_intfem(pfem ppf, const pintegration_method &ppi)
  {
    static dal::FONC_TABLE<intfem, intfem> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<intfem, intfem>(); isinit = true;
    }
    return tab->add(intfem(ppf, ppi));
  }

  typedef bgeot::ref_mesh_point_ind_ct ref_mesh_dof_ind_ct;

  const dal::bit_vector &mesh_fem::convex_on_boundary(size_type b) const
  {
    static dal::bit_vector *without;
    static bool isinit = false;
    if (!isinit) {
      without = new dal::bit_vector(); isinit = true;
    }
    return (valid_boundaries[b]) ?  boundaries[b].cvindex : *without;
  }

  const dal::bit_vector &mesh_fem::faces_of_convex_on_boundary(size_type c,
							 size_type b) const
  {
    static dal::bit_vector *without;
    static bool isinit = false;
    if (!isinit) {
      without = new dal::bit_vector(); isinit = true;
    }
    return (valid_boundaries[b]) ? boundaries[b].faces_of_convex(c)
      : *without;
  }

  base_vector mesh_fem::normal_of_face_of_convex(size_type ic, short_type f, const base_node &pt) const
  {
    size_type N  = linked_mesh().dim();
    bgeot::pgeometric_trans pgt = linked_mesh().trans_of_convex(ic);
    base_matrix S(N,pgt->nb_points());
    
    for (size_type i=0; i < pgt->nb_points(); i++) {
      std::copy(linked_mesh().points_of_convex(ic)[i].begin(), 
		linked_mesh().points_of_convex(ic)[i].end(),
		S.begin()+i*N);
    }

    return compute_normal(S, f, pgt, pt);
  }

  void mesh_fem::sup_boundaries_of_convex(size_type c)
  {
    dal::bit_vector::iterator it = valid_boundaries.begin(),
      ite = valid_boundaries.end();
    for (; it != ite; ++it)
      if (*it) boundaries[it.index()].sup_convex(c);
  }

  void mesh_fem::swap_boundaries_convex(size_type c1, size_type c2)
  {
    dal::bit_vector::iterator it = valid_boundaries.begin(),
      ite = valid_boundaries.end();
    for (; it != ite; ++it)
      if (*it) boundaries[it.index()].swap_convex(c1, c2);
  }

  dal::bit_vector mesh_fem::dof_on_boundary(size_type b) const
  {
    dal::bit_vector res;
    if (valid_boundaries[b])
    {
      dal::bit_vector::const_iterator it = boundaries[b].cvindex.begin(),
	ite = boundaries[b].cvindex.end();
      for (; it != ite; ++it)
	if (*it)
	{
	  // cout << "scanning convex " << it.index() << endl;
	  // cout << " dofs : " << res << " ]] " << endl;
	dal::bit_vector::const_iterator
	  itf = boundaries[b].faces_of_convex(it.index()).begin(),
	  itfe = boundaries[b].faces_of_convex(it.index()).end();
	for (; itf != itfe; ++itf)
	  if (*itf)
	  {
	    // cout << "scanning face " <<  itf.index() << endl;
	    // cout << " dofs : " << res << " ]] " << endl;
	    size_type nbb // a refaire ... mais bug ...
	      = structure_of_convex(it.index())->nb_points_of_face(itf.index());
	    for (size_type i = 0; i < nbb; ++i)
	      res.add(ind_points_of_face_of_convex(it.index(), itf.index())[i]);
	    
// 	             bgeot::ind_ref_mesh_point_ind_ct::const_iterator
// 	     	    itp=ind_points_of_face_of_convex(it.index(), itf.index()).begin(),
// 	     	    itpe=ind_points_of_face_of_convex(it.index(), itf.index()).end();
// 	     	  for (; itp != itpe; ++itp)
// 	     	    { cout << "et hop " << endl; cout << " ajout de " << *itp << endl;  res.add(*itp); }
	  }
	}
    }
    return res;
  }

  void mesh_fem::receipt(const MESH_CLEAR &) { clear(); }
  void mesh_fem::receipt(const MESH_SUP_CONVEX &m)
  { 
    if (fe_convex[m.icv])
      { fe_convex[m.icv] = false; dof_enumeration_made = false; }
    sup_boundaries_of_convex(m.icv);
  }
  void mesh_fem::receipt(const MESH_SWAP_CONVEX &m)
  { 
    fe_convex.swap(m.icv1, m.icv2);
    f_elems.swap(m.icv1, m.icv2);
    swap_boundaries_convex(m.icv1, m.icv2);
  }
  void mesh_fem::receipt(const MESH_REFINE_CONVEX &m)
  { 
    // ajouter la strategie au rafinement / derafinement
    assert(false);
  }
  void mesh_fem::receipt(const MESH_UNREFINE_CONVEX &m)
  { 
    // ajouter la strategie au rafinement / derafinement
    assert(false);
  }
   
  void mesh_fem::set_finite_element(size_type cv, pintfem pif)
  {
    if (pif == NULL || pif->pf == NULL )
    {
      if (fe_convex.is_in(cv))
      { fe_convex.sup(cv); dof_enumeration_made = false; }
    }
    else
    {
      #ifdef __GETFEM_VERIFY
      assert(_linked_mesh->structure_of_convex(cv)->basic_structure() 
	     == pif->pf->basic_structure());
      #endif
      if (!fe_convex.is_in(cv) || f_elems[cv] != pif)
      { fe_convex.add(cv); f_elems[cv] = pif; dof_enumeration_made = false; }
    }
  }

  void mesh_fem::set_finite_element(const dal::bit_vector &cvs, pfem ppf,
			      const pintegration_method &ppi)
  { 
    dal::bit_vector::const_iterator it = cvs.begin(), ite = cvs.end();
    pintfem pif =  give_intfem(ppf, ppi);
    for ( ; it != ite; ++it) if (*it) set_finite_element(it.index(), pif);
  }

  base_node mesh_fem::point_of_dof(size_type cv, size_type i) const
  {
    pfem pf = f_elems[cv]->pf;
    bgeot::pgeometric_trans pgt = _linked_mesh->trans_of_convex(cv);
    const base_node *pt = &(pf->node_of_dof(i));
    base_node P(_linked_mesh->dim()); P.fill(0.0);

    size_type k = pgt->nb_points();
    for (size_type l = 0; l < k; ++l)
    {
      P.addmul(pgt->poly_vector()[l].eval(pt->begin()),
	       _linked_mesh->points_of_convex(cv)[l]);
      
    }
    return P;
  }

  base_node mesh_fem::point_of_dof(size_type d) const
  { return point_of_dof(points_tab[d].first, points_tab[d].ind_in_first); }

  struct _dof_comp
  { 
    dal::approx_less<scalar_type> comp;
    int operator()(const fem_dof& m, const fem_dof& n) const
    { 
      int d = dal::lexicographical_compare(m.P.begin(), m.P.end(),
					   n.P.begin(), n.P.end(), comp);
      if (d != 0) return d;
      return dof_description_compare(m.pnd, n.pnd);
    }
    _dof_comp(double e = 1.0E-10) : comp(e) { }
  };

  void mesh_fem::enumerate_dof(void)
  {
    dal::bit_vector nn = fe_convex;
    std::queue<int> pile;
    dal::dynamic_tree_sorted<fem_dof, _dof_comp, 10> dof_sort;

    size_type cv;
    std::vector<size_type> tab;
    fem_dof fd;

    size_type count = 0;
    cv = nn.take_first();
    bgeot::mesh_structure::clear();

    while (cv != ST_NIL)
    {
      /* ajout des voisins dans la pile.                                  */

      // cout << "nbd = " <<  nb_dof_of_element(cv) << endl; getchar();
      size_type nbp = _linked_mesh->nb_points_of_convex(cv);
     
      for (size_type i = 0; i < nbp; i++)
      {
	size_type ip = _linked_mesh->ind_points_of_convex(cv)[i];
	bgeot::mesh_convex_ind_ct::const_iterator 
	  it = _linked_mesh->convex_to_point(ip).begin(),
	  ite =  _linked_mesh->convex_to_point(ip).end();
	for ( ; it != ite; ++it)
	  if (nn.is_in(*it)) { nn.sup(*it); pile.push(*it); }
      }
      
      size_type nbd = nb_dof_of_element(cv);
      tab.resize(nbd);
      for (size_type i = 0; i < nbd; i++)
      {
	fd.P = point_of_dof(cv, i); // optimisable ...
	// cout << "point of dof : " << fd.P << endl;
	fd.pnd = f_elems[cv]->pf->dof_types()[i];
	size_type j;
	if (dof_linkable(fd.pnd)) {
	  j = dof_sort.add_norepeat(fd);
	} else
	  j = dof_sort.add(fd);
	count++;
	tab[i] = j;
      }

      size_type k = add_convex(
	              f_elems[cv]->pf->structure(),
		      tab.begin());
      if (k != cv) bgeot::mesh_structure::swap_convex(k, cv);
      
      if (pile.empty()) cv = nn.take_first();
                 else { cv = pile.front(); pile.pop(); }
    }
    
    dof_enumeration_made = true;
    nb_total_dof = dof_sort.card();
    
  }

/*     void mesh_fem::enumerate_dof(void) A finir ...  */
/*   { */
     /* L'algorithme deparcours est de type avancement par front           */
     /* Le tri des ddl prend pas mal de memoire ... faire autrement ?      */
     /* on peut deja ne pas stocker les noeuds interieurs.                 */
/*     dal::int_set nn = fe_convex; */
/*     std::queue<int> pile; */
/*     dal::dynamic_tree_sorted<_dof, _dof_comp, 10 > dof_sort; */
/*     dal::dynamic_array<int> front1, front2; */

/*     int *tab = new int[10000]; */
/*     int *degree = new int[nn.last_true() + 1]; */
/*     int cv, nbp, i, ip, ncv, nbd, fd, degree_min, first_cv; */
/*     _dof P; */

     /* 1-) On cherche le convexe de degree minimal (moins de voisins).    */

/*     degree_min = -1; */
/*     dal::int_set nm = nn; */
/*     for (cv << nm; cv >= 0; cv << nm)  il faut un iterateur sur nn.     */
/*     { */
/*       degree[cv] = 0; */
/*       nbp = _linked_mesh->nb_points_of_convex(cv); */
/*       for (i = 0; i < nbp; i++) */
/*       { */
/* 	ip = _linked_mesh->ind_point_of_convex(cv, i); */
/* 	nbd = _linked_mesh->convex_to_point(ip, tab); */
/* 	assert(nbd < 10000); */
/* 	for (i = 0; i < nbd; i++) */
/* 	{ */
/* 	  ncv = tab[i]; */
/* 	  if (nn.is_in(ncv)) degree[cv]++;  cette strategie compte les  */
 	  /* voisins avec un poids correspondant au nombre de points      */
 	  /* communs.                                                     */
/* 	} */
/*       } */
/*       if (degree_min == -1 && degree_min > degree[cv]) */
/* 	{ degree_min = degree[cv]; first_cv = cv; } */
/*     } */

     /* 2-) construction du premier front.                                 */
/*     if (degree_min != -1) */
/*     { front1[0] = first_cv; front1[1] = -2; nn[first_cv] = false; } */
/*     else */
/*       front1[0] = -2; */

     /* 3-) iterations sur les fronts.                                    */
    
/*     dof_list.clear(); */
/*     while (front1[0] != -2) */
/*     { */
/*       i = 0;  + boucle sur les fronts.        */
/*       j = 0; */
/*       while (front1[i] >= 0) */
/*       { */
/* 	cv = front1[i++]; */
/* 	nbp = _linked_mesh->nb_points_of_convex(cv); */
/* 	nbd = _linked_mesh->convex_to_point(ip, tab); */
/* 	assert(nbd < 10000); extraire une fonction qui ajoute les ddl d'un cv*/
/* 	for (i = 0; i < nbd; i++) */
/* 	{ */
/* 	  ncv = tab[i]; */
/* 	  if (nn.is_in(ncv)) { nn.sup(ncv); front2[j++] = ncv;  + tri  } */
/* 	} */
/*       } */
/*       if (j > 0 && front2[j-1] >= 0) front2[j++] = -1; */
/*     } */
/*   } */

  void mesh_fem::clear(void)
  {
    fe_convex.clear(); dof_enumeration_made = false;
    bgeot::mesh_structure::clear();
    boundaries.clear();
    valid_boundaries.clear();
  }

  mesh_fem::mesh_fem(getfem_mesh &me)
  {
    _linked_mesh = &me;
 
    add_sender(me.lmsg_sender(), *this,
	   lmsg::mask(MESH_CLEAR()) | lmsg::mask(MESH_SUP_CONVEX()) |
	   lmsg::mask(MESH_SWAP_CONVEX()) | lmsg::mask(MESH_REFINE_CONVEX()) |
           lmsg::mask(MESH_UNREFINE_CONVEX()));
  }

  mesh_fem::~mesh_fem(void)
  {
    sup_sender(_linked_mesh->lmsg_sender());
  }

}  /* end of namespace getfem.                                             */


