/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Elements Methods (getfem)            */
/* File    :  getfem_mesh_fem.C : finite element methods on convex meshes. */
/*     									   */
/*                                                                         */
/* Date : December 21, 1999.                                               */
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


#include <queue>
#include <getfem_mesh_fem.h>

namespace getfem
{

  void boundary_description::add_elt(size_type c, short_type f) { 
    faces[  (cvindex[c]) ? cv_in.search(c) : cv_in.add(c)  ].add(f); 
    cvindex.add(c);
  }
  
  void boundary_description::sup_elt(size_type c, short_type f) {
    if (cvindex[c]) { 
      size_type i = cv_in.search(c);
      faces[i].sup(f);
      if (faces[i].card() == 0) cvindex.sup(c);
    }
  }

  void boundary_description::sup_convex(size_type c) {
    if (cvindex[c]) { 
      size_type i = cv_in.search(c); faces[i].clear();
      cvindex.sup(c); cv_in.sup(i);
    }
  }
  
  const dal::bit_vector &boundary_description::faces_of_convex(size_type c)
    const {
    static dal::bit_vector *without;
    static bool isinit = false;
    if (!isinit) {
      without = new dal::bit_vector(); isinit = true;
    }
    return (cvindex[c]) ? faces[cv_in.search(c)] : *without;
  }

  void boundary_description::swap_convex(size_type c1, size_type c2) {
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
    if (pi < l.pi) return true;
    return false;
  }
  
  pintfem give_intfem(pfem ppf, const pintegration_method ppi) {
    static dal::FONC_TABLE<intfem, intfem> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<intfem, intfem>(); isinit = true;
    }
    if (ppf->basic_structure() != ppi->structure())
      DAL_THROW(internal_error, 
		"Incompatibility between fem and integration method");
    return tab->add(intfem(ppf, ppi));
  }
  
  const dal::bit_vector &mesh_fem::convex_on_boundary(size_type b) const {
    static dal::bit_vector *without;
    static bool isinit = false;
    if (!isinit) {
      without = new dal::bit_vector(); isinit = true;
    }
    return (valid_boundaries[b]) ?  boundaries[b].cvindex : *without;
  }
  
  const dal::bit_vector &mesh_fem::faces_of_convex_on_boundary(size_type c,
							  size_type b) const {
    static dal::bit_vector *without;
    static bool isinit = false;
    if (!isinit) {
      without = new dal::bit_vector(); isinit = true;
    }
    return (valid_boundaries[b]) ? boundaries[b].faces_of_convex(c)
      : *without;
  }
  
  void mesh_fem::sup_boundaries_of_convex(size_type c) {
    dal::bit_vector::iterator it = valid_boundaries.begin(),
      ite = valid_boundaries.end();
    for (; it != ite; ++it)
      if (*it) boundaries[it.index()].sup_convex(c);
  }
  
  void mesh_fem::swap_boundaries_convex(size_type c1, size_type c2) {
    dal::bit_vector::iterator it = valid_boundaries.begin(),
      ite = valid_boundaries.end();
    for (; it != ite; ++it)
      if (*it) boundaries[it.index()].swap_convex(c1, c2);
  }
  
  dal::bit_vector mesh_fem::dof_on_boundary(size_type b) const {
    if (!dof_enumeration_made) this->enumerate_dof();
    dal::bit_vector res;
    if (valid_boundaries[b]) {
      dal::bit_vector::const_iterator it = boundaries[b].cvindex.begin(),
	ite = boundaries[b].cvindex.end();
      for (; it != ite; ++it)
	if (*it) {
	  dal::bit_vector::const_iterator
	    itf = boundaries[b].faces_of_convex(it.index()).begin(),
	    itfe = boundaries[b].faces_of_convex(it.index()).end();
	  for (; itf != itfe; ++itf)
	    if (*itf) {
	      size_type nbb
		= dof_structure.structure_of_convex(it.index())->
		nb_points_of_face(itf.index());
	      for (size_type i = 0; i < nbb; ++i) {
		size_type n = Qdim / fem_of_element(it.index())->target_dim();
		for (size_type ll = 0; ll < n; ++ll)
		  res.add(n * dof_structure.ind_points_of_face_of_convex(it.index(),
									 itf.index())[i] + ll);
	      }
	    }
	}
    }
    return res;
  }
  
  void mesh_fem::receipt(const MESH_CLEAR &) { clear(); }
  void mesh_fem::receipt(const MESH_DELETE &) {
    clear(); is_valid = false;
    linked_mesh().lmsg_sender().send(MESH_FEM_DELETE((void *)(this)));
    sup_sender(_linked_mesh->lmsg_sender());
  }
  void mesh_fem::receipt(const MESH_SUP_CONVEX &m) { 
    if (fe_convex[m.icv])
      { fe_convex[m.icv] = false; dof_enumeration_made = false; }
    sup_boundaries_of_convex(m.icv);
  }
  void mesh_fem::receipt(const MESH_SWAP_CONVEX &m) { 
    fe_convex.swap(m.icv1, m.icv2);
    f_elems.swap(m.icv1, m.icv2);
    swap_boundaries_convex(m.icv1, m.icv2);
  }
  void mesh_fem::receipt(const MESH_REFINE_CONVEX &) { 
    // ajouter la strategie au rafinement / derafinement
    DAL_THROW(internal_error, "internal error");
  }
  void mesh_fem::receipt(const MESH_UNREFINE_CONVEX &) { 
    // ajouter la strategie au rafinement / derafinement
    DAL_THROW(internal_error, "internal error");
  }
  void mesh_fem::receipt(const MESH_FEM_TOUCH &m) {
    if (m.ptr == (void *)(this)) { 
      dof_enumeration_made = false;
      linked_mesh().lmsg_sender().send(MESH_FEM_CHANGE((void *)(this)));
    }
  }
   
  void mesh_fem::set_finite_element(size_type cv, pintfem pif) {
    if (pif == NULL || pif->pf == NULL ) {
      if (fe_convex.is_in(cv)) {
	fe_convex.sup(cv);
	dof_enumeration_made = false;
	linked_mesh().lmsg_sender().send(MESH_FEM_CHANGE((void *)(this)));
      }
    }
    else {
      if (_linked_mesh->structure_of_convex(cv)->basic_structure() 
	  != pif->pf->basic_structure() || 
	  (pif->pf->target_dim() != Qdim && pif->pf->target_dim() != 1))
	DAL_THROW(internal_error,
		  "Incompatibility between fem and mesh element");
      if (!fe_convex.is_in(cv) || f_elems[cv] != pif) {
	fe_convex.add(cv);
	f_elems[cv] = pif;
	dof_enumeration_made = false;  
	linked_mesh().lmsg_sender().send(MESH_FEM_CHANGE((void *)(this)));
      }
    }
  }

  void mesh_fem::set_finite_element(const dal::bit_vector &cvs, pfem ppf,
			      const pintegration_method ppi) { 
    dal::bit_vector::const_iterator it = cvs.begin(), ite = cvs.end();
    pintfem pif =  give_intfem(ppf, ppi);
    for ( ; it != ite; ++it) if (*it) set_finite_element(it.index(), pif);
  }
  
  base_node mesh_fem::point_of_dof(size_type cv, size_type i) const {
    pfem pf = f_elems[cv]->pf;
    return linked_mesh().trans_of_convex(cv)->transform
      (pf->node_of_dof(i * pf->target_dim() / Qdim),
       linked_mesh().points_of_convex(cv));
  }

  base_node mesh_fem::point_of_dof(size_type d) const { 
    if (!dof_enumeration_made) enumerate_dof();
    return point_of_dof(first_convex_of_dof(d),
			ind_in_first_convex_of_dof(d));
  }

  struct _dof_comp { 
    dal::approx_less<scalar_type> comp;
    int operator()(const fem_dof& m, const fem_dof& n) const { 
      int d = dal::lexicographical_compare(m.P.begin(), m.P.end(),
					   n.P.begin(), n.P.end(), comp);
      if (d != 0) return d;
      return dof_description_compare(m.pnd, n.pnd);
    }
    _dof_comp(double e = 1.0E-10) : comp(e) { }
  };
  

  size_type mesh_fem::first_convex_of_dof(size_type d) const {
    for (size_type i = d; i != d - Qdim && i != size_type(-1); --i) {
      size_type j = dof_structure.first_convex_of_point(i);
      if (j != size_type(-1)) return j;
    }
    return size_type(-1);
  }

  size_type mesh_fem::ind_in_first_convex_of_dof(size_type d) const {
    for (size_type i = d; i != d - Qdim && i != size_type(-1); --i) {
      size_type j = dof_structure.first_convex_of_point(i);
      if (j != size_type(-1))
	return dof_structure.ind_in_first_convex_of_point(i);
    }
    return size_type(-1);
  }

  void mesh_fem::enumerate_dof(void) const {
    dal::bit_vector nn = fe_convex;
    std::queue<int> pile;
    dal::dynamic_tree_sorted<fem_dof, _dof_comp, 10> dof_sort;

    size_type cv;
    std::vector<size_type> tab;
    fem_dof fd;

    // cout << "\n\nEntering enumerate_dof\n\n\n\n\n";

    cv = nn.take_first();
    dof_structure.clear();

    while (cv != ST_NIL) {
      /* ajout des voisins dans la pile.                                  */

      size_type nbp = _linked_mesh->nb_points_of_convex(cv);

      for (size_type i = 0; i < nbp; i++) {
	size_type ip = _linked_mesh->ind_points_of_convex(cv)[i];
	bgeot::mesh_convex_ind_ct::const_iterator 
	  it = _linked_mesh->convex_to_point(ip).begin(),
	  ite =  _linked_mesh->convex_to_point(ip).end();
	for ( ; it != ite; ++it)
	  if (nn.is_in(*it)) { nn.sup(*it); pile.push(*it); }
      }

      size_type nbd = nb_dof_of_element(cv);
      pfem pf = fem_of_element(cv);
      pdof_description andof = already_numerate_dof(pf->dim());
      tab.resize(nbd);
      for (size_type i = 0; i < nbd; i++) {
	fd.P = point_of_dof(cv, i); // optimisable ...
	fd.pnd = pf->dof_types()[i];
	size_type j = 0, j_old = 0;
	if (fd.pnd == andof) {
	  // cout << "detecting a specialdof\n";
	  j = pf->index_of_already_numerate_dof(cv, i);
	  if (dof_sort.index_valid(j)) {
	    if (dof_sort[j].pnd != andof)
	      DAL_THROW(internal_error,
	      "Conflict between an already numerate dof and an existing dof.");
	  }
	  else
	    dof_sort.add_to_index(j, fd);
	}
	else {
	  if (pf->target_dim() == 1 && Qdim != 1) {

	    pdof_description paux = fd.pnd;
	    
	    for (size_type k = 0; k < Qdim; ++k) {
	      fd.pnd = to_coord_dof(paux, k);
	      
	      if (dof_linkable(fd.pnd))
		j = dof_sort.add_norepeat(fd);
	      else
		j = dof_sort.add(fd);
	      if (k == 0)
		tab[i] = j;
	      else if (j != j_old + 1)
		dof_sort.swap(j, j_old+1);
	      j_old = j;
	    }
	    
	  }
	  else {
	    if (dof_linkable(fd.pnd))
	    j = dof_sort.add_norepeat(fd);
	  else
	    j = dof_sort.add(fd);
	  }
	  tab[i] = j;
	}

      }
      
      dof_structure.add_convex_noverif(pf->structure(), tab.begin(), cv);
      
      if (pile.empty()) cv = nn.take_first();
                 else { cv = pile.front(); pile.pop(); }
    }
    
    dof_enumeration_made = true;
    nb_total_dof = (dof_sort.card() == 0) ? 0 : dof_sort.index().last_true()+1;
  }

/*     void mesh_fem::enumerate_dof(void) A finir ...  */
/*   { */
     /* L'algorithme de parcours est de type avancement par front          */
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
/* 	ssert(nbd < 10000); */
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
/* 	ssert(nbd < 10000); extraire une fonction qui ajoute les ddl d'un cv*/
/* 	for (i = 0; i < nbd; i++) */
/* 	{ */
/* 	  ncv = tab[i]; */
/* 	  if (nn.is_in(ncv)) { nn.sup(ncv); front2[j++] = ncv;  + tri  } */
/* 	} */
/*       } */
/*       if (j > 0 && front2[j-1] >= 0) front2[j++] = -1; */
/*     } */
/*   } */

  void mesh_fem::clear(void) {
    fe_convex.clear();
    dof_enumeration_made = false;
    linked_mesh().lmsg_sender().send(MESH_FEM_CHANGE((void *)(this)));
    dof_structure.clear();
    boundaries.clear();
    valid_boundaries.clear();
  }

  mesh_fem::mesh_fem(getfem_mesh &me, dim_type Q) : Qdim(Q) {
    _linked_mesh = &me;
 
    add_sender(me.lmsg_sender(), *this,
	   lmsg::mask(MESH_CLEAR()) | lmsg::mask(MESH_SUP_CONVEX()) |
	   lmsg::mask(MESH_SWAP_CONVEX()) | lmsg::mask(MESH_REFINE_CONVEX()) |
           lmsg::mask(MESH_UNREFINE_CONVEX()) | lmsg::mask(MESH_FEM_TOUCH()) |
	       lmsg::mask(MESH_DELETE()));
    is_valid = true;
  }

  mesh_fem::~mesh_fem() {
    if (is_valid) {
      linked_mesh().lmsg_sender().send(MESH_FEM_DELETE((void *)(this)));
    }
  }


  void mesh_fem::read_from_file(std::istream &ist) {
    dal::bit_vector npt;
    dal::dynamic_array<double> tmpv;
    char tmp[1024];

    ist.precision(16);
    clear();
    ist.seekg(0);
    ftool::read_untill(ist, "BEGIN MESH_FEM");

    while (true)
    {
      ftool::get_token(ist, tmp, 1023);

      if (strcmp(tmp, "END")==0) {
	break;
      } else if (strcmp(tmp, "CONVEX")==0) {
	ftool::get_token(ist, tmp, 1023);
	size_type ic = atoi(tmp);
	if (!linked_mesh().convex_index().is_in(ic)) {
	  DAL_THROW(failure_error, "Convex " << ic << " does not exist, are you sure "
		    "that the mesh attached to this object is right one ?");
	}
	
	ftool::get_token(ist, tmp, 1023);
	getfem::pfem fem = getfem::fem_descriptor(tmp);
	if (!fem) DAL_THROW(failure_error, "could not create the FEM '" << tmp << "'");

	ftool::get_token(ist, tmp, 1023);
	getfem::pintegration_method pfi = getfem::int_method_descriptor(tmp);
	if (!pfi) DAL_THROW(failure_error, "could not create the integration method '" << tmp << "'");
	
	dal::bit_vector bv; bv.add(ic);
	set_finite_element(bv, fem, pfi);
      } else if (strcmp(tmp, "BEGIN")==0) {
	ftool::get_token(ist, tmp, 1023);
	if (!strcmp(tmp, "BOUNDARY")) {
	  ftool::get_token(ist, tmp, 1023);
	  size_type bnum = atoi(tmp);
	  while (true) {
	    ftool::get_token(ist, tmp, 1023);
	    if (strcmp(tmp, "END")!=0) {
	      //	      cerr << "tmp = '" << tmp << "'" << endl;
	      size_type ic = atoi(tmp);
	      char *sf = strchr(tmp, '/');
	      if (sf) {
		size_type f = atoi(sf+1);
		add_boundary_elt(bnum, ic, f);
	      } else DAL_THROW(failure_error, "Syntax error in boundary " << bnum);
	    } else break;
	  }
	  ftool::get_token(ist, tmp, 1023);
	  ftool::get_token(ist, tmp, 1023);
	} else DAL_THROW(failure_error, "Syntax error in file at token" << tmp);
      } else if (strcmp(tmp, "QDIM")==0) {
	ftool::get_token(ist, tmp, 1023);
	int q = atoi(tmp); if (q <= 0 || q > 250) DAL_THROW(failure_error, "invalid qdim: "<<q);
	set_qdim(q);
      } else {
	DAL_THROW(failure_error, "Syntax error2 in file at token " << tmp);
      }
    }
  }

  void mesh_fem::read_from_file(const std::string &name)
  { 
    std::ifstream o(name.c_str());
    if (!o) DAL_THROW(file_not_found_error,
		      "Mesh_fem file '" << name << "' does not exist");
    read_from_file(o);
    o.close();
  }

  void mesh_fem::write_to_file(std::ostream &ost) const
  {
    ost << endl << "BEGIN MESH_FEM" << endl << endl;
    ost << "QDIM " << get_qdim() << endl;
    dal::bit_vector bv = convex_index();
    size_type cv;
    for (cv << bv; cv != size_type(-1); cv << bv) {
      ost << " CONVEX " << cv;
      ost << " " << name_of_fem(fem_of_element(cv));
      ost << " " << name_of_int_method(int_method_of_element(cv));
      ost << endl;
    }
    bv = get_valid_boundaries();

    dal::bit_vector blst = get_valid_boundaries();
    size_type bnum;
    for (bnum << blst; bnum != size_type(-1); bnum << blst) {
      ost << " BEGIN BOUNDARY " << bnum;
      dal::bit_vector cvlst = boundaries[bnum].cvindex;
      size_type cnt = 0;
      for (cv << cvlst; cv != size_type(-1); cv << cvlst) {
	dal::bit_vector cvflst = boundaries[bnum].faces_of_convex(cv);
	size_type f;
	for (f << cvflst; f != size_type(-1); f << cvflst, ++cnt) {
	  if ((cnt % 10) == 0) ost << endl << " ";
	  ost << " " << cv << "/" << f;
	}
      }
      ost << endl << " END BOUNDARY " << bnum << endl;
    }
    ost << "END MESH_FEM" << endl;
  }

  void mesh_fem::write_to_file(const std::string &name) const
  {
    std::ofstream o(name.c_str());
    if (!o)
      DAL_THROW(failure_error, "impossible to open file '" << name << "'");
    o << "% GETFEM MESH_FEM FILE " << endl;
    o << "% GETFEM VERSION " << __GETFEM_VERSION << "."
      << __GETFEM_REVISION << endl << endl << endl;
    write_to_file(o);
    o.close();
  }

}  /* end of namespace getfem.                                             */


