// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 1999-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================


#include <queue>
#include <dal_singleton.h>
#include <getfem_mesh_fem.h>

namespace getfem {
  
  dal::bit_vector mesh_fem::dof_on_region(const mesh_region &b) const {
    if (!dof_enumeration_made) this->enumerate_dof();
    dal::bit_vector res;
    for (getfem::mr_visitor v(b,linked_mesh()); !v.finished(); ++v) {
      size_type cv = v.cv();
      if (convex_index().is_in(cv)) {
	if (v.is_face()) {
	  size_type f = v.f();
	  size_type nbb =
	    dof_structure.structure_of_convex(cv)->nb_points_of_face(f);
	  for (size_type i = 0; i < nbb; ++i) {
	    size_type n = Qdim / fem_of_element(cv)->target_dim();
	    for (size_type ll = 0; ll < n; ++ll)
	      res.add(dof_structure.ind_points_of_face_of_convex(cv,f)[i]+ll);
	  }
	} else {
	  size_type nbb =
	    dof_structure.structure_of_convex(cv)->nb_points();
	  for (size_type i = 0; i < nbb; ++i) {
	    size_type n = Qdim / fem_of_element(cv)->target_dim();
	    for (size_type ll = 0; ll < n; ++ll)
	      res.add(dof_structure.ind_points_of_convex(cv)[i]+ll);
	  }
	}
      }
    }
    return res;
  }
  
  void mesh_fem::receipt(const MESH_CLEAR &) { clear(); }
  void mesh_fem::receipt(const MESH_DELETE &) { clear(); }
  void mesh_fem::receipt(const MESH_SUP_CONVEX &m) { 
    if (fe_convex[m.icv])
      { fe_convex[m.icv] = false; dof_enumeration_made = false; }
  }
  void mesh_fem::receipt(const MESH_ADD_CONVEX &m) {
    if (auto_add_elt_K != size_type(-1)) {
      pfem pf = getfem::classical_fem(linked_mesh().trans_of_convex(m.icv), 
				      auto_add_elt_K);
      set_finite_element(m.icv, pf);
    }
  }
  void mesh_fem::receipt(const MESH_SWAP_CONVEX &m) { 
    fe_convex.swap(m.icv1, m.icv2);
    f_elems.swap(m.icv1, m.icv2);
  }
  void mesh_fem::receipt(const MESH_REFINE_CONVEX &m) { 
    if (m.is_refine) {
      if (fe_convex[m.icv])
	for (size_type i = 0; i < m.sub_cv_list.size(); ++i) {
	  f_elems[m.sub_cv_list[i]] = f_elems[m.icv];
	  fe_convex.add(m.sub_cv_list[i]);
	}
    }
    else if (fe_convex[m.sub_cv_list[0]]) {
      f_elems[m.icv] = f_elems[m.sub_cv_list[0]];
      fe_convex.add(m.icv);
    }
  }
   
  void mesh_fem::set_finite_element(size_type cv, pfem pf) {
    if (pf == 0) {
      if (fe_convex.is_in(cv)) {
	fe_convex.sup(cv);
	dof_enumeration_made = false;
	touch();
      }
    }
    else {
      if (linked_mesh_->structure_of_convex(cv)->basic_structure() 
	  != pf->basic_structure(cv))
	DAL_THROW(dal::failure_error,
		  "Incompatibility between fem " << name_of_fem(pf) << 
		  " and mesh element " <<
		  name_of_geometric_trans(linked_mesh_->trans_of_convex(cv)));
      if ((Qdim % pf->target_dim()) != 0 && pf->target_dim() != 1)
	DAL_THROW(dal::failure_error,
		  "Incompatibility between Qdim=" << int(Qdim) << " and target_dim " <<
		  int(pf->target_dim()) << " of " << name_of_fem(pf));
      if (!fe_convex.is_in(cv) || f_elems[cv] != pf) {
	fe_convex.add(cv);
	f_elems[cv] = pf;
	dof_enumeration_made = false;  
	touch();
      }
    }
  }

  void mesh_fem::set_finite_element(const dal::bit_vector &cvs, pfem ppf) { 
    for (dal::bv_visitor cv(cvs); !cv.finished(); ++cv)
      set_finite_element(cv, ppf);
  }

  void mesh_fem::set_finite_element(pfem ppf)
  { set_finite_element(linked_mesh().convex_index(), ppf); }
  
  void mesh_fem::set_classical_finite_element(const dal::bit_vector &cvs, 
					      dim_type fem_degree) {
    for (dal::bv_visitor cv(cvs); !cv.finished(); ++cv) {
      pfem pf = getfem::classical_fem(linked_mesh().trans_of_convex(cv), 
				      fem_degree);
      set_finite_element(cv, pf);
    }
  }

  void mesh_fem::set_classical_finite_element(dim_type fem_degree)
  { set_classical_finite_element(linked_mesh().convex_index(), fem_degree); }

  void mesh_fem::set_classical_discontinuous_finite_element
  (const dal::bit_vector &cvs, dim_type fem_degree, scalar_type alpha) {
    for (dal::bv_visitor cv(cvs); !cv.finished(); ++cv) {
      pfem pf = getfem::classical_discontinuous_fem
	(linked_mesh().trans_of_convex(cv), fem_degree, alpha);
      set_finite_element(cv, pf);
    }
  }

  void mesh_fem::set_classical_discontinuous_finite_element
  (dim_type fem_degree, scalar_type alpha) { 
    set_classical_discontinuous_finite_element(linked_mesh().convex_index(),
					       fem_degree,alpha);
  }

  base_node mesh_fem::point_of_dof(size_type cv, size_type i) const {
    pfem pf = f_elems[cv];
    return linked_mesh().trans_of_convex(cv)->transform
      (pf->node_of_dof(cv, i * pf->target_dim() / Qdim),
       linked_mesh().points_of_convex(cv));
  }

  base_node mesh_fem::point_of_dof(size_type d) const {
    if (!dof_enumeration_made) enumerate_dof();
    for (size_type i = d; i != d - Qdim && i != size_type(-1); --i) {
      size_type cv = dof_structure.first_convex_of_point(i);
      if (cv != size_type(-1)) {
	pfem pf = f_elems[cv];
	return linked_mesh().trans_of_convex(cv)->transform
	  (pf->node_of_dof(cv, dof_structure.ind_in_convex_of_point(cv, i)),
	   linked_mesh().points_of_convex(cv));
      }
    }
    DAL_THROW(failure_error, "Inexistent dof");
  }

  dim_type mesh_fem::dof_qdim(size_type d) const {
    if (!dof_enumeration_made) enumerate_dof();
    size_type cv = first_convex_of_dof(d);
    if (cv == size_type(-1)) DAL_THROW(failure_error, "Inexistent dof");
    size_type tdim = f_elems[cv]->target_dim();
    return dof_structure.ind_in_convex_of_point(cv, d) % (Qdim / tdim);
  }

  size_type mesh_fem::first_convex_of_dof(size_type d) const {
    for (size_type i = d; i != d - Qdim && i != size_type(-1); --i) {
      size_type j = dof_structure.first_convex_of_point(i);
      if (j != size_type(-1)) return j;
    }
    return size_type(-1);
  }

//   size_type mesh_fem::ind_in_convex_of_dof(size_type d) const {
//     for (size_type i = d; i != d - Qdim && i != size_type(-1); --i) {
//       size_type j = dof_structure.first_convex_of_point(i);
//       if (j != size_type(-1))
// 	return (dof_structure.ind_in_first_convex_of_point(i)
// 		* Qdim / f_elems[j]->target_dim());
//     }
//     return size_type(-1);
//   }

  const mesh::ind_cv_ct &mesh_fem::convex_to_dof(size_type d) const {
    for (size_type i = d; i != d - Qdim && i != size_type(-1); --i) {
      size_type j = dof_structure.first_convex_of_point(i);
      if (j != size_type(-1)) return dof_structure.convex_to_point(i);
    }
    DAL_THROW(failure_error, "Inexistent dof");
  }

  struct dof_comp_ { 
    dal::approx_less<scalar_type> comp;
    int operator()(const fem_dof& m, const fem_dof& n) const { 
      int d = dal::lexicographical_compare(m.P.begin(), m.P.end(),
					   n.P.begin(), n.P.end(), comp);
      if (d != 0) return (d < 0);
      if (m.part == n.part)
	return dof_description_compare(m.pnd, n.pnd) < 0;
      else if (m.part < n.part) return true;
      else if (m.part > n.part) return false;
    }
    dof_comp_(double e = 1.0E-10) : comp(e) { }
  };

  typedef std::map<fem_dof, size_type, dof_comp_> dof_sort_type;

  // static double enumerate_dof_time = 0;

  /// Enumeration of dofs
  void mesh_fem::enumerate_dof(void) const {
    if (fe_convex.card() == 0) {
      dof_enumeration_made = true;
      nb_total_dof = 0;
      return;
    }
    const std::vector<size_type> &cmk = linked_mesh().cuthill_mckee_ordering();
    // double t = dal::uclock_sec();

    //cout << "enumerate_dof\n";

    dof_sort_type dof_sort;
    dal::bit_vector encountered_global_dof;
    dal::dynamic_array<size_type> ind_global_dof;
    std::vector<size_type> tab;
    size_type nbdof = 0;

    size_type cv = fe_convex.first_true();
    fem_dof fd;

    dof_structure.clear();
    encountered_global_dof.clear();

    bgeot::pstored_point_tab pspt_old = 0;
    bgeot::pgeometric_trans pgt_old = 0;
    bgeot::pgeotrans_precomp pgp = 0;
    for (size_type icv=0; icv < cmk.size(); ++icv) {
      cv = cmk[icv];
      if (!fe_convex.is_in(cv)) continue;
      pfem pf = fem_of_element(cv);
      bgeot::pgeometric_trans pgt = linked_mesh().trans_of_convex(cv);
      bgeot::pstored_point_tab pspt = pf->node_tab(cv);
      if (pgt != pgt_old || pspt != pspt_old)
	pgp = bgeot::geotrans_precomp(pgt, pspt);
      pgt_old = pgt; pspt_old = pspt;
      
      size_type nbd = pf->nb_dof(cv); 
      pdof_description andof = global_dof(pf->dim());
      tab.resize(nbd);

      for (size_type i = 0; i < nbd; i++) {
	fd.P.resize(linked_mesh().dim()); 
	fd.pnd = pf->dof_types()[i];
	fd.part = get_dof_partition(cv);
	if (fd.pnd == andof) {
	  size_type num = pf->index_of_global_dof(cv, i);
	  if (!(encountered_global_dof[num])) {
	    ind_global_dof[num] = nbdof;
	    nbdof += Qdim / pf->target_dim();
	    encountered_global_dof[num] = true;
	  }
	  tab[i] = ind_global_dof[num];
	} else if (!dof_linkable(fd.pnd)) {
	  tab[i] = nbdof;
	  nbdof += Qdim / pf->target_dim();
	} else {
	  pgp->transform(linked_mesh().points_of_convex(cv), i, fd.P);

	  std::pair<dof_sort_type::iterator, bool> pa = dof_sort.insert(std::make_pair(fd, nbdof));
	  if (pa.second) {
	    tab[i] = nbdof;
	    nbdof += Qdim / pf->target_dim();
	  }
	  else { tab[i] = pa.first->second; }
	}
      }
      dof_structure.add_convex_noverif(pf->structure(cv), tab.begin(), cv);
    }
    
    dof_enumeration_made = true;
    nb_total_dof = nbdof;
    
//    enumerate_dof_time += dal::uclock_sec() - t;
//    cerr << "enumerate_dof_time: " << enumerate_dof_time << " sec [nbd=" << nbdof << "]\n";
  }





  void mesh_fem::clear(void) {
    fe_convex.clear();
    dof_enumeration_made = false;
    touch();
    dof_structure.clear();
  }

  mesh_fem::mesh_fem(const mesh &me, dim_type Q)
    : dof_enumeration_made(false), auto_add_elt_K(size_type(-1)), 
      Qdim(Q), QdimM(1), QdimN(1) {
    linked_mesh_ = &me;
    this->add_dependency(me);
    add_sender(me.lmsg_sender(), *this,
	       lmsg::mask(MESH_CLEAR::ID) | lmsg::mask(MESH_SUP_CONVEX::ID) |
	       lmsg::mask(MESH_SWAP_CONVEX::ID) | lmsg::mask(MESH_DELETE::ID) |
	       lmsg::mask(MESH_ADD_CONVEX::ID)|
	       lmsg::mask(MESH_REFINE_CONVEX::ID));
  }

  mesh_fem::~mesh_fem() {}

  void mesh_fem::read_from_file(std::istream &ist) {
    dal::bit_vector npt;
    dal::dynamic_array<double> tmpv;
    std::string tmp, tmp2;
    bool dof_read = false;
    ist.precision(16);
    clear();
    ist.seekg(0);ist.clear();
    ftool::read_until(ist, "BEGIN MESH_FEM");

    while (true) {
      ist >> std::ws; ftool::get_token(ist, tmp);
      if (ftool::casecmp(tmp, "END")==0) {
	break;
      } else if (ftool::casecmp(tmp, "CONVEX")==0) {
	ftool::get_token(ist, tmp);
	size_type ic = atoi(tmp.c_str());
	if (!linked_mesh().convex_index().is_in(ic)) {
	  DAL_THROW(failure_error, "Convex " << ic <<
		    " does not exist, are you sure "
		    "that the mesh attached to this object is right one ?");
	}
	
	int rgt = ftool::get_token(ist, tmp);
	if (rgt != 3) { // for backward compatibility
	  char c; ist.get(c);
	  while (!isspace(c)) { tmp.push_back(c); ist.get(c); }
	}
	getfem::pfem fem = getfem::fem_descriptor(tmp);
	if (!fem) DAL_THROW(failure_error, "could not create the FEM '" 
			    << tmp << "'");
	set_finite_element(ic, fem);
      } else if (ftool::casecmp(tmp, "BEGIN")==0) {
	ftool::get_token(ist, tmp);
	if (!ftool::casecmp(tmp, "DOF_ENUMERATION")) {
	  dal::bit_vector doflst;
	  dof_structure.clear(); dof_enumeration_made = false; touch();
	  while (true) {
	    ftool::get_token(ist, tmp);
	    if (ftool::casecmp(tmp, "END")==0) { 
	      break;
	    }
	    ftool::get_token(ist, tmp2);

	    size_type ic = atoi(tmp.c_str());
	    std::vector<size_type> tab;
	    if (convex_index().is_in(ic) && tmp.size() &&
		isdigit(tmp[0]) && tmp2 == ":") { // && tmp[tmp.size()-1] == ':') {
	      tab.resize(nb_dof_of_element(ic));
	      for (size_type i=0; i < fem_of_element(ic)->nb_dof(ic); i++) {
		ist >> tab[i];
		for (size_type q=0; q < size_type(get_qdim())
		       / fem_of_element(ic)->target_dim(); ++q)
		  doflst.add(tab[i]+q);
	      }
	      dof_structure.add_convex_noverif
		(fem_of_element(ic)->structure(ic), tab.begin(), ic);
	    } else DAL_THROW(failure_error, "Missing convex or wrong number "
			     << "in dof enumeration: '" 
			     << tmp << "' [pos="
			     << std::streamoff(ist.tellg())<<"]");
	    /*ftool::get_token(ist, tmp);
	      cerr << " tok: '" << tmp << "'\n";*/
	  } 
	  dof_read = true;
	  this->dof_enumeration_made = true;
	  touch();
	  this->nb_total_dof = doflst.card();
	  ist >> ftool::skip("DOF_ENUMERATION");
	} else if (tmp.size())
	  DAL_THROW(failure_error, "Syntax error in file at token '"
		    << tmp << "' [pos=" << std::streamoff(ist.tellg())
							  << "]");
      } else if (ftool::casecmp(tmp, "QDIM")==0) {
	if (dof_read)
	  DAL_THROW(failure_error, "Can't change QDIM after dof enumeration");
	ftool::get_token(ist, tmp);
	int q = atoi(tmp.c_str());
	if (q <= 0 || q > 250) DAL_THROW(failure_error, "invalid qdim: "<<q);
	set_qdim(q);
      } else if (tmp.size()) {
	DAL_THROW(failure_error, "Unexpected token '" << tmp <<
		  "' [pos=" << std::streamoff(ist.tellg()) << "]");
      } else if (ist.eof()) {
	DAL_THROW(failure_error, "Unexpected end of stream "
		  << "(missing BEGIN MESH_FEM/END MESH_FEM ?)");	
      }
    }
  }

  void mesh_fem::read_from_file(const std::string &name)
  { 
    std::ifstream o(name.c_str());
    if (!o) DAL_THROW(file_not_found_error,
		      "Mesh_fem file '" << name << "' does not exist");
    read_from_file(o);
  }

  void mesh_fem::write_to_file(std::ostream &ost) const
  {
    ost << '\n' << "BEGIN MESH_FEM" << '\n' << '\n';
    ost << "QDIM " << size_type(get_qdim()) << '\n';
    for (dal::bv_visitor cv(convex_index()); !cv.finished(); ++cv) {
      ost << " CONVEX " << cv;
      ost << " \'" << name_of_fem(fem_of_element(cv));
      ost << "\'\n";
    }

    ost << " BEGIN DOF_ENUMERATION " << '\n';
    for (dal::bv_visitor cv(convex_index()); !cv.finished(); ++cv) {
      ost << "  " << cv << ": ";
      ind_dof_ct::const_iterator it = ind_dof_of_element(cv).begin();
      while (it != ind_dof_of_element(cv).end()) {
	ost << " " << *it;
	// skip repeated dofs for "pseudo" vector elements
	for (size_type i=0;
	     i < size_type(get_qdim())/fem_of_element(cv)->target_dim(); ++i)
	  ++it;
      }
      ost << '\n';
    }
    ost << " END DOF_ENUMERATION " << '\n';
    ost << "END MESH_FEM" << '\n';
  }

  void mesh_fem::write_to_file(const std::string &name, bool with_mesh) const
  {
    std::ofstream o(name.c_str());
    if (!o)
      DAL_THROW(failure_error, "impossible to open file '" << name << "'");
    o << "% GETFEM MESH_FEM FILE " << '\n';
    o << "% GETFEM VERSION " << GETFEM_VERSION << '\n' << '\n' << '\n';
    if (with_mesh) linked_mesh().write_to_file(o);
    write_to_file(o);
  }

  struct mf__key_ {
    const mesh *pmesh;
    dim_type order;
    mf__key_(const mesh &msh, dim_type o) : pmesh(&msh),order(o) {}
    bool operator <(const mf__key_ &a) const {
    if (pmesh < a.pmesh) return true; else
      if (a.pmesh < pmesh) return false; else
	if (order < a.order) return true; else return false;
    }
  };


  class classical_mesh_fem_pool : public mesh_receiver {

    typedef const mesh_fem * pmesh_fem;
    typedef std::map<mf__key_, pmesh_fem> mesh_fem_tab;

    mesh_fem_tab mfs;

    // void receipt(const MESH_CLEAR &) {}
    void receipt(const MESH_DELETE &M) {
      for (mesh_fem_tab::iterator it = mfs.begin(); it != mfs.end(); ) {
	mesh_fem_tab::iterator it2 = it; it2 ++;
	if (it->first.pmesh == M.msh) mfs.erase(it);
	it = it2;
      }
    }
//     void receipt(const MESH_ADD_CONVEX &m) {}
//     void receipt(const MESH_SUP_CONVEX &m) {}
//     void receipt(const MESH_SWAP_CONVEX &m) {}

  public :

    const mesh_fem &operator()(const mesh &msh, dim_type o) {

      mf__key_ key(msh, o);
      mesh_fem_tab::iterator it = mfs.find(key);
      assert(it == mfs.end() || it->second->is_context_valid());
      
      if (it == mfs.end()) {
	// the list of mesh pointers should be sorted ...
	for (mesh_fem_tab::iterator itt = mfs.begin(); itt != mfs.end(); ++itt)
	  if (itt->first.pmesh == &msh) goto nothing_to_do;
	add_sender(msh.lmsg_sender(), *this, lmsg::mask(MESH_DELETE::ID));
	
      nothing_to_do :
	
	mesh_fem *pmf = new mesh_fem(msh);
	pmf->set_auto_add(o);
	pmf->set_classical_finite_element(o);
	return *(mfs[key] = pmf);
      }
      else return *(it->second);
    }
    
  };

  const mesh_fem &classical_mesh_fem(const mesh &msh,
				     dim_type order) {
    classical_mesh_fem_pool &pool
      = dal::singleton<classical_mesh_fem_pool>::instance();
    return pool(msh, order);
  }

  struct dummy_mesh_fem_ {
    mesh m;
    mesh_fem mf;
    dummy_mesh_fem_() : mf(m) {}
  };
  const mesh_fem &dummy_mesh_fem(void) { return dal::singleton<dummy_mesh_fem_>::instance().mf; }




}  /* end of namespace getfem.                                             */


