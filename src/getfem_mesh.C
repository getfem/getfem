/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Elements Methods (getfem)            */
/* File    :  getfem_mesh.C : meshes for computations.                     */
/*     									   */
/*                                                                         */
/* Date : November 05, 1999.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 1999-2003  Yves Renard.                                   */
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

#include <gmm_condition_number.h>
#include <getfem_mesh.h>
#include <bgeot_precomp.h>
namespace getfem
{
  getfem_mesh::getfem_mesh(dim_type NN) {
    dimension = NN; eps_p = 1.0E-10;
    pts.comparator() = dal::lexicographical_less<base_node,
                dal::approx_less<base_node::value_type> >(eps_p);
  }

  size_type getfem_mesh::add_point(const base_node &pt) {
    if (dimension == dim_type(-1)) dimension = pt.size();
    if (pt.size() != dimension)
      throw dimension_error("getfem_mesh::add_point : dimensions mismatch");
    
    bool present;
    size_type i = pts.add_norepeat(pt, false, &present);
    if (!present) lmsg_sender().send(MESH_ADD_POINT(i));
    return i;
  }

  void getfem_mesh::sup_point(size_type i) {
    if (!point_is_valid(i))
    { pts.sup(i); lmsg_sender().send(MESH_SUP_POINT(i)); }
  }

  void getfem_mesh::swap_points(size_type i, size_type j) {
    if (i != j) {
      bgeot::mesh<base_node>::swap_points(i,j);
      lmsg_sender().send(MESH_SWAP_POINT(i, j));
    }	 
  }

  void getfem_mesh::optimize_structure() {
    size_type i, j;
    j = nb_convex();
    for (i = 0; i < j; i++)
      if (!convex_tab.index_valid(i))
	swap_convex(i, convex_tab.ind_last());
    if (points().size())
      for (i = 0, j = points().size()-1; 
	   i < j && j != ST_NIL; ++i, --j) {
	while (i < j && j != ST_NIL && points().index()[i]) ++i;
	while (i < j && j != ST_NIL && !(points().index()[j])) --j;
	if (i < j && j != ST_NIL ) swap_points(i, j);
      }
  }

  void getfem_mesh::translation(base_small_vector V)
  {
    for (dal::bv_visitor i(points().index()); !i.finished(); ++i)
      points()[i] += V;
    points().resort();
  }

  void getfem_mesh::transformation(base_matrix M)
  {
    base_small_vector w(M.nrows());
    for (dal::bv_visitor i(points().index()); !i.finished(); ++i) {
      w = points()[i]; gmm::mult(M,w,points()[i]);
    }
    points().resort();
  }

  void getfem_mesh::clear(void) { 
    bgeot::mesh<base_node>::clear();
    gtab.clear(); trans_exists.clear();
    lmsg_sender().send(MESH_CLEAR());
  }

  size_type getfem_mesh::add_segment(size_type a, size_type b) { 
    size_type ipt[2]; ipt[0] = a; ipt[1] = b;
    return add_convex(bgeot::simplex_geotrans(1, 1), &(ipt[0]));
  }
  
  size_type getfem_mesh::add_triangle(size_type a, 
				      size_type b, size_type c) {
    size_type ipt[3]; ipt[0] = a; ipt[1] = b; ipt[2] = c;
    return add_simplex(2, &(ipt[0]));
  }
  
  size_type getfem_mesh::add_triangle_by_points
    (const base_node &p1, const base_node &p2, const base_node &p3)
  { return add_triangle(add_point(p1), add_point(p2), add_point(p3)); }
  
  size_type getfem_mesh::add_tetrahedron(size_type a, size_type b,
					 size_type c, size_type d) {
    size_type ipt[4]; ipt[0] = a; ipt[1] = b; ipt[2] = c; ipt[3] = d;
    return add_simplex(3, &(ipt[0]));
  }

  size_type getfem_mesh::add_tetrahedron_by_points(const base_node &p1,
						   const base_node &p2,
						   const base_node &p3,
						   const base_node &p4) {
    return add_tetrahedron(add_point(p1), add_point(p2),
			   add_point(p3), add_point(p4));
  }

  void getfem_mesh::sup_convex(size_type ic) {
    bgeot::mesh<base_node>::sup_convex(ic);
    trans_exists[ic] = false;
    lmsg_sender().send(MESH_SUP_CONVEX(ic));
  }

  void getfem_mesh::swap_convex(size_type i, size_type j) {
    if (i != j) {
      bgeot::mesh<base_node>::swap_convex(i,j);
      trans_exists.swap(i, j);
      gtab.swap(i,j);
      lmsg_sender().send(MESH_SWAP_CONVEX(i, j));
    }
  }

  base_small_vector getfem_mesh::normal_of_face_of_convex(size_type ic, short_type f,
							  const base_node &pt) const {
    bgeot::pgeometric_trans pgt = trans_of_convex(ic);
    base_matrix G(dim(),pgt->nb_points());
    vectors_to_base_matrix(G,points_of_convex(ic));
    bgeot::geotrans_interpolation_context c(trans_of_convex(ic), pt, G);
    return bgeot::compute_normal(c, f);
  }


  base_small_vector getfem_mesh::normal_of_face_of_convex(size_type ic, short_type f,
							  size_type n) const {
    bgeot::pgeometric_trans pgt = trans_of_convex(ic);
    bgeot::pgeotrans_precomp pgp = bgeot::geotrans_precomp(pgt, &pgt->geometric_nodes());
    base_matrix G(dim(),pgt->nb_points());
    vectors_to_base_matrix(G,points_of_convex(ic));
    bgeot::geotrans_interpolation_context c(pgp,pgt->structure()->ind_points_of_face(f)[n], G);
    return bgeot::compute_normal(c, f);
  }


  scalar_type  getfem_mesh::convex_quality_estimate(size_type ic) const { 
    base_matrix G;
    bgeot::vectors_to_base_matrix(G, points_of_convex(ic));
    return getfem::convex_quality_estimate(trans_of_convex(ic), G);
  }

  scalar_type  getfem_mesh::convex_radius_estimate(size_type ic) const { 
    base_matrix G;
    bgeot::vectors_to_base_matrix(G, points_of_convex(ic));
    return getfem::convex_radius_estimate(trans_of_convex(ic), G);
  }
  
  scalar_type getfem_mesh::minimal_convex_radius_estimate() const {
    if (convex_index().empty()) return 1;
    scalar_type r = convex_radius_estimate(convex_index().first_true());
    for (dal::bv_visitor cv(convex_index()); !cv.finished(); ++cv) {
      r = std::min(r, convex_radius_estimate(cv));
    }
    return r;
  }

  void getfem_mesh::copy_from(const getfem_mesh& m) {
    clear();
    bgeot::mesh<base_node>::operator=(m);
    eps_p = m.eps_p;
    gtab = m.gtab;
    trans_exists = m.trans_exists; 
    /* --- TODO --- : send an appropriate message to the mesh_fem using this mesh */
  }

  struct mesh_convex_structure_loc
  {
    bgeot::pgeometric_trans cstruct;
    size_type pts;
  };

  void getfem_mesh::read_from_file(std::istream &ist) {
   
    dal::bit_vector npt;
    dal::dynamic_array<double> tmpv;
    char tmp[1024];
    bool te = false, please_get = true;

    ist.precision(16);
    clear();
    ist.seekg(0);ist.clear();
    ftool::read_until(ist, "BEGIN POINTS LIST");

    while (!te)
    {
      if (please_get) ftool::get_token(ist, tmp, 1023); else please_get = true;

      if (!strcmp(tmp, "END"))
      { te = true; }
      else if (!strcmp(tmp, "POINT"))
      {
	ftool::get_token(ist, tmp, 1023);
        size_type ip = atoi(tmp);
        dim_type d = 0;
	if (npt.is_in(ip))
	  DAL_THROW(failure_error, 
		    "Two points with the same index. loading aborted.");
	npt.add(ip);
	ftool::get_token(ist, tmp, 1023);
	while (isdigit(tmp[0]) || tmp[0] == '-' || tmp[0] == '+'
	                       || tmp[0] == '.')
	{ tmpv[d++] = atof(tmp); ftool::get_token(ist, tmp, 1023); }
	please_get = false;
	if (dimension == dim_type(-1)) dimension = d;
	else if (dimension != d)
	  DAL_THROW(failure_error, "Points of different dimensions.");
	base_node v(d);
	for (size_type i = 0; i < d; i++) v[i] = tmpv[i];
	size_type ipl = add_point(v);
	if (ip != ipl) {
	  if (npt.is_in(ipl))
	    DAL_THROW(failure_error, 
		    "Two points with the same coords. loading aborted.");
	  swap_points(ip, ipl);
	}
      } else if (strlen(tmp)) {
	DAL_THROW(failure_error, "Syntax error in file, at token '" << tmp
		  << "', pos=" << std::streamoff(ist.tellg()));
      } else if (ist.eof()) {
	DAL_THROW(failure_error, "Unexpected end of stream");	
      }
    }

    bool tend = false;
    dal::dynamic_alloc<size_type> cv_pt;
    dal::dynamic_array<mesh_convex_structure_loc> cv;
    dal::bit_vector ncv;
    
    ist.seekg(0);
    if (!ftool::read_until(ist, "BEGIN MESH STRUCTURE DESCRIPTION"))
      DAL_THROW(failure_error, "This seems not to be a mesh file");

    while (!tend) {
      tend = !ftool::get_token(ist, tmp, 1023);
      if (!strcmp(tmp, "END"))
      { tend = true; }
      else if (!strcmp(tmp, "CONVEX"))
      {
        size_type ic;
	ftool::get_token(ist, tmp, 1023);
        ic = dal::abs(atoi(tmp));
	if (ncv.is_in(ic)) 
	  DAL_THROW(failure_error,
		    "Negative or repeated index, loading aborted.");
	ncv.add(ic);
	ftool::get_token(ist, tmp, 1023);
	
	bgeot::pgeometric_trans pgt = bgeot::geometric_trans_descriptor(tmp);

       
	size_type nb = pgt->nb_points();

	cv[ic].cstruct = pgt;
	cv[ic].pts = cv_pt.alloc(nb);
	for (size_type i = 0; i < nb; i++)
	{
	  ftool::get_token(ist, tmp, 1023);	  
	  cv_pt[cv[ic].pts+i] = dal::abs(atoi(tmp));
	}
      }
      else if (strlen(tmp)) {
	DAL_THROW(failure_error, "Syntax error reading a mesh file "
		  " at pos " << std::streamoff(ist.tellg())
		  << "(expecting 'CONVEX' or 'END', found '" << tmp << "')"); 
      } else if (ist.eof()) {
	DAL_THROW(failure_error, "Unexpected end of stream "
		  << "(missing BEGIN MESH/END MESH ?)");
      }
    }

    for (dal::bv_visitor ic(ncv); !ic.finished(); ++ic) {
      size_type i = add_convex(cv[ic].cstruct, cv_pt.begin() + cv[ic].pts);
      if (i != ic) swap_convex(i, ic);
    }

    lmsg_sender().send(MESH_READ_FROM_FILE(ist));
  }

  void getfem_mesh::read_from_file(const std::string &name) { 
    std::ifstream o(name.c_str());
    if (!o) DAL_THROW(file_not_found_error,
		      "Mesh file '" << name << "' does not exist");
    read_from_file(o); 
    o.close();
  }

  template<class ITER>
    static void write_tab_to_file_(std::ostream &ost, ITER b, ITER e)
  { for ( ; b != e; ++b) ost << "  " << *b; }

  template<class ITER>
    static void write_convex_to_file_(const getfem_mesh &ms,
				      std::ostream &ost,
				      ITER b, ITER e) {
    for ( ; b != e; ++b) {
      size_type i = b.index();
      ost << "CONVEX " << i << "    "
	  << bgeot::name_of_geometric_trans(ms.trans_of_convex(i)).c_str()
	  << "    ";
      write_tab_to_file_(ost, ms.ind_points_of_convex(i).begin(),
			 ms.ind_points_of_convex(i).end()  );
      ost << '\n';
    }
  }

  template<class ITER> static void write_point_to_file_(std::ostream &ost,
						  ITER b, ITER e)
  { for ( ; b != e; ++b) ost << "  " << *b; ost << '\n'; }

  void getfem_mesh::write_to_file(std::ostream &ost) const {
    ost.precision(16);
    ost << '\n' << "BEGIN POINTS LIST" << '\n' << '\n';
    bgeot::mesh_point_st_ct::const_iterator b = point_structures().begin();
    bgeot::mesh_point_st_ct::const_iterator e = point_structures().end();
    for (size_type i = 0; b != e; ++b, ++i)
      if ( (*b).is_valid() ) {
	ost << "  POINT  " << i;
	write_point_to_file_(ost, points()[i].begin(), points()[i].end());
      }
    ost << '\n' << "END POINTS LIST" << '\n' << '\n' << '\n';
    
    ost << '\n' << "BEGIN MESH STRUCTURE DESCRIPTION" << '\n' << '\n';
    write_convex_to_file_(*this, ost, convex_tab.tas_begin(),
			              convex_tab.tas_end());
    ost << '\n' << "END MESH STRUCTURE DESCRIPTION" << '\n';
    lmsg_sender().send(MESH_WRITE_TO_FILE(ost));
  }

  void getfem_mesh::write_to_file(const std::string &name) const {
    std::ofstream o(name.c_str());
    if (!o)
      DAL_THROW(failure_error, "impossible to open file '" << name << "'");
    o << "% GETFEM MESH FILE " << '\n';
    o << "% GETFEM VERSION " << GETFEM_VERSION << '\n' << '\n' << '\n';
    write_to_file(o);
    o.close();
  }

  size_type getfem_mesh::memsize() const {
    return bgeot::mesh<base_node>::memsize() + 
      (pts.index().last_true()+1)*dim()*sizeof(scalar_type)+
      sizeof(getfem_mesh) - sizeof(bgeot::mesh<base_node>)
      +trans_exists.memsize() + gtab.memsize();
  }
    
  /* TODO : use the geotrans from an "equilateral" reference element to the real element 
     check if the sign of the determinants does change (=> very very bad quality of the convex)
  */
  scalar_type convex_quality_estimate(bgeot::pgeometric_trans pgt,
				      const base_matrix& G) {
    static bgeot::pgeometric_trans pgt_old = 0;
    static bgeot::pgeotrans_precomp pgp = 0;
    if (pgt_old != pgt) {
      pgt_old=pgt; pgp=bgeot::geotrans_precomp(pgt, &pgt->convex_ref()->points());
    }

    size_type n = (pgt->is_linear()) ? 1 : pgt->nb_points();
    scalar_type q = 1;
    base_matrix K(pgp->grad(0).ncols(),G.nrows());
    for (size_type ip=0; ip < n; ++ip) {
      gmm::mult(gmm::transposed(pgp->grad(ip)), gmm::transposed(G), K);
      q = std::max(q, gmm::condition_number(K));
    }
    return 1./q;
  }

  scalar_type convex_radius_estimate(bgeot::pgeometric_trans pgt,
				     const base_matrix& G) {
    static bgeot::pgeometric_trans pgt_old = 0;
    static bgeot::pgeotrans_precomp pgp = 0;
    if (pgt_old != pgt) {
      pgt_old=pgt; pgp=bgeot::geotrans_precomp(pgt, &pgt->convex_ref()->points());
    }
    size_type N = G.nrows();
    size_type n = (pgt->is_linear()) ? 1 : pgt->nb_points();
    scalar_type q = 0;
    base_matrix K(pgp->grad(0).ncols(),G.nrows());
    for (size_type ip=0; ip < n; ++ip) {
      gmm::mult(gmm::transposed(pgp->grad(ip)), gmm::transposed(G), K);
      scalar_type emax, emin; gmm::condition_number(K,emax,emin);
      q = std::max(q, emax);
    }
    return q * sqrt(scalar_type(N)) / scalar_type(N);
  }

  /* extract faces of convexes which are not shared
      + convexes whose dimension is smaller that m.dim()
  */
  void
  outer_faces_of_mesh(const getfem::getfem_mesh &m, const dal::bit_vector& cvlst, convex_face_ct& flist) {
    for (dal::bv_visitor ic(cvlst); !ic.finished(); ++ic) {
      if (m.structure_of_convex(ic)->dim() == m.dim()) {
	for (size_type f = 0; f < m.structure_of_convex(ic)->nb_faces(); f++) {
	  if (bgeot::neighbour_of_convex(m,ic,f).empty()) {
	    flist.push_back(convex_face(ic,f));
	  }
	}
      } else {
	flist.push_back(convex_face(ic,size_type(-1)));
      }
    }
  }

}  /* end of namespace getfem.                                             */
