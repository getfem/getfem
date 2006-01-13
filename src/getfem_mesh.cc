// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Elements Methods (getfem)
// File    : getfem_mesh.cc : meshes for computations.
//           
// Date    : November 05, 1999.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1999-2006 Yves Renard
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
#include <gmm_condition_number.h>
#include <getfem_mesh.h>

namespace getfem {

  void mesh::sup_convex_from_regions(size_type c) {
    for (dal::bv_visitor i(valid_cvf_sets); !i.finished(); ++i)
      cvf_sets[i].sup(c);
    touch();
  }

  void mesh::swap_convex_in_regions(size_type c1, size_type c2) {
    for (dal::bv_visitor i(valid_cvf_sets); !i.finished(); ++i)
      cvf_sets[i].swap_convex(c1, c2);
    touch();
  }

  mesh::mesh(dim_type NN) {
#if GETFEM_PARA_LEVEL > 1
    modified = true;
#endif
    cuthill_mckee_uptodate = false;
    dimension = NN; eps_p = 1.0E-10;
    pts.comparator() = dal::lexicographical_less<base_node,
                dal::approx_less<base_node::value_type> >(eps_p);
  }

#if GETFEM_PARA_LEVEL > 1

    void mesh::compute_mpi_region(void) const {

      int size, rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);

      if (size < 2) { mpi_region = mesh_region::all_convexes(); }
      else {

	int ne = int(nb_convex());
	std::vector<int> xadj(ne+1), adjncy, numelt(ne), npart(ne);
	
	int j = 0, k = 0;
	for (dal::bv_visitor ic(convex_index()); !ic.finished(); ++ic, ++j) {
	  numelt[j] = ic;
	  xadj[j] = k;
	  bgeot::mesh_structure::ind_set s(neighbours_of_convex(ic));
	  for (bgeot::mesh_structure::ind_set::iterator it = s.begin();
	       it != s.end(); ++it) { xadjncy.push_back(*it); ++k; }  
	}
	xadj[j] = k;

	int wgtflag = 0, edgecut, numflag = 0, options[5] = {0,0,0,0,0};

	METIS_PartGraphKway(&ne, &(xadj[0]), &(adjncy[0]), 0, 0, &wgtflag,
			    &numflag, &size, options, &edgecut, &(npart[0]));

	for (size_type i = 0; i < size_type(ne); ++i)
	  if (nparts[i] == rank) mpi_region.add(numelt[i]);
      }
      modified = false;
      valid_sub_regions.clear();
    }

  void mesh::compute_mpi_sub_region(size_type n) const {
    if (valid_cvf_sets.is_in(n)) {
      mpi_sub_region[n] = mesh_region::intersection(cvf_sets[n], mpi_region);
    }
    else
      mpi_sub_region[n] = mesh_region();
    valid_sub_regions.add(n);
  }

    
#endif


  size_type mesh::add_point(const base_node &pt, bool norepeat) {
    if (dimension == dim_type(-1)) dimension = pt.size();
    if (pt.size() != dimension)
      throw dimension_error("mesh::add_point : dimensions mismatch");
    
    if (norepeat) {
      bool present;
      size_type i = pts.add_norepeat(pt, false, &present);
      return i;
    }
    else return pts.add(pt);
  }

  void mesh::optimize_structure() {
    size_type i, j;
    j = nb_convex();
    for (i = 0; i < j; i++)
      if (!convex_tab.index_valid(i))
	swap_convex(i, convex_tab.ind_last());
    if (pts.size())
      for (i = 0, j = pts.size()-1; 
	   i < j && j != ST_NIL; ++i, --j) {
	while (i < j && j != ST_NIL && pts.index()[i]) ++i;
	while (i < j && j != ST_NIL && !(pts.index()[j])) --j;
	if (i < j && j != ST_NIL ) swap_points(i, j);
      }
  }

  const std::vector<size_type> &mesh::cuthill_mckee_ordering() const {
    if (!cuthill_mckee_uptodate) {
      bgeot::cuthill_mckee_on_convexes(*this, cmk_order);
      cuthill_mckee_uptodate = true;
    }
    return cmk_order;
  }

  void mesh::translation(base_small_vector V)
  {
    for (dal::bv_visitor i(pts.index()); !i.finished(); ++i)
      pts[i] += V;
    pts.resort();
  }

  void mesh::transformation(base_matrix M)
  {
    base_small_vector w(M.nrows());
    if (gmm::mat_nrows(M) == 0 || gmm::mat_ncols(M) != dim()) 
      DAL_THROW(dal::dimension_error, "invalid dimensions for the transformation matrix");
    for (dal::bv_visitor i(pts.index()); !i.finished(); ++i) {
      w = pts[i]; gmm::resize(pts[i], gmm::mat_nrows(M)); gmm::mult(M,w,pts[i]);
    }
    dimension = gmm::mat_nrows(M);
    pts.resort();
  }

  void mesh::bounding_box(base_node& Pmin, base_node& Pmax) const {
    bool is_first = true;
    Pmin.clear(); Pmax.clear(); 
    for (dal::bv_visitor i(pts.index()); !i.finished(); ++i) {
      if (is_first) { Pmin = Pmax = pts[i]; is_first = false; }
      else for (unsigned j=0; j < dim(); ++j) {
	Pmin[j] = std::min(Pmin[j], pts[i][j]);
	Pmax[j] = std::max(Pmax[j], pts[i][j]);
      }
    }
  }

  void mesh::clear(void) {
    dimension = dim_type(-1);
    mesh_structure::clear();
    pts.clear();
    gtab.clear(); trans_exists.clear();
    cvf_sets.clear(); valid_cvf_sets.clear();
    lmsg_sender().send(MESH_CLEAR()); touch();
  }

  size_type mesh::add_segment(size_type a, size_type b) { 
    size_type ipt[2]; ipt[0] = a; ipt[1] = b;
    return add_convex(bgeot::simplex_geotrans(1, 1), &(ipt[0]));
  }
  
  size_type mesh::add_triangle(size_type a, 
				      size_type b, size_type c) {
    size_type ipt[3]; ipt[0] = a; ipt[1] = b; ipt[2] = c;
    return add_simplex(2, &(ipt[0]));
  }
  
  size_type mesh::add_triangle_by_points
    (const base_node &p1, const base_node &p2, const base_node &p3)
  { return add_triangle(add_point(p1), add_point(p2), add_point(p3)); }
  
  size_type mesh::add_tetrahedron(size_type a, size_type b,
					 size_type c, size_type d) {
    size_type ipt[4]; ipt[0] = a; ipt[1] = b; ipt[2] = c; ipt[3] = d;
    return add_simplex(3, &(ipt[0]));
  }

  size_type mesh::add_tetrahedron_by_points(const base_node &p1,
						   const base_node &p2,
						   const base_node &p3,
						   const base_node &p4) {
    return add_tetrahedron(add_point(p1), add_point(p2),
			   add_point(p3), add_point(p4));
  }

  void mesh::sup_convex(size_type ic) {
    bgeot::mesh_structure::sup_convex(ic);
    trans_exists[ic] = false;
    sup_convex_from_regions(ic);
    lmsg_sender().send(MESH_SUP_CONVEX(ic)); touch();
  }

  void mesh::swap_convex(size_type i, size_type j) {
    if (i != j) {
      bgeot::mesh_structure::swap_convex(i,j);
      trans_exists.swap(i, j);
      gtab.swap(i,j);
      swap_convex_in_regions(i, j);
      lmsg_sender().send(MESH_SWAP_CONVEX(i, j)); touch();
    }
  }

  base_small_vector mesh::normal_of_face_of_convex(size_type ic, short_type f,
							  const base_node &pt) const {
    bgeot::pgeometric_trans pgt = trans_of_convex(ic);
    base_matrix G(dim(),pgt->nb_points());
    vectors_to_base_matrix(G,points_of_convex(ic));
    bgeot::geotrans_interpolation_context c(trans_of_convex(ic), pt, G);
    return bgeot::compute_normal(c, f);
  }


  base_small_vector mesh::normal_of_face_of_convex(size_type ic, short_type f,
							  size_type n) const {
    bgeot::pgeometric_trans pgt = trans_of_convex(ic);
    bgeot::pgeotrans_precomp pgp
      = bgeot::geotrans_precomp(pgt, &(pgt->geometric_nodes()));
    base_matrix G(dim(),pgt->nb_points());
    vectors_to_base_matrix(G,points_of_convex(ic));
    bgeot::geotrans_interpolation_context c(pgp,pgt->structure()->ind_points_of_face(f)[n], G);
    return bgeot::compute_normal(c, f);
  }

  base_matrix mesh::local_basis_of_face_of_convex(size_type ic, short_type f,
							 const base_node &pt) const {
    bgeot::pgeometric_trans pgt = trans_of_convex(ic);
    base_matrix G(dim(),pgt->nb_points());
    vectors_to_base_matrix(G,points_of_convex(ic));
    bgeot::geotrans_interpolation_context c(trans_of_convex(ic), pt, G);
    return bgeot::compute_local_basis(c, f);
  }

  base_matrix mesh::local_basis_of_face_of_convex(size_type ic, short_type f,
							 size_type n) const {
    bgeot::pgeometric_trans pgt = trans_of_convex(ic);
    bgeot::pgeotrans_precomp pgp = bgeot::geotrans_precomp(pgt, &pgt->geometric_nodes());
    base_matrix G(dim(),pgt->nb_points());
    vectors_to_base_matrix(G,points_of_convex(ic));
    bgeot::geotrans_interpolation_context c(pgp,pgt->structure()->ind_points_of_face(f)[n], G);
    return bgeot::compute_local_basis(c, f);
  }


  scalar_type  mesh::convex_quality_estimate(size_type ic) const { 
    base_matrix G;
    bgeot::vectors_to_base_matrix(G, points_of_convex(ic));
    return getfem::convex_quality_estimate(trans_of_convex(ic), G);
  }

  scalar_type  mesh::convex_radius_estimate(size_type ic) const { 
    base_matrix G;
    bgeot::vectors_to_base_matrix(G, points_of_convex(ic));
    return getfem::convex_radius_estimate(trans_of_convex(ic), G);
  }
  
  scalar_type mesh::minimal_convex_radius_estimate() const {
    if (convex_index().empty()) return 1;
    scalar_type r = convex_radius_estimate(convex_index().first_true());
    for (dal::bv_visitor cv(convex_index()); !cv.finished(); ++cv) {
      r = std::min(r, convex_radius_estimate(cv));
    }
    return r;
  }

  void mesh::copy_from(const mesh& m) {
    clear();
    bgeot::mesh_structure::operator=(m);
    eps_p = m.eps_p;
    gtab = m.gtab;
    trans_exists = m.trans_exists; 
    dimension = m.dimension;
    pts = m.pts;
    for (dal::bv_visitor i(convex_index()); !i.finished(); ++i)
      lmsg_sender().send(MESH_ADD_CONVEX(i));
  }

  struct mesh_convex_structure_loc {
    bgeot::pgeometric_trans cstruct;
    std::vector<size_type> pts;
  };

  void mesh::read_from_file(std::istream &ist) {
   
    dal::bit_vector npt;
    dal::dynamic_array<double> tmpv;
    std::string tmp;
    bool te = false, please_get = true;

    ist.precision(16);
    clear();
    ist.seekg(0);ist.clear();
    ftool::read_until(ist, "BEGIN POINTS LIST");

    while (!te) {
      if (please_get) ftool::get_token(ist, tmp); else please_get = true;

      if (!ftool::casecmp(tmp, "END"))
      { te = true; }
      else if (!ftool::casecmp(tmp, "POINT")) {
	ftool::get_token(ist, tmp);
        size_type ip = atoi(tmp.c_str());
        dim_type d = 0;
	if (npt.is_in(ip))
	  DAL_THROW(failure_error, 
		    "Two points with the same index. loading aborted.");
	npt.add(ip);
	ftool::get_token(ist, tmp);
	while (isdigit(tmp[0]) || tmp[0] == '-' || tmp[0] == '+'
	                       || tmp[0] == '.')
	  { tmpv[d++] = atof(tmp.c_str()); ftool::get_token(ist, tmp); }
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
		      "Two points [#" << ip << " and #" << ipl << "] with the same coords. loading aborted.");
	  swap_points(ip, ipl);
	}
      } else if (tmp.size()) {
	DAL_THROW(failure_error, "Syntax error in file, at token '" << tmp
		  << "', pos=" << std::streamoff(ist.tellg()));
      } else if (ist.eof()) {
	DAL_THROW(failure_error, "Unexpected end of stream");	
      }
    }

    bool tend = false;
    dal::dynamic_array<mesh_convex_structure_loc> cv;
    dal::bit_vector ncv;
    
    ist.seekg(0);
    if (!ftool::read_until(ist, "BEGIN MESH STRUCTURE DESCRIPTION"))
      DAL_THROW(failure_error, "This seems not to be a mesh file");

    while (!tend) {
      tend = !ftool::get_token(ist, tmp);
      if (!ftool::casecmp(tmp, "END"))
      { tend = true; }
      else if (!ftool::casecmp(tmp, "CONVEX"))
      {
        size_type ic;
	ftool::get_token(ist, tmp);
        ic = gmm::abs(atoi(tmp.c_str()));
	if (ncv.is_in(ic)) 
	  DAL_THROW(failure_error,
		    "Negative or repeated index, loading aborted.");
	ncv.add(ic);

	int rgt = ftool::get_token(ist, tmp);
	if (rgt != 3) { // for backward compatibility with version 1.7
	  char c; ist.get(c);
	  while (!isspace(c)) { tmp.push_back(c); ist.get(c); }
	}
	
	bgeot::pgeometric_trans pgt = bgeot::geometric_trans_descriptor(tmp);

       
	size_type nb = pgt->nb_points();

	cv[ic].cstruct = pgt;
	cv[ic].pts.resize(nb);
	for (size_type i = 0; i < nb; i++) {
	  ftool::get_token(ist, tmp);	  
	  cv[ic].pts[i] = gmm::abs(atoi(tmp.c_str()));
	}
      }
      else if (tmp.size()) {
	DAL_THROW(failure_error, "Syntax error reading a mesh file "
		  " at pos " << std::streamoff(ist.tellg())
		  << "(expecting 'CONVEX' or 'END', found '" << tmp << "')"); 
      } else if (ist.eof()) {
	DAL_THROW(failure_error, "Unexpected end of stream "
		  << "(missing BEGIN MESH/END MESH ?)");
      }
    }
    ist >> ftool::skip("MESH STRUCTURE DESCRIPTION");

    for (dal::bv_visitor ic(ncv); !ic.finished(); ++ic) {
      size_type i = add_convex(cv[ic].cstruct, cv[ic].pts.begin());
      if (i != ic) swap_convex(i, ic);
    }

    tend = false;
    while (!tend) {
      tend = !ftool::get_token(ist, tmp);
      // bool error = false;
      if (ftool::casecmp(tmp, "BEGIN")==0) {
	ftool::get_token(ist, tmp);
	if (ftool::casecmp(tmp, "BOUNDARY")==0 ||
	    ftool::casecmp(tmp, "REGION")==0) {
	  ftool::get_token(ist, tmp);
	  size_type bnum = atoi(tmp.c_str());
	  ftool::get_token(ist, tmp);
	  while (true) {
	    if (ftool::casecmp(tmp, "END")!=0) {
	      size_type ic = atoi(tmp.c_str());
	      ftool::get_token(ist, tmp);
	      if (tmp[0] == '/') {
		ftool::get_token(ist, tmp);
		if (!ftool::casecmp(tmp, "END")) break;
		size_type f = atoi(tmp.c_str());
		region(bnum).add(ic, f);
		ftool::get_token(ist, tmp);
	      }
	      else {
		region(bnum).add(ic);
		if (!ftool::casecmp(tmp, "END")) break;
	      }
	    } else break;
	  }
	  ftool::get_token(ist, tmp);
	  ftool::get_token(ist, tmp);
	} else tend = true; /*else DAL_THROW(failure_error, "Syntax error in file at token '"
			 << tmp << "' [pos=" << std::streamoff(ist.tellg())
			 << "]");*/
      } else tend=true;
    }
  }

  void mesh::read_from_file(const std::string &name) { 
    std::ifstream o(name.c_str());
    if (!o) DAL_THROW(file_not_found_error,
		      "Mesh file '" << name << "' does not exist");
    read_from_file(o); 
    o.close();
  }

  template<class ITER>
    void write_tab_to_file_(std::ostream &ost, const ITER& b_, const ITER& e)
  { for (ITER b(b_) ; b != e; ++b) ost << "  " << *b; }

  template<class ITER>
    static void write_convex_to_file_(const mesh &ms,
				      std::ostream &ost,
				      ITER b, ITER e) {
    for ( ; b != e; ++b) {
      size_type i = b.index();
      ost << "CONVEX " << i << "    \'"
	  << bgeot::name_of_geometric_trans(ms.trans_of_convex(i)).c_str()
	  << "\'    ";
      write_tab_to_file_(ost, ms.ind_points_of_convex(i).begin(),
			 ms.ind_points_of_convex(i).end()  );
      ost << '\n';
    }
  }

  template<class ITER> static void write_point_to_file_(std::ostream &ost,
						  ITER b, ITER e)
  { for ( ; b != e; ++b) ost << "  " << *b; ost << '\n'; }

  void mesh::write_to_file(std::ostream &ost) const {
    ost.precision(16);
    ost << '\n' << "BEGIN POINTS LIST" << '\n' << '\n';
    for (size_type i = 0; i < points_tab.size(); ++i)
      if (is_point_valid(i) ) {
	ost << "  POINT  " << i;
	write_point_to_file_(ost, pts[i].begin(), pts[i].end());
      }
    ost << '\n' << "END POINTS LIST" << '\n' << '\n' << '\n';
    
    ost << '\n' << "BEGIN MESH STRUCTURE DESCRIPTION" << '\n' << '\n';
    write_convex_to_file_(*this, ost, convex_tab.tas_begin(),
			              convex_tab.tas_end());
    ost << '\n' << "END MESH STRUCTURE DESCRIPTION" << '\n';

    for (dal::bv_visitor bnum(valid_cvf_sets); !bnum.finished(); ++bnum) {
      /*if (region(bnum).is_boundary()) {
	ost << " BEGIN BOUNDARY " << bnum;	
	size_type cnt = 0;
	for (dal::bv_visitor cv(cvf_sets[bnum].index()); !cv.finished();
	     ++cv) {
	  for (size_type f = 0; f < MAX_FACES_PER_CV; ++f)
	    if (cvf_sets[bnum].faces_of_convex(cv)[f]) {
	      if ((cnt++ % 10) == 0) ost << '\n' << " ";
	      ost << " " << cv << "/" << f;
	  }
	}
	ost << '\n' << " END BOUNDARY " << bnum << '\n';
      }
      else {
	ost << " BEGIN CONVEX SET " << bnum;
	size_type cnt = 0;
	for (dal::bv_visitor cv(cvf_sets[bnum].index()); !cv.finished();
	     ++cv, ++cnt) {
	  if ((cnt % 20) == 0) ost << '\n' << " ";
	  ost << " " << cv;
	}
	ost << '\n' << " END CONVEX SET " << bnum << '\n';
      }
      */
      ost << "BEGIN REGION " << bnum << "\n" << region(bnum) << "\n"
	  << "END REGION " << bnum << "\n";
    }
  }

  void mesh::write_to_file(const std::string &name) const {
    std::ofstream o(name.c_str());
    if (!o)
      DAL_THROW(failure_error, "impossible to write to file '" << name << "'");
    o << "% GETFEM MESH FILE " << '\n';
    o << "% GETFEM VERSION " << GETFEM_VERSION << '\n' << '\n' << '\n';
    write_to_file(o);
    o.close();
  }

  size_type mesh::memsize(void) const {
    return bgeot::mesh_structure::memsize() - sizeof(bgeot::mesh_structure)
      + pts.memsize() + (pts.index().last_true()+1)*dim()*sizeof(scalar_type)
      + sizeof(mesh) + trans_exists.memsize() + gtab.memsize()
      + cvf_sets.memsize() + valid_cvf_sets.memsize();
  }

  struct equilateral_to_GT_PK_grad_aux : public std::vector<base_matrix> {};
  static const base_matrix &equilateral_to_GT_PK_grad(dim_type N) {
    std::vector<base_matrix> &pbm = dal::singleton<equilateral_to_GT_PK_grad_aux >::instance();
    if (N > pbm.size()) pbm.resize(N);
    if (pbm[N-1].empty()) {
      bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(N,1);
      base_matrix Gr(N,N);
      base_matrix G(N,N+1); 
      vectors_to_base_matrix(G, bgeot::equilateral_simplex_of_reference(N)->points());
      //cout << "equilateral_simplex_of_reference: G[" << N << "]=" << G << "\n";
      gmm::mult(G, bgeot::geotrans_precomp(pgt, &pgt->convex_ref()->points())->grad(0), Gr);
      gmm::lu_inverse(Gr);
      pbm[N-1].swap(Gr);
    }
    return pbm[N-1];
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
    size_type N = G.nrows(), P = pgt->structure()->dim();
    base_matrix K(N,P);
    for (size_type ip=0; ip < n; ++ip) {
      gmm::mult(G, pgp->grad(ip), K);
      /* TODO : this is an ugly fix for simplexes only.. there should be a transformation
         of any pgt to the equivalent equilateral pgt (for prisms etc) */
      if (pgt->structure()->basic_structure() == bgeot::simplex_structure(P)) 
        gmm::mult(base_matrix(K),equilateral_to_GT_PK_grad(P),K);
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
  outer_faces_of_mesh(const getfem::mesh &m, const dal::bit_vector& cvlst, convex_face_ct& flist) {
    for (dal::bv_visitor ic(cvlst); !ic.finished(); ++ic) {
      if (m.structure_of_convex(ic)->dim() == m.dim()) {
	for (size_type f = 0; f < m.structure_of_convex(ic)->nb_faces(); f++) {
	  if (m.neighbour_of_convex(ic,f) == size_type(-1)) {
	    flist.push_back(convex_face(ic,f));
	  }
	}
      } else {
	flist.push_back(convex_face(ic,size_type(-1)));
      }
    }
  }

  void  outer_faces_of_mesh(const mesh &m, 
			    const mesh_region &cvlst,
			    mesh_region &flist) {
    cvlst.error_if_not_convexes();
    for (mr_visitor i(cvlst); !i.finished(); ++i) {
      if (m.structure_of_convex(i.cv())->dim() == m.dim()) {
	for (size_type f = 0; f < m.structure_of_convex(i.cv())->nb_faces();
	     f++) {
	  if (!m.is_convex_having_neighbour(i.cv(),f)) {
	    flist.add(i.cv(),f);
	  }
	}
      }
      else flist.add(i.cv());
    }
  }

  void extrude(const mesh& in, mesh& out, unsigned nb_layers) {
    unsigned dim = in.dim();
    base_node pt(dim+1);
    out.clear();
    size_type nbpt = in.points().index().last()+1;
    if (nbpt != in.points().index().card()) 
      DAL_THROW(dal::failure_error, "please optimize the mesh before using it as a base for prismatic mesh");
    for (size_type i = 0; i < nbpt; ++i) {
      std::copy(in.points()[i].begin(), in.points()[i].end(),pt.begin());
      pt[dim] = 0.0;
      for (size_type j = 0; j <= nb_layers; ++j, pt[dim] += 1.0 / nb_layers)
	out.add_point(pt);
    }
  
    std::vector<size_type> tab;
    for (dal::bv_visitor cv(in.convex_index()); !cv.finished(); ++cv) {
      size_type nbp = in.nb_points_of_convex(cv);
      tab.resize(2*nbp);
      for (size_type j = 0; j < nb_layers; ++j) {
	for (size_type k = 0; k < nbp; ++k)
	  tab[k] = (nb_layers+1)*in.ind_points_of_convex(cv)[k] + j;
	for (size_type k = 0; k < nbp; ++k)
	  tab[k+nbp] = (nb_layers+1)*in.ind_points_of_convex(cv)[k] + j + 1;
	bgeot::pgeometric_trans pgt = 
	  bgeot::product_geotrans(in.trans_of_convex(cv), bgeot::simplex_geotrans(1,1));      
	out.add_convex(pgt, tab.begin());
      }
    }
  }

}  /* end of namespace getfem.                                             */
