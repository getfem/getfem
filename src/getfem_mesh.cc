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
    Bank_info = 0;
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
	
#if GETFEM_PARA_LEVEL > 1
	double t_ref = MPI_Wtime();
#endif

	int j = 0, k = 0;
	for (dal::bv_visitor ic(convex_index()); !ic.finished(); ++ic, ++j) {
	  numelt[j] = ic;
	  xadj[j] = k;
	  bgeot::mesh_structure::ind_set s(neighbours_of_convex(ic));
	  for (bgeot::mesh_structure::ind_set::iterator it = s.begin();
	       it != s.end(); ++it) { adjncy.push_back(*it); ++k; }  
	}
	xadj[j] = k;

	int wgtflag = 0, edgecut, numflag = 0, options[5] = {0,0,0,0,0};




	METIS_PartGraphKway(&ne, &(xadj[0]), &(adjncy[0]), 0, 0, &wgtflag,
			    &numflag, &size, options, &edgecut, &(npart[0]));
	
	for (size_type i = 0; i < size_type(ne); ++i)
	  if (npart[i] == rank) mpi_region.add(numelt[i]);

#if GETFEM_PARA_LEVEL > 1
    cout << "Partition time "<< MPI_Wtime()-t_ref << endl;
#endif
	
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

  void mesh::intersect_with_mpi_region(mesh_region &rg) const {
    if (rg.id() == mesh_region::all_convexes().id()) { 
      rg = get_mpi_region(); 
    } else if (int(rg.id()) >= 0) { 
      rg = get_mpi_sub_region(rg.id()); 
    } else 
      rg = mesh_region::intersection(rg, get_mpi_region());    
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

  void mesh::translation(base_small_vector V) {
    for (dal::bv_visitor i(pts.index()); !i.finished(); ++i)
      pts[i] += V;
    pts.resort();
  }

  void mesh::transformation(base_matrix M) {
    base_small_vector w(M.nrows());
    if (gmm::mat_nrows(M) == 0 || gmm::mat_ncols(M) != dim()) 
      DAL_THROW(dal::dimension_error,
		"invalid dimensions for the transformation matrix");
    for (dal::bv_visitor i(pts.index()); !i.finished(); ++i) {
      w = pts[i];
      gmm::resize(pts[i], gmm::mat_nrows(M));
      gmm::mult(M,w,pts[i]);
    }
    dimension = gmm::mat_nrows(M);
    if (Bank_info) { delete Bank_info; Bank_info = 0; }
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
    if (Bank_info) { delete Bank_info; Bank_info = 0; }
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

  size_type mesh::add_tetrahedron_by_points
  (const base_node &p1, const base_node &p2,
   const base_node &p3, const base_node &p4) {
    return add_tetrahedron(add_point(p1), add_point(p2),
			   add_point(p3), add_point(p4));
  }

  void mesh::sup_convex(size_type ic, bool sup_points) {
    static std::vector<size_type> ipt;
    if (sup_points) {
      const ind_cv_ct &ct = ind_points_of_convex(ic);
      ipt.assign(ct.begin(), ct.end()); 
    }
    bgeot::mesh_structure::sup_convex(ic);
    if (sup_points)
      for (size_type ip = 0; ip < ipt.size(); ++ip) sup_point(ipt[ip]);
    trans_exists.sup(ic);
    sup_convex_from_regions(ic);
    if (Bank_info) Bank_sup_convex_from_green(ic);
    lmsg_sender().send(MESH_SUP_CONVEX(ic)); touch();
  }

  void mesh::swap_convex(size_type i, size_type j) {
    if (i != j) {
      bgeot::mesh_structure::swap_convex(i,j);
      trans_exists.swap(i, j);
      gtab.swap(i,j);
      swap_convex_in_regions(i, j);
      if (Bank_info) Bank_swap_convex(i,j);
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
    bgeot::geotrans_interpolation_context
      c(pgp,pgt->structure()->ind_points_of_face(f)[n], G);
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
    bgeot::pgeotrans_precomp pgp
      = bgeot::geotrans_precomp(pgt, &pgt->geometric_nodes());
    base_matrix G(dim(),pgt->nb_points());
    vectors_to_base_matrix(G,points_of_convex(ic));
    bgeot::geotrans_interpolation_context
      c(pgp,pgt->structure()->ind_points_of_face(f)[n], G);
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
    if (Bank_info)
      delete Bank_info;      
    if (m.Bank_info) {
      Bank_info = new Bank_info_struct;
      *Bank_info = *(m.Bank_info);
    }
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
		      "Two points [#" << ip << " and #" << ipl
		      << "] with the same coords. loading aborted.");
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
      else if (!ftool::casecmp(tmp, "CONVEX")) {
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
	} else tend = true;
	/*else DAL_THROW(failure_error, "Syntax error in file at token '"
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
    std::vector<base_matrix>
      &pbm = dal::singleton<equilateral_to_GT_PK_grad_aux >::instance();
    if (N > pbm.size()) pbm.resize(N);
    if (pbm[N-1].empty()) {
      bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(N,1);
      base_matrix Gr(N,N);
      base_matrix G(N,N+1); 
      vectors_to_base_matrix
	(G, bgeot::equilateral_simplex_of_reference(N)->points());
      gmm::mult(G, bgeot::geotrans_precomp
		(pgt, &pgt->convex_ref()->points())->grad(0), Gr);
      gmm::lu_inverse(Gr);
      pbm[N-1].swap(Gr);
    }
    return pbm[N-1];
  }
    
  /* TODO : use the geotrans from an "equilateral" reference element to
     the real element 
     check if the sign of the determinants does change
     (=> very very bad quality of the convex)
  */
  scalar_type convex_quality_estimate(bgeot::pgeometric_trans pgt,
				      const base_matrix& G) {
    static bgeot::pgeometric_trans pgt_old = 0;
    static bgeot::pgeotrans_precomp pgp = 0;
    if (pgt_old != pgt) {
      pgt_old=pgt;
      pgp=bgeot::geotrans_precomp(pgt, &pgt->convex_ref()->points());
    }

    size_type n = (pgt->is_linear()) ? 1 : pgt->nb_points();
    scalar_type q = 1;
    size_type N = G.nrows(), P = pgt->structure()->dim();
    base_matrix K(N,P);
    for (size_type ip=0; ip < n; ++ip) {
      gmm::mult(G, pgp->grad(ip), K);
      /* TODO : this is an ugly fix for simplexes only.. there should be
	 a transformation of any pgt to the equivalent equilateral pgt
	 (for prisms etc) */
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
      pgt_old=pgt;
      pgp=bgeot::geotrans_precomp(pgt, &pgt->convex_ref()->points());
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
  outer_faces_of_mesh(const getfem::mesh &m, const dal::bit_vector& cvlst,
		      convex_face_ct& flist) {
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
      DAL_THROW(dal::failure_error, "please optimize the mesh before using "
		"it as a base for prismatic mesh");
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
	  bgeot::product_geotrans(in.trans_of_convex(cv),
				  bgeot::simplex_geotrans(1,1));      
	out.add_convex(pgt, tab.begin());
      }
    }
  }
}


// 
// Bank refinement
//

namespace bgeot {
  size_type refinement_simplexe_tab(size_type n, size_type **tab);
}

namespace getfem {
  
  bool mesh::edge::operator <(const edge &e) const {
    if (i1 < e.i1) return true; if (i1 > e.i1) return false;
    if (i2 < e.i2) return true; return false;
  }
  
  void mesh::Bank_sup_convex_from_green(size_type i) {
    if (Bank_info && Bank_info->is_green_simplex.is_in(i)) {
      size_type igs = Bank_info->num_green_simplex[i];
      green_simplex &gs = Bank_info->green_simplices[igs];
      for (size_type j = 0; j < gs.sub_simplices.size(); ++j) {
	Bank_info->num_green_simplex.erase(gs.sub_simplices[j]);
	Bank_info->is_green_simplex.sup(gs.sub_simplices[j]);
      }
      Bank_info->green_simplices.sup(igs);
    }
  }

  void mesh::Bank_swap_convex(size_type i, size_type j) {
    if (Bank_info) {
      Bank_info->is_green_simplex.swap(i, j);
      std::map<size_type, size_type>::iterator
	iti = Bank_info->num_green_simplex.find(i);
      std::map<size_type, size_type>::iterator
	ite = Bank_info->num_green_simplex.end();
      std::map<size_type, size_type>::iterator
	itj = Bank_info->num_green_simplex.find(j);
      size_type numi(0), numj(0);
      if (iti != ite)
	{ numi = iti->second; Bank_info->num_green_simplex.erase(i); }
      if (itj != ite)
	{ numj = itj->second; Bank_info->num_green_simplex.erase(j); }
      if (iti != ite) { 
	Bank_info->num_green_simplex[j] = numi;
	green_simplex &gs = Bank_info->green_simplices[numi];
	for (size_type k = 0; k < gs.sub_simplices.size(); ++k)
	  if (gs.sub_simplices[k] == i) gs.sub_simplices[k] = j;
	  else if (gs.sub_simplices[k] == j) gs.sub_simplices[k] = i;
      }
      if (itj != ite) {
	Bank_info->num_green_simplex[i] = numj;
	if (iti == ite || numi != numj) {
	  green_simplex &gs = Bank_info->green_simplices[numj];
	  for (size_type k = 0; k < gs.sub_simplices.size(); ++k)
	    if (gs.sub_simplices[k] == i) gs.sub_simplices[k] = j;
	    else if (gs.sub_simplices[k] == j) gs.sub_simplices[k] = i;
	}
      }
    }
  }

  void mesh::Bank_build_first_mesh(mesh &m, size_type n) {
    bgeot::pconvex_ref pcr = bgeot::simplex_of_reference(n, 2);
    m.clear(); 
    for (size_type ip = 0; ip < pcr->nb_points(); ++ip)
      m.add_point(pcr->points()[ip]);
    size_type *tab;
    size_type nbc = bgeot::refinement_simplexe_tab(n, &tab);
    for (size_type ic = 0; ic < nbc; ++ic, tab+=(n+1))
      m.add_simplex(n, tab);
  }

  void mesh::Bank_basic_refine_convex(size_type i) {
    bgeot::pgeometric_trans pgt = trans_of_convex(i);
    size_type n = pgt->basic_structure()->dim();

    static bgeot::pgeometric_trans pgt1 = 0;
    static mesh mesh2;
    static bgeot::pstored_point_tab pspt = 0;
    static bgeot::pgeotrans_precomp pgp = 0;
    static std::vector<size_type> ipt, ipt2, icl;

    if (pgt != pgt1) {
      pgt1 = pgt;
      mesh mesh1;
      Bank_build_first_mesh(mesh1, n);
      
      mesh2.clear();
      ipt.resize(pgt->nb_points());
      
      for (size_type ic = 0; ic < mesh1.nb_convex(); ++ic) {
	
	bgeot::pgeometric_trans pgt2 = mesh1.trans_of_convex(ic);
	for (size_type ip = 0; ip < pgt->nb_points(); ++ip)
	  ipt[ip] =
	    mesh2.add_point(pgt2->transform(pgt->geometric_nodes()[ip],
					    mesh1.points_of_convex(ic)));
	mesh2.add_convex(pgt, &ipt[0]);
      }

      pspt = bgeot::store_point_tab(mesh2.points());
      pgp = bgeot::geotrans_precomp(pgt, pspt);
    }
    
    base_node pt(n);
    ipt.resize(pspt->size());
    for (size_type ip = 0; ip < pspt->size(); ++ip) {
      pgp->transform(points_of_convex(i), ip, pt);
      ipt[ip] = add_point(pt);
    }

    ipt2.resize(n+1); icl.resize(0);
    for (size_type ic = 0; ic < mesh2.nb_convex(); ++ic) {
      for (size_type j = 0; j <= n; ++j)
	ipt2[j] = ipt[mesh2.ind_points_of_convex(ic)[j]];
      icl.push_back(add_convex(pgt, ipt2.begin()));
    }

    lmsg_sender().send(MESH_REFINE_CONVEX(i, icl, true));
    sup_convex(i, true);
  }

  std::vector<size_type> &mesh::Bank_loc_ind_of_pgt
  (bgeot::pgeometric_trans pgt) {
    static std::vector<size_type> loc_ind;
    static bgeot::pgeometric_trans pgt0 = 0;
    if (pgt != pgt0) {
      pgt0 = pgt;
      loc_ind.resize(0);
      for (size_type ip = 0; ip < pgt->nb_points(); ++ip) {
	scalar_type ni=gmm::vect_norminf(pgt->geometric_nodes()[ip]);
	if (gmm::abs(ni) < 1e-12 || gmm::abs(ni - 1.0) < 1e-12)
	  loc_ind.push_back(ip);
      }
    }
    return loc_ind;
  }

  void mesh::Bank_refine_normal_convex(size_type i) {
    bgeot::pgeometric_trans pgt = trans_of_convex(i);
    size_type n = pgt->basic_structure()->dim();
    if (pgt->basic_structure() != bgeot::simplex_structure(n))
      DAL_THROW(failure_error, "Sorry, refinement is only working "
		"with simplices.");
    
    std::vector<size_type> &loc_ind = Bank_loc_ind_of_pgt(pgt);
    
    for (size_type ip1 = 0; ip1 < n; ++ip1)
      for (size_type ip2 = ip1+1; ip2 <= n; ++ip2) {
	size_type i1 = ind_points_of_convex(i)[loc_ind[ip1]];
	size_type i2 = ind_points_of_convex(i)[loc_ind[ip2]];
	Bank_info->edges.insert(edge(i1, i2));
      }
    Bank_basic_refine_convex(i);
  }

  void mesh::Bank_test_and_refine_convex(size_type i,
					 dal::bit_vector &b) {
    if (Bank_info->is_green_simplex[i]) {
      size_type igs = Bank_info->num_green_simplex[i];
      green_simplex &gs = Bank_info->green_simplices[igs];
      size_type icc = add_convex_by_points(gs.pgt, gs.cv.points().begin());
      for (size_type ic = 0; ic < gs.sub_simplices.size(); ++ic) {
	sup_convex(gs.sub_simplices[ic], true);
	b.sup(gs.sub_simplices[ic]);
      }
      lmsg_sender().send(MESH_REFINE_CONVEX(i, gs.sub_simplices, false));
      Bank_sup_convex_from_green(i);
      Bank_refine_normal_convex(icc);
    }
    else
      Bank_refine_normal_convex(i);
  }

  void mesh::Bank_build_green_simplexes(size_type ic,
					std::vector<size_type> &ipt) {
    size_type igs = Bank_info->green_simplices.add(green_simplex());
    green_simplex &gs(Bank_info->green_simplices[igs]);
    std::vector<base_node> pt_tab(nb_points_of_convex(ic));
    ref_mesh_pt_ct ptab = points_of_convex(ic);
    pt_tab.assign(ptab.begin(), ptab.end());
    gs.cv = bgeot::convex<base_node>(structure_of_convex(ic), pt_tab);

    bgeot::pgeometric_trans pgt = gs.pgt = trans_of_convex(ic);

    size_type d = ipt.size() - 1, n = structure_of_convex(ic)->dim();
    static size_type d0 = 0;
    static mesh mesh1;
    if (d0 != d) {
      d0 = d;
      Bank_build_first_mesh(mesh1, d);    
    }
    
    std::vector<size_type> &loc_ind = Bank_loc_ind_of_pgt(pgt);
    const bgeot::mesh_structure::ind_cv_ct &ct = ind_points_of_convex(ic);

    std::vector<size_type> ipt_loc(ipt.size()), ipt_other;
    for (size_type ip = 0; ip < loc_ind.size(); ++ip) {
      bool found = false;
      for (size_type i = 0; i < ipt.size(); ++i)
	if (ct[loc_ind[ip]] == ipt[i])
	  { ipt_loc[i] = ip; found = true; break; }
      if (!found) ipt_other.push_back(ip);
    }
    
    mesh mesh2;
    bgeot::pconvex_ref pcr = bgeot::simplex_of_reference(n);
    size_type ic0 = mesh2.add_simplex_by_points(n, pcr->points().begin());
    size_type ic1 = mesh2.add_simplex(d, ipt_loc.begin());
    bgeot::pgeometric_trans pgt1 = mesh2.trans_of_convex(ic1);

    bgeot::pstored_point_tab pspt = bgeot::store_point_tab(mesh1.points());
    bgeot::pgeotrans_precomp pgp = bgeot::geotrans_precomp(pgt1, pspt);
  
    std::vector<size_type> ipt1(pspt->size());
    base_node pt(n);
    for (size_type i = 0; i < pspt->size(); ++i) {
      pgp->transform(mesh2.points_of_convex(ic1), i, pt);
      ipt1[i] = mesh2.add_point(pt);
    }
    mesh2.sup_convex(ic1);
    
    std::vector<size_type> ipt2(n+1);
    for (size_type i = 0; i < mesh1.nb_convex(); ++i) {
      for (size_type j = 0; j <= d; ++j)
	ipt2[j] = ipt1[mesh1.ind_points_of_convex(i)[j]];
      for (size_type j = d+1; j <= n; ++j)
	ipt2[j] = ipt_other[j-d-1];
      mesh2.add_simplex(n, ipt2.begin());
    }
    mesh2.sup_convex(ic0);
    

    mesh mesh3;    
    ipt1.resize(pgt->nb_points());  
    for (dal::bv_visitor i(mesh2.convex_index()); !i.finished(); ++i) {
      bgeot::pgeometric_trans pgt2 = mesh2.trans_of_convex(i);
      for (size_type ip = 0; ip < pgt->nb_points(); ++ip)
	ipt1[ip] =
	  mesh3.add_point(pgt2->transform(pgt->geometric_nodes()[ip],
					  mesh2.points_of_convex(i)));
      mesh3.add_convex(pgt, ipt1.begin());
    }
        
    
    pspt = bgeot::store_point_tab(mesh3.points());
    pgp = bgeot::geotrans_precomp(pgt, pspt);
    
    ipt1.resize(pspt->size());
    for (size_type ip = 0; ip < pspt->size(); ++ip) {
      pgp->transform(points_of_convex(ic), ip, pt);
      ipt1[ip] = add_point(pt);
    }
    
    ipt2.resize(n+1);
    for (size_type icc = 0; icc < mesh3.nb_convex(); ++icc) {
      for (size_type j = 0; j <= n; ++j)
	ipt2[j] = ipt1[mesh3.ind_points_of_convex(icc)[j]];
      size_type i = add_convex(pgt, ipt2.begin());
      gs.sub_simplices.push_back(i);
      Bank_info->is_green_simplex.add(i);
      Bank_info->num_green_simplex[i] = igs;
    }

    for (size_type ip1 = 0; ip1 < ipt.size(); ++ip1)
      for (size_type ip2 = ip1+1; ip2 < ipt.size(); ++ip2)
	Bank_info->edges.insert(edge(ipt[ip1], ipt[ip2]));
    
    lmsg_sender().send(MESH_REFINE_CONVEX(ic, gs.sub_simplices, true));
    sup_convex(ic, true);
  }

  void mesh::Bank_refine(dal::bit_vector b) {

    if (Bank_info == 0) Bank_info = new Bank_info_struct;

    Bank_info->edges.clear();
    while (b.card() > 0) Bank_test_and_refine_convex(b.take_first(), b);

    std::vector<size_type> ipt;
    edge_set marked_convexes;
    while (Bank_info->edges.size()) {
      marked_convexes.clear();
      b = convex_index();
      edge_set::const_iterator it = Bank_info->edges.begin();
      edge_set::const_iterator ite = Bank_info->edges.end(), it2=it;
      ++it2;
      for (; it != ite; it = it2) {
	convex_with_edge(it->i1, it->i2, ipt);
	if (ipt.size() == 0) Bank_info->edges.erase(it);
	else {
	  for (size_type ic = 0; ic < ipt.size(); ++ic)
	    marked_convexes.insert(edge(ipt[ic], it->i1, it->i2));
	}
	if (it2 != ite) ++it2;
      }
      
      it = marked_convexes.begin(); ite = marked_convexes.end();
      while (it != ite) {
	size_type ic = it->i0;
	ipt.resize(0);
	while (it != ite && it->i0 == ic) {
	  bool found1 = false, found2 = false;
	  for (size_type j = 0; j < ipt.size(); ++j)
	    if (ipt[j] == it->i1) found1 = true;
	    else if (ipt[j] == it->i2) found2 = true;
	  if (!found1) ipt.push_back(it->i1);
	  if (!found2) ipt.push_back(it->i2);
	  ++it;
	}
	if (b.is_in(ic)) {
	  if (ipt.size() > structure_of_convex(ic)->dim()
	      || Bank_info->is_green_simplex[ic])
	    Bank_test_and_refine_convex(ic, b);
	  else
	    Bank_build_green_simplexes(ic, ipt);
	}
      }
    }
    Bank_info->edges.clear();
  }






}  /* end of namespace getfem.                                             */
