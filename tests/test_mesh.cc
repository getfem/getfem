/*===========================================================================

 Copyright (C) 2002-2020 Yves Renard.

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/
#include "getfem/getfem_regular_meshes.h"
#include "getfem/bgeot_poly_composite.h"
#include "getfem/bgeot_comma_init.h"
#include "getfem/getfem_export.h"
#include "getfem/bgeot_node_tab.h"
using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;
using getfem::size_type;
using getfem::base_node;
using getfem::base_small_vector;


void export_mesh(getfem::mesh &m, const std::string &name) {
  getfem::mesh_fem mf(m);
  mf.set_classical_finite_element(m.convex_index(), 0);
  std::vector<double> U(mf.nb_basic_dof());
  for (size_type i = 0; i < mf.nb_basic_dof(); ++i)
    U[i] = double(mf.first_convex_of_basic_dof(i));

  getfem::vtk_export exp(name + ".vtk", false);
  exp.exporting(mf); 
  exp.write_point_data(mf, U, "mesh");
  cout << "export done, you can view the data file with (for example)\n"
    "mayavi -d " << name << ".vtk -m Outline -m BandedSurfaceMap\n";
}

typedef base_node POINT;
typedef base_small_vector VECT;


void test_conforming(getfem::mesh &m) {
  size_type dim = m.dim();
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(m, border_faces);
  size_type nb_faces = 0;
  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i, nb_faces++) {
    for (size_type ip = 0; ip < dim; ++ip) {
      bool onbound = false;
      const POINT &pt = m.points_of_face_of_convex(i.cv(), i.f())[ip];
      for (size_type nc = 0; nc < dim; ++nc)
	if (gmm::abs(pt[nc]) < 1E-10 || gmm::abs(pt[nc] - 1.0) < 1E-10)
	  onbound = true;
      if (!onbound)
	GMM_ASSERT1(false, "Mesh is not conforming " << pt);
    }
  }
  cout << "nb faces of mesh of dim " << dim << " : " << nb_faces << endl;

}


void test_refinable(unsigned dim, unsigned degree) {

  int Nsubdiv = 2;
  std::vector<size_type> nsubdiv(dim, Nsubdiv);
  getfem::mesh m;
  getfem::regular_unit_mesh
    (m, nsubdiv, bgeot::simplex_geotrans(dim, bgeot::short_type(degree)), false);
  dal::bit_vector b; b.add(0);
  
  cout << "\nrefine mesh in dimension " << dim << " and degree "
       << degree << endl;

  // m.write_to_file(cout); getchar();

  // if (m.dim() < 4) export_mesh(m, "test_mesh_bank_ref0");
  m.Bank_refine(b);
  // if (m.dim() < 4) export_mesh(m, "test_mesh_bank_ref1");
  m.Bank_refine(m.convex_index());
  // if (m.dim() < 4) export_mesh(m, "test_mesh_bank_ref2");
  m.Bank_refine(b);
  // if (m.dim() < 4) export_mesh(m, "test_mesh_bank_ref3");
  m.Bank_refine(m.convex_index());
  // if (m.dim() < 4) export_mesh(m, "test_mesh_bank_ref4");
  m.Bank_refine(b);
  // if (m.dim() < 4) export_mesh(m, "test_mesh_bank_ref5");
  //  m.Bank_refine(m.convex_index());
  // if (m.dim() < 4) export_mesh(m, " test_mesh_bank_ref2");
  test_conforming(m);
}

void test_mesh_matching(size_type dim) {
  
  getfem::mesh m;

  int Nsubdiv = 2;
  std::vector<VECT> vects(dim);
  std::vector<int> iref(dim);
  POINT pt1(dim);
  for (size_type i = 0; i < dim; i++)
  { vects[i] = VECT(dim); vects[i][i] = 1.0 / Nsubdiv; iref[i] = Nsubdiv; }
  

  getfem::parallelepiped_regular_simplex_mesh(m, bgeot::dim_type(dim), pt1, vects.begin(),
					      iref.begin());

  test_conforming(m);
}


void test_mesh(getfem::mesh &m) {
    
  POINT pt1, pt2, pt3;
  if (pt1.size() == 0) pt1 = POINT(3);

  gmm::clear(pt1); pt2 = pt1; pt3 = pt1;
  
  pt2[0] = 1E-13;
  pt3[0] = 1.0;
  size_t i1 = m.add_point(pt1);
  size_t i3 = m.add_point(pt3);
  size_t i2 = m.add_point(pt2);
  
  int dim = m.dim();

  cout << "dimension of mesh " << dim << endl;
  cout << "Nomber of points : " << m.nb_points() << endl;
  assert(m.nb_points() == 2);
  m.sup_point(i1);

  cout << "Nomber of points : " << m.nb_points() << endl;
  assert(m.nb_points() == 1);

  m.add_segment(i2, i3);

  std::vector<POINT> pts(3);
  std::fill(pts.begin(), pts.end(), pt1);
  
  (pts[0])[0] = 2;
  (pts[1])[1] = 2;
  (pts[2])[2] = 2;

  bgeot::convex<POINT, std::vector<POINT> >
    cv(bgeot::simplex_structure(2), pts);

  cout << "point 0 of convex " << (cv.points())[0] << endl;
  cout << "point 1 of convex " << (cv.points())[1] << endl;
  cout << "point 2 of convex " << (cv.points())[2] << endl;

  size_t ic2 = m.add_convex_by_points(bgeot::simplex_geotrans(2,1),pts.begin());

  size_t ic3 = m.add_convex_by_points(bgeot::simplex_geotrans(2,1),pts.begin());

  assert(ic2 == ic3);

  cout << "second point of convex 0 : ";
  cout << (m.convex(0).points())[1] << endl;

  cout << "second point of convex 1 : ";
  cout << (m.convex(1).points())[1] << endl;

  std::vector<VECT> vects(dim);
  std::vector<int> iref(dim);
  gmm::clear(pt1);
  for (int i = 0; i < dim; i++)
    { vects[i] = VECT(dim); gmm::clear(vects[i]);; vects[i][i] = 1.0; iref[i] = 3; }
  

  getfem::parallelepiped_regular_simplex_mesh(m, bgeot::dim_type(dim),pt1, vects.begin(), iref.begin());
  m.write_to_file("test_mesh.mesh");
  getfem::mesh m3; m3.read_from_file("test_mesh.mesh");
  m.copy_from(m3);

  cout << "0 before set_exists(3): " << m.has_region(3) << "\n";
  m.region(3).is_in(87);
  cout << "0 after  set_exists(3): " << m.has_region(3) << "\n";
  cout << "1 before set_exists(3): " << m.has_region(3) << "\n";
  m.region(3).add(0);
  cout << "1 after  set_exists(3): " << m.has_region(3) << "\n";

  m.region(3).add(0);
  m.region(3).add(3,1);
  m.region(3).add(2);
  m.region(4).add(0);
  m.region(4).add(5);
  cout << "m.region(3)=" << m.region(3) << "\n";
  assert(m.region(3).is_in(0));
  m.sup_convex(0);
  m.sup_convex(1);
  assert(!m.region(3).is_in(0));
  m.optimize_structure(false);
  m.write_to_file("test_mesh.mesh");

  getfem::mesh m2; m2.read_from_file("test_mesh.mesh");
  assert(m.convex_index().card() == m2.convex_index().card());
  cout << "m2.regions_index=" << m2.regions_index() << "\n";
  assert(m2.regions_index().is_in(3));
  assert(m2.regions_index().is_in(4));
  assert(m2.region(4).index().card() == 1);
  assert(m2.region(4).is_in(5));
  assert(!m2.region(4).is_in(0,0));
  assert(m2.region(3).index().card() == 2);
  assert(m2.region(3).is_in(3,1));
  assert(m2.region(3).is_in(2));
  assert(!m2.region(3).is_in(0));
  //m.write_to_file(cout);
}

void
print_mesh_structure(const bgeot::mesh_structure *ms) {
  cout << "nb_pts=" << ms->nb_max_points() << ", nb_cvs="
       << ms->nb_convex() << endl;
  dal::bit_vector bv = ms->convex_index();
  size_type cv;
  for (cv << bv; cv != size_type(-1); cv << bv) {
    cout << "cv " << cv << " : dim=" << ms->structure_of_convex(cv)->dim() << ", nfaces=" << ms->structure_of_convex(cv)->nb_faces() << ", pts = {";
    for (size_type i=0; i < ms->nb_points_of_convex(cv); ++i) {
      if (i != 0) cout << ",";
      cout << ms->ind_points_of_convex(cv)[i];
    }
    cout << "}\n";
  }
}

void 
test_convex_simplif(void) {
  bgeot::pconvex_ref sr2 = bgeot::simplex_of_reference(2,1);
  bgeot::pconvex_ref sr3 = bgeot::simplex_of_reference(3,2);
  bgeot::pconvex_ref pr = bgeot::parallelepiped_of_reference(3);
  
  cout << "sr2->points() = faces:" << sr2->structure()->nb_faces() << " ; ";
  for (size_type i=0; i < sr2->points().size(); ++i) { if (i) cout << ","; cout << sr2->points()[i]; }
  cout << endl;

  cout << "sr3->points() = faces:" << sr3->structure()->nb_faces() << " ; ";
  for (size_type i=0; i < sr3->points().size(); ++i) { if (i) cout << ","; cout << sr3->points()[i]; }
  cout << endl;

  cout << "sr3=" << sr3 << endl;
  bgeot::pconvex_ref t=basic_convex_ref(sr3);
  cout << "sr3->basic_convex_ref()=" << t << endl;

  cout << "sr3->basic_convex_ref() = faces:" << sr3->structure()->nb_faces() << " ; ";
  for (size_type i=0; i < basic_convex_ref(sr3)->points().size(); ++i) { if (i) cout << ","; cout << basic_convex_ref(sr3)->points()[i]; }
  cout << endl;

  const bgeot::mesh_structure *msr2 = basic_convex_ref(sr2)->simplexified_convex();
  print_mesh_structure(msr2);
  const bgeot::mesh_structure *msr3 = basic_convex_ref(sr3)->simplexified_convex();
  print_mesh_structure(msr3);
  const bgeot::mesh_structure *psr = basic_convex_ref(pr)->simplexified_convex();
  print_mesh_structure(psr);
  const bgeot::basic_mesh *msrr2 = bgeot::refined_simplex_mesh_for_convex(sr2,2);
  getfem::mesh(*msrr2).write_to_file(cout);
  const bgeot::basic_mesh *msrr3 = bgeot::refined_simplex_mesh_for_convex(sr3,3);
  getfem::mesh(*msrr3).write_to_file(cout);
  const bgeot::basic_mesh *mprr = bgeot::refined_simplex_mesh_for_convex(pr,2);
  getfem::mesh(*mprr).write_to_file(cout);
}

void test_convex_quality(getfem::scalar_type dx, getfem::scalar_type dy) {
  getfem::mesh m;
  getfem::base_node A,B,C;
  getfem::scalar_type h = sqrt(3.0);
  B = {0,0};
  C = {2,0};
  A = {0,h};
  cout << "quality of triangles, and their radius estimates "
       << " (should decrease): ";
  for (size_type i=0; i < 10; ++i) {    
    size_type cv = m.add_triangle_by_points(A,B,C);
    cout << "Q=" << m.convex_quality_estimate(cv) << "\t"
	 << "R=" << m.convex_radius_estimate(cv) << endl;
    A[0] += dx; A[1] += dy;
  }
}

void test_region() {
  getfem::mesh_region a,b,r;
  a.add(4);
  a.add(9);
  a.add(2);
  a.add(5);
  a.add(3,7);
  a.add(3,3);
  b.add(2);
  b.add(3,2);
  b.add(3,7);
  b.add(9,1);
  b.add(9,5);
  b.add(8);
  r = getfem::mesh_region::intersection(a,b);
  cout << "a=" << a << "\nb=" << b << "a inter b=" << r << "\n";
}

void test_convex_ref() {
  for (bgeot::short_type k=1; k <= 2; ++k) {
    bgeot::pconvex_ref cvr  = bgeot::simplex_of_reference(1,k);
    base_node P(1); P[0] = .5;
    assert(gmm::abs(cvr->is_in(P)+.5) < 1e-6);
    P[0] = -.1;  assert(gmm::abs(cvr->is_in(P)-.1) < 1e-6);
    P[0] = 1.1;  assert(gmm::abs(cvr->is_in(P)-.1) < 1e-6);
    
    cvr  = bgeot::simplex_of_reference(2,k);
    assert(gmm::abs(cvr->is_in(base_node(0,0))) < 1e-6);
    assert(gmm::abs(cvr->is_in(base_node(0.5,0.5))) < 1e-6);
    assert(cvr->is_in(base_node(0.25,0.25)) < -.2);
    assert(cvr->is_in(base_node(0.85,0.85)) > .2);
    assert(cvr->is_in(base_node(-0.25,0.05)) > .2);
    assert(cvr->is_in(base_node(0.05, -0.25)) > .2);

    cvr = bgeot::parallelepiped_of_reference(3);
    assert(gmm::abs(cvr->is_in(base_node(.5,.5,.5))+.5) < 1e-6);
    assert(gmm::abs(cvr->is_in(base_node(-.5,-.5,-.5))-.5) < 1e-6);
  }
}



//struct basic_mesh_point_comparator2
//  : public std::binary_function<base_node, base_node, int> {
//  double eps;
//  std::vector<double> v;
//
//  int operator()(const base_node &x, const base_node &y) const {
//    double a = gmm::vect_sp(x, v), b =  gmm::vect_sp(y, v);
//    if (a < b - eps) return -1; else if (a > b + eps) return +1; else return 0;
//  }
//
//  basic_mesh_point_comparator2(unsigned dim_ = 3, double e = double(10000)
//			       *gmm::default_tol(double()))
//    : eps(e), v(dim_) {
//    gmm::fill_random(v);
//    gmm::scale(v, 1.0/gmm::vect_norm2(v));
//  }
//};


void test_mesh_building(int dim, int Nsubdiv) {

  double exectime = gmm::uclock_sec();

  std::vector<size_type> nsubdiv(dim, Nsubdiv);
  getfem::mesh m;
  bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(dim, 1);
  getfem::regular_unit_mesh(m, nsubdiv, pgt,false);

  cout << "Time to create mesh : "
       << gmm::uclock_sec() - exectime << endl;

//   // Test of a simple sort on one "random" component
//   exectime = gmm::uclock_sec();
//   typedef basic_mesh_point_comparator2 pt_comp2;
//   typedef dal::dynamic_tree_sorted<base_node, pt_comp2> PT_TAB2;
  
//   PT_TAB2 pttab2 = PT_TAB2(pt_comp2(dim));

//   for (dal::bv_visitor i(m.convex_index()); !i.finished(); ++i) {
//     for (int j = 0; j <= dim; ++j)
//       pttab2.add_norepeat(m.points_of_convex(i)[j]);
//   }
//   cout << "Time to insert points in structure with random component: "
//        << gmm::uclock_sec() - exectime << endl;
//   GMM_ASSERT1(pttab2.card() == m.nb_points(),
// 	      "Problem in identifying points " << pttab2.card()
// 	      << " : " << m.nb_points());
  

//   // Test on the new structure
//   exectime = gmm::uclock_sec();
  
//   bgeot::node_tab pttab3;

//   for (dal::bv_visitor i(m.convex_index()); !i.finished(); ++i) {
//     for (int j = 0; j <= dim; ++j)
//       pttab3.add_node(m.points_of_convex(i)[j]);
//   }
//   cout << "Time to insert points in new structure: "
//        << gmm::uclock_sec() - exectime << endl;
//   GMM_ASSERT1(pttab3.card() == m.nb_points(),
// 	      "Problem in identifying points " << pttab3.card()
// 	      << " : " << m.nb_points());
  


//   // Test of a simple sort on one "random" componentwith Bkdtree
//   exectime = gmm::uclock_sec();
//   std::vector<base_node> T0;
//   bgeot::kdtree trees[30];
//   bgeot::kdtree_tab_type ipts;
  
//   for (dal::bv_visitor ic(m.convex_index()); !ic.finished(); ++ic) {
//     for (int jp = 0; jp <= dim; ++jp) {

//       // insert a point.
//       T0.push_back(m.points_of_convex(ic)[jp]);

//       if (T0.size() >= 4) {
// 	size_type i = 0;
// 	for (; i < 30; ++i)
// 	  if (trees[i].nb_points() == 0) break;

// 	// cout << "building tree " << i << endl;
	
// 	for (size_type j = 0; j < 4; ++j)
// 	  trees[i].add_point(T0[j]);
// 	T0.resize(0);

// 	for (size_type j = 0; j < i; ++j) {
// 	  for (size_type k = 0; k < trees[j].nb_points(); ++k)
// 	    trees[i].add_point(trees[j].points()[k].n);
// 	  trees[j].clear();
// 	}
// 	trees[i].points_in_box(ipts, base_node(dim), base_node(dim));
//       }
      

//     }
//   }

//   cout << "Time to insert points in bkdtree: "
//        << gmm::uclock_sec() - exectime << endl;



  getfem::mesh_fem mf1(m);
  mf1.set_finite_element( getfem::classical_fem(pgt, 1) );
  cout << "nb points = " << m.nb_points()
       << "nb dof = " << mf1.nb_dof() <<  endl;
  GMM_ASSERT1(m.nb_points() == mf1.nb_dof(),
	      "Problem in identifying dofs");

  // shake a little bit the mesh
  getfem::mesh::PT_TAB &pts = m.points();
  for (size_type i = 0; i < m.nb_points(); ++i)
    for (int k = 0; k < dim; ++k)
      pts[i][k] += gmm::random(double()) * 1e-10;
  getfem::mesh_fem mf2(m);
  mf2.set_finite_element( getfem::classical_fem(pgt, 1) );
  cout << "nb points = " << m.nb_points()
       << "nb dof = " << mf2.nb_dof() <<  endl;
  GMM_ASSERT1(m.nb_points() == mf2.nb_dof(),
	      "Problem in identifying dofs");

}


void test_search_point() {
  const char *s = "BEGIN POINTS LIST\n"
    "  POINT  1  -4  6  2\n"
    "  POINT  2  0  6  0\n"
    "  POINT  3  0  2  0\n"
    "  POINT  4  -2  6  2\n"
    "  POINT  5  0  4  0\n"
    "  POINT  6  -1.5  4.5  0.5\n"
    "  POINT  7  1  2  0\n"
    "  POINT  8  1.5  1.5  0\n"
    "  POINT  9  5  5  0\n"
    "  POINT  10  2  1  0\n"
    "  POINT  11  6  3  0\n"
    "  POINT  12  2  0  0\n"
    "  POINT  13  6  0  0\n"
    "  POINT  14  2  4  0\n"
    "  POINT  15  4  2  0\n"
    "  POINT  17  3  6  0\n"
    "  POINT  18  2  -2  2\n"
    "  POINT  19  2  -2  -2\n"
    "  POINT  20  6  -2  2\n"
    "  POINT  21  6  -2  -2\n"
    "  POINT  22  2  -1  1\n"
    "  POINT  23  2  -2.5  0\n"
    "  POINT  24  2  -1  -1\n"
    "  POINT  25  6  -1  1\n"
    "  POINT  26  6  -2.5  0\n"
    "  POINT  27  6  -1  -1\n"
    "  POINT  28  -1  6  -1\n"
    "  POINT  29  -1  2  -1\n"
    "  POINT  30  1  6  -2\n"
    "  POINT  31  1  2  -2\n"
    "  POINT  32  0  6  -3\n"
    "  POINT  33  0  2  -3\n"
    "  POINT  34  2  -5  -2\n"
    "  POINT  35  2  -4  0\n"
    "  POINT  36  4  -5  2\n"
    "  POINT  37  6  -5  -2\n"
    "  POINT  38  6  -5  0\n"
    "  POINT  46  4  4  0\n"
    "  POINT  49  6  -5  2\n"
    "\n"
    "END POINTS LIST\n"
    "\n"
    "\n"
    "\n"
    "BEGIN MESH STRUCTURE DESCRIPTION\n"
    "\n"
    "CONVEX 0    'GT_PK(2,2)'      1  4  2  6  5  3\n"
    "CONVEX 1    'GT_QK(2,1)'      2  17  3  7\n"
    "CONVEX 2    'GT_QK(2,2)'      7  8  10  14  46  15  17  9  11\n"
    "CONVEX 3    'GT_QK(2,1)'      10  12  11  13\n"
    "CONVEX 4    'GT_PRODUCT(GT_PK(2,2),GT_PK(1,1))'      12  22  18  24  23  19  13  25  20  27  26  21\n"
    "CONVEX 5    'GT_PRODUCT(GT_PK(1,1),GT_PK(1,3))'      2  3  28  29  30  31  32  33\n"
    "CONVEX 8    'GT_PRODUCT(GT_PK(1,2),GT_PRODUCT(GT_PK(1,1),GT_PK(1,1)))'      19  23  18  21  26  20  34  35  36  37  38  49\n"
    "\n"
    "END MESH STRUCTURE DESCRIPTION\n";
  std::stringstream ss(s);
  getfem::mesh m; m.read_from_file(ss);
  cout << "read " << m.nb_points() << " points and " << m.convex_index().card() << " convexes\n";
  dal::bit_vector pid = m.points().index();
  cout << "point index: " << pid << "\n";
  for (dal::bv_visitor ii(pid); !ii.finished(); ++ii) {
    base_node P = m.points()[ii];
    size_type j = m.search_point(P);
    cerr << "search point " << ii << ": " << P << " -> " << j << "\n";
    assert(j == ii);
  }
  base_node P(m.dim()); gmm::fill_random(P);
  assert(m.search_point(P) == size_type(-1));
}



void test_incomplete_Q2(void) {
  // By Yao Koutsawa <yao.koutsawa@tudor.lu> 2012-12-10
  bgeot::pgeometric_trans pgt = bgeot::geometric_trans_descriptor("GT_Q2_INCOMPLETE(2)");
  
  const char *s =
    "BEGIN POINTS LIST\n"
    "  POINT  1  -4  6  2\n"
    "  POINT  2  0  6  0\n"
    "  POINT  3  0  2  0\n"
    "  POINT  4  -2  6  2\n"
    "  POINT  5  0  4  0\n"
    "  POINT  6  -1.5  4.5  0.5\n"
    "  POINT  7  1  2  0\n"
    "  POINT  8  1.5  1.5  0\n"
    "  POINT  9  5  5  0\n"
    "  POINT  10  2  1  0\n"
    "  POINT  11  6  3  0\n"
    "  POINT  12  2  0  0\n"
    "  POINT  13  6  0  0\n"
    "  POINT  14  2  4  0\n"
    "  POINT  15  4  2  0\n"
    "  POINT  17  3  6  0\n"
    "  POINT  18  2  -2  2\n"
    "  POINT  19  2  -2  -2\n"
    "  POINT  20  6  -2  2\n"
    "  POINT  21  6  -2  -2\n"
    "  POINT  22  2  -1  1\n"
    "  POINT  23  2  -2.5  0\n"
    "  POINT  24  2  -1  -1\n"
    "  POINT  25  6  -1  1\n"
    "  POINT  26  6  -2.5  0\n"
    "  POINT  27  6  -1  -1\n"
    "  POINT  28  -1  6  -1\n"
    "  POINT  29  -1  2  -1\n"
    "  POINT  30  1  6  -2\n"
    "  POINT  31  1  2  -2\n"
    "  POINT  32  0  6  -3\n"
    "  POINT  33  0  2  -3\n"
    "  POINT  34  2  -5  -2\n"
    "  POINT  35  2  -4  0\n"
    "  POINT  36  4  -5  2\n"
    "  POINT  37  6  -5  -2\n"
    "  POINT  38  6  -5  0\n"
    "  POINT  46  4  4  0\n"
    "  POINT  49  6  -5  2\n"
    "  POINT 100  7 0 0 \n"
    "  POINT 200  8 0 0 \n"
    "  POINT 300  8 1 0\n"
    "  POINT 400  7 1 0\n"
    "  POINT 500  7 0 1\n"
    "  POINT 600  8 0 1\n"
    "  POINT 700  8 1 1\n"
    "  POINT 800  7 1 1\n"
    "  POINT 900  7.5 0 0 \n"
    "  POINT 1000  8 0.5 0 \n"
    "  POINT 1100  7.5 1 0\n"
    "  POINT 1200  7 0.5 0\n"
    "  POINT 1300  7.5 0 1\n"
    "  POINT 1400  8 0.5 1\n"
    "  POINT 1500  7.5 1 1\n"
    "  POINT 1600  7 0.5 1\n"
    "  POINT 1700  7 0 0.5\n"
    "  POINT 1800  8 0 0.5\n"
    "  POINT 1900  8 1 0.5\n"
    "  POINT 2000  7 1 0.5\n"
    "\n"
    "END POINTS LIST\n"
    "\n"
    "\n"
    "\n"
    "BEGIN MESH STRUCTURE DESCRIPTION\n"
    "\n"
    "CONVEX 0    'GT_PK(2,2)'      1  4  2  6  5  3\n"
    "CONVEX 1    'GT_QK(2,1)'      2  17  3  7\n"
    "CONVEX 2    'GT_Q2_INCOMPLETE(2)'  7  8  10  14  15  17  9  11\n"
    "CONVEX 3    'GT_QK(2,1)'      10  12  11  13\n"
    "CONVEX 4    'GT_PRODUCT(GT_PK(2,2),GT_PK(1,1))'      12  22  18  24  23  19  13  25  20  27  26  21\n"
    "CONVEX 5    'GT_PRODUCT(GT_PK(1,1),GT_PK(1,3))'      2  3  28  29  30  31  32  33\n"
    "CONVEX 6    'GT_Q2_INCOMPLETE(3)'  100 900 200 1100 1300 500 1700 600 1000 1200 1800 1900 400 1400 300 1600 1500 800 2000 700\n"
    "\n"
    "END MESH STRUCTURE DESCRIPTION\n";
  
  std::stringstream ss(s);
  getfem::mesh m;
  m.read_from_file(ss);
  m.write_to_file("Q2_incomplete.msh");
  getfem::pos_export exp("Q2_incomplete.pos");
  exp.write(m,"mesh");
}












int main(void) {

  test_mesh_building(2, 100); 

  getfem::mesh m1;
  test_convex_ref();
  test_convex_simplif();
  test_mesh(m1);
  test_convex_simplif();
  test_convex_quality(0.2,0);
  test_convex_quality(-0.2,0);
  test_convex_quality(-0.01,-0.2);
  test_region();

  test_search_point();
  
  for (size_type d = 1; d <= 4 /* 6 */; ++d)
    test_mesh_matching(d);
  
  test_refinable(2, 1);
  test_refinable(2, 2);
  test_refinable(2, 3);
  test_refinable(3, 1);
  test_refinable(3, 2);
  test_refinable(3, 3);

  test_incomplete_Q2();
  
  return 0;
}




