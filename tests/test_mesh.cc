// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2002-2007 Yves Renard.
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
#include "getfem/getfem_regular_meshes.h"
#include "getfem/bgeot_poly_composite.h"
#include "getfem/bgeot_comma_init.h"
#include "getfem/getfem_export.h"


using getfem::size_type;
using getfem::base_node;
using getfem::base_small_vector;


void export_mesh(getfem::mesh &m, const std::string &name) {
  getfem::mesh_fem mf(m);
  mf.set_classical_finite_element(m.convex_index(), 0);
  std::vector<double> U(mf.nb_dof());
  for (size_type i = 0; i < mf.nb_dof(); ++i)
    U[i] = double(mf.first_convex_of_dof(i));

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
	DAL_THROW(gmm::failure_error, "Mesh is not conforming " << pt);
    }
  }
  cout << "nb faces of mesh of dim " << dim << " : " << nb_faces << endl;

}


void test_refinable(unsigned dim, unsigned degree) {

  int Nsubdiv = 2;
  std::vector<size_type> nsubdiv(dim, Nsubdiv);
  getfem::mesh m;
  getfem::regular_unit_mesh
    (m, nsubdiv, bgeot::simplex_geotrans(dim, degree), false);
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
  

  getfem::parallelepiped_regular_simplex_mesh(m, dim, pt1, vects.begin(),
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
  

  getfem::parallelepiped_regular_simplex_mesh(m, dim,pt1, vects.begin(), iref.begin());
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
  m.optimize_structure();
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
  bgeot::pconvex_ref t=sr3->basic_convex_ref();
  cout << "sr3->basic_convex_ref()=" << t << endl;

  cout << "sr3->basic_convex_ref() = faces:" << sr3->structure()->nb_faces() << " ; ";
  for (size_type i=0; i < sr3->basic_convex_ref()->points().size(); ++i) { if (i) cout << ","; cout << sr3->basic_convex_ref()->points()[i]; }
  cout << endl;

  const bgeot::mesh_structure *msr2 = sr2->basic_convex_ref()->simplexified_convex();
  print_mesh_structure(msr2);
  const bgeot::mesh_structure *msr3 = sr3->basic_convex_ref()->simplexified_convex();
  print_mesh_structure(msr3);
  const bgeot::mesh_structure *psr = pr->basic_convex_ref()->simplexified_convex();
  print_mesh_structure(psr);
  const bgeot::basic_mesh *msrr2 = bgeot::refined_simplex_mesh_for_convex(sr2,2);
  getfem::mesh(*msrr2).write_to_file(cout);
  const bgeot::basic_mesh *msrr3 = bgeot::refined_simplex_mesh_for_convex(sr3,3);
  getfem::mesh(*msrr3).write_to_file(cout);
  const bgeot::basic_mesh *mprr = bgeot::refined_simplex_mesh_for_convex(pr,2);
  getfem::mesh(*mprr).write_to_file(cout);
}

using bgeot::sc;

void test_convex_quality(getfem::scalar_type dx, getfem::scalar_type dy) {
  getfem::mesh m;
  getfem::base_node A,B,C;
  getfem::scalar_type h = sqrt(3.0);
  sc(B)+=0,0;
  sc(C)+=2,0;
  sc(A)+=0,h;
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
  for (unsigned k=1; k <= 2; ++k) {
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

int main(void) {
  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.
  
  try {
    getfem::mesh m1;
    test_convex_ref();
    test_convex_simplif();
    test_mesh(m1);
    test_convex_simplif();
    test_convex_quality(0.2,0);
    test_convex_quality(-0.2,0);
    test_convex_quality(-0.01,-0.2);
    test_region();

    for (size_type d = 1; d <= 4 /* 6 */; ++d)
      test_mesh_matching(d);

    test_refinable(2, 1);
    test_refinable(2, 2);
    test_refinable(2, 3);
    test_refinable(3, 1);
    test_refinable(3, 2);
    test_refinable(3, 3);
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0;
}




