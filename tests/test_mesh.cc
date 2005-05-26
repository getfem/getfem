/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
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
#include <dal_std.h>
#include <stdio.h>
#include <getfem_regular_meshes.h>
#include <getfem_poly_composite.h>
#include <bgeot_comma_init.h>
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif

using getfem::size_type;

template<class MESH> void test_mesh(MESH &m) {
  
  typedef typename MESH::point_type POINT;
  typedef typename MESH::vector_type VECT;
  
  POINT pt1, pt2, pt3;
  if (pt1.size() == 0) pt1 = POINT(3);

  pt1.fill(0.0); pt2 = pt1; pt3 = pt1;
  
  pt2[0] = 1.0E-13;
  pt3[0] = 1.0;
  size_t i1 = m.add_point(pt1);
  size_t i2 = m.add_point(pt2);
  size_t i3 = m.add_point(pt3);
  
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
  pt1.fill(0.0);
  for (int i = 0; i < dim; i++)
  { vects[i] = VECT(dim); vects[i].fill(0.0); vects[i][i] = 1.0; iref[i] = 3; }
  

  getfem::parallelepiped_regular_simplex_mesh(m, dim,pt1, vects.begin(), iref.begin());
	    

  m.add_convex_to_set(3,0);

  m.region(3).add(0);
  m.region(3).add(3);
  m.region(3).add(2);

  m.region(4).add(0);
  m.sup_convex(0);
  m.sup_convex(1);
  m.optimize_structure();
  m.write_to_file("test_mesh.mesh");
  m.write_to_file(cout);
}

void
print_mesh_structure(const bgeot::mesh_structure *ms) {
  cout << "nb_pts=" <<ms->point_structures().size() << ", nb_cvs=" << ms->nb_convex() << endl;
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
  const getfem::getfem_mesh *msrr2 = getfem::refined_simplex_mesh_for_convex(sr2,2);
  msrr2->write_to_file(cout);
  const getfem::getfem_mesh *msrr3 = getfem::refined_simplex_mesh_for_convex(sr3,3);
  msrr3->write_to_file(cout);
  const getfem::getfem_mesh *mprr = getfem::refined_simplex_mesh_for_convex(pr,2);
  mprr->write_to_file(cout);
}

using bgeot::sc;

void test_convex_quality(getfem::scalar_type dx, getfem::scalar_type dy) {
  getfem::getfem_mesh m;
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

class myexc : public dal::exception_callback {
  void callback(const std::string& s)
  { cerr << "exception launched: " << s << std::endl; }
};

int main(void)
{
#ifdef GETFEM_HAVE_FEENABLEEXCEPT /* trap SIGFPE */
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
  dal::set_exception_callback(new myexc);
  try {
    cout << "sizeof(size_type)=" << sizeof(size_type) << endl;
  getfem::getfem_mesh m1;

  test_convex_simplif();
  test_mesh(m1);
  test_convex_simplif();
  test_convex_quality(0.2,0);
  test_convex_quality(-0.2,0);
  test_convex_quality(-0.01,-0.2);
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0;
}




