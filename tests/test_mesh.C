#include <dal_std.h>
#include <stdio.h>
#include <getfem_regular_meshes.h>


typedef bgeot::fsvector<double, 3> fsvect3;
typedef bgeot::PT<fsvect3> fspoint3;
typedef bgeot::fsvector<double, 4> fsvect4;
typedef bgeot::PT<fsvect4> fspoint4;
typedef bgeot::vsvector<double> vsvect;
typedef bgeot::PT<vsvect> vspoint;


template<class MESH> void test_mesh(MESH &m)
{
  
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

  size_t ic1 = m.add_segment(i2, i3);

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

  size_t ic2 = m.add_convex_by_points(bgeot::simplex_trans(2,1), pts.begin());

  size_t ic3 = m.add_convex_by_points(bgeot::simplex_trans(2,1), pts.begin());

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
  

  parallelepiped_regular_simplex_mesh(m, dim,pt1, vects.begin(), iref.begin());
	    

  m.write_to_file("test_mesh.mesh");
  m.write_to_file(cout);

}



int main(void)
{
  try {
  getfem::getfem_mesh m1;

  test_mesh(m1);
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0;
}




