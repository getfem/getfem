#include <bgeot_comma_init.h>
#include <getfem_mesh_slice.h>

using getfem::size_type;
namespace getfem {
  std::ostream& operator<<(std::ostream& o, const bgeot::mesh_structure& m) {
    o << "mesh_structure: nb_pts=" << m.point_structures().size() << ", nb_cvs=" << m.convex_index().card() << endl;
    dal::bit_vector bv = m.convex_index();
    for (size_type cv = bv.take_first(); cv != size_type(-1); cv << bv) {
      o << "convex " << cv << ": "; 
      std::copy(m.ind_points_of_convex(cv).begin(), m.ind_points_of_convex(cv).end(), std::ostream_iterator<size_type>(o, ","));
      o << endl;
    }
    return o;
  }

#if 0  
  std::ostream& operator<<(std::ostream& o, const mesh_slice& m) {
    std::vector<size_type> e;
    o << "mesh_slice, containing " << m.nb_convex() << " convexes\n";
    for (size_type ic = 0; ic < m.nb_convex(); ++ic) {
      o << "slice convex #" << ic << " (original = " << m.convex_num(ic) << ")\n";
      for (size_type i = 0; i < m.nodes(ic).size(); ++i) {
        o << "node " << i << ": " << m.nodes(ic)[i].pt << ", ref=" << m.nodes(ic)[i].pt_ref << " flist=" << m.nodes(ic)[i].faces << endl;
      }
      for (size_type i = 0; i < m.simplexes(ic).size(); ++i) {
        o << "simplex " << i << ", inodes=";
        for (size_type j=0;j< m.simplexes(ic)[i].dim()+1;++j)
          o << m.simplexes(ic)[i].inodes[j] << " ";
        o << endl;
      }
      /*m.edges(ic,e);
      o << "edges: "; for (size_type i=0; i < e.size()/2; ++i) o << e[2*i] << "-" << e[2*i+1] << " ";
      o << endl;
      */
    }
    return o;
  }
#endif
}

int 
main() {
  getfem::getfem_mesh m;
  getfem::base_node A; bgeot::sc(A)=0,0;
  getfem::base_node B; bgeot::sc(B)=1,0;
  getfem::base_node C; bgeot::sc(C)=0,2;
  getfem::base_node D; bgeot::sc(D)=1,1;
  m.add_triangle_by_points(A,B,C);
  m.add_triangle_by_points(B,C,D);

  getfem::base_node x0; bgeot::sc(x0) = .4,0;
  getfem::base_node n0;  bgeot::sc(n0) = 1,0;  
  getfem::base_node x1; bgeot::sc(x1) = 0,0.1;
  getfem::base_node n1;  bgeot::sc(n1) = 1,-1;  

  getfem::convex_face_ct cvlst; 
  cvlst.push_back(getfem::convex_face(0,size_type(-1)));
  cvlst.push_back(getfem::convex_face(1,size_type(-1)));
  getfem::slicer_half_space sl0(x0,n0,false);
  getfem::mesh_slice msl0(m);
  msl0.build(&sl0,40,cvlst);
  //cout << msl0 << endl;

  getfem::slicer_half_space sl1(x1,n1,false);
  getfem::slicer_intersect(&sl0,&sl1);
  getfem::mesh_slice msl1(m);
  //cout << "memory 1: " << msl1.memsize() << " bytes\n";
  msl1.build(&sl1,40,cvlst);
  /*cout << msl1 << endl;
  cout << "memory 0: " << msl0.memsize() << " bytes\n";
  cout << "memory 1: " << msl1.memsize() << " bytes\n";*/
  return 0;
}
