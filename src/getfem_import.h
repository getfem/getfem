/* -*- c++ -*- (enables emacs c++ mode)                                    */
#ifndef __GETFEM_IMPORT_H
#define __GETFEM_IMPORT_H
#include <iostream>
#include <iomanip>
#include <fstream>

#include <getfem_mesh_fem.h>

namespace getfem {
  inline void import_gmsh_msh_file(const std::string& filename, getfem_mesh& m) {
    m.clear();
    cerr << "m.dim()=" << int(m.dim()) << endl;
    try {
      std::ifstream f(filename.c_str());
      if (!f.good()) DAL_THROW(failure_error, "can't open file " << filename);
      /* throw exceptions when an error occurs */
      f.exceptions(std::ifstream::badbit | std::ifstream::failbit);
      
      /* read the node list */
      ftool::read_untill(f, "$NOD");
      //std::string s;
      //std::getline(f,s);
      std::cerr << "4.3=" << 4.3 << ", 2.5=" << 2.5 << std::endl;
      size_type nb_node; 
      f >> nb_node;
      cerr << "reading nodes..[nb=" << nb_node << "]\n";
      dal::dynamic_tree_sorted<size_type> gmsh_node_2_getfem_node;
      for (size_type node_cnt=0; node_cnt < nb_node; ++node_cnt) {
	size_type node_id;
        base_node n(3); n[0]=0.0; n[1]=0.1; n[2]=1e30;
	cerr << "node " << node_cnt << std::flush;
	f >> node_id >> n[0] >> n[1] >> n[2];
	cerr << " id=" << node_id << " n=" << n[0] << "," << n[1] << "," << n[2] << " m.dim() =" << int(m.dim()) << endl;
	gmsh_node_2_getfem_node.add_to_index(m.add_point(n), node_id);
      }
      ftool::read_untill(f, "$ENDNOD");
      
      /* read the convexes */
      cerr << "reading convexes..\n";
      ftool::read_untill(f, "$ELM");
      size_type nb_cv;
      f >> nb_cv;
      //std::vector<size_type> gmsh_cv_2_getfem_cv(nb_cv, size_type(-1));
      std::vector<size_type> cv_nodes;
      for (size_type cv=0; cv < nb_cv; ++cv) {
	size_type cv_id, cv_type, cv_region, cv_dummy, cv_nb_nodes;
	f >> cv_id >> cv_type >> cv_region >> cv_dummy >> cv_nb_nodes;
	cerr << "convex " << cv << " id=" << cv_id << ", type=" << cv_type << ", nb_nodes=" << cv_nb_nodes << endl;
	cv_id--; /* numbering starts at 1 */
	cv_nodes.resize(cv_nb_nodes);
	for (size_type i=0; i < cv_nb_nodes; ++i) {
	  size_type j;
	  f >> j;
	  cv_nodes[i] = gmsh_node_2_getfem_node.search(j);
	  if (cv_nodes[i] == size_type(-1)) 
	    DAL_THROW(failure_error, "Invalid node ID " << j 
		      << " in gmsh convex " << cv_id);
	}
	switch (cv_type) {
	case 1: { /* LINE */
	  cerr << "add line: cv_nodes = " << cv_nodes << endl;
	  cerr << "pt0 = " << m.points()[0] << ", pt1=" << m.points()[1] << endl;
	  m.add_segment(cv_nodes[0], cv_nodes[1]);
	} break;
	case 2: { /* TRIANGLE */
	  cerr << "and triangle: cv_nodes = " << cv_nodes << endl;
	  cerr << "pt0 = " << m.points()[0] << ", pt1=" << m.points()[1] << ", pt2=" << m.points()[2] << endl;
	  m.add_triangle(cv_nodes[0], cv_nodes[1], cv_nodes[2]);
	} break;
	case 3: { /* QUADRANGLE */
	  m.add_parallelepiped(2,cv_nodes.begin());
	} break;
	case 4: { /* TETRAHEDRON */
	  m.add_simplex(3,cv_nodes.begin());
	} break;
	case 5: { /* HEXAHEDRON */
	  m.add_parallelepiped(3,cv_nodes.begin());
	} break;
	case 6: { /* PRISM */
	  m.add_prism(3,cv_nodes.begin());
	} break;
	case 7: { /* PYRAMID */
	  DAL_THROW(failure_error, "sorry pyramidal convexes not done for the moment..");
	} break;
	case 15: { /* POINT */
	  DAL_WARNING(2, "ignoring point element");
	} break;
	default: { /* UNKNOWN .. */
	  DAL_THROW(failure_error, "the gmsh element type " << cv_type << "is unknown..");
	} break;
	}	
      }
    }
    catch (std::ios_base::failure& exc) {
      m.clear();
      DAL_THROW(dal::failure_error, "error while reading gmsh file \"" << filename << "\" : " << exc.what());
    }
  }
}
#endif /* __GETFEM_IMPORT_H  */
