/* -*- c++ -*- (enables emacs c++ mode)                                    */
#ifndef __GETFEM_IMPORT_H
#define __GETFEM_IMPORT_H
#include <iostream>
#include <iomanip>
#include <fstream>

#include <getfem_mesh_fem.h>

namespace getfem {
  /* mesh file from gmsh [http://www.geuz.org/gmsh/] 
     structure: $NOD list_of_nodes $ENDNOD $ELT list_of_elt $ENDELT
  */
  inline void import_gmsh_msh_file(const std::string& filename, getfem_mesh& m) {
    m.clear();
    try {
      std::ifstream f(filename.c_str());
      if (!f.good()) DAL_THROW(failure_error, "can't open file " << filename);
      /* throw exceptions when an error occurs */
      f.exceptions(std::ifstream::badbit | std::ifstream::failbit);
      
      /* read the node list */
      ftool::read_untill(f, "$NOD");

      size_type nb_node; 
      f >> nb_node;
      cerr << "reading nodes..[nb=" << nb_node << "]\n";
      dal::dynamic_tree_sorted<size_type> msh_node_2_getfem_node;
      for (size_type node_cnt=0; node_cnt < nb_node; ++node_cnt) {
	size_type node_id;
        base_node n(3); n[0]=0.0; n[1]=0.1; n[2]=1e30;
	cerr << "node " << node_cnt << std::flush;
	f >> node_id >> n[0] >> n[1] >> n[2];
	cerr << " id=" << node_id << " n=" << n[0] << "," << n[1] << "," << n[2] << " m.dim() =" << int(m.dim()) << endl;
	msh_node_2_getfem_node.add_to_index(m.add_point(n), node_id);
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
	  cv_nodes[i] = msh_node_2_getfem_node.search(j);
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
      f.close();
    }
    catch (std::ios_base::failure& exc) {
      m.clear();
      DAL_THROW(dal::failure_error, "error while reading gmsh file \"" << filename << "\" : " << exc.what());
    }
  }

  /* mesh file from GiD [http://gid.cimne.upc.es/]

    supports linear and quadratic elements (quadrilaterals, use 9(or 27)-noded elements)
  */
  inline void import_gid_msh_file(const std::string& filename, getfem_mesh& m) {
    m.clear();
    try {
      std::ifstream f(filename.c_str());
      if (!f.good()) DAL_THROW(failure_error, "can't open file " << filename);
      /* throw exceptions when an error occurs */
      f.exceptions(std::ifstream::badbit | std::ifstream::failbit);
      
      /* read the node list */
      size_type dim;
      enum { LIN,TRI,QUAD,TETR, PRISM, HEX,BADELTYPE } eltype=BADELTYPE;
      size_type nnode = 0; 
      dal::dynamic_tree_sorted<size_type> msh_node_2_getfem_node;
      std::vector<size_type> cv_nodes, getfem_cv_nodes;
      do {
	cerr << "hello..." << int(f.tellg()) << "\n";
	if (!f.eof()) f >> std::ws;
	if (f.eof() || !ftool::read_untill(f, "MESH")) break;
	cerr << "hello2..." << int(f.tellg()) << "\n";
	std::string selemtype;
	f >> ftool::skip("DIMENSION") >> dim 
	  >> ftool::skip("ELEMTYPE") >> std::ws 
	  >> selemtype 
	  >> ftool::skip("NNODE") >> nnode;
	cerr << "elem type = " << "'" << selemtype << "', nnode=" << nnode << endl;

// 	ftool::read_untill(f, "DIMENSION");
// 	cerr << "read dim " << int(f.tellg()) << "\n";
// 	f >> dim; if (dim < 1 || dim > 3) DAL_THROW(dal::failure_error, "wrong mesh dimension");
// 	cerr << "dim = " << dim << endl;
// 	ftool::read_untill(f, "ELEMTYPE");
// 	char selemtype[512]; ftool::get_token(f, selemtype, 512);
// cerr << "before  coordinates " << int(f.tellg()) << "\n";
// 	f >> nnode;
// 	cerr << "elem type = " << "'" << selemtype << "', nnode" << nnode << endl;
	if (ftool::casecmp(selemtype, "linear")==0) { eltype = LIN;  }
	else if (ftool::casecmp(selemtype, "triangle")==0) { eltype = TRI; }
	else if (ftool::casecmp(selemtype, "quadrilateral")==0) { eltype = QUAD; }
	else if (ftool::casecmp(selemtype, "tetrahedra")==0) { eltype = TETR; }
	else if (ftool::casecmp(selemtype, "prisma")==0) { eltype = PRISM; }
	else if (ftool::casecmp(selemtype, "hexahedra")==0) { eltype = HEX; }
	else DAL_THROW(dal::failure_error, "unknown element type '"<< selemtype << "'");
	cerr << "before  coordinates " << int(f.tellg()) << "\n";
	assert(!f.eof());
	//ftool::read_untill(f,"COORDINATES");
	f >> ftool::skip("COORDINATES");
	do {
	  cerr << "reading coordinates " << int(f.tellg()) << "\n";
	  std::string ls;
	  f >> std::ws;
	  std::getline(f,ls);
	  if (ftool::casecmp(ls, "END COORDINATES", 15)==0) break;
	  std::stringstream s(ls); 
	  base_node n(dim);
	  size_type id;
	  s >> id;
	  for (size_type i=0; i < dim; ++i) s >> n[i];
	  cerr << "ppoint " << id << ", " << n << endl;
	  msh_node_2_getfem_node.add_to_index(m.add_point(n), id);	  
	} while (true);
	
	ftool::read_untill(f, "ELEMENTS");
	bgeot::pgeometric_trans pgt = NULL;
	std::vector<size_type> order(nnode); // ordre de GiD cf http://gid.cimne.upc.es/support/gid_11.subst#SEC160
	for (size_type i=0; i < nnode; ++i) order[i]=i;
	cerr << "reading elements " << int(f.tellg()) << "\n";
	switch (eltype) {
	case LIN: {
	  if (nnode == 2) pgt = bgeot::simplex_geotrans(1,1);
	  else if (nnode == 3) {
	    pgt = bgeot::simplex_geotrans(1,2); // A VERIFIER TOUT CA
	    std::swap(order[1],order[2]);
	  }
	} break;
	case TRI: {
	  if (nnode == 3) pgt = bgeot::simplex_geotrans(2,1);
	  else if (nnode == 6) {
	    static size_type lorder[6] = {0,3,1,5,4,2};
	    pgt = bgeot::simplex_geotrans(2,2);
	    std::copy(lorder,lorder+nnode,order.begin());
	  }
	} break;
	case QUAD: {
	  if (nnode == 4) pgt = bgeot::parallelepiped_geotrans(2,1);
	  else if (nnode == 9) {
	    static size_type lorder[9] = {0,4,1, 7,8,5, 3,6,2};
	    pgt = bgeot::parallelepiped_geotrans(2,2);
	    std::copy(lorder,lorder+nnode,order.begin());
	  }
	} break;
	case TETR: {
	  if (nnode == 4) pgt = bgeot::simplex_geotrans(3,1);
	  else if (nnode == 10) {
	    static size_type lorder[10] = {0,4,1, 7,8, 3, 5,6, 9, 2};
	    pgt = bgeot::simplex_geotrans(3,2);
	    std::copy(lorder,lorder+nnode,order.begin());
	  }
	} break;
	case PRISM: {
	  if (nnode == 6) pgt = bgeot::prism_geotrans(3,1);
	} break;
	case HEX: {
	  if (nnode == 8) pgt = bgeot::parallelepiped_geotrans(3,1);
	  else if (nnode == 27) {
	    static size_type lorder[27] = {0,8,1, 12,21,13, 4,16,5,
					   11,20,9, 22,26,24, 19,25,17,
					   3,10,2, 15,23,14, 7,18,6};
	    pgt = bgeot::parallelepiped_geotrans(3,2);
	    std::copy(lorder,lorder+nnode,order.begin());
	  }
	} break;
	default: {
	  DAL_INTERNAL_ERROR("");
	} break;
	}
	if (pgt == NULL) DAL_THROW(dal::failure_error, "unknown element type " << selemtype 
				   << " with " << nnode << "nodes");
	do {
	  std::string ls;
	  f >> std::ws;
	  std::getline(f,ls);
	  if (ftool::casecmp(ls, "END ELEMENTS", 12)==0) break;
	  cerr << "line: '" << ls << "'" << endl;
	  std::stringstream s(ls); 
	  size_type cv_id;
	  s >> cv_id;
	  cerr << "cv_id=" << cv_id << endl;
	  cv_nodes.resize(nnode);
	  for (size_type i=0; i < nnode; ++i) {
	    size_type j;
	    s >> j;
	    cv_nodes[i] = msh_node_2_getfem_node.search(j);
	    if (cv_nodes[i] == size_type(-1)) 
	      DAL_THROW(failure_error, "Invalid node ID " << j 
			<< " in GiD mesh convex num " << cv_id);
	  }
	  getfem_cv_nodes.resize(nnode);
	  for (size_type i=0; i < nnode; ++i) {
	    getfem_cv_nodes[i] = cv_nodes[order[i]];
	  }
	  cerr << "elt " << cv_id << ", " << getfem_cv_nodes << endl;

	  // envisager la "simplification" quand on a une transfo non 
	  // lineaire mais que la destination est lineaire 
	  m.add_convex(pgt, getfem_cv_nodes.begin());
	} while (true);
      } while (!f.eof());
      f.close();
    }
    catch (dal::failure_error& exc) {
      m.clear();
      throw exc;
    }
    catch (std::ios_base::failure& exc) {
      m.clear();
      DAL_THROW(dal::failure_error, "error while reading GiD mesh file \"" << filename << "\" : " << exc.what());
    }
  }

}
#endif /* __GETFEM_IMPORT_H  */
