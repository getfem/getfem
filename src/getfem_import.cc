// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_import.cc : misc. imports.
//           
// Date    : October 03, 2003.
// Author  : Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2000-2005 Julien Pommier
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

#include <iostream>
#include <iomanip>
#include <fstream>

#include <getfem_mesh.h>
#include <getfem_import.h>

namespace getfem {

  /* mesh file from gmsh [http://www.geuz.org/gmsh/] 
     structure: $NOD list_of_nodes $ENDNOD $ELT list_of_elt $ENDELT
  */
  static void import_gmsh_msh_file(std::ifstream& f, getfem_mesh& m) {
    /* read the node list */
    ftool::read_until(f, "$NOD");
    
    size_type nb_node; 
    f >> nb_node;
    //cerr << "reading nodes..[nb=" << nb_node << "]\n";
    dal::dynamic_tree_sorted<size_type> msh_node_2_getfem_node;
    for (size_type node_cnt=0; node_cnt < nb_node; ++node_cnt) {
      size_type node_id;
      base_node n(3); n[0]=0.0; n[1]=0.1; n[2]=1e30;
      f >> node_id >> n[0] >> n[1] >> n[2];
      msh_node_2_getfem_node.add_to_index(m.add_point(n), node_id);
    }
    ftool::read_until(f, "$ENDNOD");
    
    /* read the convexes */
    ftool::read_until(f, "$ELM");
    size_type nb_cv;
    f >> nb_cv;
    std::vector<size_type> cv_nodes;
    for (size_type cv=0; cv < nb_cv; ++cv) {
      size_type cv_id, cv_type, cv_region, cv_dummy, cv_nb_nodes;
      f >> cv_id >> cv_type >> cv_region >> cv_dummy >> cv_nb_nodes;
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
	m.add_segment(cv_nodes[0], cv_nodes[1]);
      } break;
      case 2: { /* TRIANGLE */
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

  /* mesh file from GiD [http://gid.cimne.upc.es/]

    supports linear and quadratic elements (quadrilaterals, use 9(or 27)-noded elements)
  */
  static void import_gid_msh_file(std::ifstream& f, getfem_mesh& m) {
    /* read the node list */
    size_type dim;
    enum { LIN,TRI,QUAD,TETR, PRISM, HEX,BADELTYPE } eltype=BADELTYPE;
    size_type nnode = 0; 
    dal::dynamic_tree_sorted<size_type> msh_node_2_getfem_node;
    std::vector<size_type> cv_nodes, getfem_cv_nodes;
    bool nodes_done = false;
    do {
      if (!f.eof()) f >> std::ws;
      if (f.eof() || !ftool::read_until(f, "MESH")) break;
      std::string selemtype;
      f >> ftool::skip("DIMENSION") >> dim 
	>> ftool::skip("ELEMTYPE") >> std::ws 
	>> selemtype 
	>> ftool::skip("NNODE") >> nnode;
      if (ftool::casecmp(selemtype, "linear")==0) { eltype = LIN;  }
      else if (ftool::casecmp(selemtype, "triangle")==0) { eltype = TRI; }
      else if (ftool::casecmp(selemtype, "quadrilateral")==0) { eltype = QUAD; }
      else if (ftool::casecmp(selemtype, "tetrahedra")==0) { eltype = TETR; }
      else if (ftool::casecmp(selemtype, "prisma")==0) { eltype = PRISM; }
      else if (ftool::casecmp(selemtype, "hexahedra")==0) { eltype = HEX; }
      else DAL_THROW(dal::failure_error, "unknown element type '"<< selemtype << "'");
      assert(!f.eof());
      f >> ftool::skip("COORDINATES");
      if (!nodes_done) {
	dal::dynamic_array<base_node> gid_nodes;
	dal::bit_vector gid_nodes_used;
	do {
	  //cerr << "reading coordinates " << std::streamoff(f.tellg()) << "\n";
	  std::string ls;
	  f >> std::ws;
	  std::getline(f,ls);
	  if (ftool::casecmp(ls, "END COORDINATES", 15)==0) break;
	  std::stringstream s; s << ls; 
	  size_type id;
	  s >> id;
	  
	  gid_nodes[id].resize(dim); gid_nodes_used.add(id);
	  for (size_type i=0; i < dim; ++i) s >> gid_nodes[id][i];
	  //cerr << "ppoint " << id << ", " << n << endl;
	} while (true);
	
	if (gid_nodes_used.card() == 0) {
	  DAL_THROW(dal::failure_error,"no nodes in the mesh!");
	}
	
	/* suppression of unused dimensions */
	std::vector<bool> direction_useless(3,true);
	base_node first_pt = gid_nodes[gid_nodes_used.first()];
	for (dal::bv_visitor ip(gid_nodes_used); !ip.finished(); ++ip) {
          for (size_type j=0; j < first_pt.size(); ++j) {
            if (direction_useless[j] && (gmm::abs(gid_nodes[ip][j]-first_pt[j]) > 1e-13))
              direction_useless[j] = false;
          }
	}
	size_type dim2=0;
	for (size_type j=0; j < dim; ++j) if (!direction_useless[j]) dim2++;
	for (dal::bv_visitor ip(gid_nodes_used); !ip.finished(); ++ip) {
          base_node n(dim2);
          for (size_type j=0, cnt=0; j < dim; ++j) if (!direction_useless[j]) n[cnt++]=gid_nodes[ip][j];
          msh_node_2_getfem_node.add_to_index(m.add_point(n), ip);
        }
      }
      
      ftool::read_until(f, "ELEMENTS");
      bgeot::pgeometric_trans pgt = NULL;
      std::vector<size_type> order(nnode); // ordre de GiD cf http://gid.cimne.upc.es/support/gid_11.subst#SEC160
      for (size_type i=0; i < nnode; ++i) order[i]=i;
      //cerr << "reading elements " << std::streamoff(f.tellg()) << "\n";
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
	else if (nnode == 6) { // validé
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
	else if (nnode == 10) { // validé
	  static size_type lorder[10] = {0,4,1, 7,8, 3, 6, 5, 9, 2};
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
	std::stringstream s; s << ls; 
	size_type cv_id;
	s >> cv_id;
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
	//cerr << "elt " << cv_id << ", " << getfem_cv_nodes << endl;
	
	// envisager la "simplification" quand on a une transfo non 
	// lineaire mais que la destination est lineaire 
	m.add_convex(pgt, getfem_cv_nodes.begin());
      } while (true);
    } while (!f.eof());
  }

  static double round_to_nth_significant_number(double x, int ndec) {
    double p = 1.;
    double s = (x < 0 ? -1 : 1);
    double pdec = pow(10.,ndec);
    if (x == 0) return 0.;
    x = gmm::abs(x);
    while (x > 1) { x /= 10.0; p*=10; }
    while (x < 0.1) { x *= 10.0; p/=10; }
    //cerr << "x=" << x << ", p=" << p << ", pdec=" << pdec << "\n";
    x = s * (floor(x * pdec + 0.5) / pdec) * p;
    return x;
  }

  /* mesh file from emc2 [http://pauillac.inria.fr/cdrom/prog/unix/emc2/eng.htm], am_fmt format

    (only triangular 2D meshes)
  */
  static void import_am_fmt_file(std::ifstream& f, getfem_mesh& m) {
    /* read the node list */
    std::vector<size_type> tri;      
    size_type nbs,nbt;
    base_node P(2);
    
    f >> nbs >> nbt; ftool::read_until(f,"\n");
    tri.resize(nbt*3);
    for (size_type i=0; i < nbt*3; ++i) f >> tri[i];
    for (size_type j=0; j < nbs; ++j) {
      f >> P[0] >> P[1]; 
      cerr.precision(16);
      P[0]=round_to_nth_significant_number(P[0],6); // force 9.999999E-1 to be converted to 1.0
      P[1]=round_to_nth_significant_number(P[1],6);
      if (m.add_point(P) != j) DAL_INTERNAL_ERROR("ouch");
    }
    for (size_type i=0; i < nbt*3; i+=3)
      m.add_triangle(tri[i]-1,tri[i+1]-1,tri[i+2]-1);
  }

  /* mesh file from emc2 [http://pauillac.inria.fr/cdrom/prog/unix/emc2/eng.htm], am_fmt format

    triangular/quadrangular 2D meshes
  */
  static void import_emc2_mesh_file(std::ifstream& f, getfem_mesh& m) {
    /* read the node list */
    std::vector<size_type> tri;      
    size_type nbs=0,nbt=0,nbq=0,dummy;
    base_node P(2);
    ftool::read_until(f,"Vertices");
    f >> nbs;
    for (size_type j=0; j < nbs; ++j) {
      f >> P[0] >> P[1] >> dummy; 
      if (m.add_point(P) != j) DAL_INTERNAL_ERROR("ouch");
    }
    while (!f.eof()) {
      size_type ip[4];
      std::string ls;
      std::getline(f,ls);
      if (ls.find("Triangles")+1) {
        f >> nbt;
        for (size_type i=0; i < nbt; ++i) {
          f >> ip[0] >> ip[1] >> ip[2] >> dummy; ip[0]--; ip[1]--; ip[2]--;
          m.add_triangle(ip[0],ip[1],ip[2]);
        }
      } else if (ls.find("Quadrangles")+1) {
        f >> nbq;
        for (size_type i=0; i < nbq; ++i) {
          f >> ip[0] >> ip[1] >> ip[2] >> ip[3] >> dummy; ip[0]--; ip[1]--; ip[2]--; ip[3]--;
          m.add_parallelepiped(2, &ip[0]);
        }
      } else if (ls.find("End")+1) break;
    }
  }

  void import_mesh(const std::string& filename, const std::string& format, getfem_mesh& m) {
    m.clear();
    try {
      std::ifstream f(filename.c_str());
      if (!f.good()) DAL_THROW(failure_error, "can't open file " << filename);
      /* throw exceptions when an error occurs */
      f.exceptions(std::ifstream::badbit | std::ifstream::failbit);
      import_mesh(f,format,m);
      f.close();
    }
    catch (dal::failure_error& exc) {
      m.clear();
      throw exc;
    }
    catch (std::ios_base::failure& exc) {
      m.clear();
      DAL_THROW(dal::failure_error, "error while importing " << format << " mesh file \"" << filename << "\" : " << exc.what());
    }
  }

  void import_mesh(std::ifstream& f, const std::string& format,
		   getfem_mesh& m) {
    if (ftool::casecmp(format,"gmsh")==0)
      import_gmsh_msh_file(f,m);
    else if (ftool::casecmp(format,"gid")==0)
      import_gid_msh_file(f,m);
    else if (ftool::casecmp(format,"am_fmt")==0)
      import_am_fmt_file(f,m);
    else if (ftool::casecmp(format,"emc2_mesh")==0)
      import_emc2_mesh_file(f,m);
    else DAL_THROW(dal::failure_error, "cannot import "
		   << format << " mesh type : unknown mesh type");
  }

  void import_mesh(const std::string& filename, getfem_mesh& mesh) {
    if (filename.compare(0,4,"gid:")==0)
      getfem::import_mesh(filename.substr(4), "gid", mesh);
    else if (filename.compare(0,5,"gmsh:") == 0) 
      getfem::import_mesh(filename.substr(5), "gmsh", mesh);
    else if (filename.compare(0,7,"am_fmt:") == 0) 
      getfem::import_mesh(filename.substr(7), "am_fmt", mesh);
    else if (filename.compare(0,10,"emc2_mesh:") == 0)
      getfem::import_mesh(filename.substr(10), "emc2_mesh", mesh);
    else mesh.read_from_file(filename);
  }

}

