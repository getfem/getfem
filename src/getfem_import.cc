/*===========================================================================
 
 Copyright (C) 2000-2012 Julien Pommier
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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

#include <iostream>
#include <iomanip>
#include <fstream>

#include "getfem/getfem_mesh.h"
#include "getfem/getfem_import.h"
#include "getfem/getfem_regular_meshes.h"

namespace getfem {

  /* mesh file from gmsh [http://www.geuz.org/gmsh/]*/

  struct gmsh_cv_info {
    unsigned id, type, region;
    bgeot::pgeometric_trans pgt;
    std::vector<size_type> nodes;
    void set_pgt() {
      switch (type) {
      case 1: { /* LINE */
        pgt = bgeot::simplex_geotrans(1,1);
      } break;
      case 2: { /* TRIANGLE */
        pgt = bgeot::simplex_geotrans(2,1);
      } break;
      case 3: { /* QUADRANGLE */
        pgt = bgeot::parallelepiped_geotrans(2,1);
      } break;
      case 4: { /* TETRAHEDRON */
        pgt = bgeot::simplex_geotrans(3,1);
      } break;
      case 5: { /* HEXAHEDRON */
        pgt = bgeot::parallelepiped_geotrans(3,1);
      } break;
      case 6: { /* PRISM */
        pgt = bgeot::prism_geotrans(3,1);
      } break;
      case 7: { /* PYRAMID */
        GMM_ASSERT1(false, "sorry pyramidal elements not yet supported.");
      } break;
      case 8: { /* 2ND ORDER LINE */
        pgt = bgeot::simplex_geotrans(1,2);
      } break;
      case 9: { /* 2ND ORDER TRIANGLE */
        pgt = bgeot::simplex_geotrans(2,2);
      } break;
	  case 11: { /* 2ND ORDER TETRAEDRON */
		pgt = bgeot::simplex_geotrans(3,2);
      } break;
      case 15: { /* POINT */
        GMM_WARNING2("ignoring point element");
      } break;
      default: { /* UNKNOWN .. */
        /* higher order elements : to be done .. */
        GMM_ASSERT1(false, "the gmsh element type "<< type <<"is unknown..");
      } break;
      }
    }

    void set_nb_nodes() {
      /* Especially for the gmsh file format version 2.*/
      switch (type) {
      case 1: { /* LINE */
        nodes.resize(2);
      } break;
      case 2: { /* TRIANGLE */
        nodes.resize(3);
      } break;
      case 3: { /* QUADRANGLE */
        nodes.resize(4);
      } break;
      case 4: { /* TETRAHEDRON */
        nodes.resize(4);
      } break;
      case 5: { /* HEXAHEDRON */
        nodes.resize(8);
      } break;
      case 6: { /* PRISM */
        nodes.resize(6);
      } break;
      case 7: { /* PYRAMID */
        GMM_ASSERT1(false,
                    "sorry pyramidal convexes not done for the moment..");
      } break;
      case 8: { /* 2ND ORDER LINE */
        nodes.resize(3);
      } break;
      case 9: { /* 2ND ORDER TRIANGLE */
        nodes.resize(6);
      } break;
	  case 11: { /*2ND ORDER TETRAHEDRON */
		nodes.resize(10);
      } break;
      case 15: { /* POINT */
        GMM_WARNING2("ignoring point element");
      } break;
      default: { /* UNKNOWN .. */
        /* higher order elements : to be done .. */
        GMM_ASSERT1(false, "the gmsh element type " << type << "is unknown..");
      } break;
      }
    }

    bool operator<(const gmsh_cv_info& other) const {
      return pgt->dim() > other.pgt->dim();
    }
  };

  /*
     Format version 1 [for gmsh version < 2.0].
     structure: $NOD list_of_nodes $ENDNOD $ELT list_of_elt $ENDELT

     Format version 2 [for gmsh version >= 2.0]. Modification of some keywords
     and more than one tag is authorized for a specific region.

     structure: $Nodes list_of_nodes $EndNodes $Elements list_of_elt
     $EndElements
  */
  static void import_gmsh_msh_file(std::ifstream& f, mesh& m, int deprecate=0)
  {
    gmm::stream_standard_locale sl(f);
    /* print general warning */
    GMM_WARNING3("  All regions must have different number!");

    /* print deprecate warning */
    if (deprecate!=0){
      GMM_WARNING4("" << endl
                << "  deprecate: " << endl
                << "   static void" << endl
                << "   import_gmsh_msh_file(std::ifstream& f,"
                << " mesh& , int version)" << endl
                << "  replace with:" << endl
                << "   static void" << endl
                << "   import_gmsh_msh_file(std::ifstream& f,"
                << " mesh&)");
    }

    /* read the version */
    int version;
    std::string header;
    f >> header;
    if (bgeot::casecmp(header,"$MeshFormat")==0)
      f >> version;
    else if (bgeot::casecmp(header,"$NOD")==0)
      version = 1;
    else
      GMM_ASSERT1(false, "can't read Gmsh format: " << header);

    /* read the node list */
    if (version == 2)
      bgeot::read_until(f, "$Nodes"); /* Format version 2 */

    size_type nb_node;
    f >> nb_node;
    //cerr << "reading nodes..[nb=" << nb_node << "]\n";
    std::map<size_type, size_type> msh_node_2_getfem_node;
    for (size_type node_cnt=0; node_cnt < nb_node; ++node_cnt) {
      size_type node_id;
      base_node n(3); n[0]=n[1]=n[2]=0.0;
      f >> node_id >> n[0] >> n[1] >> n[2];
      msh_node_2_getfem_node[node_id] = m.add_point(n);
    }

    if (version == 2)
      bgeot::read_until(f, "$Endnodes"); /* Format version 2 */
    else
      bgeot::read_until(f, "$ENDNOD");

    /* read the convexes */
    if (version == 2)
      bgeot::read_until(f, "$Elements"); /* Format version 2 */
    else
      bgeot::read_until(f, "$ELM");

    size_type nb_cv;
    f >> nb_cv;
    std::vector<gmsh_cv_info> cvlst; cvlst.reserve(nb_cv);
    for (size_type cv=0; cv < nb_cv; ++cv) {
      unsigned id, type, region;
      unsigned dummy, cv_nb_nodes;

      if (version == 2) { /* Format version 2 */
        unsigned nbtags, mesh_part;
        f >> id >> type >> nbtags;
        if (nbtags == 0 || nbtags > 3)
          GMM_ASSERT1(false, "Number of tags " << nbtags
                      << " is not managed.");

        f >> region;
        if (nbtags > 1) f >> dummy;
        if (nbtags > 2) f >> mesh_part;
      }
      else
        f >> id >> type >> region >> dummy >> cv_nb_nodes;

      id--; /* gmsh numbering starts at 1 */
      if (type == 15) { // just a node
        // get rid of the rest of the line
        f >> dummy;
      } else {
        cvlst.push_back(gmsh_cv_info());
        gmsh_cv_info &ci = cvlst.back();
        ci.id = id; ci.type = type; ci.region = region;

        if (version == 2) { /* For version 2 */
          ci.set_nb_nodes();
          cv_nb_nodes = unsigned(ci.nodes.size());
        }
        else
          ci.nodes.resize(cv_nb_nodes);

        // cout << "cv_nb_nodes = " << cv_nb_nodes << endl;

        for (size_type i=0; i < cv_nb_nodes; ++i) {
          size_type j;
          f >> j;
	  std::map<size_type, size_type>::iterator
	    it = msh_node_2_getfem_node.find(j);
          GMM_ASSERT1(it != msh_node_2_getfem_node.end(),
		      "Invalid node ID " << j << " in gmsh convex "
		      << (ci.id + 1));
	  ci.nodes[i] = it->second;
        }
        ci.set_pgt();
        // Reordering nodes for certain elements (should be completed ?)
        switch(ci.type) {
        case 3 : std::swap(ci.nodes[2], ci.nodes[3]); break;
        case 5 : { /* First order hexaedron */
          std::vector<size_type> tmp_nodes(8);
          tmp_nodes[0] = ci.nodes[0]; 
          tmp_nodes[1] = ci.nodes[1];
          tmp_nodes[2] = ci.nodes[3]; 
          tmp_nodes[3] = ci.nodes[2];
          tmp_nodes[4] = ci.nodes[4]; 
          tmp_nodes[5] = ci.nodes[5];
          tmp_nodes[6] = ci.nodes[7]; 
          tmp_nodes[7] = ci.nodes[6];
          ci.nodes[0] = tmp_nodes[0], ci.nodes[1] = tmp_nodes[1],
            ci.nodes[2] = tmp_nodes[2];
          ci.nodes[3] = tmp_nodes[3], ci.nodes[4] = tmp_nodes[4],
            ci.nodes[5] = tmp_nodes[5];
          ci.nodes[6] = tmp_nodes[6], ci.nodes[7] = tmp_nodes[7];
        }
          break;
        case 8 : { /* Second order line */
          std::vector<size_type> tmp_nodes(3);
          tmp_nodes[0] = ci.nodes[0]; tmp_nodes[1] = ci.nodes[2];
          tmp_nodes[2] = ci.nodes[1];

          ci.nodes[0] = tmp_nodes[0]; ci.nodes[1] = tmp_nodes[1];
          ci.nodes[2] = tmp_nodes[2];
        }
          break;
		 case 11: { /* Second order tetrahedron */
			std::vector<bgeot::size_type> tmp_nodes(10);
			tmp_nodes[0] = ci.nodes[0], tmp_nodes[1] = ci.nodes[4],
		    tmp_nodes[2] = ci.nodes[1];
			tmp_nodes[3] = ci.nodes[6], tmp_nodes[4] = ci.nodes[5],
		    tmp_nodes[5] = ci.nodes[2];
			tmp_nodes[6] = ci.nodes[7], tmp_nodes[7] = ci.nodes[9],
		    tmp_nodes[8] = ci.nodes[8],
		    tmp_nodes[9] = ci.nodes[3];

			ci.nodes[0] = tmp_nodes[0], ci.nodes[1] = tmp_nodes[1],
		    ci.nodes[2] = tmp_nodes[2];
			ci.nodes[3] = tmp_nodes[3], ci.nodes[4] = tmp_nodes[4],
		    ci.nodes[5] = tmp_nodes[5];
			ci.nodes[6] = tmp_nodes[6], ci.nodes[7] = tmp_nodes[7],
		    ci.nodes[8] = tmp_nodes[8],
		    ci.nodes[9] = tmp_nodes[9];
		  }
		 break;
        case 9 : /* Second order triangle */
          std::vector<size_type> tmp_nodes(6);
          tmp_nodes[0] = ci.nodes[0], tmp_nodes[1] = ci.nodes[3],
            tmp_nodes[2] = ci.nodes[1];
          tmp_nodes[3] = ci.nodes[5], tmp_nodes[4] = ci.nodes[4],
            tmp_nodes[5] = ci.nodes[2];

          ci.nodes[0] = tmp_nodes[0], ci.nodes[1] = tmp_nodes[1],
            ci.nodes[2] = tmp_nodes[2];
          ci.nodes[3] = tmp_nodes[3], ci.nodes[4] = tmp_nodes[4],
            ci.nodes[5] = tmp_nodes[5];
          break;
		  
        }
      }
    }
    nb_cv = cvlst.size();
    if (cvlst.size()) {
      std::sort(cvlst.begin(), cvlst.end());
      unsigned N = cvlst.front().pgt->dim();
      for (size_type cv=0; cv < nb_cv; ++cv) {
        bool cvok = false;
        gmsh_cv_info &ci = cvlst[cv];
        //cout << "importing cv dim=" << int(ci.pgt->dim()) << " N=" << N
        //     << " region: " << ci.region << "\n";
        if (ci.pgt->dim() == N) {
          size_type ic = m.add_convex(ci.pgt, ci.nodes.begin()); cvok = true;
          m.region(ci.region).add(ic);
        } else if (ci.pgt->dim() == N-1) {
          bgeot::mesh_structure::ind_cv_ct ct =
            m.convex_to_point(ci.nodes[0]);
          for (bgeot::mesh_structure::ind_cv_ct::const_iterator
                 it = ct.begin(); it != ct.end(); ++it) {
            for (unsigned face=0;
                 face < m.structure_of_convex(*it)->nb_faces(); ++face) {
              if (m.is_convex_face_having_points(*it,face,
                                                 short_type(ci.nodes.size()),
                                                 ci.nodes.begin())) {
                m.region(ci.region).add(*it,face);
                cvok = true;
              }
            }
          }
          if (!cvok)
            GMM_WARNING2("face not found ... " << endl);
        }
        if (!cvok)
          GMM_WARNING2("gmsh import ignored a convex of type "
                       << bgeot::name_of_geometric_trans(ci.pgt));
      }
    }
    maybe_remove_last_dimension(m);
  }

  /* mesh file from GiD [http://gid.cimne.upc.es/]

  supports linear and quadratic elements (quadrilaterals, use 9(or 27)-noded elements)
  */
  static void import_gid_msh_file(std::ifstream& f, mesh& m) {
    gmm::stream_standard_locale sl(f);
    /* read the node list */
    size_type dim;
    enum { LIN,TRI,QUAD,TETR, PRISM, HEX,BADELTYPE } eltype=BADELTYPE;
    size_type nnode = 0;
    std::map<size_type, size_type> msh_node_2_getfem_node;
    std::vector<size_type> cv_nodes, getfem_cv_nodes;
    bool nodes_done = false;
    do {
      if (!f.eof()) f >> std::ws;
      if (f.eof() || !bgeot::read_until(f, "MESH")) break;
      std::string selemtype;
      f >> bgeot::skip("DIMENSION") >> dim
        >> bgeot::skip("ELEMTYPE") >> std::ws
        >> selemtype
        >> bgeot::skip("NNODE") >> nnode;
      if (bgeot::casecmp(selemtype, "linear")==0) { eltype = LIN;  }
      else if (bgeot::casecmp(selemtype, "triangle")==0) { eltype = TRI; }
      else if (bgeot::casecmp(selemtype, "quadrilateral")==0) { eltype = QUAD; }
      else if (bgeot::casecmp(selemtype, "tetrahedra")==0) { eltype = TETR; }
      else if (bgeot::casecmp(selemtype, "prisma")==0) { eltype = PRISM; }
      else if (bgeot::casecmp(selemtype, "hexahedra")==0) { eltype = HEX; }
      else GMM_ASSERT1(false, "unknown element type '"<< selemtype << "'");
      GMM_ASSERT1(!f.eof(), "File ended before coordinates");
      f >> bgeot::skip("COORDINATES");
      if (!nodes_done) {
        dal::dynamic_array<base_node> gid_nodes;
        dal::bit_vector gid_nodes_used;
        do {
          //cerr << "reading coordinates " << std::streamoff(f.tellg()) << "\n";
          std::string ls;
          f >> std::ws;
          std::getline(f,ls);
          if (bgeot::casecmp(ls, "END COORDINATES", 15)==0) break;
          std::stringstream s; s << ls;
          size_type id;
          s >> id;

          gid_nodes[id].resize(dim); gid_nodes_used.add(id);
          for (size_type i=0; i < dim; ++i) s >> gid_nodes[id][i];
          //cerr << "ppoint " << id << ", " << n << endl;
        } while (true);

        GMM_ASSERT1(gid_nodes_used.card() != 0, "no nodes in the mesh!");

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
          msh_node_2_getfem_node[ip] = m.add_point(n);
        }
      }

      bgeot::read_until(f, "ELEMENTS");
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
      default: GMM_ASSERT1(false, ""); break;
      }
      GMM_ASSERT1(pgt != NULL, "unknown element type " << selemtype
                  << " with " << nnode << "nodes");
      do {
        std::string ls;
        f >> std::ws;
        std::getline(f,ls);
        if (bgeot::casecmp(ls, "END ELEMENTS", 12)==0) break;
        std::stringstream s; s << ls;
        size_type cv_id;
        s >> cv_id;
        cv_nodes.resize(nnode);
        for (size_type i=0; i < nnode; ++i) {
          size_type j;
          s >> j;
	  std::map<size_type, size_type>::iterator
	    it = msh_node_2_getfem_node.find(j);
          GMM_ASSERT1(it != msh_node_2_getfem_node.end(),
		      "Invalid node ID " << j << " in GiD convex " << cv_id);
	  cv_nodes[i] = it->second;
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
    double pdec = pow(10.,double(ndec));
    if (x == 0) return 0.;
    x = gmm::abs(x);
    while (x > 1) { x /= 10.0; p*=10; }
    while (x < 0.1) { x *= 10.0; p/=10; }
    //cerr << "x=" << x << ", p=" << p << ", pdec=" << pdec << "\n";
    x = s * (floor(x * pdec + 0.5) / pdec) * p;
    return x;
  }


  /* mesh file from noboite [http://www.distene.com/fr/corp/newsroom16.html] */
  static void import_noboite_msh_file(std::ifstream& f, mesh& m) {
    
    using namespace std;
    gmm::stream_standard_locale sl(f);
    
    ofstream fichier_GiD("noboite_to_GiD.gid",    ios::out | ios::trunc );  //déclaration du flux et ouverture du fichier
    
    fichier_GiD << "MESH    dimension 3 ElemType Tetrahedra  Nnode 4"<<endl;
    
    
    int i,NE,NP,ligne_debut_NP;
    
    /*NE: nombre d'elements (premier nombre du fichier .noboite)
      NP: nombre de points (deuxieme nombre du fichier .noboite)
      ligne_debut_NP: la ligne commence la liste des coordonnees des points dans le 		fichier .noboite*/
    
    f >> NE>>NP;
    ligne_debut_NP=NE*4/6+3;
    
    //passer 3 premiers lignes du fichier .noboite (la liste des elements commence a la 		quatrieme ligne)
    string contenu;
    for (i=1; i<=ligne_debut_NP; i++){
      getline(f, contenu);
    }
    

    /*-----------------------------------------------------------------------
      Lire les coordonnees dans .noboite
      -----------------------------------------------------------------------*/
    fichier_GiD << "Coordinates" <<endl;
    
    for (i=1; i<=NP; i++){
      float coor_x,coor_y,coor_z;
      
      fichier_GiD << i<<" ";
      
      f>>coor_x >>coor_y >>coor_z;
      fichier_GiD<<coor_x<<" " <<coor_y <<" "<<coor_z <<endl;
      
    }
    
    fichier_GiD << "end coordinates" <<endl<<endl;
    
    /*-----------------------------------------------------------------------
      Lire les elements dans .noboite et ecrire au . gid
      ------------------------------------------------------------------------*/
    
    //revenir au debut du fichier . noboite, puis passer les trois premiere lignes
    f.seekg(0, ios::beg);
    for (i=1; i<=3; i++){
      getline(f, contenu);
    }
    
    
    fichier_GiD << "Elements" <<endl;
    
    for (i=1; i<=NE; i++){
      float elem_1,elem_2,elem_3,elem_4;
      
      fichier_GiD << i<<" ";
      f>>elem_1>>elem_2>>elem_3>>elem_4;
      fichier_GiD<<elem_1<<" " <<elem_2 <<" "<<elem_3<<" "<<elem_4<<" 1"<<endl;
      
    }
    fichier_GiD << "end elements" <<endl<<endl;
    
    
    
    if(fichier_GiD)  // si l'ouverture a réussi
      {
	// instructions
	fichier_GiD.close();  // on referme le fichier
      }
    else  // sinon
      cerr << "Erreur à l'ouverture !" << endl;
    
    if(f)  // si l'ouverture a réussi
      {
	// instructions
	f.close();  // on referme le fichier
      }
    else  // sinon
      cerr << "Erreur à l'ouverture !" << endl;
    
    // appeler sunroutine import_gid_msh_file
    //import_mesh(const std::string& "noboite_to_GiD.gid", mesh& msh)
    ifstream fichier1_GiD("noboite_to_GiD.gid", ios::in);
    import_gid_msh_file(fichier1_GiD, m);
    
    //      return 0;
  }
  
  /* mesh file from emc2 [http://pauillac.inria.fr/cdrom/prog/unix/emc2/eng.htm], am_fmt format

  (only triangular 2D meshes)
  */
  static void import_am_fmt_file(std::ifstream& f, mesh& m) {
    gmm::stream_standard_locale sl(f);
    /* read the node list */
    std::vector<size_type> tri;
    size_type nbs,nbt;
    base_node P(2);
    f >> nbs >> nbt; bgeot::read_until(f,"\n");
    tri.resize(nbt*3);
    for (size_type i=0; i < nbt*3; ++i) f >> tri[i];
    for (size_type j=0; j < nbs; ++j) {
      f >> P[0] >> P[1];
      cerr.precision(16);
      P[0]=round_to_nth_significant_number(P[0],6); // force 9.999999E-1 to be converted to 1.0
      P[1]=round_to_nth_significant_number(P[1],6);
      size_type jj = m.add_point(P);
      GMM_ASSERT1(jj == j, "ouch");
    }
    for (size_type i=0; i < nbt*3; i+=3)
      m.add_triangle(tri[i]-1,tri[i+1]-1,tri[i+2]-1);
  }

  /* mesh file from emc2 [http://pauillac.inria.fr/cdrom/prog/unix/emc2/eng.htm], am_fmt format

  triangular/quadrangular 2D meshes
  */
  static void import_emc2_mesh_file(std::ifstream& f, mesh& m) {
    gmm::stream_standard_locale sl(f);
    /* read the node list */
    std::vector<size_type> tri;
    size_type nbs=0,nbt=0,nbq=0,dummy;
    base_node P(2);
    bgeot::read_until(f,"Vertices");
    f >> nbs;
    for (size_type j=0; j < nbs; ++j) {
      f >> P[0] >> P[1] >> dummy;
      size_type jj = m.add_point(P);
      GMM_ASSERT1(jj == j, "ouch");
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

  void import_mesh(const std::string& filename, const std::string& format,
                   mesh& m) {
    m.clear();
    try {

      if (bgeot::casecmp(format,"structured")==0)
        { regular_mesh(m, filename); return; }

      std::ifstream f(filename.c_str());
      GMM_ASSERT1(f.good(), "can't open file " << filename);
      /* throw exceptions when an error occurs */
      f.exceptions(std::ifstream::badbit | std::ifstream::failbit);
      import_mesh(f,format,m);
      f.close();
    }
    catch (failure_error& exc) {
      m.clear();
      throw exc;
    }
    catch (std::ios_base::failure& exc) {
      m.clear();
      GMM_ASSERT1(false, "error while importing " << format
                  << " mesh file \"" << filename << "\" : " << exc.what());
    }
  }

  void import_mesh(std::ifstream& f, const std::string& format,
                   mesh& m) {
    if (bgeot::casecmp(format,"gmsh")==0)
      import_gmsh_msh_file(f,m);
    else if (bgeot::casecmp(format,"gmshv2")==0)/* deprecate */
      import_gmsh_msh_file(f,m,2);
    else if (bgeot::casecmp(format,"gid")==0)
      import_gid_msh_file(f,m);
    else if (bgeot::casecmp(format,"noboite")==0)
      import_noboite_msh_file(f,m);
    else if (bgeot::casecmp(format,"am_fmt")==0)
      import_am_fmt_file(f,m);
    else if (bgeot::casecmp(format,"emc2_mesh")==0)
      import_emc2_mesh_file(f,m);
    else GMM_ASSERT1(false, "cannot import "
                     << format << " mesh type : unknown mesh type");
  }

  void import_mesh(const std::string& filename, mesh& msh) {
    if (filename.compare(0,4,"gid:")==0)
      getfem::import_mesh(filename.substr(4), "gid", msh);
    else if (filename.compare(0,8,"noboite:") == 0)
      getfem::import_mesh(filename.substr(8), "noboite", msh);
    else if (filename.compare(0,5,"gmsh:") == 0)
      getfem::import_mesh(filename.substr(5), "gmsh", msh);
    else if (filename.compare(0,7,"gmshv2:") == 0)
      getfem::import_mesh(filename.substr(7), "gmshv2", msh);
    else if (filename.compare(0,7,"am_fmt:") == 0)
      getfem::import_mesh(filename.substr(7), "am_fmt", msh);
    else if (filename.compare(0,10,"emc2_mesh:") == 0)
      getfem::import_mesh(filename.substr(10), "emc2_mesh", msh);
    else if (filename.compare(0,11,"structured:") == 0)
      getfem::import_mesh(filename.substr(11), "structured", msh);
    else msh.read_from_file(filename);
  }

  void maybe_remove_last_dimension(mesh &m) {
    bool is_flat = true;
    unsigned N = m.dim(); if (N < 1) return;
    for (dal::bv_visitor i(m.points().index()); !i.finished(); ++i)
      if (m.points()[i][N-1] != 0) is_flat = 0;
    if (is_flat) {
      base_matrix M(N-1,N);
      for (unsigned i=0; i < N-1; ++i) M(i,i) = 1;
      m.transformation(M);
    }
  }

}
