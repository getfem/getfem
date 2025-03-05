/*===========================================================================

 Copyright (C) 2000-2020 Julien Pommier

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
        pgt = bgeot::pyramid_QK_geotrans(1);
      } break;
      case 8: { /* 2ND ORDER LINE */
        pgt = bgeot::simplex_geotrans(1,2);
      } break;
      case 9: { /* 2ND ORDER TRIANGLE */
        pgt = bgeot::simplex_geotrans(2,2);
      } break;
      case 10: { /* 2ND ORDER QUADRANGLE */
        pgt = bgeot::parallelepiped_geotrans(2,2);
      } break;
      case 11: { /* 2ND ORDER TETRAHEDRON (10-NODE) */
        pgt = bgeot::simplex_geotrans(3,2);
      } break;
      case 12: { /* 2ND ORDER HEXAHEDRON (27-NODE) */
        pgt = bgeot::parallelepiped_geotrans(3,2);
      } break;
      case 14: { /* 2ND ORDER PYRAMID (14-NODE) */
        pgt = bgeot::pyramid_QK_geotrans(2);
      } break;
      case 15: { /* POINT */
        GMM_WARNING2("ignoring point element");
      } break;
      case 16: { /* INCOMPLETE 2ND ORDER QUADRANGLE (8-NODE) */
        pgt = bgeot::Q2_incomplete_geotrans(2);
      } break;
      case 17: { /* INCOMPLETE 2ND ORDER HEXAHEDRON (20-NODE) */
        pgt = bgeot::Q2_incomplete_geotrans(3);
      } break;
      case 19: { /* INCOMPLETE 2ND ORDER PYRAMID (13-NODE) */
        pgt = bgeot::pyramid_Q2_incomplete_geotrans();
      } break;
      case 26: { /* 3RD ORDER LINE */
        pgt = bgeot::simplex_geotrans(1,3);
      } break;
      case 21: { /* 3RD ORDER TRIANGLE */
        pgt = bgeot::simplex_geotrans(2,3);
      } break;
      case 23: { /* 4TH ORDER TRIANGLE */
        pgt = bgeot::simplex_geotrans(2, 4);
      } break;
      case 27: { /* 4TH ORDER LINE */
        pgt = bgeot::simplex_geotrans(1, 4);
      } break;
      default: { /* UNKNOWN .. */
        /* higher order elements : to be done .. */
        GMM_ASSERT1(false, "gmsh element type " << type << " is unknown.");
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
        nodes.resize(5);
      } break;
      case 8: { /* 2ND ORDER LINE */
        nodes.resize(3);
      } break;
      case 9: { /* 2ND ORDER TRIANGLE */
        nodes.resize(6);
      } break;
      case 10: { /* 2ND ORDER QUADRANGLE */
        nodes.resize(9);
      } break;
      case 11: { /* 2ND ORDER TETRAHEDRON (10-NODE) */
        nodes.resize(10);
      } break;
      case 12: { /* 2ND ORDER HEXAHEDRON (27-NODE) */
        nodes.resize(27);
      } break;
      case 14: { /* 2ND ORDER PYRAMID (14-NODE) */
        nodes.resize(14);
      } break;
      case 15: { /* POINT */
        nodes.resize(1);
      } break;
      case 16: { /* INCOMPLETE 2ND ORDER QUADRANGLE (8-NODE) */
        nodes.resize(8);
      } break;
      case 17: { /* INCOMPLETE 2ND ORDER HEXAHEDRON (20-NODE) */
        nodes.resize(20);
      } break;
      case 19: { /* INCOMPLETE 2ND ORDER PYRAMID (13-NODE) */
        nodes.resize(13);
      } break;
      case 26: { /* 3RD ORDER LINE */
        nodes.resize(4);
      } break;
      case 21: { /* 3RD ORDER TRIANGLE */
        nodes.resize(10);
      } break;
      case 23: { /* 4TH ORDER TRIANGLE */
        nodes.resize(15);
      } break;
      case 27: { /* 4TH ORDER LINE */
        nodes.resize(5);
      } break;
      default: { /* UNKNOWN .. */
        /* higher order elements : to be done .. */
        GMM_ASSERT1(false, "the gmsh element type " << type << " is unknown..");
      } break;
      }
    }

    bool operator<(const gmsh_cv_info& other) const {
      unsigned this_dim = (type == 15) ? 0 : pgt->dim();
      unsigned other_dim = (other.type == 15) ? 0 : other.pgt->dim();
      if (this_dim == other_dim) return region < other.region;
      return this_dim > other_dim;
    }
  };

  std::map<std::string, size_type> read_region_names_from_gmsh_mesh_file(std::istream& f)
  {
    std::map<std::string, size_type> region_map;
    bgeot::read_until(f, "$PhysicalNames");
    size_type nb_regions;
    f >> nb_regions;
    size_type rt,ri;
    std::string region_name;
    for (size_type region_cnt=0; region_cnt < nb_regions; ++ region_cnt) {
      f >> rt >> ri;
      std::getline(f, region_name);
      /* trim the string to the quote character front and back*/
      size_t pos = region_name.find_first_of("\"");
      if (pos != region_name.npos) {
        region_name.erase(0, pos+1);
        pos = region_name.find_last_of("\"");
        region_name.erase(pos);
      }
      region_map[region_name] = ri;
    }

    return region_map;
  }

  /*
     Format version 1 [for gmsh version < 2.0].
     structure: $NOD list_of_nodes $ENDNOD $ELT list_of_elt $ENDELT

     Format version 2 [for gmsh version >= 2.0]. Modification of some keywords
     and more than one tag is authorized for a specific region.

     structure: $Nodes list_of_nodes $EndNodes $Elements list_of_elt
     $EndElements

     Lower dimensions elements in the regions of lower_dim_convex_rg will
     be imported as independant convexes.

     If add_all_element_type is set to true, elements with lower dimension
     than highest dimension and that are not part of other element's face will
     be imported as independent elements.

     for gmsh and gid meshes, the mesh nodes are always 3D, so for a 2D mesh
     if remove_last_dimension == true the z-component of nodes will be removed
  */
  static void import_gmsh_mesh_file
  (std::istream& f, mesh& m, int deprecate=0,
   std::map<std::string, size_type> *region_map=NULL,
   std::set<size_type> *lower_dim_convex_rg=NULL,
   bool add_all_element_type = false,
   bool remove_last_dimension = true,
   std::map<size_type, std::set<size_type>> *nodal_map = NULL,
   bool remove_duplicated_nodes = true) {
    gmm::stream_standard_locale sl(f);
    // /* print general warning */
    // GMM_WARNING3("  All regions must have different number!");

    /* print deprecate warning */
    if (deprecate!=0){
      GMM_WARNING4("" << endl
                << "  deprecate: " << endl
                << "   static void" << endl
                << "   import_gmsh_mesh_file(std::istream& f,"
                << " mesh& , int version)" << endl
                << "  replace with:" << endl
                << "   static void" << endl
                << "   import_gmsh_mesh_file(std::istream& f,"
                << " mesh&)");
    }

    /* read the version */
    double version;
    std::string header;
    f >> header;
    if (bgeot::casecmp(header,"$MeshFormat")==0)
      f >> version;
    else if (bgeot::casecmp(header,"$NOD")==0)
      version = 1;
    else
      GMM_ASSERT1(false, "can't read Gmsh format: " << header);

    /* read the region names */
    if (region_map != NULL) {
      if (version >= 2.) {
        *region_map = read_region_names_from_gmsh_mesh_file(f);
      }
    }
    /* read the node list */
    if (version >= 2.)
      bgeot::read_until(f, "$Nodes"); /* Format versions 2 and 4 */

    size_type nb_block, nb_node, dummy;
    std::string dummy2;
    // cout << "version = " << version << endl;
    if (version >= 4.05) {
      f >> nb_block >> nb_node; bgeot::read_until(f, "\n");
    } else if (version >= 4.) {
      f >> nb_block >> nb_node;
    } else {
      nb_block = 1;
      f >> nb_node;
    }

    // cerr << "reading nodes..[nb=" << nb_node << "]\n";
    std::map<size_type, size_type> msh_node_2_getfem_node;
     std::vector<size_type> inds(nb_node);
    for (size_type block=0; block < nb_block; ++block) {
      if (version >= 4.)
        f >> dummy >> dummy >> dummy >> nb_node;
      // cout << "nb_nodes = " << nb_node << endl;

      inds.resize(nb_node);
      if (version >= 4.05) {
        for (size_type node_cnt=0; node_cnt < nb_node; ++node_cnt)
          f >> inds[node_cnt];
      }

      for (size_type node_cnt=0; node_cnt < nb_node; ++node_cnt) {
        size_type node_id;
        base_node n{0,0,0};
        if (version < 4.05) f >> node_id; else node_id = inds[node_cnt];

        f >> n[0] >> n[1] >> n[2];
        msh_node_2_getfem_node[node_id]
          = m.add_point(n, remove_duplicated_nodes ? 0. : -1.);
      }
    }

    if (version >= 2.)
      bgeot::read_until(f, "$Endnodes"); /* Format versions 2 and 4 */
    else
      bgeot::read_until(f, "$ENDNOD");

    /* read the elements */
    if (version >= 2.)
      bgeot::read_until(f, "$Elements"); /* Format versions 2 and 4 */
    else
      bgeot::read_until(f, "$ELM");

    size_type nb_cv;
    if (version >= 4.05) {
      f >> nb_block >> nb_cv; bgeot::read_until(f, "\n");
    } else if (version >= 4.) { /* Format version 4 */
      f >> nb_block >> nb_cv;
    } else {
      nb_block = 1;
      f >> nb_cv;
    }
    // cout << "nb_bloc = " << nb_block << " nb_cv = " << nb_cv << endl;

    std::vector<gmsh_cv_info> cvlst; cvlst.reserve(nb_cv);
    dal::bit_vector reg;
    for (size_type block=0; block < nb_block; ++block) {
      unsigned dimr, type, region;
      if (version >= 4.) { /* Format version 4 */
        f >> dimr >> region >> type >> nb_cv;
        if (reg.is_in(region)) {
          GMM_WARNING2("Two regions share the same number, "
                       "the region numbering is modified");
          while (reg.is_in(region)) region += 5;
        }
        reg.add(region);
      }
      for (size_type cv=0; cv < nb_cv; ++cv) {

        cvlst.push_back(gmsh_cv_info());
        gmsh_cv_info &ci = cvlst.back();
        f >> ci.id;
        ci.id--; /* gmsh numbering starts at 1 */

        unsigned cv_nb_nodes;
        if (version >= 2.) { /* For versions 2 and 4 */
          if (int(version) == 2) { /* Format version 2 */
            unsigned nbtags;
            f >> type >> nbtags;
            GMM_ASSERT1(nbtags > 0 && nbtags <= 3,
                        "Number of tags " << nbtags << " is not managed.");
            f >> region;
            if (nbtags > 1) f >> dummy;
            if (nbtags > 2) f >> dummy;
          }
          ci.type = type;
          ci.set_nb_nodes();
          cv_nb_nodes = unsigned(ci.nodes.size());
        } else if (int(version) == 1) {
          f >> type >> region >> dummy >> cv_nb_nodes;
          ci.type = type;
          ci.nodes.resize(cv_nb_nodes);
        }
        ci.region = region;

        // cout << "cv_nb_nodes = " << cv_nb_nodes << endl;

        for (size_type i=0; i < cv_nb_nodes; ++i) {
          size_type j;
          f >> j;
          const auto it = msh_node_2_getfem_node.find(j);
          GMM_ASSERT1(it != msh_node_2_getfem_node.end(),
                      "Invalid node ID " << j << " in gmsh element "
                      << (ci.id + 1));
          ci.nodes[i] = it->second;
        }
        if (ci.type != 15)
          ci.set_pgt();
        // Reordering nodes for certain elements (should be completed ?)
        // http://www.geuz.org/gmsh/doc/texinfo/gmsh.html#Node-ordering
        std::vector<size_type> tmp_nodes(ci.nodes);
        switch(ci.type) {
        case 3 : {
          ci.nodes[2] = tmp_nodes[3];
          ci.nodes[3] = tmp_nodes[2];
        } break;
        case 5 : { /* First order hexaedron */
          //ci.nodes[0] = tmp_nodes[0];
          //ci.nodes[1] = tmp_nodes[1];
          ci.nodes[2] = tmp_nodes[3];
          ci.nodes[3] = tmp_nodes[2];
          //ci.nodes[4] = tmp_nodes[4];
          //ci.nodes[5] = tmp_nodes[5];
          ci.nodes[6] = tmp_nodes[7];
          ci.nodes[7] = tmp_nodes[6];
        } break;
        case 7 : { /* first order pyramid */
          //ci.nodes[0] = tmp_nodes[0];
          ci.nodes[1] = tmp_nodes[2];
          ci.nodes[2] = tmp_nodes[1];
          // ci.nodes[3] = tmp_nodes[3];
          // ci.nodes[4] = tmp_nodes[4];
        } break;
        case 8 : { /* Second order line */
          //ci.nodes[0] = tmp_nodes[0];
          ci.nodes[1] = tmp_nodes[2];
          ci.nodes[2] = tmp_nodes[1];
        } break;
        case 9 : { /* Second order triangle */
          //ci.nodes[0] = tmp_nodes[0];
          ci.nodes[1] = tmp_nodes[3];
          ci.nodes[2] = tmp_nodes[1];
          ci.nodes[3] = tmp_nodes[5];
          //ci.nodes[4] = tmp_nodes[4];
          ci.nodes[5] = tmp_nodes[2];
        } break;
        case 10 : { /* Second order quadrangle */
          //ci.nodes[0] = tmp_nodes[0];
          ci.nodes[1] = tmp_nodes[4];
          ci.nodes[2] = tmp_nodes[1];
          ci.nodes[3] = tmp_nodes[7];
          ci.nodes[4] = tmp_nodes[8];
          //ci.nodes[5] = tmp_nodes[5];
          ci.nodes[6] = tmp_nodes[3];
          ci.nodes[7] = tmp_nodes[6];
          ci.nodes[8] = tmp_nodes[2];
        } break;
        case 11: { /* Second order tetrahedron */
          //ci.nodes[0] = tmp_nodes[0];
          ci.nodes[1] = tmp_nodes[4];
          ci.nodes[2] = tmp_nodes[1];
          ci.nodes[3] = tmp_nodes[6];
          ci.nodes[4] = tmp_nodes[5];
          ci.nodes[5] = tmp_nodes[2];
          ci.nodes[6] = tmp_nodes[7];
          ci.nodes[7] = tmp_nodes[9];
          //ci.nodes[8] = tmp_nodes[8];
          ci.nodes[9] = tmp_nodes[3];
        } break;
        case 12: { /* Second order hexahedron */
          //ci.nodes[0] = tmp_nodes[0];
          ci.nodes[1] = tmp_nodes[8];
          ci.nodes[2] = tmp_nodes[1];
          ci.nodes[3] = tmp_nodes[9];
          ci.nodes[4] = tmp_nodes[20];
          ci.nodes[5] = tmp_nodes[11];
          ci.nodes[6] = tmp_nodes[3];
          ci.nodes[7] = tmp_nodes[13];
          ci.nodes[8] = tmp_nodes[2];
          ci.nodes[9] = tmp_nodes[10];
          ci.nodes[10] = tmp_nodes[21];
          ci.nodes[11] = tmp_nodes[12];
          ci.nodes[12] = tmp_nodes[22];
          ci.nodes[13] = tmp_nodes[26];
          ci.nodes[14] = tmp_nodes[23];
          //ci.nodes[15] = tmp_nodes[15];
          ci.nodes[16] = tmp_nodes[24];
          ci.nodes[17] = tmp_nodes[14];
          ci.nodes[18] = tmp_nodes[4];
          ci.nodes[19] = tmp_nodes[16];
          ci.nodes[20] = tmp_nodes[5];
          ci.nodes[21] = tmp_nodes[17];
          ci.nodes[22] = tmp_nodes[25];
          ci.nodes[23] = tmp_nodes[18];
          ci.nodes[24] = tmp_nodes[7];
          ci.nodes[25] = tmp_nodes[19];
          ci.nodes[26] = tmp_nodes[6];
        } break;
        case 14: { /* Second order pyramid */
          //ci.nodes[0] = tmp_nodes[0];
          ci.nodes[1] = tmp_nodes[5];
          ci.nodes[2] = tmp_nodes[1];
          ci.nodes[3] = tmp_nodes[6];
          ci.nodes[4] = tmp_nodes[13];
          ci.nodes[5] = tmp_nodes[8];
          ci.nodes[6] = tmp_nodes[3];
          ci.nodes[7] = tmp_nodes[10];
          ci.nodes[8] = tmp_nodes[2];
          ci.nodes[9] = tmp_nodes[7];
          ci.nodes[10] = tmp_nodes[9];
          ci.nodes[11] = tmp_nodes[12];
          ci.nodes[12] = tmp_nodes[11];
          ci.nodes[13] = tmp_nodes[4];
        } break;
        case 16 : { /* Incomplete second order quadrangle */
          //ci.nodes[0] = tmp_nodes[0];
          ci.nodes[1] = tmp_nodes[4];
          ci.nodes[2] = tmp_nodes[1];
          ci.nodes[3] = tmp_nodes[7];
          ci.nodes[4] = tmp_nodes[5];
          ci.nodes[5] = tmp_nodes[3];
          ci.nodes[6] = tmp_nodes[6];
          ci.nodes[7] = tmp_nodes[2];
        } break;
        case 17: { /* Incomplete second order hexahedron */
          //ci.nodes[0] = tmp_nodes[0];
          ci.nodes[1] = tmp_nodes[8];
          ci.nodes[2] = tmp_nodes[1];
          ci.nodes[3] = tmp_nodes[9];
          ci.nodes[4] = tmp_nodes[11];
          ci.nodes[5] = tmp_nodes[3];
          ci.nodes[6] = tmp_nodes[13];
          ci.nodes[7] = tmp_nodes[2];
          ci.nodes[8] = tmp_nodes[10];
          ci.nodes[9] = tmp_nodes[12];
          ci.nodes[10] = tmp_nodes[15];
          ci.nodes[11] = tmp_nodes[14];
          ci.nodes[12] = tmp_nodes[4];
          ci.nodes[13] = tmp_nodes[16];
          ci.nodes[14] = tmp_nodes[5];
          ci.nodes[15] = tmp_nodes[17];
          ci.nodes[16] = tmp_nodes[18];
          ci.nodes[17] = tmp_nodes[7];
          ci.nodes[18] = tmp_nodes[19];
          ci.nodes[19] = tmp_nodes[6];
        } break;
        case 19: { /* Incomplete second order pyramid */
          //ci.nodes[0] = tmp_nodes[0];
          ci.nodes[1] = tmp_nodes[5];
          ci.nodes[2] = tmp_nodes[1];
          ci.nodes[3] = tmp_nodes[6];
          ci.nodes[4] = tmp_nodes[8];
          ci.nodes[5] = tmp_nodes[3];
          ci.nodes[6] = tmp_nodes[10];
          ci.nodes[7] = tmp_nodes[2];
          ci.nodes[8] = tmp_nodes[7];
          //ci.nodes[9] = tmp_nodes[9];
          ci.nodes[10] = tmp_nodes[12];
          //ci.nodes[11] = tmp_nodes[11];
          ci.nodes[12] = tmp_nodes[4];
        } break;
        case 26 : { /* Third order line */
          //ci.nodes[0] = tmp_nodes[0];
          ci.nodes[1] = tmp_nodes[2];
          ci.nodes[2] = tmp_nodes[3];
          ci.nodes[3] = tmp_nodes[1];
        } break;
        case 21 : { /* Third order triangle */
          //ci.nodes[0] = tmp_nodes[0];
          ci.nodes[1] = tmp_nodes[3];
          ci.nodes[2] = tmp_nodes[4];
          ci.nodes[3] = tmp_nodes[1];
          ci.nodes[4] = tmp_nodes[8];
          ci.nodes[5] = tmp_nodes[9];
          ci.nodes[6] = tmp_nodes[5];
          //ci.nodes[7] = tmp_nodes[7];
          ci.nodes[8] = tmp_nodes[6];
          ci.nodes[9] = tmp_nodes[2];
        } break;
        case 23: { /* Fourth order triangle */
        //ci.nodes[0]  = tmp_nodes[0];
          ci.nodes[1]  = tmp_nodes[3];
          ci.nodes[2]  = tmp_nodes[4];
          ci.nodes[3]  = tmp_nodes[5];
          ci.nodes[4]  = tmp_nodes[1];
          ci.nodes[5]  = tmp_nodes[11];
          ci.nodes[6]  = tmp_nodes[12];
          ci.nodes[7]  = tmp_nodes[13];
          ci.nodes[8]  = tmp_nodes[6];
          ci.nodes[9]  = tmp_nodes[10];
          ci.nodes[10] = tmp_nodes[14];
          ci.nodes[11] = tmp_nodes[7];
          ci.nodes[12] = tmp_nodes[9];
          ci.nodes[13] = tmp_nodes[8];
          ci.nodes[14] = tmp_nodes[2];
        } break;
        case 27: { /* Fourth order line */
        //ci.nodes[0]  = tmp_nodes[0];
          ci.nodes[1]  = tmp_nodes[2];
          ci.nodes[2]  = tmp_nodes[3];
          ci.nodes[3]  = tmp_nodes[4];
          ci.nodes[4]  = tmp_nodes[1];
        } break;
        }
      }
    }

    nb_cv = cvlst.size();
    if (cvlst.size()) {
      std::sort(cvlst.begin(), cvlst.end());
      if (cvlst.front().type == 15) {
        GMM_WARNING2("Only nodes defined in the mesh! No elements are added.");
        return;
      }

      unsigned N = cvlst.front().pgt->dim();
      for (size_type cv=0; cv < nb_cv; ++cv) {
        bool cvok = false;
        gmsh_cv_info &ci = cvlst[cv];
        bool is_node = (ci.type == 15);
        unsigned ci_dim = (is_node) ? 0 : ci.pgt->dim();
        //  cout << "importing cv dim=" << ci_dim << " N=" << N
        //       << " region: " << ci.region << " type: " << ci.type << "\n";

        //main convex import
        if (ci_dim == N) {
          size_type ic = m.add_convex(ci.pgt, ci.nodes.begin());
          cvok = true;
          m.region(ci.region).add(ic);

        //convexes with lower dimensions
        }
        else {
          //convex that lies within the regions of lower_dim_convex_rg
          //is imported explicitly as a convex.
          if (lower_dim_convex_rg != NULL &&
              lower_dim_convex_rg->find(ci.region) != lower_dim_convex_rg->end()
              && !is_node) {
              size_type ic = m.add_convex(ci.pgt, ci.nodes.begin());
              cvok = true; m.region(ci.region).add(ic);
          }
          //find if the convex is part of a face of higher dimension convex
          else{
            bgeot::mesh_structure::ind_cv_ct ct=m.convex_to_point(ci.nodes[0]);
            for (bgeot::mesh_structure::ind_cv_ct::const_iterator
                   it = ct.begin(); it != ct.end(); ++it) {
              if (m.structure_of_convex(*it)->dim() == ci_dim + 1) {
                for (short_type face=0;
                     face < m.structure_of_convex(*it)->nb_faces(); ++face) {
                  if (m.is_convex_face_having_points(*it, face,
                                                    short_type(ci.nodes.size()),
                                                    ci.nodes.begin())) {
                    m.region(ci.region).add(*it,face);
                    cvok = true;
                  }
                }
              }
            }
            if (is_node && (nodal_map != NULL)) {
              for (auto i : ci.nodes) (*nodal_map)[ci.region].insert(i);
            }
            // if the convex is not part of the face of others
            if (!cvok) {
              if (is_node) {
                if (nodal_map == NULL){
                  GMM_WARNING2("gmsh import ignored a node id: "
                               << ci.id << " region :" << ci.region <<
                               " point is not added explicitly as an element.");
                }
              }
              else if (add_all_element_type) {
                size_type ic = m.add_convex(ci.pgt, ci.nodes.begin());
                m.region(ci.region).add(ic);
                cvok = true;
              } else {
                GMM_WARNING2("gmsh import ignored an element of type "
                             << bgeot::name_of_geometric_trans(ci.pgt) <<
                    " as it does not belong to the face of another element");
              }
            }
          }
        }
      }
    }
    if (remove_last_dimension) maybe_remove_last_dimension(m);
  }

  /* mesh file from GiD [http://gid.cimne.upc.es/]

  supports linear and quadratic elements (quadrilaterals, use 9(or 27)-noded elements)
  */
  static void import_gid_mesh_file(std::istream& f, mesh& m) {
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
                      "Invalid node ID " << j << " in GiD element " << cv_id);
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

  /* mesh file from ANSYS

  Supports solid structural elements stored with cdwrite in blocked format.
  Use the following command in ANSYS for exporting the mesh:

  cdwrite,db,filename,cdb
  */
  static void import_cdb_mesh_file(std::istream& f, mesh& m,
                                   size_type imat_filt=size_type(-1)) {

    std::map<size_type, size_type> cdb_node_2_getfem_node;
    std::vector<size_type> getfem_cv_nodes;
    std::vector<std::string> elt_types;
    std::vector<size_type> elt_cnt;
    std::vector<dal::bit_vector> regions;

    size_type pos, pos2;
    std::string line;
    while (true) {
      std::getline(f,line);
      pos = line.find_first_not_of(" ");
      if (bgeot::casecmp(line.substr(pos,2),"ET") == 0) {
        size_type itype;
        std::string type_name;
        pos = line.find_first_of(",")+1;
        pos2 = line.find_first_of(",", pos);
        itype = std::stol(line.substr(pos, pos2-pos));
        pos = line.find_first_not_of(" ,\n\r\t", pos2);
        pos2 = line.find_first_of(" ,\n\r\t", pos);
        type_name = line.substr(pos, pos2-pos);
        bool only_digits
          = (type_name.find_first_not_of("0123456789") == std::string::npos);
        for (auto&& c : type_name) c = char(std::toupper(c));

        if (elt_types.size() < itype+1)
          elt_types.resize(itype+1);

        elt_types[itype] = "";
        if (only_digits) {
          size_type type_num = std::stol(type_name);
          if (type_num == 42 || type_num == 82 ||
              type_num == 182 || type_num == 183)
            elt_types[itype] = "PLANE";
          else if (type_num == 45 || type_num == 73 || type_num == 87 ||
                   type_num == 90 || type_num == 92 || type_num == 95 ||
                   type_num == 162 || type_num == 185 || type_num == 186 ||
                   type_num == 187 || type_num == 191)
            elt_types[itype] = "SOLID";
          else if (type_num == 89)
            elt_types[itype] = "VISCO";
        }
        elt_types[itype].append(type_name);
      }
      else if (bgeot::casecmp(line.substr(pos,5),"KEYOP") == 0) {
        size_type itype, knum, keyval;
        pos = line.find_first_of(",")+1;
        pos2 = line.find_first_of(",", pos);
        itype = std::stol(line.substr(pos, pos2-pos));
        pos = pos2+1;
        pos2 = line.find_first_of(",", pos);
        knum = std::stol(line.substr(pos, pos2-pos));
        keyval = std::stol(line.substr(pos2+1));
        if (knum == 1 && itype < elt_types.size() &&
            elt_types[itype].size() == 7 &&
            bgeot::casecmp(elt_types[itype].substr(0,7),"MESH200") == 0) {
          std::stringstream ss;
          ss << "MESH200_" << keyval;
          elt_types[itype] = ss.str();
        }
      }
      else if (bgeot::casecmp(line.substr(pos,6),"NBLOCK") == 0)
        break;
      else if (f.eof())
        return;
    }
    elt_cnt.resize(elt_types.size());

    //(3i8,6e20.13)
    size_type fields1, fieldwidth1, fields2, fieldwidth2; // 3,8,6,20
    { // "%8lu%*8u%*8u%20lf%20lf%20lf"
      std::string fortran_fmt; // "(%lu%*[i]%lu,%lu%*[e,E]%lu.%*u)"
      std::getline(f,fortran_fmt);
      pos = fortran_fmt.find_first_of("(")+1;
      pos2 = fortran_fmt.find_first_of("iI", pos);
      fields1 = std::stol(fortran_fmt.substr(pos, pos2-pos));
      pos = pos2+1;
      pos2 = fortran_fmt.find_first_of(",", pos);
      fieldwidth1 = std::stol(fortran_fmt.substr(pos, pos2-pos));
      pos = pos2+1;
      pos2 = fortran_fmt.find_first_of("eE", pos);
      fields2 = std::stol(fortran_fmt.substr(pos, pos2-pos));
      pos = pos2+1;
      pos2 = fortran_fmt.find_first_of(".", pos);
      fieldwidth2 = std::stol(fortran_fmt.substr(pos, pos2-pos));
      GMM_ASSERT1(fields1 >= 1 && fields2 >= 3 ,
                  "Ansys mesh import routine requires NBLOCK entries with at least "
                  "1 integer field and 3 float number fields");
    }

    base_node pt(3);
    for (size_type i=0; i < size_type(-1); ++i) {
      size_type nodeid;
      std::getline(f,line);
      if (line.compare(0,1,"N") == 0 || line.compare(0,1,"!") == 0)
        break;
      //       1       0       0-3.0000000000000E+00 2.0000000000000E+00 1.0000000000000E+00
      nodeid = std::stol(line.substr(0, fieldwidth1));
      pos = fields1*fieldwidth1;
      for (size_type j=0; j < 3; ++j, pos += fieldwidth2)
        if (pos < line.length())
          pt[j] = std::stod(line.substr(pos, fieldwidth2));
        else
          pt[j] = 0;

      cdb_node_2_getfem_node[nodeid] = m.add_point(pt, -1.);
    }

    while (bgeot::casecmp(line.substr(0,6),"EBLOCK") != 0) {
      if (f.eof())
        return;
      std::getline(f,line);
    }


    //(19i8)
    size_type fieldsno, fieldwidth; // 19,8
    { // "%8lu%8lu%8lu%8lu%8lu%8lu%8lu%8lu"
      std::string fortran_fmt;
      std::getline(f,fortran_fmt);

      pos = fortran_fmt.find_first_of("(")+1;
      pos2 = fortran_fmt.find_first_of("iI", pos);
      fieldsno = std::stol(fortran_fmt.substr(pos, pos2-pos));
      pos = pos2+1;
      pos2 = fortran_fmt.find_first_of(")\n", pos);
      fieldwidth = std::stol(fortran_fmt.substr(pos, pos2-pos));
      GMM_ASSERT1(fieldsno == 19, "Ansys mesh import routine requires EBLOCK "
                                  "entries with 19 fields");
    }

    size_type II,JJ,KK,LL,MM,NN,OO,PP,QQ,RR,SS,TT,UU,VV,WW,XX,YY,ZZ,AA,BB;
    for (size_type i=0; i < size_type(-1); ++i) {
      GMM_ASSERT1(!f.eof(), "File ended before all elements could be read");
      size_type imat, itype, nodesno(0);
      std::getline(f,line);
      {
        long int ii = std::stol(line.substr(0,fieldwidth));
        if (ii < 0)
          break;
        else
          imat = size_type(ii);
      }
      itype = std::stol(line.substr(fieldwidth,fieldwidth));
      nodesno = std::stol(line.substr(8*fieldwidth,fieldwidth));
      line = line.substr(11*fieldwidth);

      if (imat_filt != size_type(-1) && imat != imat_filt) { // skip current element
        if (nodesno > 8)
          std::getline(f,line);
        continue;
      }

      if (nodesno >= 1) II = std::stol(line.substr(0,fieldwidth));
      if (nodesno >= 2) JJ = std::stol(line.substr(1*fieldwidth,fieldwidth));
      if (nodesno >= 3) KK = std::stol(line.substr(2*fieldwidth,fieldwidth));
      if (nodesno >= 4) LL = std::stol(line.substr(3*fieldwidth,fieldwidth));
      if (nodesno >= 5) MM = std::stol(line.substr(4*fieldwidth,fieldwidth));
      if (nodesno >= 6) NN = std::stol(line.substr(5*fieldwidth,fieldwidth));
      if (nodesno >= 7) OO = std::stol(line.substr(6*fieldwidth,fieldwidth));
      if (nodesno >= 8) PP = std::stol(line.substr(7*fieldwidth,fieldwidth));
      if (nodesno >= 9) {
        std::getline(f,line);
        if (nodesno >= 9) QQ = std::stol(line.substr(0,fieldwidth));
        if (nodesno >= 10) RR = std::stol(line.substr(1*fieldwidth,fieldwidth));
        if (nodesno >= 11) SS = std::stol(line.substr(2*fieldwidth,fieldwidth));
        if (nodesno >= 12) TT = std::stol(line.substr(3*fieldwidth,fieldwidth));
        if (nodesno >= 13) UU = std::stol(line.substr(4*fieldwidth,fieldwidth));
        if (nodesno >= 14) VV = std::stol(line.substr(5*fieldwidth,fieldwidth));
        if (nodesno >= 15) WW = std::stol(line.substr(6*fieldwidth,fieldwidth));
        if (nodesno >= 16) XX = std::stol(line.substr(7*fieldwidth,fieldwidth));
        if (nodesno >= 17) YY = std::stol(line.substr(8*fieldwidth,fieldwidth));
        if (nodesno >= 18) ZZ = std::stol(line.substr(9*fieldwidth,fieldwidth));
        if (nodesno >= 19) AA = std::stol(line.substr(10*fieldwidth,fieldwidth));
        if (nodesno >= 20) BB = std::stol(line.substr(11*fieldwidth,fieldwidth));
      }

      if (imat+1 > regions.size())
        regions.resize(imat+1);

      if (nodesno == 3) {
        // assume MESH200_4 (3-node triangular)
        std::string eltname("MESH200_4");
        if (elt_types.size() > itype && elt_types[itype].size() > 0)
          eltname = elt_types[itype];

        if (eltname.compare("MESH200_4") == 0) {
          getfem_cv_nodes.resize(3);
          getfem_cv_nodes[0] = cdb_node_2_getfem_node[II];
          getfem_cv_nodes[1] = cdb_node_2_getfem_node[JJ];
          getfem_cv_nodes[2] = cdb_node_2_getfem_node[KK];
          regions[imat].add(m.add_convex(bgeot::simplex_geotrans(2,1),
                                         getfem_cv_nodes.begin()));
          if (itype < elt_cnt.size())
            elt_cnt[itype] += 1;
        } else {
          // TODO MESH200_2, MESH200_3, MESH200_4
          GMM_WARNING2("Ignoring ANSYS element " << eltname
                       << ". Import not supported yet.");
        }
      }
      else if (nodesno == 4) {

        // assume MESH200_6 (4-node quadrilateral)
        std::string eltname("MESH200_6");
        if (elt_types.size() > itype && elt_types[itype].size() > 0)
          eltname = elt_types[itype];

        if (eltname.compare("MESH200_6") == 0 ||
            eltname.compare("PLANE42") == 0 ||
            eltname.compare("PLANE182") == 0) {
          getfem_cv_nodes.resize(4);
          getfem_cv_nodes[0] = cdb_node_2_getfem_node[II];
          getfem_cv_nodes[1] = cdb_node_2_getfem_node[JJ];
          getfem_cv_nodes[2] = cdb_node_2_getfem_node[LL];
          getfem_cv_nodes[3] = cdb_node_2_getfem_node[KK];
          regions[imat].add(m.add_convex(bgeot::parallelepiped_geotrans(2,1),
                                         getfem_cv_nodes.begin()));
          if (itype < elt_cnt.size())
            elt_cnt[itype] += 1;
        }
        else if (eltname.compare("MESH200_8") == 0 ||
                 eltname.compare("SOLID72") == 0) {
          getfem_cv_nodes.resize(4);
          getfem_cv_nodes[0] = cdb_node_2_getfem_node[II];
          getfem_cv_nodes[1] = cdb_node_2_getfem_node[JJ];
          getfem_cv_nodes[2] = cdb_node_2_getfem_node[KK];
          getfem_cv_nodes[3] = cdb_node_2_getfem_node[LL];
          regions[imat].add(m.add_convex(bgeot::simplex_geotrans(3,1),
                                         getfem_cv_nodes.begin()));
          if (itype < elt_cnt.size())
            elt_cnt[itype] += 1;
        }
      }
      else if (nodesno == 6) {
        // assume MESH200_5 (6-node triangular)
        std::string eltname("MESH200_5");
        if (elt_types.size() > itype && elt_types[itype].size() > 0)
          eltname = elt_types[itype];
        if (eltname.compare("MESH200_5") == 0 ||
            eltname.compare("PLANE183") == 0) {
          getfem_cv_nodes.resize(6);
          getfem_cv_nodes[0] = cdb_node_2_getfem_node[II];
          getfem_cv_nodes[1] = cdb_node_2_getfem_node[LL];
          getfem_cv_nodes[2] = cdb_node_2_getfem_node[JJ];
          getfem_cv_nodes[3] = cdb_node_2_getfem_node[NN];
          getfem_cv_nodes[4] = cdb_node_2_getfem_node[MM];
          getfem_cv_nodes[5] = cdb_node_2_getfem_node[KK];
          regions[imat].add(m.add_convex(bgeot::simplex_geotrans(2,2),
                                         getfem_cv_nodes.begin()));
          if (itype < elt_cnt.size())
            elt_cnt[itype] += 1;
        }
      }
      else if (nodesno == 8) {

        // assume MESH200_10
        std::string eltname("MESH200_10");
        if (elt_types.size() > itype && elt_types[itype].size() > 0)
          eltname = elt_types[itype];

        if (eltname.compare("MESH200_7") == 0 ||
            eltname.compare("PLANE82") == 0 ||
            eltname.compare("PLANE183") == 0) {
          if (KK == LL && KK == OO) { // 6-node triangular
            getfem_cv_nodes.resize(6);
            getfem_cv_nodes[0] = cdb_node_2_getfem_node[II];
            getfem_cv_nodes[1] = cdb_node_2_getfem_node[MM];
            getfem_cv_nodes[2] = cdb_node_2_getfem_node[JJ];
            getfem_cv_nodes[3] = cdb_node_2_getfem_node[PP];
            getfem_cv_nodes[4] = cdb_node_2_getfem_node[NN];
            getfem_cv_nodes[5] = cdb_node_2_getfem_node[KK];
            regions[imat].add(m.add_convex(bgeot::simplex_geotrans(2,2),
                                           getfem_cv_nodes.begin()));
          } else {
            getfem_cv_nodes.resize(8);
            getfem_cv_nodes[0] = cdb_node_2_getfem_node[II];
            getfem_cv_nodes[1] = cdb_node_2_getfem_node[MM];
            getfem_cv_nodes[2] = cdb_node_2_getfem_node[JJ];
            getfem_cv_nodes[3] = cdb_node_2_getfem_node[PP];
            getfem_cv_nodes[4] = cdb_node_2_getfem_node[NN];
            getfem_cv_nodes[5] = cdb_node_2_getfem_node[LL];
            getfem_cv_nodes[6] = cdb_node_2_getfem_node[OO];
            getfem_cv_nodes[7] = cdb_node_2_getfem_node[KK];
            regions[imat].add(m.add_convex(bgeot::Q2_incomplete_geotrans(2),
                                           getfem_cv_nodes.begin()));
          }
          if (itype < elt_cnt.size())
            elt_cnt[itype] += 1;
        }
        else if (eltname.compare("MESH200_10") == 0 ||
                 eltname.compare("SOLID45") == 0 ||
                 eltname.compare("SOLID185") == 0) {
          if (KK == LL && OO == PP) {
            if (MM == NN && NN == OO) { // 4-node tetrahedral
              getfem_cv_nodes.resize(4);
              getfem_cv_nodes[0] = cdb_node_2_getfem_node[II];
              getfem_cv_nodes[1] = cdb_node_2_getfem_node[JJ];
              getfem_cv_nodes[2] = cdb_node_2_getfem_node[KK];
              getfem_cv_nodes[3] = cdb_node_2_getfem_node[MM];
              regions[imat].add(m.add_convex(bgeot::simplex_geotrans(3,1),
                                             getfem_cv_nodes.begin()));
            }
            else { // 6-node prism
              getfem_cv_nodes.resize(6);
              getfem_cv_nodes[0] = cdb_node_2_getfem_node[II];
              getfem_cv_nodes[1] = cdb_node_2_getfem_node[JJ];
              getfem_cv_nodes[2] = cdb_node_2_getfem_node[KK];
              getfem_cv_nodes[3] = cdb_node_2_getfem_node[MM];
              getfem_cv_nodes[4] = cdb_node_2_getfem_node[NN];
              getfem_cv_nodes[5] = cdb_node_2_getfem_node[OO];
              regions[imat].add(m.add_convex(bgeot::prism_geotrans(3,1),
                                             getfem_cv_nodes.begin()));
            }
          }
          else { // 8-node hexahedral
            getfem_cv_nodes.resize(8);
            getfem_cv_nodes[0] = cdb_node_2_getfem_node[II];
            getfem_cv_nodes[1] = cdb_node_2_getfem_node[JJ];
            getfem_cv_nodes[2] = cdb_node_2_getfem_node[LL];
            getfem_cv_nodes[3] = cdb_node_2_getfem_node[KK];
            getfem_cv_nodes[4] = cdb_node_2_getfem_node[MM];
            getfem_cv_nodes[5] = cdb_node_2_getfem_node[NN];
            getfem_cv_nodes[6] = cdb_node_2_getfem_node[PP];
            getfem_cv_nodes[7] = cdb_node_2_getfem_node[OO];
            regions[imat].add(m.add_convex(bgeot::parallelepiped_geotrans(3,1),
                                           getfem_cv_nodes.begin()));
          }
          if (itype < elt_cnt.size())
            elt_cnt[itype] += 1;
        }
      }
      else if (nodesno == 10) {

        // assume MESH200_9 (10-node tetrahedral)
        std::string eltname("MESH200_9");
        if (elt_types.size() > itype && elt_types[itype].size() > 0)
          eltname = elt_types[itype];

        if (eltname.compare("MESH200_9") == 0 ||
            eltname.compare("SOLID87") == 0 ||
            eltname.compare("SOLID92") == 0 ||
            eltname.compare("SOLID162") == 0 ||
            eltname.compare("SOLID187") == 0) {
          getfem_cv_nodes.resize(10);
          getfem_cv_nodes[0] = cdb_node_2_getfem_node[II];
          getfem_cv_nodes[1] = cdb_node_2_getfem_node[MM];
          getfem_cv_nodes[2] = cdb_node_2_getfem_node[JJ];
          getfem_cv_nodes[3] = cdb_node_2_getfem_node[OO];
          getfem_cv_nodes[4] = cdb_node_2_getfem_node[NN];
          getfem_cv_nodes[5] = cdb_node_2_getfem_node[KK];
          getfem_cv_nodes[6] = cdb_node_2_getfem_node[PP];
          getfem_cv_nodes[7] = cdb_node_2_getfem_node[QQ];
          getfem_cv_nodes[8] = cdb_node_2_getfem_node[RR];
          getfem_cv_nodes[9] = cdb_node_2_getfem_node[LL];
          regions[imat].add(m.add_convex(bgeot::simplex_geotrans(3,2),
                                         getfem_cv_nodes.begin()));
          if (itype < elt_cnt.size())
            elt_cnt[itype] += 1;
        }
      }
      else if (nodesno == 20) { //  # assume SOLID186/SOLID95

        // assume MESH200_11 (20-node hexahedral)
        std::string eltname("MESH200_11");
        if (elt_types.size() > itype && elt_types[itype].size() > 0)
          eltname = elt_types[itype];

        if (eltname.compare("MESH200_11") == 0 ||
            eltname.compare("VISCO89") == 0 ||
            eltname.compare("SOLID90") == 0 ||
            eltname.compare("SOLID95") == 0 ||
            eltname.compare("SOLID186") == 0 ||
            eltname.compare("SOLID191") == 0) {
          if (KK == LL && MM == NN && NN == OO && OO == PP) { // assume 10-node tetrahedral
            getfem_cv_nodes.resize(10);
            getfem_cv_nodes[0] = cdb_node_2_getfem_node[II];
            getfem_cv_nodes[1] = cdb_node_2_getfem_node[QQ];
            getfem_cv_nodes[2] = cdb_node_2_getfem_node[JJ];
            getfem_cv_nodes[3] = cdb_node_2_getfem_node[TT];
            getfem_cv_nodes[4] = cdb_node_2_getfem_node[RR];
            getfem_cv_nodes[5] = cdb_node_2_getfem_node[KK];
            getfem_cv_nodes[6] = cdb_node_2_getfem_node[YY];
            getfem_cv_nodes[7] = cdb_node_2_getfem_node[ZZ];
            getfem_cv_nodes[8] = cdb_node_2_getfem_node[AA];
            getfem_cv_nodes[9] = cdb_node_2_getfem_node[MM];
            regions[imat].add(m.add_convex(bgeot::simplex_geotrans(3,2),
                                           getfem_cv_nodes.begin()));
            if (itype < elt_cnt.size())
              elt_cnt[itype] += 1;
          } else if (MM == NN && NN == OO && OO == PP) { // assume 13-node pyramid
            getfem_cv_nodes.resize(13);
            getfem_cv_nodes[0] = cdb_node_2_getfem_node[II];
            getfem_cv_nodes[1] = cdb_node_2_getfem_node[QQ];
            getfem_cv_nodes[2] = cdb_node_2_getfem_node[JJ];
            getfem_cv_nodes[3] = cdb_node_2_getfem_node[TT];
            getfem_cv_nodes[4] = cdb_node_2_getfem_node[RR];
            getfem_cv_nodes[5] = cdb_node_2_getfem_node[LL];
            getfem_cv_nodes[6] = cdb_node_2_getfem_node[SS];
            getfem_cv_nodes[7] = cdb_node_2_getfem_node[KK];
            getfem_cv_nodes[8] = cdb_node_2_getfem_node[YY];
            getfem_cv_nodes[9] = cdb_node_2_getfem_node[ZZ];
            getfem_cv_nodes[10] = cdb_node_2_getfem_node[BB];
            getfem_cv_nodes[11] = cdb_node_2_getfem_node[AA];
            getfem_cv_nodes[12] = cdb_node_2_getfem_node[MM];
            regions[imat].add(m.add_convex(bgeot::pyramid_Q2_incomplete_geotrans(),
                                           getfem_cv_nodes.begin()));
            if (itype < elt_cnt.size())
              elt_cnt[itype] += 1;

          } else if (KK == LL && OO == PP) { // assume 15-node prism
            getfem_cv_nodes.resize(15);
            getfem_cv_nodes[0]  = cdb_node_2_getfem_node[II];
            getfem_cv_nodes[1]  = cdb_node_2_getfem_node[QQ];
            getfem_cv_nodes[2]  = cdb_node_2_getfem_node[JJ];
            getfem_cv_nodes[3]  = cdb_node_2_getfem_node[TT];
            getfem_cv_nodes[4]  = cdb_node_2_getfem_node[RR];
            getfem_cv_nodes[5]  = cdb_node_2_getfem_node[LL];
            getfem_cv_nodes[6]  = cdb_node_2_getfem_node[YY];
            getfem_cv_nodes[7]  = cdb_node_2_getfem_node[ZZ];
            getfem_cv_nodes[8]  = cdb_node_2_getfem_node[AA];
            getfem_cv_nodes[9]  = cdb_node_2_getfem_node[MM];
            getfem_cv_nodes[10] = cdb_node_2_getfem_node[UU];
            getfem_cv_nodes[11] = cdb_node_2_getfem_node[NN];
            getfem_cv_nodes[12] = cdb_node_2_getfem_node[XX];
            getfem_cv_nodes[13] = cdb_node_2_getfem_node[VV];
            getfem_cv_nodes[14] = cdb_node_2_getfem_node[OO];
            regions[imat].add(m.add_convex
                              (bgeot::prism_incomplete_P2_geotrans(),
                               getfem_cv_nodes.begin()));
            if (itype < elt_cnt.size())
              elt_cnt[itype] += 1;
          } else {
            getfem_cv_nodes.resize(20);
            getfem_cv_nodes[0] = cdb_node_2_getfem_node[II];
            getfem_cv_nodes[1] = cdb_node_2_getfem_node[QQ];
            getfem_cv_nodes[2] = cdb_node_2_getfem_node[JJ];
            getfem_cv_nodes[3] = cdb_node_2_getfem_node[TT];
            getfem_cv_nodes[4] = cdb_node_2_getfem_node[RR];
            getfem_cv_nodes[5] = cdb_node_2_getfem_node[LL];
            getfem_cv_nodes[6] = cdb_node_2_getfem_node[SS];
            getfem_cv_nodes[7] = cdb_node_2_getfem_node[KK];
            getfem_cv_nodes[8] = cdb_node_2_getfem_node[YY];
            getfem_cv_nodes[9] = cdb_node_2_getfem_node[ZZ];
            getfem_cv_nodes[10] = cdb_node_2_getfem_node[BB];
            getfem_cv_nodes[11] = cdb_node_2_getfem_node[AA];
            getfem_cv_nodes[12] = cdb_node_2_getfem_node[MM];
            getfem_cv_nodes[13] = cdb_node_2_getfem_node[UU];
            getfem_cv_nodes[14] = cdb_node_2_getfem_node[NN];
            getfem_cv_nodes[15] = cdb_node_2_getfem_node[XX];
            getfem_cv_nodes[16] = cdb_node_2_getfem_node[VV];
            getfem_cv_nodes[17] = cdb_node_2_getfem_node[PP];
            getfem_cv_nodes[18] = cdb_node_2_getfem_node[WW];
            getfem_cv_nodes[19] = cdb_node_2_getfem_node[OO];
            regions[imat].add(m.add_convex(bgeot::Q2_incomplete_geotrans(3),
                                           getfem_cv_nodes.begin()));
            if (itype < elt_cnt.size())
              elt_cnt[itype] += 1;
          }
        }
      }
    }

    int nonempty_regions=0;
    for (size_type i=0; i < regions.size(); ++i)
      if (regions[i].card() > 0)
        ++nonempty_regions;

    if (nonempty_regions > 1)
      for (size_type i=0; i < regions.size(); ++i)
        if (regions[i].card() > 0)
          m.region(i).add(regions[i]);

    for (size_type i=1; i < elt_types.size(); ++i)
      if (elt_cnt[i] > 0)
        cout << "Imported " << elt_cnt[i] << " " << elt_types[i] << " elements." << endl;
    cout << "Imported " << m.convex_index().card() << " elements in total." << endl;

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
  static void import_noboite_mesh_file(std::istream& f, mesh& m) {

    using namespace std;
    gmm::stream_standard_locale sl(f);

    ofstream fichier_GiD("noboite_to_GiD.gid",    ios::out | ios::trunc );  //déclaration du flux et ouverture du fichier

    fichier_GiD << "MESH    dimension 3 ElemType Tetrahedra  Nnode 4"<<endl;


    int i,NE,NP,ligne_debut_NP;

    /*NE: nombre d'elements (premier nombre du fichier .noboite)
      NP: nombre de points (deuxieme nombre du fichier .noboite)
      ligne_debut_NP: la ligne commence la liste des coordonnees des points dans le
                      fichier .noboite*/

    f >> NE>>NP;
    ligne_debut_NP=NE*4/6+3;

    //passer 3 premiers lignes du fichier .noboite (la liste des elements commence a la
    //quatrieme ligne)
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
        //f.close();  // on referme le fichier
      }
    else  // sinon
      cerr << "Erreur à l'ouverture !" << endl;

    // appeler sunroutine import_gid_mesh_file
    //import_mesh(const std::string& "noboite_to_GiD.gid", mesh& msh)
    ifstream fichier1_GiD("noboite_to_GiD.gid", ios::in);
    import_gid_mesh_file(fichier1_GiD, m);

    //      return 0;
  }

  /* mesh file from emc2 [http://pauillac.inria.fr/cdrom/prog/unix/emc2/eng.htm], am_fmt format

  (only triangular 2D meshes)
  */
  static void import_am_fmt_mesh_file(std::istream& f, mesh& m) {
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
  static void import_emc2_mesh_file(std::istream& f, mesh& m) {
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
      else if (bgeot::casecmp(format,"structured_ball")==0)
        { regular_ball_mesh(m, filename); return; }
      else if (bgeot::casecmp(format,"structured_ball_shell")==0)
        { regular_ball_shell_mesh(m, filename); return; }

      std::ifstream f(filename.c_str());
      GMM_ASSERT1(f.good(), "can't open file " << filename);
      /* throw exceptions when an error occurs */
      f.exceptions(std::ifstream::badbit | std::ifstream::failbit);
      import_mesh(f, format, m);
      f.close();
    }
    catch (std::logic_error& exc) {
      m.clear();
      throw exc;
    }
    catch (std::ios_base::failure& exc) {
      m.clear();
      GMM_ASSERT1(false, "error while importing " << format
                  << " mesh file \"" << filename << "\" : " << exc.what());
    }
    catch (std::runtime_error& exc) {
      m.clear();
      throw exc;
    }
  }

  void import_mesh_gmsh(std::istream& f, mesh &m,
                        std::map<std::string, size_type> &region_map,
                        bool remove_last_dimension,
                        std::map<size_type, std::set<size_type>> *nodal_map,
                        bool remove_duplicated_nodes)
  {
    import_gmsh_mesh_file(f, m, 0, &region_map, nullptr, false, remove_last_dimension, nodal_map,
                          remove_duplicated_nodes);
  }

  void import_mesh_gmsh(std::istream& f, mesh& m,
                        bool add_all_element_type,
                        std::set<size_type> *lower_dim_convex_rg,
                        std::map<std::string, size_type> *region_map,
                        bool remove_last_dimension,
                        std::map<size_type, std::set<size_type>> *nodal_map,
                        bool remove_duplicated_nodes)
  {
    import_gmsh_mesh_file(f, m, 0, region_map, lower_dim_convex_rg, add_all_element_type,
                          remove_last_dimension, nodal_map, remove_duplicated_nodes);
  }

  void import_mesh_gmsh(const std::string& filename, mesh& m,
                        bool add_all_element_type,
                        std::set<size_type> *lower_dim_convex_rg,
                        std::map<std::string, size_type> *region_map,
                        bool remove_last_dimension,
                        std::map<size_type, std::set<size_type>> *nodal_map,
                        bool remove_duplicated_nodes)
  {
    m.clear();
    try {
      std::ifstream f(filename.c_str());
      GMM_ASSERT1(f.good(), "can't open file " << filename);
      /* throw exceptions when an error occurs */
      f.exceptions(std::ifstream::badbit | std::ifstream::failbit);
      import_gmsh_mesh_file(f, m, 0, region_map, lower_dim_convex_rg, add_all_element_type,
                            remove_last_dimension, nodal_map, remove_duplicated_nodes);
      f.close();
    }
    catch (std::logic_error& exc) {
      m.clear();
      throw exc;
    }
    catch (std::ios_base::failure& exc) {
      m.clear();
      GMM_ASSERT1(false, "error while importing " << "gmsh"
                  << " mesh file \"" << filename << "\" : " << exc.what());
    }
    catch (std::runtime_error& exc) {
      m.clear();
      throw exc;
    }
  }

  void import_mesh_gmsh(const std::string& filename,
    mesh& m, std::map<std::string, size_type> &region_map,
    bool remove_last_dimension,
    std::map<size_type, std::set<size_type>> *nodal_map,
    bool remove_duplicated_nodes) {
    import_mesh_gmsh(filename, m, false, NULL, &region_map, remove_last_dimension, nodal_map,
                     remove_duplicated_nodes);
  }

  void import_mesh(std::istream& f, const std::string& format,
                   mesh& m) {
    if (bgeot::casecmp(format,"gmsh")==0)
      import_gmsh_mesh_file(f,m);
    else if (bgeot::casecmp(format,"gmsh_with_lower_dim_elt")==0)
      import_gmsh_mesh_file(f,m,0,NULL,NULL,true);
    else if (bgeot::casecmp(format,"gmshv2")==0)/* deprecate */
      import_gmsh_mesh_file(f,m,2);
    else if (bgeot::casecmp(format,"gid")==0)
      import_gid_mesh_file(f,m);
    else if (bgeot::casecmp(format,"noboite")==0)
      import_noboite_mesh_file(f,m);
    else if (bgeot::casecmp(format,"am_fmt")==0)
      import_am_fmt_mesh_file(f,m);
    else if (bgeot::casecmp(format,"emc2_mesh")==0)
      import_emc2_mesh_file(f,m);
    else if (bgeot::casecmp(format,"cdb")==0)
      import_cdb_mesh_file(f,m);
    else if (bgeot::casecmp(format.substr(0,4),"cdb:")==0) {
      size_type imat(-1);
      bool success(true);
      try {
        size_t sz;
        imat = std::stol(format.substr(4), &sz);
        success = (sz == format.substr(4).size() && imat != size_type(-1));
      } catch (const std::invalid_argument&) {
        success = false;
      } catch (const std::out_of_range&) {
        success = false;
      }
      if (success)
        import_cdb_mesh_file(f,m,imat);
      else GMM_ASSERT1(false, "cannot import "
                       << format << " mesh type : wrong cdb mesh type input");
    }
    else GMM_ASSERT1(false, "cannot import "
                     << format << " mesh type : unknown mesh type");
  }

  void import_mesh(const std::string& filename, mesh& msh) {
    size_type pos = filename.find_last_of(":");
    if (pos != std::string::npos)
      getfem::import_mesh(filename.substr(pos+1), filename.substr(0,pos), msh);
    else
      msh.read_from_file(filename);
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
