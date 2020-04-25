/*===========================================================================

 Copyright (C) 2020-2020 Tetsuo Koyama

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
#include "getfem/getfem_export.h"

using bgeot::base_node;

int main(void) {

  getfem::mesh m0;

  base_node org(1);
  org[0] = 0.0;

  bgeot::base_small_vector h(1);
  h[0] = 1;
  std::vector<bgeot::base_small_vector> vect(1);
  vect[0] = bgeot::base_small_vector(h);

  std::vector<int> ref(1);
  ref[0] = 1;

  getfem::parallelepiped_regular_simplex_mesh(m0, 1, org, vect.begin(), ref.begin());

  getfem::vtk_export vtk_exp("m0.vtk", true);
  vtk_exp.exporting(m0);
  vtk_exp.write_mesh();

  getfem::vtu_export vtu_exp("m0.vtu", true);
  vtu_exp.exporting(m0);
  vtu_exp.write_mesh();

  return 0;
}

