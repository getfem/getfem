/*===========================================================================

 Copyright (C) 2004-2016 Yves Renard
 Copyright (C) 2016      Konstantinos Poulios

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

#include <getfem/getfem_mesh_fem_global_function.h>

namespace getfem {

  void mesh_fem_global_function::set_functions
  (pglobal_function f, const mesh_im &mim) { 
    std::vector<pglobal_function> funcs(1);
    funcs[0]=f;
    set_functions(funcs, mim);
  }

  void mesh_fem_global_function::set_functions
  (pglobal_function f1, pglobal_function f2, const mesh_im &mim) { 
    std::vector<pglobal_function> funcs(2);
    funcs[0]=f1;
    funcs[1] = f2;
    set_functions(funcs, mim);
  }

  void mesh_fem_global_function::set_functions
  (const std::vector<pglobal_function>& funcs, const mesh_im &mim) {
    GMM_ASSERT1(linked_mesh_ != 0, "Mesh fem need to be initialized with"
                                   " a mesh first.");
    clear();
    if (&mim == &dummy_mesh_im())
      fem_ = getfem::new_fem_global_function(funcs, linked_mesh());
    else {
      GMM_ASSERT1(&(mim.linked_mesh()) == linked_mesh_,
                  "The provided mesh_im has to be linked to the same mesh"
                  " as this mesh_fem.");
      fem_ = getfem::new_fem_global_function(funcs, mim);
    }
    set_finite_element(fem_);
  }

  // mesh_fem_global_function::size_type memsize() const;

  void mesh_fem_global_function::clear() {
    mesh_fem::clear();
    if (fem_.get()) {
      getfem::del_fem_global_function(fem_);
      fem_.reset();
    }
  }

}

/* end of namespace getfem  */
