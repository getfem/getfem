/*===========================================================================

 Copyright (C) 2004-2022 Yves Renard
 Copyright (C) 2016-2022 Konstantinos Poulios

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

#include <getfem/getfem_mesh_fem_global_function.h>

namespace getfem {


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


  void define_bspline_basis_functions_for_mesh_fem
  (mesh_fem_global_function &mf,
   size_type NX, size_type NY, size_type order,
   const mesh_im &mim) {

    base_node Pmin, Pmax;
    mf.linked_mesh().bounding_box(Pmin, Pmax);
    scalar_type x0=Pmin[0], y0=Pmin[1], x1=Pmax[0], y1=Pmax[1];

    scalar_type xmin, xmax, ymin, ymax;
    size_type xtype, ytype;
    base_node pt(2);

    std::vector<pglobal_function> funcs((NX+order-1)*(NY+order-1));
    for (size_type i=0; i < NX+order-1; ++i) {
      if (i < order-1) {
        xmin = x0;
        xmax = x0+scalar_type(i+1)*(x1-x0)/scalar_type(NX);
        xtype = i+1;
        pt[0] = (i == 0) ? xmin : (xmin+(xmax-xmin)/3);
      } else if (i >= NX) {
        xmin = x1;
        xmax = x1+(scalar_type(i-NX)-scalar_type(order-1))*(x1-x0)/scalar_type(NX);
        xtype = NX-i+order-1;
        pt[0] = (i == NX+1) ? xmin : (xmin+(xmax-xmin)/3);
      } else {
        xmin = x0+scalar_type(i-order+1)*(x1-x0)/scalar_type(NX);
        xmax = x0+scalar_type(i+1)*(x1-x0)/scalar_type(NX);
        xtype = order;
        pt[0] = (xmin+xmax)/2;
      }
      for (size_type j=0; j < NY+order-1; ++j) {
        if (j < order-1) {
          ymin = y0;
          ymax = y0+scalar_type(j+1)*(y1-y0)/scalar_type(NY);
          ytype = j+1;
          pt[1] = (j == 0) ? ymin : (ymin+(ymax-ymin)/3);
        } else if (j >= NY) {
          ymin = y1;
          ymax = y1+(scalar_type(j-NY)-scalar_type(order-1))*(y1-y0)/scalar_type(NY);
          ytype = NY-j+order-1;
          pt[1] = (j == NY+1) ? ymin : (ymin+(ymax-ymin)/3);
        } else {
          ymin = y0+scalar_type(j-order+1)*(y1-y0)/scalar_type(NY);
          ymax = y0+scalar_type(j+1)*(y1-y0)/scalar_type(NY);
          ytype = order;
          pt[1] = (ymin+ymax)/2;
        }
        funcs[i*(NY+order-1)+j] = global_function_bspline
                                  (xmin, xmax, ymin, ymax, order, xtype, ytype);
      }
    }

    mf.set_functions(funcs, mim);
  }

}

/* end of namespace getfem  */
