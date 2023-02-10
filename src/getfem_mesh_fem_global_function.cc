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


// examples of bspline basis functions assigned to 8 elements in 1D

// n=8,k=3, free-free  --> n-k+1 + 2*(k-1) = n+k-1 = 8+3-1 = 10
//   1 2 3 4 5 6 7 8 |
//1  *               |
//2  * *             |
//3  * * *           |
//4    * * *         |
//5      * * *       |
//6        * * *     |
//7          * * *   |
//8            * * * |
//9              * * |
//10               * |

// n=8,k=4, free-free  --> n-k+1 + 2*(k-1) = n+k-1 = 8+4-1 = 11
//   1 2 3 4 5 6 7 8
//1  *               |
//2  * *             |
//3  * * *           |
//4  * * * *         |
//5    * * * *       |
//6      * * * *     |
//7        * * * *   |
//8          * * * * |
//9            * * * |
//10             * * |
//11               * |

// n=8,k=3, periodic  --> n-k+1 + k-1 = n
//   1 2 3 4 5 6 7 8
//1  * * *           |
//2    * * *         |
//3      * * *       |
//4        * * *     |
//5          * * *   |
//6            * * * |
//7  *           * * |
//8  * *           * |

// n=8,k=4, periodic
//   1 2 3 4 5 6 7 8
//1  * * * *         |
//2    * * * *       |
//3      * * * *     |
//4        * * * *   |
//5          * * * * |
//6  *         * * * |
//7  * *         * * |
//8  * * *         * |

// n=8,k=3, symmetry-symmetry  --> n-k+1 + 2*floor(k/2) = 6 + 2 = 8
//   1 2 3 4 5 6 7 8
//1  + *             |
//2  * * *           |
//3    * * *         |
//4      * * *       |
//5        * * *     |
//6          * * *   |
//7            * * * |
//8              * + |

// n=8,k=4, symmetry-symmetry  --> n-k+1 + 2*floor(k/2) = 5 + 4 = 9
//   1 2 3 4 5 6 7 8 |
//1  * *             |
//2  + * *           |
//3  * * * *         |
//4    * * * *       |
//5      * * * *     |
//6        * * * *   |
//7          * * * * |
//8            * * + |
//9              * * |

// n=8,k=5, symmetry-symmetry  --> n-k+1 + 2*floor(k/2) = 4 + 4 = 8
//   1 2 3 4 5 6 7 8 |
//1  + + *           |
//2  + * * *         |
//3  * * * * *       |
//4    * * * * *     |
//5      * * * * *   |
//6        * * * * * |
//7          * * * + |
//8            * + + |

// n=8,k=6, symmetry-symmetry  --> n-k+1 + 2*floor(k/2) = 3 + 6 = 9
//   1 2 3 4 5 6 7 8 |
//1  * * *           |
//2  + + * *         |
//3  + * * * *       |
//4  * * * * * *     |
//5    * * * * * *   |
//6      * * * * * * |
//7        * * * * + |
//8          * * + + |
//9            * * * |

// n=8,k=3, free-symmetry  --> n-k+1 + k-1 + floor(k/2) = 6 + 2 + 1 = 9
//   1 2 3 4 5 6 7 8 |
//1  *               |
//2  + *             |
//3  * * *           |
//4    * * *         |
//5      * * *       |
//6        * * *     |
//7          * * *   |
//8            * * * |
//9              * + |

  void params_for_uniform_1d_bspline_basis_functions
  (scalar_type x0, scalar_type x1, size_type N, size_type order,
   bspline_boundary bc_low, bspline_boundary bc_high,
   std::vector<scalar_type> &xmin, std::vector<scalar_type> &xmax,
   std::vector<scalar_type> &xshift, std::vector<size_type> &xtype) {

    if (bc_low == bspline_boundary::PERIODIC ||
        bc_high == bspline_boundary::PERIODIC)
      GMM_ASSERT1(bc_low == bc_high,
                  "Periodic BC has to be assigned to both matching sides");
    const scalar_type dx = (x1-x0)/scalar_type(N);
    size_type n_low, n_mid, n_high;
    n_low = (bc_low == bspline_boundary::PERIODIC) ? 0 :
            (bc_low == bspline_boundary::SYMMETRY  ? order/2 :
                                      /* FREE */     order-1);
    n_high = (bc_high == bspline_boundary::PERIODIC) ? order-1 :
             (bc_high == bspline_boundary::SYMMETRY  ? order/2 :
                                        /* FREE */     order-1);
    n_mid = N - order + 1;
    size_type n = n_low + n_mid + n_high; // number of basis functions

    xmin.resize(n);
    xmax.resize(n);
    xshift.resize(n);
    xtype.resize(n);
    for (size_type i=0; i < n; ++i) {
      xshift[i] = 0.;
      if (bc_low == bspline_boundary::FREE && i < n_low) {
        xtype[i] = i+1;
        xmin[i] = x0;
        xmax[i] = xmin[i] + scalar_type(xtype[i])*dx;
      } else if (bc_high == bspline_boundary::FREE && i >= n_low+n_mid) {
        xtype[i] = n-i; // safe unsigned
        xmin[i] = x1;
        xmax[i] = xmin[i] - scalar_type(xtype[i])*dx; // yes, xmax < xmin
      } else if (bc_low == bspline_boundary::SYMMETRY && i < n_low) {
        xtype[i] = order;
        xmin[i] = x0 - scalar_type(n_low-i)*dx;
        xmax[i] = xmin[i] + scalar_type(xtype[i])*dx;
        xshift[i] = -(xmin[i]+xmax[i]-2*x0); // this is 0 for already symmetric basis functions
      } else if (bc_high == bspline_boundary::SYMMETRY && i >= n_low+n_mid) {
        xtype[i] = order;
        xmin[i] = x0 + scalar_type(i-n_low)*dx; // safe unsigned
        xmax[i] = xmin[i] + scalar_type(xtype[i])*dx;
        xshift[i] = 2*x1-xmin[i]-xmax[i]; // this is 0 for already symmetric basis functions
      } else { // mid functions for periodic, free-free or free-symmetry or symmetry-free
        GMM_ASSERT1(i >= n_low, "Internal error");
        xtype[i] = order;
        xmin[i] = x0 + scalar_type(i-n_low)*dx; // safe unsigned
        xmax[i] = xmin[i] + scalar_type(xtype[i])*dx;
      }
//if (order==5) // && bc_low == bspline_boundary::SYMMETRY && bc_high == bspline_boundary::FREE)
//std::cout<<i<<":"<<xmin[i]<<","<<xmax[i]<<std::endl;

      if (bc_low == bspline_boundary::PERIODIC && xmax[i] > x1)
        xshift[i] = -(x1-x0); // this will apply to the last order-1 functions
    }
  }

  void define_uniform_bspline_basis_functions_for_mesh_fem
  (mesh_fem_global_function &mf, size_type NX, size_type order,
   bspline_boundary bcX_low, bspline_boundary bcX_high,
   const mesh_im &mim) {

    GMM_ASSERT1(mf.linked_mesh().dim() == 1,
                "This function expects a mesh_fem defined in 1d");

    base_node Pmin, Pmax;
    mf.linked_mesh().bounding_box(Pmin, Pmax);
    const scalar_type x0=Pmin[0], x1=Pmax[0];

    std::vector<scalar_type> xmin, xmax, xshift;
    std::vector<size_type> xtype;
    params_for_uniform_1d_bspline_basis_functions
      (x0, x1, NX, order, bcX_low, bcX_high, // input
       xmin, xmax, xshift, xtype);           // output

    std::vector<pglobal_function> funcs(0);
    for (size_type i=0; i < xtype.size(); ++i) {
      if (gmm::abs(xshift[i]) < 1e-10)
        funcs.push_back(global_function_bspline
                        (xmin[i], xmax[i], order, xtype[i]));
      else {
        std::vector<pglobal_function> sum;
        sum.push_back(global_function_bspline
                      (xmin[i], xmax[i], order, xtype[i]));
        sum.push_back(global_function_bspline
                      (xmin[i]+xshift[i], xmax[i]+xshift[i],
                       order, xtype[i]));
        funcs.push_back(std::make_shared<getfem::global_function_sum>(sum));
      }
    }
    mf.set_functions(funcs, mim);
  }

  void define_uniform_bspline_basis_functions_for_mesh_fem
  (mesh_fem_global_function &mf,
   size_type NX, size_type NY, size_type order,
   bspline_boundary bcX_low, bspline_boundary bcY_low,
   bspline_boundary bcX_high, bspline_boundary bcY_high,
   const mesh_im &mim) {

    GMM_ASSERT1(mf.linked_mesh().dim() == 2,
                "This function expects a mesh_fem defined in 2d");

    base_node Pmin, Pmax;
    mf.linked_mesh().bounding_box(Pmin, Pmax);
    const scalar_type x0=Pmin[0], x1=Pmax[0],
                      y0=Pmin[1], y1=Pmax[1];

    std::vector<scalar_type> xmin, xmax, xshift;
    std::vector<size_type> xtype;
    params_for_uniform_1d_bspline_basis_functions
      (x0, x1, NX, order, bcX_low, bcX_high, // input
       xmin, xmax, xshift, xtype);           // output
    std::vector<scalar_type> ymin, ymax, yshift;
    std::vector<size_type> ytype;
    params_for_uniform_1d_bspline_basis_functions
      (y0, y1, NY, order, bcY_low, bcY_high, // input
       ymin, ymax, yshift, ytype);           // output

    std::vector<pglobal_function> funcs(0);
    for (size_type i=0; i < xtype.size(); ++i) {
      for (size_type j=0; j < ytype.size(); ++j) {
        if (gmm::abs(xshift[i]) < 1e-10 &&
            gmm::abs(yshift[j]) < 1e-10)
          funcs.push_back(global_function_bspline
                          (xmin[i], xmax[i], ymin[j], ymax[j],
                           order, xtype[i], ytype[j]));
        else {
          std::vector<pglobal_function> sum;
          sum.push_back(global_function_bspline
                        (xmin[i], xmax[i], ymin[j], ymax[j],
                         order, xtype[i], ytype[j]));
          if (gmm::abs(xshift[i]) >= 1e-10)
            sum.push_back(global_function_bspline
                          (xmin[i]+xshift[i], xmax[i]+xshift[i],
                           ymin[j], ymax[j],
                           order, xtype[i], ytype[j]));
          if (gmm::abs(yshift[j]) >= 1e-10) {
            sum.push_back(global_function_bspline
                          (xmin[i], xmax[i],
                           ymin[j]+yshift[j], ymax[j]+yshift[j],
                           order, xtype[i], ytype[j]));
            if (gmm::abs(xshift[i]) >= 1e-10)
              sum.push_back(global_function_bspline
                            (xmin[i]+xshift[i], xmax[i]+xshift[i],
                             ymin[j]+yshift[j], ymax[j]+yshift[j],
                             order, xtype[i], ytype[j]));
          }
          funcs.push_back(std::make_shared<getfem::global_function_sum>(sum));
        }
      }
    }
    mf.set_functions(funcs, mim);
  }

}

/* end of namespace getfem  */
