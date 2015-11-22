/*===========================================================================

 Copyright (C) 2007-2015 Yves Renard, Julien Pommier.

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

#include "getfem/getfem_mesher.h"

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;
using getfem::base_node;
using getfem::scalar_type;

int main(int argc, char **argv) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  try {
    getfem::mesh m;
    getfem::scalar_type h = .3;
    int K=1;
    int opt = 0;
    int max_iter = 400;
    int prefind = 1;
    if (argc > 1) { opt = atoi(argv[1]); }
    if (argc > 2) { h = atof(argv[2]); cout << "h = " << h << "\n"; }
    if (argc > 3) { K = atoi(argv[3]); }
    if (argc > 4) { max_iter = atoi(argv[4]); }
    if (argc > 5) { prefind = atoi(argv[5]); }
    assert(K>0); assert(h>0.); 
    std::vector<getfem::base_node> fixed;

    getfem::pmesher_signed_distance
      D0  = getfem::new_mesher_ball(base_node(0.,0.),1.),
      B1  = getfem::new_mesher_ball(base_node(1.5,0.),2.),     
      B2  = getfem::new_mesher_ball(base_node(-1.5,0.),2.),
      D2  = getfem::new_mesher_union(B1, B2),
      D3  = getfem::new_mesher_rectangle(base_node(0.,0.,0.),
					  base_node(1.,1.,1.)),
      D4  = getfem::new_mesher_cylinder(base_node(0.,0.,0.),
					base_node(0.,1.,1.), 2.0, 1.0),
      D5  = getfem::new_mesher_ball(base_node(0.,0.,0.),1.),
      B3  = getfem::new_mesher_ball(base_node(1.5,0.,0.),2.),
      B4  = getfem::new_mesher_ball(base_node(-1.5,0.,0.),2.),
      D6  = getfem::new_mesher_union(B3, B4),
      B5  = getfem::new_mesher_ball(base_node(0.,1.,0.),2.),
      D7  = getfem::new_mesher_union(D6, B5),
      B6  = getfem::new_mesher_ball(base_node(0.,1.),2.),
      D8  = getfem::new_mesher_union(D2, B6),
      P1  = getfem::new_mesher_half_space(base_node(0,.1,.1),
					  base_node(.1,-.03,1.)),
      D9  = getfem::new_mesher_intersection(D7, P1),
      D10 = getfem::new_mesher_ball(base_node(0.,20.), 20.),
      P2  = getfem::new_mesher_half_space(base_node(0.,.1,.1), 
					  base_node(.1,-.03,.5)),
      D11 = getfem::new_mesher_intersection(D5, P2),
      B7  = getfem::new_mesher_ball(base_node(-1.,0.,0.),1.5),
      B8  = getfem::new_mesher_ball(base_node(1.,0.,0.),1.5),
      I1  = getfem::new_mesher_intersection(B7, B8),
      B9  = getfem::new_mesher_ball(base_node(0,-1.,.5),1.5),
      B10 = getfem::new_mesher_ball(base_node(0.,1.,.5),1.5),
      I2  = getfem::new_mesher_intersection(B9, B10),
      D12 = getfem::new_mesher_union(I1, I2),
      B11 = getfem::new_mesher_ball(base_node(0.,0.),2.),
      B12 = getfem::new_mesher_ball(base_node(1.,0.),2),
      D13 = getfem::new_mesher_setminus(B11, B12),
      B13 = getfem::new_mesher_ball(base_node(0.,0.,0.),2.),
      B14 = getfem::new_mesher_ball(base_node(1.,0.,0.),2),
      D14 = getfem::new_mesher_setminus(B13, B14),
      D15 = getfem::new_mesher_torus(2.0, 0.5),
      C1  = getfem::new_mesher_cylinder(base_node(0.5,0.5,-0.5),
					base_node(0.,0.,1.), 2.0, 0.2),
      C2  = getfem::new_mesher_cylinder(base_node(0.5,0.5,-0.5),
					base_node(0.,0.,1.), 2.0, 0.65),
      S1  = getfem::new_mesher_setminus(D3, C1),
      D16 = getfem::new_mesher_intersection(S1, C2),
      B15 = getfem::new_mesher_ball(base_node(0.,0.,0.),1.2),
      B16 = getfem::new_mesher_ball(base_node(-2.,-2.,0.),0.7),
      B17 = getfem::new_mesher_ball(base_node(2.,-2.,0.),0.7),
      B18 = getfem::new_mesher_ball(base_node(2.,2.,0.),0.7),
      B19 = getfem::new_mesher_ball(base_node(-2.,2.,0.),0.7),
      C3  = getfem::new_mesher_cylinder(base_node(-2.,-2.,0.),
					base_node(1.,0.,0.), 4.0, 0.3),
      C4  = getfem::new_mesher_cylinder(base_node(2.,-2.,0.),
					base_node(0.,1.,0.), 4.0, 0.3),
      C5  = getfem::new_mesher_cylinder(base_node(2.,2.,0.),
					base_node(-1.,0.,0.), 4.0, 0.3),
      C6  = getfem::new_mesher_cylinder(base_node(-2.,2.,0.),
					base_node(0.,-1.,0.), 4.0, 0.3),
      C7  = getfem::new_mesher_cylinder(base_node(0.,-2.,0.),
					base_node(0.,1.,0.), 4.0, 0.3),
      C8  = getfem::new_mesher_cylinder(base_node(-2.,0.,0.),
					base_node(1.,0.,0.), 4.0, 0.3),
      C9  = getfem::new_mesher_cylinder(base_node(0.,0.,-1.5),
					base_node(0.,0.,1.), 3.0, 0.3),
      U1  = getfem::new_mesher_union(B15, B16, B17, B18, B19,
				     C3, C4, C5, C6, C7, C8),
      D17 = getfem::new_mesher_setminus(U1, C9),
      D18 = getfem::new_mesher_rectangle(base_node(0.,0.,0.),
					 base_node(.1,.3,1.)),
      C10 = getfem::new_mesher_cylinder(base_node(0.5,0.5,1.2),
					base_node(0.,0.,-1.), 3.0, 0.1),
      C11 = getfem::new_mesher_cylinder(base_node(0.5,0.5,0.),
					base_node(0.,0.,-1.), 1.0, 0.3),
      P3  = getfem::new_mesher_rectangle(base_node(0.,0.,0.),
					 base_node(1.,1.,1.)),
      P4  = getfem::new_mesher_rectangle(base_node(0.3,0.3,1.0),
					base_node(0.6,0.6,1.2)),
      U2  = getfem::new_mesher_union(P3, C11),
      M1  = getfem::new_mesher_setminus(U2, C10),
      D19 = getfem::new_mesher_setminus(M1, P4),
      D22 = getfem::new_mesher_ball(base_node(0.,0.,20.),20.),
      D23_1 = getfem::new_mesher_rectangle(base_node(0.,1.),
					   base_node(20.0, 10.0)),
      D23_2 = getfem::new_mesher_rectangle(base_node(1.,0.),
					   base_node(19.0, 1.0)),
      D23_3 = getfem::new_mesher_ball(base_node(1.0,1.0),1.0),
      D23_4 = getfem::new_mesher_ball(base_node(19.0,1.0),1.0),
      D23 = getfem::new_mesher_union(D23_1, D23_2, D23_3, D23_4);


    getfem::pmesher_signed_distance dist = D0;
    switch (opt) {
    case 0: dist = D0; break; /* disk */
    case 1: {
      getfem::scalar_type z = sqrt(4 - 1.5*1.5);
      fixed.push_back(base_node(0,z));
      fixed.push_back(base_node(0,-z));
    }
    case 2: dist = D2; break; /* union of 2 disks */
    case 3: dist = D3; break; /* cube */
    case 4: dist = D4; break; /* cylinder */
    case 5: dist = D5; break; /* ball */
    case 6: dist = D6; break; /* union of 2 balls */
    case 7: dist = D7; break; /* union of 3 balls */
    case 8: dist = D8; break; /* union of 3 disks */
    case 9: dist = D9; break; /* union of 3 half balls */
    case 10: dist = D10;      /* disk r=20 with fixed points */
      fixed.push_back(base_node(0.,0.));
      fixed.push_back(base_node(0.,20.));
      fixed.push_back(base_node(-20.,20.));
      fixed.push_back(base_node(20.,20.));
      fixed.push_back(base_node(0.,40.));
      dist = D10; 
      break;
    case 11: dist = D11; break; /* half-balls */
    case 12: dist = D12; break; /* UFO */
    case 13: dist = D13; break; /* moon */
    case 14: dist = D14; break; /* subtraction of two balls */
    case 15: dist = D15; break; /* torus */
    case 16: dist = D16; break; /* cube with a hole */
    case 17: dist = D17; break; /* space station */
    case 18: dist = D18; break;
    case 19: dist = D19; break;
    case 20: {                   /* ladder */
      scalar_type H = 20; 
      unsigned nb_step = 3;
      dist = getfem::new_mesher_union
	(getfem::new_mesher_rectangle(base_node(-3, -1.2, 0),
				      base_node(-1.8, 1.2, H)), 
	 getfem::new_mesher_rectangle(base_node(1.8, -1.2, 0),
				      base_node(3, 1.2, H)));
      for (unsigned i=0; i < nb_step; ++i) {
	scalar_type z = H/nb_step/2 + i*(H/nb_step);
	getfem::pmesher_signed_distance d = 
	  getfem::new_mesher_cylinder(base_node(-2.4, 0, z),
				      base_node(1, 0, 0), 4.8, .8);
	dist = getfem::new_mesher_union(dist, d);
      }
    } break;
    case 21: { /* bar with holes */
      scalar_type H=20;
      unsigned nb_holes = 4;
      dist = getfem::new_mesher_rectangle(base_node(-2, -1, 0), 
					  base_node(+2, +1, H));
      for (unsigned i=0; i < nb_holes; ++i) {
	scalar_type z = H/nb_holes/2 + i*(H/nb_holes);
	getfem::pmesher_signed_distance d = 
	  getfem::new_mesher_cylinder(base_node(0, -20, z), 
				      base_node(0, 1, 0), 40, .7);
	dist = getfem::new_mesher_setminus(dist, d);
      }
    } break;
    case 22: dist = D22; break; /* ball for test of radius 20. */
    case 23:                   /* rectangle with rounded corners */
      fixed.push_back(base_node(0.,1.));
      fixed.push_back(base_node(1.,0.));
      fixed.push_back(base_node(19.,0.));
      fixed.push_back(base_node(20.,1.));
      dist = D23;
      break;
    }
    getfem::build_mesh(m, dist, h, fixed, K, 2, max_iter, prefind);
    cout << "You can view the result with"
	 << "\n mayavi -d totoq.vtk -m BandedSurfaceMap\n";
  }
  GMM_STANDARD_CATCH_ERROR;
  return 0;
}
