/*===========================================================================
 
 Copyright (C) 2007-2012 Yves Renard, Julien Pommier.
 
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

    getfem::mesher_ball D0(base_node(0.,0.),1.);

    getfem::mesher_ball B1(base_node(1.5,0.),2.);
    getfem::mesher_ball B2(base_node(-1.5,0.),2.);
    getfem::mesher_union D2(B1, B2);

    getfem::mesher_rectangle D3(base_node(0.,0.,0.), base_node(1.,1.,1.));

    getfem::mesher_cylinder D4(base_node(0.,0.,0.), base_node(0.,1.,1.),
			       2.0, 1.0);

    getfem::mesher_ball D5(base_node(0.,0.,0.),1.);

    getfem::mesher_ball B3(base_node(1.5,0.,0.),2.);
    getfem::mesher_ball B4(base_node(-1.5,0.,0.),2.);
    getfem::mesher_union D6(B3, B4);

    getfem::mesher_ball B5(base_node(0,1.,0.),2.);
    getfem::mesher_union D7(D6, B5);

    getfem::mesher_ball B6(getfem::base_node(0,1.),2.);
    getfem::mesher_union D8(D2, B6);

    getfem::mesher_half_space P1(getfem::base_node(0,.1,.1),
				 getfem::base_node(.1,-.03,1.));
    getfem::mesher_intersection D9(D7, P1);

    getfem::mesher_ball D10(getfem::base_node(0.,20.),20.);

    getfem::mesher_half_space P2(getfem::base_node(0,.1,.1), 
				 getfem::base_node(.1,-.03,.5));
    getfem::mesher_intersection D11(D5, P2);

    getfem::mesher_ball B7(getfem::base_node(-1.,0.,0.),1.5);
    getfem::mesher_ball B8(getfem::base_node(1.,0.,0.),1.5);
    getfem::mesher_intersection I1(B7, B8);
    getfem::mesher_ball B9(getfem::base_node(0,-1.,.5),1.5);
    getfem::mesher_ball B10(getfem::base_node(0.,1.,.5),1.5);
    getfem::mesher_intersection I2(B9, B10);
    getfem::mesher_union D12(I1, I2);
    
    getfem::mesher_ball B11(getfem::base_node(0.,0.),2.);
    getfem::mesher_ball B12(getfem::base_node(1.,0.),2);
    getfem::mesher_setminus D13(B11, B12);

    getfem::mesher_ball B13(getfem::base_node(0.,0.,0.),2.);
    getfem::mesher_ball B14(getfem::base_node(1.,0.,0.),2);
    getfem::mesher_setminus D14(B13, B14);
    
    getfem::mesher_torus D15(2.0, 0.5);

    getfem::mesher_cylinder C1(base_node(0.5,0.5,-0.5), base_node(0.,0.,1.),
			       2.0, 0.2);

    getfem::mesher_cylinder C2(base_node(0.5,0.5,-0.5), base_node(0.,0.,1.),
			       2.0, 0.65);
    getfem::mesher_setminus S1(D3, C1);
    getfem::mesher_intersection D16(S1, C2);


    getfem::mesher_ball B15(getfem::base_node(0.,0.,0.),1.2);
    getfem::mesher_ball B16(getfem::base_node(-2.,-2.,0.),0.7);
    getfem::mesher_ball B17(getfem::base_node(2.,-2.,0.),0.7);
    getfem::mesher_ball B18(getfem::base_node(2.,2.,0.),0.7);
    getfem::mesher_ball B19(getfem::base_node(-2.,2.,0.),0.7);
    getfem::mesher_cylinder C3(base_node(-2.,-2.,0.), base_node(1.,0.,0.),
			       4.0, 0.3);
    getfem::mesher_cylinder C4(base_node(2.,-2.,0.), base_node(0.,1.,0.),
			       4.0, 0.3);
    getfem::mesher_cylinder C5(base_node(2.,2.,0.), base_node(-1.,0.,0.),
			       4.0, 0.3);
    getfem::mesher_cylinder C6(base_node(-2.,2.,0.), base_node(0.,-1.,0.),
			       4.0, 0.3);
    getfem::mesher_cylinder C7(base_node(0.,-2.,0.), base_node(0.,1.,0.),
			       4.0, 0.3);
    getfem::mesher_cylinder C8(base_node(-2.,0.,0.), base_node(1.,0.,0.),
			       4.0, 0.3);
    getfem::mesher_cylinder C9(base_node(0.,0.,-1.5), base_node(0.,0.,1.),
			       3.0, 0.3);
    
    getfem::mesher_union U1(B15, B16, B17, B18, B19, C3, C4, C5, C6, C7, C8);
    getfem::mesher_setminus D17(U1, C9);
    getfem::mesher_rectangle D18(base_node(0.,0.,0.), base_node(.1,.3,1.));



    getfem::mesher_cylinder C10(base_node(0.5,0.5,1.2), base_node(0.,0.,-1.),
			       3.0, 0.1);
    getfem::mesher_cylinder C11(base_node(0.5,0.5,0.), base_node(0.,0.,-1.),
			       1.0, 0.3);
    getfem::mesher_rectangle P3(base_node(0.,0.,0.), base_node(1.,1.,1.));
    getfem::mesher_rectangle P4(base_node(0.3,0.3,1.0), base_node(0.6,0.6,1.2));
    getfem::mesher_union U2(P3, C11);
    getfem::mesher_setminus M1(U2, C10);
    getfem::mesher_setminus D19(M1, P4);

    getfem::mesher_ball D22(base_node(0.,0.,20.),20.);

    getfem::mesher_rectangle D23_1(base_node(0.,1.), base_node(20.0, 10.0));
    getfem::mesher_rectangle D23_2(base_node(1.,0.), base_node(19.0, 1.0));
    getfem::mesher_ball D23_3(getfem::base_node(1.0,1.0),1.0);
    getfem::mesher_ball D23_4(getfem::base_node(19.0,1.0),1.0);
    getfem::mesher_union D23(D23_1, D23_2, D23_3, D23_4);


    getfem::mesher_signed_distance *dist = &D0;
    switch (opt) {
    case 0: dist = &D0; break; /* disk */
    case 1: {
      getfem::scalar_type z = sqrt(4 - 1.5*1.5);
      fixed.push_back(base_node(0,z));
      fixed.push_back(base_node(0,-z));
    }
    case 2: dist = &D2; break; /* union of 2 disks */
    case 3: dist = &D3; break; /* cube */
    case 4: dist = &D4; break; /* cylinder */
    case 5: dist = &D5; break; /* ball */
    case 6: dist = &D6; break; /* union of 2 balls */
    case 7: dist = &D7; break; /* union of 3 balls */
    case 8: dist = &D8; break; /* union of 3 disks */
    case 9: dist = &D9; break; /* union of 3 half balls */
    case 10: dist = &D10;      /* disk r=20 with fixed points */
      fixed.push_back(getfem::base_node(0.,0.));
      fixed.push_back(getfem::base_node(0.,20.));
      fixed.push_back(getfem::base_node(-20.,20.));
      fixed.push_back(getfem::base_node(20.,20.));
      fixed.push_back(getfem::base_node(0.,40.));
      dist = &D10; 
      break;
    case 11: dist = &D11; break; /* half-balls */
    case 12: dist = &D12; break; /* UFO */
    case 13: dist = &D13; break; /* moon */
    case 14: dist = &D14; break; /* subtraction of two balls */
    case 15: dist = &D15; break; /* torus */
    case 16: dist = &D16; break; /* cube with a hole */
    case 17: dist = &D17; break; /* space station */
    case 18: dist = &D18; break;
    case 19: dist = &D19; break;
    case 20: {                   /* ladder */
      scalar_type H = 20; 
      unsigned nb_step = 3;
      dist = new getfem::mesher_union
	(*new getfem::mesher_rectangle(base_node(-3, -1.2, 0), base_node(-1.8, 1.2, H)), 
	 *new getfem::mesher_rectangle(base_node(1.8, -1.2, 0), base_node(3, 1.2, H)));
      for (unsigned i=0; i < nb_step; ++i) {
	scalar_type z = H/nb_step/2 + i*(H/nb_step);
	getfem::mesher_signed_distance *d = 
	  new getfem::mesher_cylinder(base_node(-2.4, 0, z), base_node(1, 0, 0), 4.8, .8);
	dist = new getfem::mesher_union(*dist, *d);
      }
    } break;
    case 21: { /* bar with holes */
      scalar_type H=20;
      unsigned nb_holes = 4;
      dist = new getfem::mesher_rectangle(base_node(-2, -1, 0), 
					  base_node(+2, +1, H));
      for (unsigned i=0; i < nb_holes; ++i) {
	scalar_type z = H/nb_holes/2 + i*(H/nb_holes);
	getfem::mesher_signed_distance *d = 
	  new getfem::mesher_cylinder(base_node(0, -20, z), 
				      base_node(0, 1, 0), 40, .7);
	dist = new getfem::mesher_setminus(*dist, *d);
      }
    } break;
    case 22: dist = &D22; break; /* ball for test of radius 20. */
    case 23:                   /* rectangle with rounded corners */
      fixed.push_back(getfem::base_node(0.,1.));
      fixed.push_back(getfem::base_node(1.,0.));
      fixed.push_back(getfem::base_node(19.,0.));
      fixed.push_back(getfem::base_node(20.,1.));
      dist = &D23;
      break;
    }
    getfem::build_mesh(m, *dist, h, fixed, K, 2, max_iter, prefind);
    cerr << "You can view the result with\n mayavi -d totoq.vtk -m BandedSurfaceMap\n";
  }
  GMM_STANDARD_CATCH_ERROR;
  return 0;
}
