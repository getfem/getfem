
#include <getfem_mesher.h>
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif


int main(int argc, char **argv) {
  dal::exception_callback_debug cb;
  dal::exception_callback::set_exception_callback(&cb);

#ifdef GETFEM_HAVE_FEENABLEEXCEPT /* trap SIGFPE */
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
  try {
    getfem::getfem_mesh m;
    getfem::scalar_type h = .3;
    int K=1;
    int opt = 0;
    if (argc > 1) { opt = atoi(argv[1]); }
    if (argc > 2) { h = atof(argv[2]); cout << "h = " << h << "\n"; }
    if (argc > 3) { K = atoi(argv[3]); }
    assert(K>0); assert(h>0.); 
    std::vector<getfem::base_node> fixed;
    switch (opt) {
    case 0: {
      getfem::build_mesh(m, getfem::mesher_ball(getfem::base_node(0.,0.),1.), 
			 h, fixed, K);
    } break;
    case 1: {
      getfem::scalar_type z = sqrt(4 - 1.5*1.5);
      fixed.push_back(getfem::base_node(0,z));
      fixed.push_back(getfem::base_node(0,-z));
    }
    case 2: {
      getfem::build_mesh(m, getfem::mesher_union
			 (getfem::mesher_ball(getfem::base_node(1.5,0.),2.),
			  getfem::mesher_ball(getfem::base_node(-1.5,0.),2.)),
			 h, fixed, K);
    } break;
    case 3: {
      getfem::build_mesh(m, getfem::mesher_rectangle
			 (getfem::base_node(0.,0.,0.),
			  getfem::base_node(1.,1.,1.)),
			 h, fixed, K);
    } break;
    case 4: {
      getfem::build_mesh(m, getfem::mesher_cylinder(), 
			 h, fixed, K);
    } break;
    case 5: {
      getfem::build_mesh
	(m, getfem::mesher_ball(getfem::base_node(0.,0.,0.),1.), 
	 h, fixed, K);
    } break;
    case 6: {
      getfem::build_mesh
	(m, getfem::mesher_union
	 (getfem::mesher_ball(getfem::base_node(1.5,0.,0.),2.),
	  getfem::mesher_ball(getfem::base_node(-1.5,0.,0.),2.)),
	 h, fixed, K);
    } break;
    case 7: { /* union of 2 balls */
      getfem::build_mesh
	(m, getfem::mesher_union
	 (getfem::mesher_union
	  (getfem::mesher_ball(getfem::base_node(1.5,0.,0.),2.),
	   getfem::mesher_ball(getfem::base_node(-1.5,0.,0.),2.)),
	  getfem::mesher_ball(getfem::base_node(0,1.,0.),2.)),
	 h, fixed, K);
    } break;
    case 8: { /* union of 2 disks */
      getfem::build_mesh
	(m, getfem::mesher_union
	 (getfem::mesher_union
	  (getfem::mesher_ball(getfem::base_node(1.5,0.),2.),
	   getfem::mesher_ball(getfem::base_node(-1.5,0.),2.)),
	  getfem::mesher_ball(getfem::base_node(0,1.),2.)),
	 h, fixed, K);
    } break;
    case 9: { /* tri-balls + half-plan */
      getfem::build_mesh
	(m, getfem::mesher_intersection
	 (getfem::mesher_union
	  (getfem::mesher_union
	   (getfem::mesher_ball(getfem::base_node(1.5,0.,0.),2.),
	    getfem::mesher_ball(getfem::base_node(-1.5,0.,0.),2.)),
	   getfem::mesher_ball(getfem::base_node(0,1.,0.),2.)),
	  getfem::mesher_half_space(getfem::base_node(0,.1,.1),
				      getfem::base_node(.1,-.03,1.))),
	 h, fixed, K);
    } break;

    case 10: { /* disk r=20 */
      fixed.push_back(getfem::base_node(0.,0.));
      fixed.push_back(getfem::base_node(0.,20.));
      fixed.push_back(getfem::base_node(-20.,20.));
      fixed.push_back(getfem::base_node(20.,20.));
      fixed.push_back(getfem::base_node(0.,40.));
      getfem::build_mesh
	(m, getfem::mesher_ball(getfem::base_node(0.,20.),20.), 
	 h, fixed, K);
    } break;
    case 11: { /* half-balls */
      getfem::build_mesh
	(m, getfem::mesher_intersection
	 (getfem::mesher_ball(getfem::base_node(0.,0.,0.),1.),
	  getfem::mesher_half_space(getfem::base_node(0,.1,.1), 
				    getfem::base_node(.1,-.03,.5))),
	 h, fixed, K);
    } break;
    case 12: { /* UFO */
      getfem::build_mesh
	(m, getfem::mesher_union
	 (getfem::mesher_intersection
	  (getfem::mesher_ball(getfem::base_node(-1.,0.,0.),1.5),
	   getfem::mesher_ball(getfem::base_node(1.,0.,0.),1.5)),
	  getfem::mesher_intersection
	  (getfem::mesher_ball(getfem::base_node(0,-1.,.5),1.5),
	   getfem::mesher_ball(getfem::base_node(0.,1.,.5),1.5))),
	 h, fixed, K);
    } break;
    case 13: { /* moon */
      getfem::build_mesh
	(m, getfem::mesher_setminus
	 (getfem::mesher_ball(getfem::base_node(0.,0.),2.),
	  getfem::mesher_ball(getfem::base_node(1.,0.),2)),
	 h, fixed, K);
    } break;
    case 14: { /* substraction of two balls */
      getfem::build_mesh
	(m, getfem::mesher_setminus
	 (getfem::mesher_ball(getfem::base_node(0.,0.,0.),2.),
	  getfem::mesher_ball(getfem::base_node(1.,0.,0.),2)),
	 h, fixed, K);
    } break;
    case 15: { /* torus */
      getfem::build_mesh(m, getfem::mesher_torus(2.0, 3.0), h, fixed, K);
    } break;
  }
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0;
}
