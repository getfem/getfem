#include <getfem_mesh_im_level_set.h>
/* try to enable the SIGFPE if something evaluates to a Not-a-number
 * of infinity during computations
 */
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */


int main(int argc, char **argv) {
  dal::exception_callback_debug cb;
  dal::exception_callback::set_exception_callback(&cb);

#ifdef GETFEM_HAVE_FEENABLEEXCEPT /* trap SIGFPE */
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
  try {
    //getfem::getfem_mesh_im_level_set_noisy();

    getfem::getfem_mesh m; m.read_from_file("meshes/disc_2D_degree3.mesh");
    getfem::mesh_fem mf(m);
    getfem::mesh_im_level_set mim(m, getfem::int_method_descriptor("IM_TRIANGLE(6)"));
    getfem::level_set ls1(m, 2), ls2(m, 2), ls3(m, 3);
    const getfem::mesh_fem &ls1mf = ls1.get_mesh_fem();
    scalar_type R=.4;
    for (unsigned i=0; i < ls1mf.nb_dof(); ++i) {
      ls1.values()[i] = gmm::vect_dist2_sqr(ls1mf.point_of_dof(i), 
					    getfem::base_node(0,0)) -R*R;
    }
    const getfem::mesh_fem &ls2mf = ls2.get_mesh_fem();
    R=.1;
    for (unsigned i=0; i < ls2mf.nb_dof(); ++i) {
      ls2.values()[i] = gmm::vect_dist2_sqr(ls2mf.point_of_dof(i), 
					    getfem::base_node(0,0.3)) -R*R;
    }
    const getfem::mesh_fem &ls3mf = ls3.get_mesh_fem();
    R=.03;
    for (unsigned i=0; i < ls3mf.nb_dof(); ++i) {
      ls3.values()[i] = gmm::vect_dist2_sqr(ls3mf.point_of_dof(i), 
					    getfem::base_node(0,0.38)) -R*R;
    }
    
    mim.add_level_set(ls1);
    mim.add_level_set(ls2);
    mim.add_level_set(ls3);
    mim.adapt();
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0;
}
