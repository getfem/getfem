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
    getfem::getfem_mesh m; m.read_from_file("disc_2D_degree3.mesh");
    getfem::mesh_fem mf(m);
    getfem::mesh_im_level_set  mim(m, getfem::int_method_descriptor("IM_TRIANGLE(6)"));
    getfem::level_set ls(m, 4);
    const getfem::mesh_fem &lsmf = ls.get_mesh_fem();
    scalar_type R=.4;
    for (unsigned i=0; i < lsmf.nb_dof(); ++i) {
      ls.values()[i] = gmm::vect_dist2(lsmf.point_of_dof(i), 
				       getfem::base_node(0,0)) -R;
    }
    mim.add_level_set(ls);
    mim.adapt();
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0;
}
