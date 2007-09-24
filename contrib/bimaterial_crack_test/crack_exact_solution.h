#include "getfem/getfem_mesh_fem_global_function.h"


/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

#define VALIDATE_XFEM

#ifdef VALIDATE_XFEM

struct crack_exact_solution_function : public getfem::abstract_xy_function {
  unsigned function_num;
  unsigned component_num; /* 0 -> x component, 1 -> y component */
  scalar_type lambda, mu;
  virtual scalar_type val(scalar_type x, scalar_type y) const;
  virtual base_small_vector grad(scalar_type x, scalar_type y) const;
  virtual base_matrix hess(scalar_type, scalar_type) const
  { GMM_ASSERT1(false, "Sorry, to be done ..."); }
  crack_exact_solution_function(unsigned fnum, 
				unsigned cnum,
				scalar_type l, scalar_type m) {
    function_num = fnum; 
    component_num = cnum;
    lambda = l; mu = m;
  }
  base_small_vector eval(const base_node &x, base_matrix *pgrad) const;
};


struct crack_exact_solution {
  getfem::mesh_fem_global_function mf;
  getfem::base_vector U;

  crack_exact_solution(getfem::mesh &me) : mf(me) {}
  
  void init(int function_num, scalar_type lambda, scalar_type mu,
	    getfem::level_set &ls);
};

inline base_small_vector sol_f(const base_node &x) {
  int N = x.size();
  base_small_vector res(N);
  return res;
}

#else

inline base_small_vector sol_f(const base_node &x) {
  int N = x.size();
  base_small_vector res(N); res[N-1] = x[N-1];
  return res;
}

#endif
