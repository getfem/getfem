#include <getfem_derivatives.h>

namespace getfem
{  
  /* ********************************************************************* */
  /*       Precomputation on geometric transformations.                    */
  /* ********************************************************************* */

  struct _ipre_geot_light
  {
    bgeot::pgeometric_trans pgt;
    pfem pf;
    bool operator < (const _ipre_geot_light &ls) const
    {
      if (pgt < ls.pgt) return true; if (pgt > ls.pgt) return false; 
      if (pf < ls.pf) return true; return false;
    }
    _ipre_geot_light(bgeot::pgeometric_trans pg, pfem pff)
    { pgt = pg; pf = pff; }
    _ipre_geot_light(void) { }
   
  };


  _geotrans_iprecomp::_geotrans_iprecomp(const _ipre_geot_light &ls)
  {
    base_poly P, Q;
    dim_type N = ls.pgt->structure()->dim();
    assert(ls.pgt->basic_structure() == ls.pf->basic_structure());
    pgt = ls.pgt; pf = ls.pf;
    pc.resize(pf->nb_dof());
    hpc.resize(pf->nb_dof());
    std::fill(pc.begin(), pc.end(),
	      base_matrix(pgt->nb_points() , N));
    std::fill(hpc.begin(), hpc.end(),
	      base_matrix(dal::sqr(N), pgt->nb_points()));
    for (size_type i = 0; i < pgt->nb_points(); ++i)
      for (dim_type n = 0; n < N; ++n) {
	P = pgt->poly_vector()[i];
	P.derivative(n);
	for (size_type j = 0; j < pf->nb_dof(); ++j)
	  pc[j](i,n) = P.eval(pf->node_of_dof(j).begin());
	for (dim_type m = 0; m <= n; ++m) {
	  Q = P; Q.derivative(m);
	  for (size_type j = 0; j < pf->nb_dof(); ++j)
	    hpc[j](m * N + n, i) = hpc[j](n * N + m, i)
	      = P.eval(pf->node_of_dof(j).begin());
	}
      }
  }
  
  pgeotrans_iprecomp geotrans_iprecomp(bgeot::pgeometric_trans pg,
				       pfem pf)
  { 
    static dal::FONC_TABLE<_ipre_geot_light, _geotrans_iprecomp> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_ipre_geot_light, _geotrans_iprecomp>();
      isinit = true;
    }
    return tab->add(_ipre_geot_light(pg, pf));
  }
}
