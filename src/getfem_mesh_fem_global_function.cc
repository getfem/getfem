// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2004-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

#include <getfem_mesh_fem_global_function.h>
#include <getfem_level_set.h>

namespace getfem {
  
  void global_function_fem::init(void) {
    is_pol = is_lag = false; es_degree = 5;
    is_equiv = real_element_defined = true;
    ntarget_dim = 1; // An extension for vectorial elements should be easy
    dim_ = cvr->structure()->dim();
    init_cvs_node();
    base_node P(dim()); P.fill(1./30.);
    for (size_type i = 0; i < functions.size(); ++i) {
      add_node(global_dof(dim()), P);
    }
    std::stringstream nm;
    nm << "GLOBAL_FEM(" << (void*)this << ")";
    debug_name_ = nm.str();
  }

  size_type global_function_fem::nb_dof(size_type) const
  { return functions.size(); }
  
  size_type global_function_fem::index_of_global_dof
  (size_type /*cv*/, size_type i) const
  { return i; }
   
  void global_function_fem::base_value(const base_node &, base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void global_function_fem::grad_base_value(const base_node &,
					 base_tensor &) const
  { DAL_THROW(internal_error, "No grad values, real only element."); }
  void global_function_fem::hess_base_value(const base_node &,
					    base_tensor &) const
  { DAL_THROW(internal_error, "No hess values, real only element."); }
  
  void global_function_fem::real_base_value(const fem_interpolation_context& c,
					    base_tensor &t, bool) const {

    mib.resize(2); mib[0] = target_dim(); mib[1] = functions.size();
    assert(target_dim() == 1);
    t.adjust_sizes(mib);
    for (size_type i=0; i < functions.size(); ++i) {
      /*cerr << "global_function_fem: real_base_value(" << c.xreal() << ")\n";
      if (c.have_G()) cerr << "G = " << c.G() << "\n"; 
      else cerr << "no G\n";*/
      t[i] = (*functions[i]).val(c);
    }
    
  } 
   
  void global_function_fem::real_grad_base_value
  (const fem_interpolation_context& c, base_tensor &t, bool) const {
    mig.resize(3); 
    mig[2] = dim(); mig[1] = target_dim(); mig[0] = functions.size();
    assert(target_dim() == 1);
    t.adjust_sizes(mig);
    base_small_vector G(dim());
    for (size_type i=0; i < functions.size(); ++i) {
      (*functions[i]).grad(c,G);
      for (unsigned j=0; j < dim(); ++j)
	t[j*functions.size() + i] = G[j];
    }
  }
  
  void global_function_fem::real_hess_base_value
  (const fem_interpolation_context &c, base_tensor &t, bool) const { 
    mih.resize(4); 
    mih[3] = mih[2] = dim(); mih[1] = target_dim(); mih[0] = functions.size();
    assert(target_dim() == 1);
    t.adjust_sizes(mih);
    base_matrix H(dim(),dim());
    for (size_type i=0; i < functions.size(); ++i) {
      (*functions[i]).hess(c,H);
      for (unsigned k=0; k < dim(); ++k)
	for (unsigned j=0; j < dim(); ++j)
	  t.at((k*dim() + j)*functions.size() + i) = H.at((k*dim() + j));
    }
  }
  
  DAL_SIMPLE_KEY(special_int_globf_fem_key, pfem);

  pfem new_global_function_fem(bgeot::pconvex_ref cvr, 
			       const std::vector<pglobal_function> &f) {
    pfem pf = new global_function_fem(cvr,f);
    dal::add_stored_object(new special_int_globf_fem_key(pf), pf);
    return pf;
  }

  void mesh_fem_global_function::receipt(const MESH_CLEAR &) { clear(); }
  void mesh_fem_global_function::receipt(const MESH_DELETE &) { clear(); }
  void mesh_fem_global_function::clear_build_methods() {
    for (std::map<bgeot::pconvex_ref,pfem>::const_iterator 
	   it = build_methods.begin(); 
	 it != build_methods.end(); ++it) 
      del_global_function_fem((*it).second);
    build_methods.clear();
  }
  void mesh_fem_global_function::clear(void) {
    mesh_fem::clear();
    clear_build_methods();
  }
  
  void mesh_fem_global_function::adapt(void) {
    clear();
    for (dal::bv_visitor cv(linked_mesh().convex_index()); 
	 !cv.finished(); ++cv) {
      bgeot::pconvex_ref cvr = 
	linked_mesh().trans_of_convex(cv)->convex_ref()->basic_convex_ref();
      
      std::map<bgeot::pconvex_ref,pfem>::iterator it = build_methods.find(cvr);
      pfem pf;
      if (it == build_methods.end()) {
	build_methods[cvr] = pf = new_global_function_fem(cvr, fun);
      } else pf = (*it).second;
      set_finite_element(cv, pf);
    }
    touch();
  }


  /* the basic 4 singular functions for 2D cracks */
  scalar_type 
  crack_singular_xy_function::val(scalar_type x, scalar_type y) const {
    scalar_type sgny = (y < 0 ? -1.0 : 1.0);
    scalar_type r = sqrt(x*x + y*y);

    if (r < 1e-10) {
      return 0;
    }

    /* ci-dessous: la valeur absolue est malheureusement necessaire,
     * sinon il peut arriver qu'on cherche sqrt(-1e-16) ...
     */
    scalar_type sin2 = sqrt(gmm::abs(.5-x/(2*r))) * sgny;
    scalar_type cos2 = sqrt(gmm::abs(.5+x/(2*r)));

    switch(l){
      case 0 : return sqrt(r)*sin2;
      case 1 : return sqrt(r)*cos2; 
      case 2 : return sin2*y/sqrt(r); 
      case 3 : return cos2*y/sqrt(r); 
      default:
	DAL_INTERNAL_ERROR("arg");
    }
  }

  base_small_vector 
  crack_singular_xy_function::grad(scalar_type x, scalar_type y) const {
    base_small_vector res(2);
    scalar_type sgny = (y < 0 ? -1.0 : 1.0);
    scalar_type r = sqrt(x*x + y*y);

    if (r < 1e-10) {
      DAL_WARNING0("Warning, point close to the singularity (r=" << r << ")");
    }

    /* ci-dessous: la valeur absolue est malheureusement necessaire,
     * sinon il peut arriver qu'on cherche sqrt(-1e-16) ...
     */
    scalar_type sin2 = sqrt(gmm::abs(.5-x/(2*r))) * sgny;
    scalar_type cos2 = sqrt(gmm::abs(.5+x/(2*r)));

    switch(l){
    case 0 :
      res[0] = -sin2/(2*sqrt(r));
      res[1] = cos2/(2*sqrt(r));
      break;
    case 1 :
      res[0] = cos2/(2*sqrt(r));  
      res[1] = sin2/(2*sqrt(r));
      break;
    case 2 :
      res[0] = cos2*(-5*cos2*cos2 + 1. + 4*(cos2*cos2*cos2*cos2))/sqrt(r);
      res[1] = sin2*(-3*cos2*cos2 + 1. + 4*(cos2*cos2*cos2*cos2))/sqrt(r);
      break;
    case 3 :
      res[0] = -cos2*cos2*sin2*(4*cos2*cos2 - 3)/sqrt(r);
      res[1] = cos2*(4*cos2*cos2*cos2*cos2 + 2 - 5*cos2*cos2)/sqrt(r);
      break;
    default: 
      DAL_INTERNAL_ERROR("oups");
    }
    return res;
  }


  scalar_type 
  cutoff_xy_function::val(scalar_type x, scalar_type y) const { 
    switch (fun) {
      case EXPONENTIAL_CUTOFF:
	return (a4>0) ? exp(-a4 * gmm::sqr(x*x+y*y)) : 1;
      case POLYNOMIAL_CUTOFF: {
	assert(r0 > r1);
	scalar_type r = gmm::sqrt(x*x+y*y);
	
	if (r <= r1) return 1;
	else if (r >= r0) return 0;
	else {
	  scalar_type c = 6./(pow(r0,3.) - pow(r1,3.) + 3*r1*r0*(r1-r0));
	  scalar_type k = -(c/6.)*(-pow(r0,3.) + 3*r1*pow(r0,2.));
	  return (c/3.)*pow(r,3.) - (c*(r0 + r1)/2.)*pow(r,2.)
	    + c*r0*r1*r + k;
	}
      }
      case POLYNOMIAL2_CUTOFF: {
	assert(r0 > r1);
	scalar_type r = gmm::sqrt(x*x+y*y);
	if (r <= r1) return scalar_type(1);
	else if (r >= r0) return scalar_type(0);
	else {
	  scalar_type c = 1./((-1./30.)*(pow(r1,5) - pow(r0,5)) 
			      + (1./6.)*(pow(r1,4)*r0 - r1*pow(r0,4)) 
			      - (1./3.)*(pow(r1,3)*pow(r0,2) - 
					 pow(r1,2)*pow(r0,3)));
	  scalar_type k = 1. - c*((-1./30.)*pow(r1,5) + 
				  (1./6.)*pow(r1,4)*r0 - 
				  (1./3.)*pow(r1,3)*pow(r0,2));
	  return c*( (-1./5.)*pow(r,5) + (1./2.)* (r1+r0)*pow(r,4) - 
		     (1./3.)*(pow(r1,2)+pow(r0,2) + 4.*r0*r1)*pow(r,3) + 
		     r0*r1*(r0+r1)*pow(r,2) - pow(r0,2)*pow(r1,2)*r) + k;
	}
      }
      default : return scalar_type(1);
    }
  }
    
  base_small_vector 
  cutoff_xy_function::grad(scalar_type x, scalar_type y) const {
    switch (fun) {
      case EXPONENTIAL_CUTOFF: {
	scalar_type r2 = x*x+y*y, ratio = -4.*exp(-a4*r2*r2)*a4*r2;
	return base_small_vector(ratio*x, ratio*y);
      } break;
      case POLYNOMIAL_CUTOFF: {
	scalar_type r = gmm::sqrt(x*x+y*y);
	scalar_type ratio = 0;
	
	if ( r > r1 && r < r0 ) {
	  scalar_type c = 6./(pow(r0,3.) - pow(r1,3.) + 3*r1*r0*(r1-r0));
	  ratio = c*(r - r0)*(r - r1);
	}
	return base_small_vector(ratio*x/r,ratio*y/r);
      } break;
      case POLYNOMIAL2_CUTOFF: {
	scalar_type r = gmm::sqrt(x*x+y*y);
	scalar_type ratio = 0;
	if (r > r1 && r < r0) {
	  scalar_type c = 1./((-1./30.)*(pow(r1,5) - pow(r0,5)) 
			      + (1./6.)*(pow(r1,4)*r0 - r1*pow(r0,4)) 
			      - (1./3.)*(pow(r1,3)*pow(r0,2) - 
					 pow(r1,2)*pow(r0,3)));
	  ratio = - c*gmm::sqr(r-r0)*gmm::sqr(r-r1);
	}
	return base_small_vector(ratio*x/r,ratio*y/r);
      } break;
      default : return base_small_vector(2);
    }
  }

  cutoff_xy_function::cutoff_xy_function(int fun_num, scalar_type r, 
					 scalar_type r1_, scalar_type r0_) 
  {
    fun = fun_num;
    r1 = r1_; r0 = r0_;
    a4 = 0;
    if (r > 0) a4 = pow(2.7/r,4.);
  }


  struct global_function_on_levelset_ : 
    public global_function, public context_dependencies {
    const level_set &ls;
    mutable mesher_level_set mls_x, mls_y;
    mutable size_type cv;

    const abstract_xy_function &fn;

    void update_mls(size_type cv_) const { 
      if (cv_ != cv) { 
	cv=cv_; 
	mls_x=ls.mls_of_convex(cv, 1);
	mls_y=ls.mls_of_convex(cv, 0);
      }
    }

    virtual scalar_type val(const fem_interpolation_context& c) const {
      update_mls(c.convex_num());
      scalar_type x = mls_x(c.xref()), y = mls_y(c.xref());
      return fn.val(x,y);
    }
    virtual void grad(const fem_interpolation_context& c,
		      base_small_vector &g) const {
      update_mls(c.convex_num());
      size_type P = c.xref().size();
      base_small_vector dx(P), dy(P), dfr(2);
      scalar_type x = mls_x.grad(c.xref(), dx);
      scalar_type y = mls_y.grad(c.xref(), dy);

      base_small_vector gfn = fn.grad(x,y);
      gmm::mult(c.B(), gfn[0]*dx + gfn[1]*dy, g);
    }
    virtual void hess(const fem_interpolation_context&, 
		      base_matrix &) const
    { DAL_THROW(dal::to_be_done_error, "hessian to be done ..."); }
    
    void update_from_context(void) const { cv =  size_type(-1); }

    global_function_on_levelset_(const level_set &ls_, 
				 const abstract_xy_function &fn_)
      : ls(ls_), fn(fn_) {
      cv = size_type(-1);
      this->add_dependency(ls);
    }

  };
  

  pglobal_function 
  global_function_on_level_set(const level_set &ls, 
			      const abstract_xy_function &fn) {
    return new global_function_on_levelset_(ls, fn);
  }


  void interpolator_on_mesh_fem::init() {
    base_node min, max;
    scalar_type EPS=1E-13;
    cv_stored = size_type(-1);
    boxtree.clear();
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
      bounding_box(min, max, mf.linked_mesh().points_of_convex(cv),
		   mf.linked_mesh().trans_of_convex(cv));
      for (unsigned k=0; k < min.size(); ++k) { min[k]-=EPS; max[k]+=EPS; }
      boxtree.add_box(min, max, cv);
    }
  }

  bool interpolator_on_mesh_fem::find_a_point(base_node pt, base_node &ptr,
					      size_type &cv) const {
    bool gt_invertible;
    if (cv_stored != size_type(-1) && gic.invert(pt, ptr, gt_invertible))
      { cv = cv_stored; if (gt_invertible) return true; }
    
    boxtree.find_boxes_at_point(pt, boxlst);
    bgeot::rtree::pbox_set::const_iterator it = boxlst.begin(),
      ite = boxlst.end();
    for (; it != ite; ++it) {
      gic = bgeot::geotrans_inv_convex
	(mf.linked_mesh().convex((*it)->id),
	 mf.linked_mesh().trans_of_convex((*it)->id));
      cv_stored = (*it)->id;
      if (gic.invert(pt, ptr, gt_invertible)) { 
	cv = (*it)->id; 
	return true; 
      }
    }
    return false;
  }

  bool interpolator_on_mesh_fem::eval(const base_node pt, base_vector &val, base_matrix &grad) const {
    base_node ptref;
    size_type cv;
    base_vector coeff;
    size_type q = mf.get_qdim(), N = mf.linked_mesh().dim();
    if (find_a_point(pt, ptref, cv)) {
      pfem pf = mf.fem_of_element(cv);
      bgeot::pgeometric_trans pgt = 
	mf.linked_mesh().trans_of_convex(cv);
      base_matrix G; 
      vectors_to_base_matrix(G, mf.linked_mesh().points_of_convex(cv));
      fem_interpolation_context fic(pgt, mf.fem_of_element(cv), ptref, G, cv);
      coeff.resize(mf.nb_dof_of_element(cv));
      gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf.ind_dof_of_element(cv))), coeff);
      val.resize(q);
      mf.fem_of_element(cv)->interpolation(fic, coeff, val, q);
      grad.resize(q, N);
      mf.fem_of_element(cv)->interpolation_grad(fic, coeff, grad, q);
      return true;
    } else return false;
  }
}

/* end of namespace getfem  */




