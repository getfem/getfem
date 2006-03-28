// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2000-2006 Yves Renard
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

#include <dal_singleton.h>
#include <getfem_mat_elem_type.h>
#include <dal_tree_sorted.h>
#include <getfem_mat_elem.h> /* for mat_elem_forget_mat_elem_type */
namespace getfem {

  bool operator < (const constituant &m, const constituant &n) {
    if (m.t < n.t) return true; if (m.t > n.t) return false;
    if (m.t == GETFEM_NONLINEAR_) {
      if (m.nlt < n.nlt) return true; if (n.nlt < m.nlt) return false;
      if (m.nl_part < n.nl_part) return true; 
    }
    else if (m.pfi < n.pfi) return true;
    return false;
  }

  struct mat_elem_type_key : virtual public dal::static_stored_object_key {
    const mat_elem_type *pmet;
  public :
    virtual bool compare(const static_stored_object_key &oo) const {
      const mat_elem_type_key &o
	= dynamic_cast<const mat_elem_type_key &>(oo);
      if (dal::lexicographical_less<mat_elem_type>()(*pmet, *(o.pmet)) < 0)
	return true;
      return false;
    }
    mat_elem_type_key(const mat_elem_type *p) : pmet(p) {}
  };

  static pmat_elem_type add_to_met_tab(const mat_elem_type &f) {
    dal::pstatic_stored_object o
      = dal::search_stored_object(mat_elem_type_key(&f));
    if (o) return dal::stored_cast<mat_elem_type>(o);
    pmat_elem_type p = new mat_elem_type(f);
    dal::add_stored_object(new mat_elem_type_key(p.get()), p,
			   dal::AUTODELETE_STATIC_OBJECT);
    for (size_type i=0; i < f.size(); ++i) {
      if (f[i].pfi) dal::add_dependency(p, f[i].pfi);
      if (f[i].t == GETFEM_NONLINEAR_ && f[i].nl_part==0)
	f[i].nlt->register_mat_elem(p);
    }
    return p;
  }

  /* on destruction, all occurences of the nonlinear term are removed
     from the mat_elem_type cache; */
  nonlinear_elem_term::~nonlinear_elem_term() {
    for (std::set<pmat_elem_type>::iterator it=melt_list.begin();
	 it != melt_list.end(); ++it)
      dal::del_stored_object(*it);
  }

  pmat_elem_type mat_elem_base(pfem pfi) {
    mat_elem_type f; f.resize(1); f[0].t = GETFEM_BASE_; f[0].pfi = pfi;
    f[0].nlt = 0;
    if (pfi->target_dim() == 1)
      { f.get_mi().resize(1); f.get_mi()[0] = 1; }
    else {
      f.get_mi().resize(2); f.get_mi()[0] = 1;
      f.get_mi()[1] = pfi->target_dim();
    }
    return add_to_met_tab(f);
  }

  pmat_elem_type mat_elem_unit_normal(void) {
    mat_elem_type f; f.resize(1); f[0].t = GETFEM_UNIT_NORMAL_;
    f[0].pfi = 0; f[0].nlt = 0;
    f.get_mi().resize(1); f.get_mi()[0] = 1;
    return add_to_met_tab(f);
  }

  pmat_elem_type mat_elem_grad_geotrans(bool inverted) {
    mat_elem_type f; f.resize(1); 
    f[0].t = (!inverted) ? GETFEM_GRAD_GEOTRANS_ : GETFEM_GRAD_GEOTRANS_INV_;
    f[0].pfi = 0; f[0].nlt = 0;
    f.get_mi().resize(2); f.get_mi()[0] = f.get_mi()[1] = 1;
    return add_to_met_tab(f);    
  }

  pmat_elem_type mat_elem_grad(pfem pfi) {
    mat_elem_type f; f.resize(1); f[0].t = GETFEM_GRAD_; f[0].pfi = pfi;
    f[0].nlt = 0;
    if (pfi->target_dim() == 1) { 
      f.get_mi().resize(2); f.get_mi()[0] = 1;
      f.get_mi()[1] = pfi->dim();
    }
    else {
      f.get_mi().resize(3); f.get_mi()[0] = 1;
      f.get_mi()[1] = pfi->target_dim();
      f.get_mi()[2] = pfi->dim();
    }
    return add_to_met_tab(f);
  }

  pmat_elem_type mat_elem_hessian(pfem pfi) {
    mat_elem_type f; f.resize(1);  f[0].t = GETFEM_HESSIAN_; f[0].pfi = pfi;
    f[0].nlt = 0;
    if (pfi->target_dim() == 1) { 
      f.get_mi().resize(2); f.get_mi()[0] = 1;
      f.get_mi()[1] = gmm::sqr(pfi->dim());
    }
    else {
      f.get_mi().resize(3); f.get_mi()[0] = 1;
      f.get_mi()[1] = pfi->target_dim();
      f.get_mi()[2] = gmm::sqr(pfi->dim());
    }
    return add_to_met_tab(f);
  }

  static pmat_elem_type mat_elem_nonlinear_(pnonlinear_elem_term nlt,
					    pfem pfi, unsigned nl_part) {
    mat_elem_type f; f.resize(1); 
    f[0].t = GETFEM_NONLINEAR_; f[0].nl_part = nl_part;
    f[0].pfi = pfi;
    f[0].nlt = nlt;
    if (nl_part) {
      f.get_mi().resize(1); f.get_mi()[0] = 1;
    } else f.get_mi() = nlt->sizes();
    return add_to_met_tab(f);
  }

  pmat_elem_type mat_elem_nonlinear(pnonlinear_elem_term nlt,
				    std::vector<pfem> pfi) {
    if (pfi.size() == 0)
      DAL_THROW(dal::dimension_error, "mat_elem_nonlinear with no pfem!");
    pmat_elem_type me = mat_elem_nonlinear_(nlt, pfi[0], 0);
    for (size_type i=1; i < pfi.size(); ++i)
      me = mat_elem_product(mat_elem_nonlinear_(nlt, pfi[i], i),me);
    return me;
  }

  pmat_elem_type mat_elem_product(pmat_elem_type a, pmat_elem_type b) {
    mat_elem_type f; f.reserve(a->size() + b->size());
    f.get_mi().reserve(a->get_mi().size() + b->get_mi().size());
    f.insert(f.end(), (*a).begin(), (*a).end());
    f.insert(f.end(), (*b).begin(), (*b).end());
    f.get_mi().insert(f.get_mi().end(), (*a).get_mi().begin(),
		      (*a).get_mi().end());
    f.get_mi().insert(f.get_mi().end(), (*b).get_mi().begin(),
		      (*b).get_mi().end());
    return add_to_met_tab(f);
  }

  pmat_elem_type mat_elem_empty() {
    return add_to_met_tab(mat_elem_type());
  }
  bgeot::multi_index mat_elem_type::sizes(size_type cv) const {
    bgeot::multi_index mii = mi;
    for (size_type i = 0, j = 0; i < size(); ++i, ++j) {
      switch ((*this)[i].t) {
	case GETFEM_BASE_ :
	  mii[j] = (*this)[i].pfi->nb_base(cv);
	  if ((*this)[i].pfi->target_dim() != 1) ++j;
	  break;
	case GETFEM_GRAD_ :
	  mii[j] = (*this)[i].pfi->nb_base(cv); ++j;
	  if ((*this)[i].pfi->target_dim() != 1) ++j;
	  break;     
	case GETFEM_HESSIAN_   :
	  mii[j] = (*this)[i].pfi->nb_base(cv); ++j;
	  if ((*this)[i].pfi->target_dim() != 1) ++j;
	  break;
	case GETFEM_UNIT_NORMAL_ :
	  break;
	case GETFEM_NONLINEAR_ :
	  if ((*this)[i].nl_part == 0)
	    { j+=(*this)[i].nlt->sizes().size(); --j; }
	  break;
	case GETFEM_GRAD_GEOTRANS_:
	case GETFEM_GRAD_GEOTRANS_INV_: 
	  break;
      }
    }
    return mii;
  }


}  /* end of namespace getfem.                                            */

