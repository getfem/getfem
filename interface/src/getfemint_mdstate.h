// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2005-2008 Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.
//
//===========================================================================

/**\file getfemint_mdbrick.h
   \brief getfem::mdbrick_abstract (and derivatives) interface
*/

#include <getfemint_std.h>
#include <getfemint_object.h>
#include <getfem/getfem_modeling.h>

namespace getfemint
{
  typedef getfem::standard_model_state         real_model_state;
  typedef getfem::standard_complex_model_state cplx_model_state;


  class getfemint_mdstate : public getfemint::getfem_object {
  private:
    std::auto_ptr<real_model_state> r_mds;
    std::auto_ptr<cplx_model_state> c_mds;
  public:
    ~getfemint_mdstate() {}
    id_type class_id() const { return MDSTATE_CLASS_ID; }
    bool is_complex() const { return c_mds.get() != 0; }
    size_type memsize() const { 
      if (!is_complex()) return memsize(*r_mds);
      else               return memsize(*c_mds); 
    }
    void clear() const {
      if (!is_complex()) clear(*r_mds);
      else               clear(*c_mds);
    }
    real_model_state &real_mdstate() 
      { if (is_complex()) THROW_INTERNAL_ERROR; return *r_mds; }
    cplx_model_state &cplx_mdstate() 
      { if (!is_complex()) THROW_INTERNAL_ERROR; return *c_mds; }

    void set(real_model_state *p) { r_mds.reset(p); }
    void set(cplx_model_state *p) { c_mds.reset(p); }
  private:
    template <typename MDSTATE> size_type memsize(MDSTATE &md) const {
      typedef typename MDSTATE::vector_type::value_type T;
      return 
	/* quite optimistic with sparse matrices! */
	gmm::nnz(md.tangent_matrix()) * (sizeof(T) + sizeof(size_type)) +
	gmm::nnz(md.constraints_matrix()) * (sizeof(T) + sizeof(size_type)) +
	(gmm::vect_size(md.state()) + gmm::vect_size(md.residual()) +
	 gmm::vect_size(md.constraints_rhs()))*sizeof(T);
    }
    template <typename MDSTATE> void clear(MDSTATE &md) const {
      gmm::clear(md.residual());
      gmm::clear(md.state());
      gmm::clear(md.tangent_matrix());
      gmm::clear(md.constraints_matrix());
      gmm::clear(md.constraints_rhs());
    }

  };

  inline bool object_is_mdstate(getfem_object *o) {
    return o->class_id() == MDSTATE_CLASS_ID;
  }

  inline getfemint_mdstate* object_to_mdstate(getfem_object *o) {
    if (object_is_mdstate(o)) return (getfemint_mdstate*)o;
    else THROW_INTERNAL_ERROR;
  }
}
