/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2005-2012 Julien Pommier.
 
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
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/

/**\file getfemint_mdbrick.h
   \brief getfem::mdbrick_abstract (and derivatives) interface
*/

#ifndef GETFEMINT_MDBRICK
#define GETFEMINT_MDBRICK

#include <getfemint_std.h>
#include <getfemint_misc.h>
#include <getfemint_object.h>
#include <getfem/getfem_modeling.h>
#include <getfem/getfem_nonlinear_elasticity.h>
#include <getfem/getfem_plasticity.h>

namespace getfemint
{
  typedef getfem::standard_model_state         real_model_state;
  typedef getfem::standard_complex_model_state cplx_model_state;


  typedef getfem::mdbrick_abstract<real_model_state> real_mdbrick_abstract;
  typedef getfem::mdbrick_abstract<cplx_model_state> cplx_mdbrick_abstract;

  typedef getfem::mdbrick_parameter<real_model_state::vector_type> real_mdbrick_parameter;
  typedef getfem::mdbrick_parameter<cplx_model_state::vector_type> cplx_mdbrick_parameter;


  class getfemint_mdbrick : public getfemint::getfem_object {
  private:
    std::auto_ptr<getfem::mdbrick_abstract_common_base> b;
    bool is_complex_;
    std::string subclass;

  public:
    /* used only by the nonlinear elasticity brick to avoid a memory leak*/
    std::auto_ptr<getfem::abstract_hyperelastic_law> hyperelastic_law;
    std::auto_ptr<getfem::abstract_constraints_projection>  plasticity_stress_projection;

  public:
    ~getfemint_mdbrick() {}
    id_type class_id() const { return MDBRICK_CLASS_ID; }
    bool is_complex() const { return is_complex_; }
    const std::string &sub_class() { return subclass; }
    size_type memsize() const;
    real_mdbrick_abstract &real_mdbrick();
    cplx_mdbrick_abstract &cplx_mdbrick();
    getfem::mdbrick_abstract_common_base &mdbrick()
    { return *b.get(); }

    template <typename T> T *cast(const char *errmsg=0) {
      T *p = dynamic_cast<T*>(b.get());
      if (!p) {
	if (errmsg == 0) { THROW_INTERNAL_ERROR; }
	else THROW_ERROR(errmsg);
      }
      return p;
    }
    template <typename T> T *cast0() {
      return dynamic_cast<T*>(b.get());
    }
    void set_brick(real_mdbrick_abstract *p, 
		   const std::string &sclass);
    void set_brick(cplx_mdbrick_abstract *p, 
		   const std::string &sclass);

    getfem::mdbrick_abstract_parameter *param(const std::string &pname);
    void set_constraints_type(getfem::constraints_type ctype);
  };

  inline bool object_is_mdbrick(getfem_object *o) {
    return o->class_id() == MDBRICK_CLASS_ID;
  }

  inline getfemint_mdbrick* object_to_mdbrick(getfem_object *o) {
    if (object_is_mdbrick(o)) return (getfemint_mdbrick*)o;
    else THROW_INTERNAL_ERROR;
  }
}

#endif 
