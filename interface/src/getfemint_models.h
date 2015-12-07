/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2009-2015 Yves Renard.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

/**\file getfemint_models.h
   \brief getfem::models interface
*/

#include <getfemint.h>
#include <getfemint_object.h>
#include <getfem/getfem_models.h>

namespace getfemint {

  class getfemint_model : public getfemint::getfem_object {
  private:
    std::unique_ptr<getfem::model> md;
  public:
    ~getfemint_model() {}
    id_type class_id() const { return MODEL_CLASS_ID; }
    bool is_complex() const { return md->is_complex() != 0; }
    size_type memsize() const { 
      /* Rough estimate. Quite optimistic with sparse matrices! */
      if (!md->is_complex()) {
	size_type szd = sizeof(double), szs = sizeof(size_type);
	return gmm::nnz(md->real_tangent_matrix()) * (szd + szs)
	  + gmm::vect_size(md->real_rhs()) *  szd * 3
	  + sizeof(getfem::model);
      }
      else {
	size_type szc = sizeof(std::complex<double>), szs = sizeof(size_type);
	return gmm::nnz(md->complex_tangent_matrix()) * (szc + szs)
	  + gmm::vect_size(md->complex_rhs()) *  szc * 3
	  + sizeof(getfem::model);
      }
    }
    void clear() const { md->clear(); }
    getfem::model &model() { return *md; }
    void set(getfem::model *p) { md.reset(p); }
  private:
    
    void clear(getfem::model &mdd) const { mdd.clear(); }

  };

  inline bool object_is_model(getfem_object *o) {
    return o->class_id() == MODEL_CLASS_ID;
  }

  inline getfemint_model* object_to_model(getfem_object *o) {
    if (object_is_model(o)) return (getfemint_model*)o;
    else THROW_INTERNAL_ERROR;
  }
}
