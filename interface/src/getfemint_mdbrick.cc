// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2006-2006 Yves Renard, Julien Pommier.
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

#include <getfemint_mdbrick.h>

namespace getfemint {
  size_type getfemint_mdbrick::memsize() const { 
    return 0; 
  }

  real_mdbrick_abstract &getfemint_mdbrick::real_mdbrick() { 
    if (is_complex()) 
      THROW_ERROR("cannot use a real-valued model brick in this context");
    return (real_mdbrick_abstract&)(*b.get()); 
  }

  cplx_mdbrick_abstract &getfemint_mdbrick::cplx_mdbrick() {
    if (!is_complex())
      THROW_ERROR("cannot use a complex-valued model brick in this context");
    return (cplx_mdbrick_abstract&)(*b.get()); 
  }
  

  void getfemint_mdbrick::set_brick(real_mdbrick_abstract *p, 
				    const std::string &sclass) { 
    b.reset(p); is_complex_ = false; subclass = sclass; 
  }

  void getfemint_mdbrick::set_brick(cplx_mdbrick_abstract *p, 
				    const std::string &sclass) {
    b.reset(p); is_complex_ = true; subclass = sclass; 
  }

  getfem::mdbrick_abstract_parameter *
  getfemint_mdbrick::param(const std::string &pname) {
    getfem::mdbrick_abstract_common_base::PARAM_MAP::iterator it = 
      mdbrick().get_parameters().find(pname);
    if (it != mdbrick().get_parameters().end()) return it->second;
    else return 0;
  }

  void getfemint_mdbrick::set_constraints_type(getfem::constraints_type ctype) {
    if (!is_complex())
      this->cast<getfem::mdbrick_constraint<real_model_state> >
	("not a constraints brick!")->set_constraints_type(ctype);
    else 
      this->cast<getfem::mdbrick_constraint<cplx_model_state> >
	("not a constraints brick!")->set_constraints_type(ctype);
  }

  
}
