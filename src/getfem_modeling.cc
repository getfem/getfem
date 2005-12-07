// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_modeling.cc : Defines model bricks to build
//           complete models.
// Date    : June 15, 2004.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2004-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//



namespace getfem {


  void mdbrick_abstract_common_base::parameters_set_uptodate(void) {
    for (PARAM_MAP::iterator it = parameters.begin(); it != parameters.end(); ++it)
      it->second->set_uptodate();
  }
  
  bool mdbrick_abstract_common_base::parameters_is_any_modified(void) const {
    for (PARAM_MAP::const_iterator it = parameters.begin(); it != parameters.end(); ++it)
      if (it->second->is_modified()) return true;
    return false;
  }


}

