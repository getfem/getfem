// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_level_set.h : Dealing with level set representation.
//           
// Date    : January 31, 2005.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1999-2005 Yves Renard, Julien Pommier
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
//========================================================================


#ifndef GETFEM_LEVEL_SET_H__
#define GETFEM_LEVEL_SET_H__

#include <dal_shared_ptr.h>
#include <getfem_mesh_fem.h>

namespace getfem {


  
  class level_set {
    
  public :
    struct mf_key {
      getfem_mesh *pmesh;
      dim_type order;
      mf_key(getfem_mesh &mesh, dim_type o) : pmesh(&mesh),order(o) {}
      bool operator <(const mf_key &a) const;
    };

  protected :
    typedef dal::shared_ptr<mesh_fem> pmesh_fem;
    pmesh_fem add_mesh_fem(getfem_mesh &mesh, dim_type o);
    void sup_mesh_fem(getfem_mesh &mesh, dim_type o);

    static std::map<mf_key, pmesh_fem> mesh_fems;
    getfem_mesh *pmesh;
    dim_type order_;
    pmesh_fem mf;
    std::vector<scalar_type> primary_, secondary_;

  public :
    
    std::vector<scalar_type> &primary(void) { return primary_; }
    const std::vector<scalar_type> &primary(void) const { return primary_; }
    std::vector<scalar_type> &secondary(void) { return secondary_; }
    const std::vector<scalar_type> &secondary(void) const {return secondary_;}

    mesh_fem &get_mesh_fem(void) { return *mf; }

    level_set(getfem_mesh &mesh, dim_type o = dim_type(1))
      : pmesh(&mesh), order_(o), mf(add_mesh_fem(mesh, o)) {
      primary_.resize(mf->nb_dof());
      secondary_.resize(mf->nb_dof());
    }

    ~level_set() { sup_mesh_fem(*pmesh, order_); }

  };
 

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_LEVEL_SET_H__  */
