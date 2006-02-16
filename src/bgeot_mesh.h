// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_mesh.h : mesh.
//           
// Date    : February 15, 2006.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2006-2006 Yves Renard
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

/**@file bgeot_mesh.h
   @brief Basic mesh definition
*/

#ifndef BGEOT_MESH_H__
#define BGEOT_MESH_H__

#include <bgeot_mesh_structure.h>
#include <bgeot_geometric_trans.h>

namespace bgeot {


  class basic_mesh :  public bgeot::mesh_structure {

  public :
    
    typedef dal::approx_less<scalar_type> val_comp;
    typedef dal::lexicographical_less<base_node, val_comp> pt_comp;
    typedef dal::dynamic_tree_sorted<base_node, pt_comp> PT_TAB;
    typedef bgeot::mesh_structure::ind_cv_ct ind_cv_ct;
    typedef bgeot::mesh_structure::ind_set ind_set;
    typedef bgeot::mesh_structure::ind_pt_face_ct ind_pt_face_ct;
    typedef bgeot::mesh_structure::point_ct point_ct;

    typedef dal::tab_ref_index_ref
    <PT_TAB::const_iterator, ind_cv_ct::const_iterator> ref_mesh_pt_ct;
    typedef dal::tab_ref_index_ref
    <PT_TAB::const_iterator, ind_pt_face_ct::const_iterator> ref_mesh_face_pt_ct;
    
    
  protected :
    
    dim_type dimension;
    PT_TAB pts;

    scalar_type eps_p; /* infinity distance under wich two points are equal. */
    dal::dynamic_array<bgeot::pgeometric_trans> gtab;
    dal::bit_vector trans_exists;
 
  public :
    basic_mesh(dim_type NN = dim_type(-1)) {
      dimension = NN; eps_p = 1.0E-10;
      pts.comparator() = dal::lexicographical_less<base_node,
	dal::approx_less<base_node::value_type> >(eps_p);
    }
  };



}

#endif /* BGEOT_MESH_H__                                         */
