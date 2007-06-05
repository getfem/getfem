// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 1999-2007 Yves Renard
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


/**
   @file getfem_inter_element.h 
   @author Yves Renard <Yves.Renard@insa-lyon.fr>
   @date June 05, 2007.
   @brief A tool to compute contributions on a face between two elements
          using the finite element methods of both the two elements.
*/

#ifndef GETFEM_INTER_ELEMENT
#define GETFEM_INTER_ELEMENT

#include "getfem_mesh_im.h"
#include "getfem_mesh_fem.h"

namespace getfem {

  papprox_integration get_approx_im_or_fail(pintegration_method pim);

  class interelt_boundary_integration_
    : virtual public dal::static_stored_object {
    
    papprox_integration pai1, pai2;
    mutable std::vector<base_node> add_points;
    mutable std::vector< std::vector<size_type> > indices;
    mutable bool warn_msg;

  public :
    
    interelt_boundary_integration_(pintegration_method pa1,
				   pintegration_method pa2);
    
    std::vector<size_type> &face_indices(size_type f1, size_type f2) const;
    const base_node &additional_point(size_type i) const
    { return add_points[i]; }

  };

  typedef boost::intrusive_ptr<const getfem::interelt_boundary_integration_>
    pinterelt_boundary_integration;

  pinterelt_boundary_integration interelt_boundary_integration
    (pintegration_method pa1, pintegration_method pa2);



  /** This object is to be derived in order to compute some contributions
      on a face between two elements, using the finite element of both
      the two elements (for instance to compute a jump).
      The compute() method has to be redefined. It is called on each point
      of the integration method defined on the considered face.
  **/
  class compute_on_inter_element {

  protected :

    const mesh_im &mim;
    const mesh_fem &mf;

    pfem pf1_old, pf2_old;
    papprox_integration pai_old, pai1_old, pai2_old;
    pfem_precomp pfp1, pfp2;
    pinterelt_boundary_integration pibi;
    base_matrix G1, G2;


    virtual void compute_on_gauss_point
    (getfem::fem_interpolation_context ctx1, getfem::pfem pf1,
     getfem::fem_interpolation_context ctx2, getfem::pfem pf2,
     getfem::papprox_integration pai1) = 0;
    
  public :

    compute_on_inter_element(const mesh_im &mmim, const mesh_fem &mmf)
      : mim(mmim), mf(mmf), pf1_old(0), pf2_old(0), pai_old(0),
	pai1_old(0), pai2_old(0), pfp1(0), pfp2(0), pibi(0) {}
    
    void compute_on_face(size_type cv, size_type f1);
    
    virtual ~compute_on_inter_element() {}
    
  };





 
}


#endif
