/*===========================================================================
 
 Copyright (C) 2013-2013 Yves Renard.
 
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
 
===========================================================================*/

#include <getfemint.h>
#include <getfemint_workspace.h>
#include <getfemint_models.h>
#include <getfemint_multi_contact_frame.h>


using namespace getfemint;

/*@GFDOC
  This object serves for describing a multi-contact situation between
  potentially several deformable bodies and eventually some rigid obstacles.
  (for more details see the Getfem++ user documentation).
@*/

void gf_multi_contact_frame(getfemint::mexargs_in& in, getfemint::mexargs_out& out) {
  getfemint_multi_contact_frame *pgs = NULL;  
  if (check_cmd("MultiContactFrame", "MultiContactFrame", in, out, 3, 9, 0, 1)) {
    
    /*@INIT S = ('.init', @tmodel md, @int N, @scalar release_distance[, @bool delaunay[, @bool self_contact[, @scalar cut_angle[, @bool use_raytrace[, @int fem_nodes_mode[, @bool ref_conf]]]]]])
    Build a new multi contact frame object linked to the model `md`.
    with `N` the space dimension (typically, 2 or 3), `release_distance` is
    the limit distance beyond which two points are not considered in
    potential contact (should be typically comparable to element sizes).
    There is several optional parameters.
    If `fem_node_mode=0` (default value), then contact is considered
    on Gauss points, `fem_node_mode=1` then contact is considered on
    Gauss points for slave surfaces and on f.e.m. nodes for master surfaces
    (in that case, the f.e.m. should be of Lagrange type) and
    `fem_node_mode=2` then contact is considered on f.e.m. nodes for
    both slave and master surfaces. if `use_delaunay` is true (default value),
    then contact detection is done calling
    `Qhull <http://www.qhull.org>`_ package to perform a Delaunay
    triangulation on potential contact points. Otherwise, contact
    detection is performed by conputing some influences boxes of the element
    of master surfaces. If `ref_conf` is true (default value : false),
    the contact detection
    is made on the reference configuration (without taking into account a
    displacement) CAUTION: not fully implemented for the moment.
    If `self_contact` is true (default value), the contact detection is
    also made
    between master surfaces and for a master surface with itself.
    The parameter `cut_angle` (default value: 0.3) is an angle in radian
    which is used
    for the simplification of unit normal cones in the case of f.e.m.
    node contact : if a contact cone has an angle less than `cut_angle`
    it is reduced to a mean unit normal to simplify the contact detection.
    if `use_raytrace` is set to true (default is false) raytracing is used
    insted of projection.
    @*/
    
    getfemint_model *md = in.pop().to_getfemint_model();
    int N = in.pop().to_integer(1, 4);
    scalar_type rd = in.pop().to_scalar();
    bool delaunay = true;
    if (in.remaining()) delaunay = in.pop().to_bool();
    bool self_contact = true;
    if (in.remaining()) self_contact = in.pop().to_bool();
    scalar_type cut_angle = 0.2;
    if (in.remaining()) cut_angle = in.pop().to_scalar();
    bool raytrace = false;
    if (in.remaining()) raytrace = in.pop().to_bool();
    int fem_node_mode = 0;
    if (in.remaining()) fem_node_mode = in.pop().to_integer(0, 2);
    bool ref_conf = false;
    if (in.remaining()) ref_conf = in.pop().to_bool();
    
    getfem::multi_contact_frame *ps
      = new getfem::multi_contact_frame(md->model(), size_type(N), rd,
                                        delaunay, self_contact, cut_angle,
                                        raytrace, fem_node_mode, ref_conf);

       pgs = getfemint_multi_contact_frame::get_from(ps);
       workspace().set_dependance(pgs, md);
  }
  out.pop().from_object_id(pgs->get_id(), MULTI_CONTACT_FRAME_CLASS_ID);
}
