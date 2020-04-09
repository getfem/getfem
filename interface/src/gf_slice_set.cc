/*===========================================================================

 Copyright (C) 2006-2020 Yves Renard, Julien Pommier.

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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

#include <map>
#include <getfemint_misc.h>
#include <getfem/getfem_mesh_slice.h>
#include <getfem/getfem_mesh_slice.h>


using namespace getfemint;

/*@GFDOC
  Edition of mesh slices.
@*/


void gf_slice_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out) {
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfem::stored_mesh_slice *sl = to_slice_object(in.pop());

  std::string cmd                  = in.pop().to_string();
  if (check_cmd(cmd, "pts", in, out, 1, 1, 0, 0)) {
    /*@SET ('pts', @dmat P)
    Replace the points of the slice.

    The new points `P` are stored in the columns the matrix. Note that
    you can use the function to apply a deformation to a slice, or to
    change the dimension of the slice (the number of rows of `P` is not
    required to be equal to SLICE:GET('dim')).@*/
    darray w = in.pop().to_darray(-1, int(sl->nb_points()));
    size_type min_dim = 0;
    for (size_type ic=0; ic < sl->nb_convex(); ++ic) {
      for (getfem::mesh_slicer::cs_simplexes_ct::const_iterator it = sl->simplexes(ic).begin();
	   it != sl->simplexes(ic).end(); ++it)
	min_dim = std::max(min_dim, it->dim());
    }
    if (w.getm() < min_dim)
      GMM_THROW(getfemint_error, "can't reduce the dimension of the slice to " <<
		w.getm() << " (it contains simplexes of dimension " << min_dim << ")");
    sl->set_dim(w.getm()); /* resize the points */
    for (size_type ic=0, cnt=0; ic < sl->nb_convex(); ++ic) {
      for (getfem::mesh_slicer::cs_nodes_ct::iterator it=sl->nodes(ic).begin();
           it != sl->nodes(ic).end(); ++it) {
        for (size_type k=0; k < sl->dim(); ++k)
          (*it).pt[k] = w[cnt++];
      }
    }
  } else bad_cmd(cmd);
}
