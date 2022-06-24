/*===========================================================================

 Copyright (C) 2014-2020 Konstantinos Poulios.

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

#include <getfemint.h>
#include <getfem/getfem_im_data.h>


using namespace getfemint;


/*@GFDOC
  This object represents data defined on a mesh_im object.
@*/

void gf_mesh_im_data(getfemint::mexargs_in& m_in, getfemint::mexargs_out& m_out) {

  if (m_in.narg() < 1 || m_in.narg() > 3) THROW_BADARG("Wrong number of input arguments");
  if (!m_out.narg_in_range(1, 1)) THROW_BADARG("Wrong number of output arguments");
  if (!is_meshim_object(m_in.front())) THROW_BADARG("Wrong type of input argument, mesh_im expected");

  /*@INIT MIMD = ('.mesh_im', @tmim mim, @int region, @ivec size)
    Build a new @tmimd object linked to a @tmim object. If `region` is
    provided, considered integration points are filtered in this region.
    `size` is a vector of integers that specifies the dimensions of the
    stored data per integration point. If not given, the scalar stored
    data are considered.
  @*/
  getfem::mesh_im *mim = to_meshim_object(m_in.pop());
  size_type rnum = size_type(-1);
  if (m_in.remaining())
    rnum = size_type(m_in.pop().to_integer());
  bgeot::multi_index tensor_size(1);
  tensor_size[0] = 1;
  if (m_in.remaining()) {
    iarray v = m_in.pop().to_iarray(-1);
    tensor_size.resize(v.size());
    std::copy(v.begin(), v.end(), tensor_size.begin());
  }
  
  auto mimd = std::make_shared<getfem::im_data>(*mim);
  mimd->set_region(rnum);
  mimd->set_tensor_size(tensor_size);
  id_type id = store_meshimdata_object(mimd);
  m_out.pop().from_object_id(id, MESHIMDATA_CLASS_ID);
}
