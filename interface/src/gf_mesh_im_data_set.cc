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
  General function for modifying mesh_im objects
  @*/

void gf_mesh_im_data_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2)
    THROW_BADARG( "Wrong number of input arguments");

  getfem::im_data *mimd = to_meshimdata_object(in.pop());
  std::string cmd = in.pop().to_string();

  if (check_cmd(cmd, "region", in, out, 1, 1, 0, 0)) {
    /*@SET ('region', @int rnum)
    Set the considered region to `rnum`.
    @*/
    size_type rnum = size_type(in.pop().to_integer());
    mimd->set_region(rnum);
  } else if (check_cmd(cmd, "tensor size", in, out, 1, 1, 0, 0)) {
    /*@SET ('tensor size', )
    Set the size of the data per integration point.
    @*/
    iarray v = in.pop().to_iarray(-1);
    bgeot::multi_index tensor_size(v.size());
    std::copy(v.begin(), v.end(), tensor_size.begin());
    mimd->set_tensor_size(tensor_size);
  } else bad_cmd(cmd);
}
