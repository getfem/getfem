/*===========================================================================

 Copyright (C) 2011-2020 Yves Renard.

 This file is a part of GetFEM

 GetFEM is free software;  you can  redistribute it  and/or modify it under
 the  terms  of the  GNU  Lesser General Public License as published by the
 Free Software Foundation;  either version 3  of  the License,  or (at your
 option) any  later  version  along with  the GCC Runtime Library Exception
 either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/
// $Id: gf_mesher_object_get.cc 3468 2010-02-24 20:12:15Z renard $

#include <getfemint.h>
#include <getfem/getfem_mesher.h>

using namespace getfemint;

/*@GFDOC
    General function for querying information about mesher_object objects.
@*/


void gf_mesher_object_get(getfemint::mexargs_in& in,
                          getfemint::mexargs_out& out) {

  if (in.narg() < 2) THROW_BADARG("Wrong number of input arguments");

  getfem::pmesher_signed_distance paf = to_mesher_object(in.pop());
  std::string init_cmd                = in.pop().to_string();
  std::string cmd                     = cmd_normalize(init_cmd);

  if (check_cmd(cmd, "char", in, out, 0, 0, 0, 1)) {
    /*@GET s = ('char')
      Output a (unique) string representation of the @tmo.

      This can be used to perform comparisons between two
      different @tmo objects.
      This function is to be completed.@*/
   GMM_ASSERT1(false, "Sorry, function to be done");
   // std::string s = ...;
   // out.pop().from_string(s.c_str());
  } else if (check_cmd(cmd, "display", in, out, 0, 0, 0, 0)) {
    /*@GET ('display')
      displays a short summary for a @tmo object.@*/
    infomsg() << "gfMesherObject object\n";
  } else
    bad_cmd(init_cmd);

}
