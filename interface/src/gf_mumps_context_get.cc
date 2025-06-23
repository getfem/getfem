/*===========================================================================

 Copyright (C) 2025-2025 Konstantinos Poulios.

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

#include <getfemint_workspace.h>
#include <getfemint_gmumps.h>

using namespace getfemint;

/*@GFDOC
  General function for querying information about a mumps_context object.
@*/

void gf_mumps_context_get(getfemint::mexargs_in& in,
                          getfemint::mexargs_out& out) {

  if (in.narg() < 2) THROW_BADARG("Wrong number of input arguments");

  const gmumps *pctx   = to_mumps_context_object(in.pop());
  std::string init_cmd = in.pop().to_string();
  std::string cmd      = cmd_normalize(init_cmd);

  if (check_cmd(cmd, "display", in, out, 0, 0, 0, 0)) {
    /*@GET ('display')
      Display a short summary for a @tmct object.@*/
    infomsg() << "gfMumpsContext object\n";
  } else if (check_cmd(cmd, "matrix", in, out, 0, 0, 0, 0)) {
    /*@GET K = ('matrix')
      The @tmct object in the scripting API cannot return
      a reference to the ija matrix stored in C++. It just
      displays information about the stored matrix.@*/
    infomsg() << "This object does not provide access to the stored matrix\n";
  } else if (check_cmd(cmd, "distributed matrix", in, out, 0, 0, 0, 0)) {
    /*@GET K = ('distributed matrix')
      The @tmct object in the scripting API cannot return
      a reference to the ija matrix stored in C++. It just
      displays information about the stored matrix.@*/
    infomsg() << "This object does not provide access to the stored matrix\n";
  } else if (check_cmd(cmd, "vector", in, out, 0, 0, 0, 1)) {
    /*@GET vec = ('vector')
      Outputs a copy of the vector stored in the @tmct object.
      Before the solution, it contains the right hand side of the system.
      After the solution, it contains the solution.@*/
    if (pctx->is_complex()) out.pop().from_dcvector(pctx->vector_c());
    else                    out.pop().from_dcvector(pctx->vector_r());
  } else if (check_cmd(cmd, "ICNTL", in, out, 1, 1, 0, 1)) {
    /*@GET VAL = ('ICNTL', @int ind)
      Output the ICNTL array entry at 1-based index `ind` stored in
      the @tmct object.

      Capital naming convention is used to imply Fortran indexing.@*/
    int ind = in.pop().to_integer();
    if (ind <= 0 || ind > 60)
      THROW_BADARG("Invalid ICNTL parameter index " << ind);
    out.pop().from_integer(pctx->ICNTL(ind));
  } else if (check_cmd(cmd, "CNTL", in, out, 1, 1, 0, 1)) {
    /*GET ('CNTL', @int ind)
      Output the CNTL array entry at 1-based index `ind` stored in
      the @tmct object.

      Capital naming convention is used to imply Fortran indexing.@*/
    int ind = in.pop().to_integer();
    if (ind <= 0 || ind > 15)
       THROW_BADARG("Invalid CNTL parameter index " << ind);
     out.pop().from_scalar(pctx->CNTL(ind));
  } else if (check_cmd(cmd, "INFO", in, out, 1, 1, 0, 1)) {
    /*GET ('INFO', @int ind)
      Output the INFO array entry at 1-based index `ind` stored in
      the @tmct object.

      Capital naming convention is used to imply Fortran indexing.@*/
    int ind = in.pop().to_integer();
    if (ind <= 0 || ind > 80)
      THROW_BADARG("Invalid INFO array index " << ind);
    out.pop().from_integer(pctx->INFO(ind));
  } else if (check_cmd(cmd, "INFOG", in, out, 1, 1, 0, 1)) {
    /*@GET VAL = ('INFOG', @int ind)
      Output the INFOG array entry at 1-based index `ind` stored in
      the @tmct object.

      Capital naming convention is used to imply Fortran indexing.@*/
    int ind = in.pop().to_integer();
    if (ind <= 0 || ind > 80)
      THROW_BADARG("Invalid INFOG array index " << ind);
    out.pop().from_integer(pctx->INFOG(ind));
  } else if (check_cmd(cmd, "RINFO", in, out, 1, 1, 0, 1)) {
    /*@GET VAL = ('RINFO', @int ind)
      Output the RINFO array entry at 1-based index `ind` stored in
      the @tmct object.

    Capital naming convention is used to imply Fortran indexing.@*/
    int ind = in.pop().to_integer();
    if (ind <= 0 || ind > 20)
      THROW_BADARG("Invalid RINFO array index " << ind);
    out.pop().from_scalar(pctx->RINFO(ind));
  } else if (check_cmd(cmd, "RINFOG", in, out, 1, 1, 0, 1)) {
    /*@GET VAL = ('RINFOG', @int ind)
      Output the RINFOG array entry at 1-based index `ind` stored in
      the @tmct object.

      Capital naming convention is used to imply Fortran indexing.@*/
    int ind = in.pop().to_integer();
    if (ind <= 0 || ind > 20)
      THROW_BADARG("Invalid RINFOG array index " << ind);
    out.pop().from_scalar(pctx->RINFOG(ind));
  } else
    bad_cmd(init_cmd);

}
