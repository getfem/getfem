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

/*GFDOC
  This object is used for setting up a MUMPS context.
*/
void gf_mumps_context(getfemint::mexargs_in& in, getfemint::mexargs_out& out) {
  if (check_cmd("MumpsContext", "MumpsContext", in, out, 0, 2, 0, 1)) {
    /*@INIT CTX = ('.init', @str sym, @str datatype)
      The argument `sym` should be 'symmetric' or 'unsymmetric' or empty.
      Default value is 'unsymmetric'.
      The argument `datatype` should be 'real' or 'complex' or empty.
      Default value is 'real'.@*/
    bool sym(false);
    if (in.remaining()) {
      if (in.front().is_string()) {
        std::string opt = in.pop().to_string();
        if (cmd_strmatch(opt, "symmetric"))
          sym = true;
        else if (cmd_strmatch(opt, "unsymmetric"))
          sym = false;
        else
          THROW_BADARG("invalid symmetry option: " << opt);
      } else
        THROW_BADARG("invalid argument to MumpsContext, expected string");
    }
    bool complex(false);
    if (in.remaining()) {
      if (in.front().is_string()) {
        std::string opt = in.pop().to_string();
        if (cmd_strmatch(opt, "complex"))
          complex = true;
        else if (cmd_strmatch(opt, "real"))
          complex = false;
        else
          THROW_BADARG("invalid data type option: " << opt);
      } else
        THROW_BADARG("invalid argument to MumpsContext, expected string");
    }

    auto pstored = std::make_shared<gmumps>(sym, complex);
    id_type id = store_mumps_context_object(pstored);
    out.pop().from_object_id(id, MUMPS_CONTEXT_CLASS_ID);
  }
}


