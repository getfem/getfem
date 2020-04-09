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

#include <getfemint_misc.h>

using namespace getfemint;

void
print_poly(bgeot::base_poly *pp) {
  bool first = true; bgeot::size_type n = 0;
  bgeot::base_poly::const_iterator it = pp->begin(), ite = pp->end();
  bgeot::power_index mi(pp->dim());

  if (it != ite && *it != 0.0)
  { mexPrintf("%g", double(*it)); first = false; ++it; ++n; ++mi; }

  for ( ; it != ite ; ++it, ++mi )
  {
    if (*it != 0.0)
    {
      if (!first) { if (*it < 0.0) mexPrintf(" - "); else mexPrintf(" + "); }
      else if (*it < 0.0) mexPrintf("-");
      if (dal::abs(*it) != 1.0) mexPrintf("%g", double(dal::abs(*it)));
      for (int j = 0; j < pp->dim(); ++j)
	if (mi[j] != 0)
	{ 
	  mexPrintf("%c", (j < 3) ? char(int('x')+ j) : char(int('x')+2-j));
	  if (mi[j]>1) mexPrintf("^%d", int(mi[j]));
	}
      first = false; ++n;
    }
  }
  if (n == 0) mexPrintf("0");
  mexPrintf("\n");
}

/*@GFDOC
  @ARGS{@tpoly P}
  Performs various operations on the polynom POLY.
@*/

void gf_poly(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG("Wrong number of input arguments");
  }
  std::string cmd = in.pop().to_string();
  getfemint_poly *pp = &(to_poly_object(in.pop())->p);

  if (check_cmd(cmd, "print", in, out, 0, 0, 0, 0)) {
    /*@FUNC ('print')
      Prints the content of P.
      @*/
    print_poly(pp->p);
  } else if (check_cmd(cmd, "product", in, out, 0, 0, 0, 0)) {
    /*@FUNC ('product')
      To be done ... !
    @*/
    mexPrintf("to be done!\n");
  } else bad_cmd(cmd);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  catch_errors(nlhs, plhs, nrhs, prhs, gf_poly, "gf_poly");
}
