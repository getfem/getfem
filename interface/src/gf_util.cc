// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2006-2006 Yves Renard, Julien Pommier.
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

#include <getfemint.h>
#include <getfemint_misc.h>
#include <gmm/gmm_inoutput.h>
#include <getfemint_gsparse.h>
#include <getfemint_gsparse_misc.h>

using namespace getfemint;

/*MLABCOM
  FUNCTION gf_util(operation, [,args])

  Performs various operations which do not fit elsewhere..

  * gf_util('save matrix', string FMT, string FILENAME, spmat A);

  Exports a sparse matrix into the file named FILENAME, using
  Harwell-Boeing (FMT='hb') or Matrix-Market (FMT='mm') formatting.

  * A = gf_util('load matrix', string FMT, string FILENAME);

  Imports a sparse matrix from a file.

  @FUNC ::UTIL('trace level')

  @FUNC ::UTIL('warning level')
MLABCOM*/

void gf_util(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG("Wrong number of input arguments");
  }
  std::string cmd = in.pop().to_string();

  if (check_cmd(cmd, "save matrix", in, out, 3, 3, 0, 0)) {
    std::string fmt = in.pop().to_string();
    int ifmt;
    if (cmd_strmatch(fmt, "hb") || cmd_strmatch(fmt, "harwell-boeing")) ifmt = 0;
    else if (cmd_strmatch(fmt, "mm") || cmd_strmatch(fmt, "matrix-market")) ifmt = 1;
    else THROW_BADARG("unknown sparse matrix file-format : " << fmt);
    std::string fname = in.pop().to_string();
    if (!in.front().is_complex()) {
      gf_real_sparse_csc_const_ref H;  in.pop().to_sparse(H);
      gmm::csc_matrix<double> cscH; gmm::copy(H,cscH);
      if (ifmt == 0) gmm::Harwell_Boeing_save(fname.c_str(), cscH);
      else           gmm::MatrixMarket_save(fname.c_str(), cscH);
    } else {
      gf_cplx_sparse_csc_const_ref H;  in.pop().to_sparse(H);
      gmm::csc_matrix<complex_type> cscH; gmm::copy(H,cscH);
      if (ifmt == 0) gmm::Harwell_Boeing_save(fname.c_str(), cscH);
      else           gmm::MatrixMarket_save(fname.c_str(), cscH);      
    }
  } else if (check_cmd(cmd, "load matrix", in, out, 2, 2, 1, 1)) {
    spmat_load(in, out, mexarg_out::USE_NATIVE_SPARSE);
  } else if (check_cmd(cmd, "trace level", in, out, 1, 1, 0, 0)) {
    /*@FUNC ::UTIL('trace level', @int level)
      Set the verbosity of some getfem++ routines (typically the messages
      printed by the model bricks), 0 means no trace message (default is 3).
      @*/
    gmm::set_traces_level(in.pop().to_integer(0, 100));
  } else if (check_cmd(cmd, "warning level", in, out, 1, 1, 0, 0)) {
    /*@FUNC ::UTIL('warning level', @int level)       
      Filter the less important warnings displayed by getfem. 0 means no
      warnings, default level is 3.
      @*/
    gmm::set_warning_level(in.pop().to_integer(0, 100));
  } else bad_cmd(cmd);
}

