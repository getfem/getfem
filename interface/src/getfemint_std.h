/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2006-2020 Yves Renard, Julien Pommier.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

#ifndef GETFEMINT_STD_H__
#define GETFEMINT_STD_H__

#include <getfem/bgeot_config.h>
#include <getfem/dal_backtrace.h>
#include <gfi_array.h>

namespace getfemint
{  
  using std::endl; using std::cout; using std::cerr;
  using std::ends; using std::cin;

  void attach_gdb();

  #define THROW_INTERNAL_ERROR { dal::dump_backtrace(); GMM_THROW(getfemint::getfemint_error, "getfem-interface: internal error\n") }

  typedef size_t size_type;
  typedef bgeot::short_type short_type;
  typedef dal::uint32_type id_type;

  class getfemint_error : public std::logic_error {
    //    std::string what_;
  public:
    getfemint_error(const std::string& what_arg) : std::logic_error (what_arg) //what_(what_arg) //
      { }
    //    const char * what() const { return what_.c_str(); }
  };


  class getfemint_interrupted : public getfemint_error {
  public:
    getfemint_interrupted() : getfemint_error("") {}
  };

  class getfemint_bad_arg : public getfemint_error {
  public:
    getfemint_bad_arg(const std::string& what_) : getfemint_error(what_) {}
      //    ~getfemint_bad_arg() throw() {} /* pourquoi g++-3.0 demande de definir cette fonction ? */
  };

  /* no callback for these exceptions */
#define THROW_ERROR(thestr) {                \
    std::stringstream msg;                   \
    msg << thestr << ends;                   \
    throw getfemint_error(msg.str());        \
  }
#define THROW_BADARG(thestr) {               \
    std::stringstream msg;                   \
    msg << thestr << ends;                   \
    throw getfemint_bad_arg(msg.str());      \
  }

#define THROW_INTERRUPTED() {		     \
    throw getfemint_interrupted();	     \
  }
  
#define GFI_WARNING(thestr) { infomsg() << "WARNING: " << thestr; }
  std::ostream& infomsg(); // defined in getfem_matlab.C

  /* see getfem_interface.C */
  struct config {
    static config *cfg;
    gfi_interface_type interface_type_;
    int base_index_; /* base indexing of arrays (matlab starts at 1, python at 0 */
    bool can_return_integer_; /* matlab < 7 is brain-damaged with respect to int32 type */
    bool has_native_sparse_; /* python has no builtin sparse matrix type */
    bool prefer_native_sparse_;
    bool has_1D_arrays_; /* true if 1D arrays do exist (for example python),
                           false if they do not existe (i.e. in matlab everything is at least a matrix) */
    const char *current_function_;
    static int base_index() { return cfg->base_index_; }
    static bool has_native_sparse() { return cfg->has_native_sparse_; }
    static bool prefer_native_sparse() { return cfg->prefer_native_sparse_; }
    static bool can_return_integer() { return cfg->can_return_integer_; }
    static bool has_1D_arrays() { return cfg->has_1D_arrays_; }
    static std::string current_function() { return std::string(cfg->current_function_); } 
    static void set_current_config(config *p) { cfg = p; }
    config(gfi_interface_type);
  };
}
#endif
