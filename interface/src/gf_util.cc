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

#include <getfemint.h>
#include <getfemint_misc.h>
#include <gmm/gmm_inoutput.h>
#include <getfemint_gsparse.h>

using namespace getfemint;

/*@GFDOC
  Performs various operations which do not fit elsewhere.
@*/




// Object for the declaration of a new sub-command.

struct sub_gf_util : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out) = 0;
};

typedef std::shared_ptr<sub_gf_util> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_util {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out)			\
      { dummy_func(in); dummy_func(out); code }				\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }






void gf_util(getfemint::mexargs_in& m_in, getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {


    /*@FUNC ('save matrix', @str FMT, @str FILENAME, @mat A)
    Exports a sparse matrix into the file named FILENAME, using
    Harwell-Boeing (FMT='hb') or Matrix-Market (FMT='mm') formatting. @*/
    sub_command
      ("save matrix", 3, 3, 0, 0,
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

       );


    /*@FUNC A = ('load matrix', @str FMT, @str FILENAME)
     Imports a sparse matrix from a file. @*/
    sub_command
      ("load matrix", 2, 2, 1, 1,
       spmat_load(in, out, mexarg_out::USE_NATIVE_SPARSE);
       );


    /*@FUNC tl = ('trace level' [, @int level])
      Set the verbosity of some GetFEM routines.

      Typically the messages printed by the model bricks, 0 means no
      trace message (default is 3). if no level is given,
      the current trace level is returned. @*/
    sub_command
      ("trace level", 0, 1, 0, 1,
       if (in.remaining())
	 gmm::set_traces_level(in.pop().to_integer(0, 100));
       else
	 out.pop().from_integer(int(gmm::traces_level::level()));
       );


    /*@FUNC tl = ('warning level', @int level)
      Filter the less important warnings displayed by getfem.

      0 means no warnings, default level is 3. if no level is given,
      the current warning level is returned. @*/
    sub_command
      ("warning level", 0, 1, 0, 1,
       if (in.remaining())
	 gmm::set_warning_level(in.pop().to_integer(0, 100));
       else
	 out.pop().from_integer(int(gmm::warning_level::level()));
       );

  }


  if (m_in.narg() < 1)  THROW_BADARG("Wrong number of input arguments");

  std::string init_cmd  = m_in.pop().to_string();
  std::string cmd       = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out);
  }
  else bad_cmd(init_cmd);

}
