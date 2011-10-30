// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2011-2011 Yves Renard.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

#include <getfemint_misc.h>
#include <getfemint_workspace.h>
#include <getfemint_mesher_object.h>
#include <getfem/getfem_mesher.h>





using namespace getfemint;

/*@GFDOC
  This object represents a geometric object to be meshed by the (very)
  experimental mesher of Getfem.
@*/


// Object for the declaration of a new sub-command.

typedef getfemint_mesher_object *pgetfemint_mesher_object;

struct sub_mesher_object : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   pgetfemint_mesher_object &pmo) = 0;
};

typedef boost::intrusive_ptr<sub_mesher_object> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_mesher_object {		       		\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       pgetfemint_mesher_object &pmo)			\
      { dummy_func(in); dummy_func(out); code }				\
    };									\
    psub_command psubc = new subc;					\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }



void gf_mesher_object(getfemint::mexargs_in& m_in,
		 getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {
    
    /*@INIT MF = ('ball', @vec center, @scalar radius)
      Represent a ball of corresponding center and radius.
      @*/
    sub_command
      ("ball", 2, 2, 0, 1,
       darray center = in.pop().to_darray();
       double radius = in.pop().to_scalar();

       getfem::base_node bncenter(gmm::vect_size(center));
       gmm::copy(center, bncenter);

       getfem::mesher_signed_distance *ball
         = new getfem::mesher_ball(bncenter, radius);

       pmo = getfemint_mesher_object::get_from(ball);
       );

    /*@INIT MF = ('intersect', @tmo object1 , @tmo object2, ...)
      Intersection of several objects.
      @*/
    sub_command
      ("intersect", 2, 100, 0, 1,

       std::vector<const getfem::mesher_signed_distance *> vd;

       const getfem::mesher_signed_distance *psd
          = in.pop().to_const_mesher_object();

       vd.push_back(psd);

       while (in.remaining()) {
	 psd = in.pop().to_const_mesher_object();
	 vd.push_back(psd);
       }
       
       getfem::mesher_signed_distance *psd2 
         = new getfem::mesher_intersection(vd);
       pmo = getfemint_mesher_object::get_from(psd2);
       );

    /*@INIT MF = ('union', @tmo object1 , @tmo object2, ...)
      Union of several objects.
      @*/
    sub_command
      ("union", 2, 100, 0, 1,

       std::vector<const getfem::mesher_signed_distance *> vd;

       const getfem::mesher_signed_distance *psd
          = in.pop().to_const_mesher_object();

       vd.push_back(psd);

       while (in.remaining()) {
	 psd = in.pop().to_const_mesher_object();
	 vd.push_back(psd);
       }
       
       getfem::mesher_signed_distance *psd2 = new getfem::mesher_union(vd);
       pmo = getfemint_mesher_object::get_from(psd2);
       );

    /*@INIT MF = ('set minus', @tmo object1 , @tmo object2)
      Geometric object being object1 minus object2.
      @*/
    sub_command
      ("set minus", 2, 100, 0, 1,

       std::vector<const getfem::mesher_signed_distance *> vd;

       const getfem::mesher_signed_distance *psd1
          = in.pop().to_const_mesher_object();
       const getfem::mesher_signed_distance *psd2
          = in.pop().to_const_mesher_object();
       
       getfem::mesher_signed_distance *psd
         = new getfem::mesher_setminus(*psd1, *psd2);
       pmo = getfemint_mesher_object::get_from(psd);
       );
  }


  if (m_in.narg() < 1) THROW_BADARG("Wrong number of input arguments");
  getfemint_mesher_object *pmo = NULL;
  
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);
  

  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out, pmo);
  }
  else bad_cmd(init_cmd);

 
  m_out.pop().from_object_id(pmo->get_id(), MESHER_OBJECT_CLASS_ID);
}
