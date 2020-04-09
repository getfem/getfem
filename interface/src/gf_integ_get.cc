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
#include <getfem/getfem_integration.h>

using namespace getfemint;

static void check_not_exact(getfem::pintegration_method im) {
  if (im->type() != getfem::IM_APPROX) THROW_ERROR("this has no meaning for exact integration methods");
}
/*@GFDOC
  General function for querying information about integration method objects.
@*/


// Object for the declaration of a new sub-command.

struct sub_gf_integ_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   const getfem::pintegration_method &im,
		   getfem::papprox_integration pai, size_type imdim) = 0;
};

typedef std::shared_ptr<sub_gf_integ_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_integ_get {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       const getfem::pintegration_method &im,		\
		       getfem::papprox_integration pai,			\
		       size_type imdim) {				\
	dummy_func(in); dummy_func(out); dummy_func(im);		\
	dummy_func(imdim); dummy_func(pai); code			\
      }	                                            			\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }                           



void gf_integ_get(getfemint::mexargs_in& m_in,
		  getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;
  
  if (subc_tab.size() == 0) {

    /*@RDATTR b = ('is_exact')
    Return 0 if the integration is an approximate one.@*/
    sub_command
      ("is_exact", 0, 0, 0, 1,
       out.pop().from_scalar(im->type() != getfem::IM_APPROX ? 1. : 0.);
       );


    /*@RDATTR d = ('dim')
    Return the dimension of the reference convex of
    the method.@*/
    sub_command
      ("dim", 0, 0, 0, 1,
       out.pop().from_scalar(double(imdim));
       );


    /*@RDATTR n = ('nbpts')
    Return the total number of integration points.

    Count the points for the volume integration, and points for
    surface integration on each face of the reference convex.

    Only for approximate methods, this has no meaning for exact
    integration methods!@*/
    sub_command
      ("nbpts", 0, 0, 0, 1,
       check_not_exact(im);
       iarray w = out.pop().create_iarray_h(1+pai->structure()->nb_faces());
       w[0] = int(pai->nb_points_on_convex());
       for (short_type i=0; i < pai->structure()->nb_faces(); ++i)
	 w[i+1] = int(pai->nb_points_on_face(i));
       );


    /*@GET Pp = ('pts')
      Return the list of integration points
      
      Only for approximate methods, this has no meaning for exact
      integration methods!@*/
    sub_command
      ("pts", 0, 0, 0, 1,
       check_not_exact(im);
       out.pop().from_vector_container(*(pai->pintegration_points()));
       );


    /*@GET Pf = ('face_pts',F)
      Return the list of integration points for a face.
      
      Only for approximate methods, this has no meaning for exact
      integration methods!@*/
    sub_command
      ("face_pts", 1, 1, 0, 1,
       check_not_exact(im);
       short_type nbf = pai->structure()->nb_faces();
       short_type f = short_type(in.pop().to_face_number(nbf));
       darray w = out.pop().create_darray(unsigned(imdim),
					  unsigned(pai->nb_points_on_face(f)));
       for (size_type j=0; j < pai->nb_points_on_face(f); ++j)
	 for (size_type i=0; i < imdim; ++i)
	   w(i,j)=pai->point_on_face(f,j)[i];
       );


    /*@GET Cp = ('coeffs')
    Returns the coefficients associated to each integration point.

    Only for approximate methods, this has no meaning for exact
    integration methods!@*/
    sub_command
      ("coeffs", 0, 0, 0, 1,
       check_not_exact(im);
       out.pop().from_dcvector
       (im->approx_method()->integration_coefficients());
       );


    /*@GET Cf = ('face_coeffs',F)
    Returns the coefficients associated to each integration of a face.

    Only for approximate methods, this has no meaning for exact
    integration methods!@*/
    sub_command
      ("face_coeffs", 1, 1, 0, 1,
       check_not_exact(im);
       short_type f =
         short_type(in.pop().to_face_number(pai->structure()->nb_faces()));
       darray w =
         out.pop().create_darray_h(unsigned(pai->nb_points_on_face(f)));
       for (size_type j=0; j < pai->nb_points_on_face(f); ++j)
	 w[j]=pai->coeff_on_face(f,j);
       );


    /*@GET s = ('char')
    Ouput a (unique) string representation of the integration method.

    This can be used to  comparisons between two different @tinteg
    objects.@*/
    sub_command
      ("char", 0, 0, 0, 1,
       std::string s = getfem::name_of_int_method(im);
       out.pop().from_string(s.c_str());
       );

    /*@GET ('display')
    displays a short summary for a @tinteg object.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       infomsg() << "gfInteg object " << getfem::name_of_int_method(im);
       if (im->type() != getfem::IM_APPROX)
	 infomsg() << "Exact method in dimension " << int(imdim) << endl;
       else
	 infomsg() << "Cubature method in dimension " << int(imdim)
		   << " with " << pai->nb_points_on_convex()
		   << " Gauss points \n";
       );
  
  }


  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");


  getfem::pintegration_method im = to_integ_object(m_in.pop());
  getfem::papprox_integration
    pai = im->type() == getfem::IM_APPROX ? im->approx_method() : 0;
  size_type imdim = 0;
  if (im->type() == getfem::IM_EXACT) imdim = im->exact_method()->dim();
  else if (im->type() == getfem::IM_APPROX) imdim = im->approx_method()->dim();
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out, im, pai, imdim);
  }
  else bad_cmd(init_cmd);

}
