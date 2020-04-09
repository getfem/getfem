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
#include <getfem/bgeot_geometric_trans.h>

using namespace getfemint;

/*@GFDOC
    General function for querying information about geometric transformations
    objects.
@*/

// Object for the declaration of a new sub-command.

struct sub_gf_geotrans : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   const bgeot::pgeometric_trans &pgt) = 0;
};

typedef std::shared_ptr<sub_gf_geotrans> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_geotrans {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       const bgeot::pgeometric_trans &pgt)		\
      { dummy_func(in); dummy_func(out); code }				\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }                           





void gf_geotrans_get(getfemint::mexargs_in& m_in,
		     getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {

    /*@RDATTR d = ('dim')
      Get the dimension of the @tgt.
      
      This is the dimension of the source space, i.e. the dimension of
      the reference convex.@*/
    sub_command
      ("dim", 0, 0, 0, 1,
       out.pop().from_scalar(pgt->dim());
       );


    /*@RDATTR b = ('is_linear')
      Return 0 if the @tgt is not linear.@*/
    sub_command
      ("is_linear", 0, 0, 0, 1,
       out.pop().from_scalar(pgt->is_linear() ? 1. : 0.);
       );


    /*@RDATTR n = ('nbpts')
      Return the number of points of the @tgt.@*/
    sub_command
      ("nbpts", 0, 0, 0, 1,
       out.pop().from_scalar(double(pgt->nb_points()));
       );


    /*@GET P = ('pts')
      Return the reference convex points of the @tgt.

      The points are stored in the columns of the output matrix.@*/
    sub_command
      ("pts", 0, 0, 0, 1,
       out.pop().from_vector_container(pgt->convex_ref()->points());
       );


    /*@GET N = ('normals')
      Get the normals for each face of the reference convex of the @tgt.

      The normals are stored in the columns of the output matrix.@*/
    sub_command
      ("normals", 0, 0, 0, 1,
       out.pop().from_vector_container(pgt->normals());
       );


    /*@GET Pt = ('transform',@mat G, @mat Pr)
      Apply the @tgt to a set of points.
      
      `G` is the set of vertices of the real convex, `Pr` is the set
      of points (in the reference convex) that are to be transformed.
      The corresponding set of points in the real convex is returned.@*/
    sub_command
      ("transform", 2, 2, 0, 1,
       getfem::base_matrix G = in.pop().to_darray(-1, -1).row_col_to_bm();
       darray pts = in.pop().to_darray(pgt->dim(), -1);
       darray w = out.pop().create_darray(unsigned(G.nrows()), pts.getn());
       for (unsigned i=0; i < pts.getn(); ++i) {
	 getfem::base_node P = pgt->transform(pts.col_to_bn(i), G);
	 for (size_type k=0; k < P.size(); ++k) w(k,i) = P[k];
       }
       );

    
    /*@GET s = ('char')
      Output a (unique) string representation of the @tgt.

      This can be used to perform comparisons between two
      different @tgt objects. @*/
    sub_command
      ("char", 0, 0, 0, 1,
       std::string s = bgeot::name_of_geometric_trans(pgt);
       out.pop().from_string(s.c_str());
       );


    /*@GET ('display')
    displays a short summary for a @tgt object.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       infomsg() << "gfGeoTrans object " << bgeot::name_of_geometric_trans(pgt)
       << " in dimension " << int(pgt->dim())
       << ", with " << pgt->nb_points() << " points \n";
       );

  }

  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  bgeot::pgeometric_trans pgt = to_geotrans_object(m_in.pop());
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out, pgt);
  }
  else bad_cmd(init_cmd);


}
