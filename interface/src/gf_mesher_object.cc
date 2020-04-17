/*===========================================================================

 Copyright (C) 2011-2020 Yves Renard.

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
#include <getfemint_workspace.h>
#include <getfem/getfem_mesher.h>


using namespace getfemint;

/*@GFDOC
  This object represents a geometric object to be meshed by the
  experimental meshing procedure of Getfem.
@*/


// Object for the declaration of a new sub-command.

struct sub_mesher_object : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   getfem::pmesher_signed_distance &psd) = 0;
};

typedef std::shared_ptr<sub_mesher_object> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_mesher_object {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       getfem::pmesher_signed_distance &psd) override	\
      { dummy_func(in); dummy_func(out); code }				\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
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
      Represents a ball of corresponding center and radius.
      @*/
    sub_command
      ("ball", 2, 2, 0, 1,
       darray center = in.pop().to_darray();
       double radius = in.pop().to_scalar();

       getfem::base_node bncenter(gmm::vect_size(center));
       gmm::copy(center, bncenter);

       psd = getfem::new_mesher_ball(bncenter, radius);
       );

    /*@INIT MF = ('half space', @vec origin, @vec normal_vector)
      Represents an half space delimited by the plane which contains the
      origin and normal to `normal_vector`. The selected part is the part
      in the direction of the normal vector. This allows to cut a geometry
      with a plane for instance to build a polygon or a polyhedron.
      @*/
    sub_command
      ("half space", 2, 2, 0, 1,
       darray origin = in.pop().to_darray();
       darray n = in.pop().to_darray();

       getfem::base_node bnorigin(gmm::vect_size(origin));
       gmm::copy(origin, bnorigin);
       getfem::base_node bnn(gmm::vect_size(n));
       gmm::copy(n, bnn);

       psd = getfem::new_mesher_half_space(bnorigin, bnn);
       );

    /*@INIT MF = ('cylinder', @vec origin, @vec n, @scalar length, @scalar radius)
      Represents a cylinder (in any dimension) of a certain radius whose axis
      is determined by the origin, a vector `n` and a certain length.
      @*/
    sub_command
      ("cylinder", 4, 4, 0, 1,
       darray origin = in.pop().to_darray();
       darray n = in.pop().to_darray();
       double length = in.pop().to_scalar();
       double radius = in.pop().to_scalar();
       
       getfem::base_node bnorigin(gmm::vect_size(origin));
       gmm::copy(origin, bnorigin);
       getfem::base_node bnn(gmm::vect_size(n));
       gmm::copy(n, bnn);

       psd = getfem::new_mesher_cylinder(bnorigin, bnn, length, radius);
       );

    /*@INIT MF = ('cone', @vec origin, @vec n, @scalar length, @scalar half_angle)
      Represents a cone (in any dimension) of a certain half-angle (in radians)
      whose axis is determined by the origin, a vector `n` and a certain length.
      @*/
    sub_command
      ("cone", 4, 4, 0, 1,
       darray origin = in.pop().to_darray();
       darray n = in.pop().to_darray();
       double length = in.pop().to_scalar();
       double half_angle = in.pop().to_scalar();
       
       getfem::base_node bnorigin(gmm::vect_size(origin));
       gmm::copy(origin, bnorigin);
       getfem::base_node bnn(gmm::vect_size(n));
       gmm::copy(n, bnn);

       psd = getfem::new_mesher_cone(bnorigin, bnn, length, half_angle);
       );

    /*@INIT MF = ('torus', @scalar R, @scalar r)
      Represents a torus in 3d of axis along the z axis with a great radius
      equal to `R` and small radius equal to `r`. For the moment, the
      possibility to change the axis is not given.
      @*/
    sub_command
      ("torus", 2, 2, 0, 1,
       double R = in.pop().to_scalar();
       double r = in.pop().to_scalar();
       psd  = getfem::new_mesher_torus(R, r);
       );

    
    /*@INIT MF = ('rectangle', @vec rmin, @vec rmax)
      Represents a rectangle (or parallelepiped in 3D) parallel to the axes.
      @*/
    sub_command
      ("rectangle", 2, 2, 0, 1,
       darray rmin = in.pop().to_darray();
       darray rmax = in.pop().to_darray();
       
       size_type N = gmm::vect_size(rmin);
       GMM_ASSERT1(N == gmm::vect_size(rmax),
		   "Extreme points should be the same lenght");

       getfem::base_node rrmin(N); getfem::base_node rrmax(N);
       gmm::copy(rmin, rrmin); gmm::copy(rmax, rrmax);

       psd = getfem::new_mesher_rectangle(rrmin, rrmax);
       );


    /*@INIT MF = ('intersect', @tmo object1 , @tmo object2, ...)
      Intersection of several objects.
      @*/
    sub_command
      ("intersect", 2, 100, 0, 1,
       std::vector<getfem::pmesher_signed_distance> vd;
       vd.push_back(to_mesher_object(in.pop()));
       while (in.remaining())
	 vd.push_back(to_mesher_object(in.pop()));
       
       psd = getfem::new_mesher_intersection(vd);
       );

    /*@INIT MF = ('union', @tmo object1 , @tmo object2, ...)
      Union of several objects.
      @*/
    sub_command
      ("union", 2, 100, 0, 1,
       std::vector<getfem::pmesher_signed_distance> vd;
       vd.push_back(to_mesher_object(in.pop()));
       while (in.remaining())
	 vd.push_back(to_mesher_object(in.pop()));
       
       psd = getfem::new_mesher_union(vd);
       );

    /*@INIT MF = ('set minus', @tmo object1 , @tmo object2)
      Geometric object being object1 minus object2.
      @*/
    sub_command
      ("set minus", 2, 100, 0, 1,
       getfem::pmesher_signed_distance psd1 = to_mesher_object(in.pop());
       getfem::pmesher_signed_distance psd2 = to_mesher_object(in.pop());
       psd = getfem::new_mesher_setminus(psd1, psd2);
       );
  }

  if (m_in.narg() < 1) THROW_BADARG("Wrong number of input arguments");
  getfem::pmesher_signed_distance psd;
  
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);
  

  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out, psd);
  }
  else bad_cmd(init_cmd);

  id_type id = store_mesher_object(psd);
  m_out.pop().from_object_id(id, MESHER_OBJECT_CLASS_ID);
}
