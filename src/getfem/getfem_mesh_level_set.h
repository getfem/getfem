// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2005-2008 Julien Pommier
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
// As a special exception, you  may use  this file  as it is a part of a free
// software  library  without  restriction.  Specifically,  if   other  files
// instantiate  templates  or  use macros or inline functions from this file,
// or  you compile this  file  and  link  it  with other files  to produce an
// executable, this file  does  not  by itself cause the resulting executable
// to be covered  by the GNU Lesser General Public License.  This   exception
// does not  however  invalidate  any  other  reasons why the executable file
// might be covered by the GNU Lesser General Public License.
//
//===========================================================================

/**@file getfem_mesh_level_set.h
   @author Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @author Yves Renard <Yves.Renard@insa-toulouse.fr>
   @date March 04, 2005.
   @brief Keep informations about a mesh crossed by level-sets. 
*/
#ifndef GETFEM_MESH_LEVEL_SET_H__
#define GETFEM_MESH_LEVEL_SET_H__

#include "getfem_integration.h"
#include "getfem_level_set.h"
#include "getfem_fem.h"


namespace getfem {
  /** @brief Keep informations about a mesh crossed by level-sets.
      Cut convexes with respect to the level sets.

      Note that the cutting won't be conformal.
  */
  class mesh_level_set : public context_dependencies {
  public:
    typedef std::string subzone;
    typedef std::set<const subzone *> zone;
    typedef std::set<const zone*> zoneset;

  protected :

    mutable std::set<subzone> allsubzones;
    mutable std::set<zone> allzones;

    dal::dynamic_array<const std::string *> zones_of_convexes;
    mesh *linked_mesh_;
    mutable bool is_adapted_;

    typedef level_set *plevel_set;
    std::vector<plevel_set> level_sets; // set of level set
    
    typedef boost::intrusive_ptr<mesh> pmesh;

    struct convex_info {
      pmesh pmsh;
      zoneset zones;
      mesh_region ls_border_faces;
      convex_info() : pmsh(0) {}
    };

    std::map<size_type, convex_info> cut_cv;

    mutable dal::bit_vector crack_tip_convexes_;

  public :
    /// Get number of level-sets referenced in this object.
    size_type nb_level_sets(void) const { return level_sets.size(); }
    plevel_set get_level_set(size_type i) const { return level_sets[i]; }
    void update_from_context(void) const { is_adapted_= false; }
    bool is_convex_cut(size_type i) const
    { return (cut_cv.find(i) != cut_cv.end()); }
    const mesh& mesh_of_convex(size_type i) const {
      if (is_convex_cut(i)) return *((cut_cv.find(i))->second.pmsh);
      GMM_ASSERT1(false, "This element is not cut !");
    }
    
    const dal::bit_vector &crack_tip_convexes() const;

    /// Gives a reference to the linked mesh of type mesh.
    mesh &linked_mesh(void) const { return *linked_mesh_; }
    void clear(void);
    
    size_type memsize() const {
      size_type res = sizeof(mesh_level_set)
	+ level_sets.size() * sizeof(plevel_set);
      for (std::map<size_type, convex_info>::const_iterator it=cut_cv.begin();
	   it != cut_cv.end(); ++it) {
	res += sizeof(convex_info)
	  + it->second.pmsh->memsize()
	  + it->second.zones.size()
	  * (level_sets.size() + sizeof(std::string *) + sizeof(std::string));
      }
      return res;
    }
    /** add a new level set. Only a reference is kept, no copy done. */
    void add_level_set(level_set &ls) {
      if (std::find(level_sets.begin(), level_sets.end(), &ls)
	  == level_sets.end()) {
	level_sets.push_back(&ls); touch();
	is_adapted_ = false;
      }
    }
    void sup_level_set(level_set &ls) {
      std::vector<plevel_set>::iterator
	it = std::find(level_sets.begin(), level_sets.end(), &ls);
      if (it != level_sets.end()) {
	level_sets.erase(it);
	is_adapted_ = false;
	touch();
      }
    }

    /** fill m with the (non-conformal) "cut" mesh. */
    void global_cut_mesh(mesh &m) const;
    /** do all the work (cut the convexes wrt the levelsets) */
    void adapt(void);
    void merge_zoneset(zoneset &zones1, const zoneset &zones2) const;
    void merge_zoneset(zoneset &zones1, const std::string &subz) const;
    const std::string &primary_zone_of_convex(size_type cv) const
    { return *(zones_of_convexes[cv]); }
    const zoneset &zoneset_of_convex(size_type cv) const {
      std::map<size_type, convex_info>::const_iterator it = cut_cv.find(cv);
      if (it != cut_cv.end()) return (*it).second.zones;
      GMM_ASSERT1(false, "You cannot call this function for uncut convexes");
    }
    
    void init_with_mesh(mesh &me);
    mesh_level_set(mesh &me);
    mesh_level_set(void);
    virtual ~mesh_level_set();
    mesh_level_set(const mesh_level_set &mls) {
      GMM_ASSERT1(linked_mesh_ == 0 && mls.linked_mesh_ == 0,
		  "Copy constructor is not allowed for mesh_level_set");
    }
    mesh_level_set & operator=(const mesh_level_set &mls) {
      GMM_ASSERT1(linked_mesh_ == 0 && mls.linked_mesh_ == 0,
		  "Copy operator is not allowed for mesh_level_set");
      return *this;
    }


  private:
    void cut_element(size_type cv, const dal::bit_vector &primary,
		     const dal::bit_vector &secondary, scalar_type radius);
    int is_not_crossed_by(size_type c, plevel_set ls, unsigned lsnum,
			  scalar_type radius);
    int sub_simplex_is_not_crossed_by(size_type cv, plevel_set ls,
				      size_type sub_cv, scalar_type radius);
    void find_zones_of_element(size_type cv, std::string &prezone,
			       scalar_type radius);

    /** For each levelset, if the convex cv is crossed, add the levelset number
	into 'prim' (and 'sec' is the levelset has a secondary part).
	zone is also filled with '0', '+', and '-'.
    */
    void find_crossing_level_set(size_type cv, 
				 dal::bit_vector &prim, 
				 dal::bit_vector &sec, std::string &zone,
				 scalar_type radius);
    void run_delaunay(std::vector<base_node> &fixed_points,
		      gmm::dense_matrix<size_type> &simplexes,
		      std::vector<dal::bit_vector> &fixed_points_constraints);
    
    void update_crack_tip_convexes();
  };

  void getfem_mesh_level_set_noisy(void);

  std::ostream &operator<<(std::ostream &os, const mesh_level_set::zone &z);
  std::ostream &operator<<(std::ostream &os, const mesh_level_set::zoneset &zs);

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MESH_LEVEL_SET_H__  */
