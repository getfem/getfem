// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_mesh_level_set.h : description of a mesh cut by a
//           number of levelsets.
//           
// Date    : March 04, 2005.
// Author  : Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//           Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2005-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================


#ifndef GETFEM_MESH_LEVEL_SET_H__
#define GETFEM_MESH_LEVEL_SET_H__

#include <getfem_integration.h>
#include <getfem_level_set.h>
#include <getfem_fem.h>


namespace getfem {

  class mesh_level_set : public getfem_mesh_receiver,
			 public context_dependencies {
  public:
    typedef std::string subzone;
    typedef std::set<const subzone *> zone;
    typedef std::set<const zone*> zoneset;
  protected :

    mutable std::set<subzone> allsubzones;

    mutable std::set<zone> allzones;

    

    dal::dynamic_array<const std::string *> zones_of_convexes;
    getfem_mesh *linked_mesh_;
    mutable bool is_valid_, is_adapted_;

    typedef level_set *plevel_set;
    std::vector<plevel_set> level_sets; // set of level set
    
    typedef boost::intrusive_ptr<getfem_mesh> pgetfem_mesh;

    struct convex_info {
      pgetfem_mesh pmesh;
      zoneset zones;
      convex_info() : pmesh(0) {}
    };

    std::map<size_type, convex_info> cut_cv;

  public :
    size_type nb_level_sets(void) const { return level_sets.size(); }
    plevel_set get_level_set(size_type i) const { return level_sets[i]; }
    bool is_valid() const { return is_valid_; }
    void update_from_context(void) const { is_adapted_= false; }
    bool is_convex_cut(size_type i) const
    { return (cut_cv.find(i) != cut_cv.end()); }
    const getfem_mesh& mesh_of_convex(size_type i) const {
      if (is_convex_cut(i)) return *((cut_cv.find(i))->second.pmesh);
      DAL_THROW(failure_error, "This element is not cut !");
    }
    
    /// Gives a reference to the linked mesh of type getfem\_mesh.
    getfem_mesh &linked_mesh(void) const { return *linked_mesh_; }
    void clear(void);
    /* explicit calls to parent class 
       for HP aCC and mipspro CC who complain about hidden functions 
       (they're right)
    */
    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_DELETE &);
    void receipt(const MESH_ADD_CONVEX &m) {getfem_mesh_receiver::receipt(m); }
    void receipt(const MESH_SUP_CONVEX &m);
    void receipt(const MESH_SWAP_CONVEX &m);
    
    size_type memsize() const {
      size_type res = sizeof(mesh_level_set)
	+ level_sets.size() * sizeof(plevel_set);
      for (std::map<size_type, convex_info>::const_iterator it=cut_cv.begin();
	   it != cut_cv.end(); ++it) {
	res += sizeof(convex_info)
	  + it->second.pmesh->memsize()
	  + it->second.zones.size()
	  * (level_sets.size() + sizeof(std::string *) + sizeof(std::string));
      }
      return res;
    }

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
    void global_cut_mesh(getfem_mesh &m) const;
    /** do all the work (cut the convexes wrt the levelsets) */
    void adapt(void);
    void merge_zoneset(zoneset &zones1, const zoneset &zones2) const;
    void merge_zoneset(zoneset &zones1, const std::string &subz) const;
    bool convex_is_cut(size_type cv) const
    { return (cut_cv.find(cv) != cut_cv.end()); }
    const std::string &primary_zone_of_convex(size_type cv) const
    { return *(zones_of_convexes[cv]); }
    const zoneset &zoneset_of_convex(size_type cv) const {
      std::map<size_type, convex_info>::const_iterator it = cut_cv.find(cv);
      if (it != cut_cv.end()) return (*it).second.zones;
      DAL_THROW(internal_error, "You cannot call this function for uncut convexes");
    }
    
    mesh_level_set(getfem_mesh &me);
    virtual ~mesh_level_set();
    

  private:
    mesh_level_set(const mesh_level_set &);
    mesh_level_set & operator=(const mesh_level_set &);
    void cut_element(size_type cv, const dal::bit_vector &primary,
		     const dal::bit_vector &secondary);     
    int is_not_crossed_by(size_type c, plevel_set ls, unsigned lsnum);
    int sub_simplex_is_not_crossed_by(size_type cv, plevel_set ls,
				      size_type sub_cv);
    void find_zones_of_element(size_type cv, std::string &prezone);

    /** For each levelset, if the convex cv is crossed, add the levelset number
	into 'prim' (and 'sec' is the levelset has a secondary part).
	zone is also filled with '0', '+', and '-'.
    */
    void find_crossing_level_set(size_type cv, 
				 dal::bit_vector &prim, 
				 dal::bit_vector &sec, std::string &zone);
    void run_delaunay(std::vector<base_node> &fixed_points,
		      gmm::dense_matrix<size_type> &simplexes,
		      std::vector<dal::bit_vector> &fixed_points_constraints);
  };

  void getfem_mesh_level_set_noisy(void);

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MESH_LEVEL_SET_H__  */
