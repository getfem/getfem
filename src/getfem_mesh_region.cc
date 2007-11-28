// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2005-2008 Yves Renard
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

#include "getfem/getfem_mesh_region.h"
#include "getfem/getfem_mesh.h"

namespace getfem {
  typedef mesh_region::face_bitset face_bitset;

  void mesh_region::touch_parent_mesh() {
    if (parent_mesh) {
      parent_mesh->touch_from_region(id_);
    }
  }

  const mesh_region& mesh_region::from_mesh(const mesh &m) const {
    if (!p.get()) {
      mesh_region *r = const_cast<mesh_region*>(this);
      if (id_ == size_type(-1)) {
	r->p.reset(new impl); 
	r->add(m.convex_index());
      } else if (id_ != size_type(-2))
	r->p = m.region(id_).p;
    }
    /* TODO? : verifier que la liste des convexes est bien inclue dans m.convex_index */
    return *this;
  }

  face_bitset mesh_region::operator[](size_t cv) const {
    map_t::const_iterator it = rp().m.find(cv);
    if (it != rp().m.end()) return (*it).second;
    else return face_bitset();
  }

  /* may be optimized .. */
  const dal::bit_vector &mesh_region::index() const {
    rp().index_.clear();
    for (map_t::const_iterator it = rp().m.begin(); 
	 it != rp().m.end(); ++it) {
      if (it->second.any()) rp().index_.add(it->first);
    }
    return rp().index_;
  }

  void mesh_region::add(const dal::bit_vector &bv) {
    for (dal::bv_visitor i(bv); !i.finished(); ++i) add(i);
  }

  void mesh_region::add(size_type cv, size_type f) {
    wp().m[cv].set(f+1,1); 
    touch_parent_mesh();
  }

  void mesh_region::sup_all(size_type cv) { 
    map_t::iterator it = wp().m.find(cv);
    if (it != wp().m.end()) {
      wp().m.erase(it);
      touch_parent_mesh();
    }
  }

  void mesh_region::sup(size_type cv, size_type f) { 
    map_t::iterator it = wp().m.find(cv);
    if (it != wp().m.end()) {
      (*it).second.set(f+1,0); 
      if ((*it).second.none()) wp().m.erase(it); 
      touch_parent_mesh();
    }
  }

  void mesh_region::clear() { 
    wp().m.clear(); wp().index_.clear(); touch_parent_mesh();
  }

  void mesh_region::clean() {
    for (map_t::iterator it = wp().m.begin(), itn; 
	 it != wp().m.end(); it = itn) {
      itn = it; ++itn;
      if (!(*it).second.any()) {
	wp().m.erase(it);
      }
    }
    touch_parent_mesh();
  }

  void mesh_region::swap_convex(size_type cv1, size_type cv2) {
    map_t::iterator it1 = wp().m.find(cv1), it2 = wp().m.find(cv2),
      ite = wp().m.end();
    face_bitset f1, f2;
    
    if (it1 != ite) f1 = it1->second;
    if (it2 != ite) f2 = it2->second;
    if (!f1.none()) wp().m[cv2] = f1;
    else if (it2 != ite) wp().m.erase(it2);
    if (!f2.none()) wp().m[cv1] = f2;
    else if (it1 != ite) wp().m.erase(it1);
    touch_parent_mesh();
  }
  
  bool mesh_region::is_in(size_type cv, size_type f) const {
    map_t::const_iterator it = rp().m.find(cv);
    if (it == rp().m.end() || f+1 >= MAX_FACES_PER_CV) return false;
    return ((*it).second)[f+1];
  }

  bool mesh_region::is_empty() const {
    return rp().m.empty();
  }

  bool mesh_region::is_only_convexes() const { 
    return is_empty() || 
      (and_mask()[0] == true && and_mask().count() == 1); 
  }

  bool mesh_region::is_only_faces() const { 
    return is_empty() || (and_mask()[0] == false); 
  }
  
  face_bitset mesh_region::faces_of_convex(size_type cv) const {
    map_t::const_iterator it = rp().m.find(cv);
    if (it != rp().m.end()) return ((*it).second) >> 1; 
    else return face_bitset();
  }

  face_bitset mesh_region::and_mask() const {
    face_bitset bs; 
    if (rp().m.empty()) return bs;
    bs.set();
    for (map_t::const_iterator it = rp().m.begin(); it != rp().m.end(); ++it)
      if ((*it).second.any())
	bs &= (*it).second;
    return bs;
  }

  size_type mesh_region::size() const {
    size_type sz=0;
    for (map_t::const_iterator it = rp().m.begin(); it != rp().m.end(); ++it)
      sz += (*it).second.count();
    return sz;
  }

  mesh_region mesh_region::intersection(const mesh_region &a, 
					const mesh_region &b) {
    
    mesh_region r;
    /* we do not allow the "all_convexes" kind of regions
       for these operations as there are not intended to be manipulated
       (they only exist to provide a default argument to the mesh_region
       parameters of assembly procedures etc. */
    GMM_ASSERT1(a.id() != all_convexes().id() &&
		b.id() != all_convexes().id(), "the 'all_convexes' regions "
		"are not supported for set operations");
    map_t::const_iterator 
      ita = a.rp().m.begin(), enda = a.rp().m.end(),
      itb = b.rp().m.begin(), endb = b.rp().m.end();
    while (ita != enda && itb != endb) {
      if (ita->first < itb->first) ++ita;
      else if (ita->first > itb->first) ++itb;
      else {
        face_bitset maska = ita->second, maskb = itb->second, bs;
        if (maska[0] && !maskb[0]) bs = maskb;
        else if (maskb[0] && !maska[0]) bs = maska;
        else bs = maska & maskb;
        if (bs.any()) r.wp().m.insert(r.wp().m.end(), std::make_pair(ita->first,bs));
        ++ita; ++itb;
      }
    }
    return r;
  }

  mesh_region mesh_region::merge(const mesh_region &a, 
				 const mesh_region &b) {
    mesh_region r;
    GMM_ASSERT1(a.id() != all_convexes().id() &&
		b.id() != all_convexes().id(), "the 'all_convexes' regions "
		"are not supported for set operations");
    r.wp() = a.rp();
    for (map_t::const_iterator it = b.rp().m.begin(); 
         it != b.rp().m.end(); ++it) {
      r.wp().m[it->first] |= it->second;
    }
    return r;
  }

  mesh_region mesh_region::substract(const mesh_region &a,
                                     const mesh_region &b) {
    mesh_region r;
    GMM_ASSERT1(a.id() != all_convexes().id() &&
		b.id() != all_convexes().id(), "the 'all_convexes' regions "
		"are not supported for set operations");
    r.wp() = a.rp();
    for (map_t::const_iterator itb = b.rp().m.begin(); 
         itb != b.rp().m.end(); ++itb) {
      size_type cv = itb->first;
      map_t::iterator it = r.wp().m.find(cv);
      if (it != r.wp().m.end()) {
        it->second &= ~(itb->second);
      }
    }
    return r;
  }

  void mesh_region::error_if_not_faces() const {
    GMM_ASSERT1(is_only_faces(), "Expecting a set of faces, not convexes");
  }
  void mesh_region::error_if_not_convexes() const {
    GMM_ASSERT1(is_only_convexes(), "Expecting a set of convexes, not faces");
  }
  void mesh_region::error_if_not_homogeneous() const {
    GMM_ASSERT1(is_only_faces() || is_only_convexes(), "Expecting a set "
		"of convexes or a set of faces, but not a mixed set");
  }

  mesh_region::visitor::visitor(const mesh_region &s, const mesh &m) : 
    cv_(size_type(-1)), f_(size_type(-1)), finished_(false) {
    s.from_mesh(m);
    init(s);
  }
  

  mesh_region::visitor::visitor(const mesh_region &s) :
    cv_(size_type(-1)), f_(size_type(-1)), finished_(false) {
    init(s);
  }

  void mesh_region::visitor::init(const mesh_region &s) {
    GMM_ASSERT1(&s != 0, "Attemps to use an invalid mesh_region "
		"(need to call 'from_mesh')");
    it = s.rp().m.begin(); ite = s.rp().m.end(); 
    next();
  }

  std::ostream & operator <<(std::ostream &os, const mesh_region &w) {
    if (w.id() == mesh_region::all_convexes().id())
      os << " ALL_CONVEXES";
    else 
      for (mr_visitor cv(w); !cv.finished(); cv.next()) {
	os << cv.cv();
	if (cv.is_face()) os << "/" << cv.f();
	os << " ";
      }
    return os;
  }
}
