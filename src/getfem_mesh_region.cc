#include <getfem_mesh_region.h>
#include <getfem_mesh.h>

namespace getfem {
  typedef mesh_region::face_bitset face_bitset;

  void mesh_region::touch_parent_mesh() {
    if (parent_mesh) {
      parent_mesh->touch_from_region(id_);
    }
  }

  const mesh_region& mesh_region::from_mesh(const getfem_mesh &m) const {
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

  void mesh_region::error_if_not_faces() const {
    if (!is_only_faces()) 
      DAL_THROW(dal::failure_error, 
		"Expecting a set of faces, not convexes");
  }
  void mesh_region::error_if_not_convexes() const {
    if (!is_only_convexes()) 
      DAL_THROW(dal::failure_error, 
		"Expecting a set of convexes, not faces");
  }
  void mesh_region::error_if_not_homogeneous() const {
    if (!is_only_faces() && !is_only_convexes()) 
      DAL_THROW(dal::failure_error, 
		"Expecting a set of convexes or a set of faces, but not a mixed set");
  }

  mesh_region::visitor::visitor(const mesh_region &s, const getfem_mesh &m) : 
    cv_(size_type(-1)), f_(size_type(-1)), finished_(false) {
    s.from_mesh(m);
    init(s);
  }
  

  mesh_region::visitor::visitor(const mesh_region &s) :
    cv_(size_type(-1)), f_(size_type(-1)), finished_(false) {
    init(s);
  }

  void mesh_region::visitor::init(const mesh_region &s) {
    if (&s == 0) 
      DAL_THROW(dal::failure_error, 
		"Attemps to use an invalid mesh_region "
		"(need to call 'from_mesh')");
    it = s.rp().m.begin(); ite = s.rp().m.end(); 
    next();
  }

  std::ostream & operator <<(std::ostream &os, const mesh_region &w) {
    for (mr_visitor cv(w); !cv.finished(); cv.next()) {
      os << cv.cv();
      if (cv.is_face()) os << "/" << cv.f();
      os << " ";
    }
    return os;
  }
}
