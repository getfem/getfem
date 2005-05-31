// -*- c++ -*- (enables emacs c++ mode)
#ifndef GETFEM_MESH_REGION
#define GETFEM_MESH_REGION

#include <map>
#include <bitset>
#include <iostream>
#include <dal_bit_vector.h>
#include <dal_shared_ptr.h>
#include <bgeot_convex_structure.h>
#include <getfem_config.h>

//#define MAX_FACES_PER_CV 63

namespace getfem {
  class getfem_mesh;

  /**
     mesh_region : the mesh_region structure is used to hold 
     a set of convexes and/or convex faces
  */
  class mesh_region {
  public:
    typedef std::bitset<MAX_FACES_PER_CV+1> face_bitset;
    typedef std::map<size_type,face_bitset> map_t;
  private:
    struct impl {
      mutable dal::bit_vector index_;
      mutable map_t m;
    };
    dal::shared_ptr<impl> p;
    size_type id_;
    getfem_mesh *parent_mesh;
    impl &wp() { return *p.get(); }
    const impl &rp() const { return *p.get(); }
    void clean();
    void touch_parent_mesh();
  public:
    mesh_region() : p(new impl), id_(size_type(-3)), parent_mesh(0) {}
    mesh_region(size_type boundid) : id_(boundid), parent_mesh(0) {}
    mesh_region(getfem_mesh& m, size_type boundid) : 
      p(new impl), id_(boundid), parent_mesh(&m) {}
    mesh_region(const dal::bit_vector &bv) : 
      p(new impl), id_(size_type(-3)), parent_mesh(0) { add(bv); }
    static mesh_region all_convexes() {
      return mesh_region(size_type(-1)); 
    }
    size_type id() const { return id_; }
    const mesh_region& from_mesh(const getfem_mesh &m) const;

    face_bitset operator[](size_t cv) const;
    const dal::bit_vector &index() const;
    void add(const dal::bit_vector &bv);
    void add(size_type cv, size_type f = size_type(-1));
    void sup(size_type cv, size_type f = size_type(-1));
    void clear();
    void swap_convex(size_type cv1, size_type cv2);
    bool is_in(size_type cv, size_type f = size_type(-1)) const;
    size_type size() const;
    size_type nb_convex() const { return rp().m.size(); }  
    bool is_empty() const;
    bool is_only_faces() const;
    bool is_boundary() const { return is_only_faces(); }
    bool is_only_convexes() const;
    face_bitset faces_of_convex(size_type cv) const;
    face_bitset and_mask() const;
    void error_if_not_faces() const;
    void error_if_not_convexes() const;
    void error_if_not_homogeneous() const;
    
    class visitor {
      mesh_region::map_t::const_iterator it,ite;
      face_bitset c;
      size_type cv_, f_;
      bool finished_;
      void init(const mesh_region &s);
    public: 
      visitor(const mesh_region &s);
      visitor(const mesh_region &s, const getfem_mesh &m);
      size_type cv() const { return cv_; }
      size_type is_face() const { return f_ != 0; }
      size_type f() const { return f_-1; }
      bool next() {
	while (c.none()) {
	  if (it == ite) { finished_=true; return false; }
	  c = (*it).second; cv_ = (*it).first; f_ = size_type(-1);
	  ++it; 
	  if (c.none()) continue;
	}
	next_face();
	return true;
      }
      bool operator++() { return next(); }
      bool finished() const { return finished_; }//it == ite && c.none(); }	
      bool next_face() {
	if (c.none()) return false;
	do { ++f_; } while (!c.test(f_));
	c.set(f_,0);
	return true;
      }
    };

    friend std::ostream & operator <<(std::ostream &os, const mesh_region &w);
  };
  
  typedef mesh_region::visitor mr_visitor;
}


#endif
