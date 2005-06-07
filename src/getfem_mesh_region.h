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
    dal::shared_ptr<impl> p;  /* the real region data */
    size_type id_;            /* used temporarily when the 
				 mesh_region(size_type) constructor is used */
    getfem_mesh *parent_mesh; /* used for mesh_region "extracted" from
				 a getfem_mesh (to provide feedback) */
    impl &wp() { return *p.get(); }
    const impl &rp() const { return *p.get(); }
    void clean();
    /** tells the owner mesh that the region is valid */
    void touch_parent_mesh();
  public:
    mesh_region() : p(new impl), id_(size_type(-3)), parent_mesh(0) {}
    /** a mesh_region can be built from a integer parameter 
	(a region number in a mesh),
	but it won't be usable until 'from_mesh(m)' has been 
	called 
	Note that these regions are read-only, this constructor is
	mostly used for backward-compatibiliy.
    */
    mesh_region(size_type boundid) : id_(boundid), parent_mesh(0) {}
    /** internal constructor. You should used m.region(id) instead. */
    mesh_region(getfem_mesh& m, size_type id__) : 
      p(new impl), id_(id__), parent_mesh(&m) {}
    /** build a mesh_region from a convex list stored in a bit_vector. */
    mesh_region(const dal::bit_vector &bv) : 
      p(new impl), id_(size_type(-3)), parent_mesh(0) { add(bv); }
    static mesh_region all_convexes() {
      return mesh_region(size_type(-1)); 
    }
    static mesh_region intersection(const mesh_region& a, const mesh_region& b);
    size_type id() const { return id_; }

    /** for regions which have been built with just a number 'id',
	from_mesh(m) sets the current region to 'm.region(id)'.  
	(works only once) 
    */
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
    /** return true if the region do contain only convex faces */
    bool is_only_faces() const;
    bool is_boundary() const { return is_only_faces(); }
    /** return true if the region do not contain any convex face */
    bool is_only_convexes() const;
    face_bitset faces_of_convex(size_type cv) const;
    face_bitset and_mask() const;
    void error_if_not_faces() const;
    void error_if_not_convexes() const;
    void error_if_not_homogeneous() const;
    
    /** "iterator" class for regions. Usage similar to bv_visitor:
	for (mr_visitor i(region); !i.finished(); ++i) {
	  ...
        }
    */
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
