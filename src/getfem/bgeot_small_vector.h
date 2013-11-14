/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2000-2012 Julien Pommier
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/

/**@file bgeot_small_vector.h
   @author  Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date January 2004.
   @brief Small (dim < 8) vectors.
*/
#ifndef BGEOT_SMALL_VECTOR_H
#define BGEOT_SMALL_VECTOR_H

#include "dal_singleton.h"
#include "bgeot_config.h"
#ifdef DEBUG_SMALL_VECTOR
# include <cassert>
# define SVEC_ASSERT(x) assert(x)
#else
# define SVEC_ASSERT(x)
#endif 

namespace bgeot {
  class block_allocator {
  public:
    typedef gmm::uint16_type uint16_type;
    typedef gmm::uint32_type node_id;
    typedef gmm::uint32_type size_type;
    /* number of objects stored in a same block, power of 2 */
    enum { p2_BLOCKSZ = 8, BLOCKSZ = 1<<p2_BLOCKSZ }; 
    enum { OBJ_SIZE_LIMIT = 129 }; /* object size limit */
    enum { MAXREF = 256 }; /* reference count limit before copying is used */
  protected:
    /* definition of a block (container of BLOCKSZ chunks) */
    struct block {
      /* effective data + reference count (stored in the BLOCKSZ first bytes) */
      unsigned char * data;
      /* keep track of unused chunks */
      uint16_type first_unused_chunk, count_unused_chunk;
      /* "pointers" for the list of free (or partially filled) blocks */
      size_type prev_unfilled, next_unfilled; 
      size_type objsz; /* size (in bytes) of the chunks stored in this block */
      block() : data(0) {}
      block(size_type objsz_) : data(0), 
				prev_unfilled(size_type(-1)), 
				next_unfilled(size_type(-1)), 
				objsz(objsz_) {}
      ~block() {} /* no cleanup of data, no copy constructor : it's on purpose
		     since the block will be moved a lot when the vector container
		     will be resized (cleanup done by ~block_allocator) */
      void init() {
	clear(); 
	data = static_cast<unsigned char*>(::operator new(BLOCKSZ*objsz + BLOCKSZ)); 
	/* first BLOCKSZ bytes are used for reference counting */
	memset(data, 0, BLOCKSZ);
	//cout << "init block&" << this << " allocated data: " << (void*)data << "\n";
      }
      void clear() { 
	//cout << "clear block&" << this << " frees data: " << (void*)data << "\n";
	if (data) { ::operator delete(data); }; 
	data = 0; first_unused_chunk = 0; count_unused_chunk = BLOCKSZ;
      }
      unsigned char& refcnt(size_type pos) { return data[pos]; }
      bool empty() const { return data == 0; }
      /* could be smarter .. */
    };
    /* container of all blocks .. a vector ensures fast access to 
       any element (better than deque) */
    std::vector<block> blocks; 
    /* pointers to free (or partially free) blocks for each object size */
    size_type first_unfilled[OBJ_SIZE_LIMIT];
  public:
    block_allocator();
    ~block_allocator();
    /* gets the data pointer for an object given its "id" */
    void * obj_data(node_id id) {
      return blocks[id/BLOCKSZ].data + BLOCKSZ + (id%BLOCKSZ)*blocks[id/BLOCKSZ].objsz;
    }
    dim_type obj_sz(node_id id) {
      return dim_type(blocks[id/BLOCKSZ].objsz);
    }
    /* reference counting */
    unsigned char& refcnt(node_id id) {
      return blocks[id/BLOCKSZ].refcnt(id%BLOCKSZ);
    }
    node_id inc_ref(node_id id) {
      if (id && ++refcnt(id) == 0) {
	--refcnt(id);
	id = duplicate(id);
      }
      return id;
    }
    void dec_ref(node_id id) {
      SVEC_ASSERT(id==0 || refcnt(id));
      if (id && --refcnt(id) == 0) {
	++refcnt(id);
	deallocate(id);
      }
    }
    void duplicate_if_aliased(node_id& id) {
      if (refcnt(id) != 1) {
	--refcnt(id);
	id = duplicate(id); SVEC_ASSERT(id == 0 || refcnt(id)==1);
      }
    }
    /* allocation of a chunk */
    node_id allocate(block_allocator::size_type n);
    /* deallocation of a chunk */
    void deallocate(node_id nid);
    void memstats();
  protected:
    /* won't work for non-POD types ... */
    node_id duplicate(node_id id) {
      node_id id2 = allocate(obj_sz(id));
      memcpy(obj_data(id2),obj_data(id),obj_sz(id));
      return id2;
    }
    void insert_block_into_unfilled(block_allocator::size_type bid);
    void remove_block_from_unfilled(block_allocator::size_type bid);
  };
  
  /* common class for all mini_vec, provides access to the common static allocator */
  struct static_block_allocator {
    /* must be a pointer ... sgi CC is not able to order correctly the 
       destructors of static variables */
    static block_allocator *palloc;
    static_block_allocator() { if (!palloc) palloc=&dal::singleton<block_allocator,1000>::instance(); } //new block_allocator(); }
  };
  
  /** container for small vectors of POD (Plain Old Data) types. Should be as fast as 
      std::vector<T> while beeing smaller and uses copy-on-write. The gain is especially
      valuable on 64 bits architectures.
  */
  template<typename T> class small_vector : public static_block_allocator {
    typedef block_allocator::node_id node_id;
    node_id id;
  public:
    typedef small_vector<T> this_type;
    typedef this_type vector_type;
    typedef T value_type;
    typedef T * pointer;
    typedef const T * const_pointer;
    typedef T& reference;
    typedef const T & const_reference;
    typedef T *iterator;
    typedef const T * const_iterator;

    reference operator[](size_type l)
    { GMM_ASSERT2(l < size(), "out of range"); return base()[l]; }
    value_type operator[](size_type l) const
    { GMM_ASSERT2(l < size(), "out of range"); return const_base()[l]; }
    value_type at(size_type l) const { return const_base()[l]; }
    iterator begin() { return base(); }
    const_iterator begin() const { return const_base(); }
    const_iterator const_begin() const { return const_base(); }
    iterator end() { return base()+size(); }
    const_iterator end() const { return const_base()+size(); }
    const_iterator const_end() const { return const_base()+size(); }
    void resize(size_type n) { 
      if (n == size()) return;
      if (n) {
	small_vector<T> other(n); SVEC_ASSERT(other.refcnt() == 1);
	memcpy(other.base(), const_base(), std::min(size(),other.size())*sizeof(value_type));
	SVEC_ASSERT(id==0 || refcnt()); 
	swap(other);
	SVEC_ASSERT(refcnt()); SVEC_ASSERT(other.id == 0 || other.refcnt());
      } else { allocator().dec_ref(id); id=0; }
    }
    const small_vector<T>& operator=(const small_vector<T>& other) { 
      /* order very important when &other == this */
      node_id id2 = allocator().inc_ref(other.id); 
      allocator().dec_ref(id); id = id2;
      SVEC_ASSERT(id == 0 || refcnt()); SVEC_ASSERT(other.id == 0 || other.refcnt());
      return *this;
    }
    void swap(small_vector<T> &v) { std::swap(id,v.id); }
    small_vector() : id(0) {}
    explicit small_vector(size_type n) : id(allocate(n)) {}
    small_vector(const small_vector<T>& v) : static_block_allocator(), id(allocator().inc_ref(v.id)) {}
    explicit small_vector(const std::vector<T>& v) : id(allocate(v.size())) { 
      std::copy(v.begin(),v.end(),begin());
    }
    ~small_vector() { 
      // in the wonderful world of static objects, the order of destruction 
      // can be really important when the memory allocator is destroyed
      // before , for ex. a global variable of type small_vector...
      // that's why there is a check on the state of the allocator..
      if (!allocator_destroyed()) 
	allocator().dec_ref(id); 
    }

    small_vector(T v1, T v2) : id(allocate(2)) 
    { begin()[0] = v1; begin()[1] = v2; }
    small_vector(T v1, T v2, T v3) : id(allocate(3)) 
    { begin()[0] = v1; begin()[1] = v2; begin()[2] = v3; }
    template<class UNOP> small_vector(const small_vector<T>& a, UNOP op) 
      : id(allocate(a.size())) { std::transform(a.begin(), a.end(), begin(), op); }
    template<class BINOP> small_vector(const small_vector<T>& a, const small_vector<T>& b, BINOP op) 
      : id(allocate(a.size())) { std::transform(a.begin(), a.end(), b.begin(), begin(), op); }
    bool empty() const { return id==0; }
    unsigned char refcnt() const { return allocator().refcnt(id); }
    dim_type size() const
    { return dim_type(allocator().obj_sz(id)/sizeof(value_type)); }
    small_vector<T> operator+(const small_vector<T>& other) const 
    { return small_vector<T>(*this,other,std::plus<T>()); }
    small_vector<T> operator-(const small_vector<T>& other) const 
    { return small_vector<T>(*this,other,std::minus<T>()); }
    small_vector<T> operator-() const 
    { return -1.*(*this); }
    small_vector<T> operator*(T v) const 
    { return small_vector<T>(*this, std::bind2nd(std::multiplies<T>(),v)); }
    small_vector<T> operator/(T v) const { return (*this)*(T(1)/v); }
    small_vector<T>& operator+=(const small_vector<T>& other) {
      const_iterator b = other.begin(); iterator it = begin(); 
      for (size_type i=0; i < size(); ++i) *it++ += *b++; 
      return *this;
    }
    small_vector<T>& addmul(T v, const small_vector<T>& other) IS_DEPRECATED;
    //{ std::transform(begin(), end(), other.begin(), begin(), std::plus<T>()); return *this; }
    small_vector<T>& operator-=(const small_vector<T>& other) { 
      const_iterator b = other.begin(); iterator it = begin();
      for (size_type i=0; i < size(); ++i) *it++ -= *b++; 
      return *this;
    }
    small_vector<T> operator*=(T v) { iterator it = begin(), ite=end(); while(it < ite) *it++ *= v; return *this; }
    small_vector<T> operator/=(T v) { return operator*=(T(1)/v); }
    bool operator<(const small_vector<T>& other) const;
    void fill(T v) { for (iterator it=begin(); it != end(); ++it) *it = v; }
    small_vector<T>& operator<<(T x) { push_back(x); return *this; }
    small_vector<T>& clear() { resize(0); return *this; }
    void push_back(T x) { resize(size()+1); begin()[size()-1] = x; }
    size_type memsize() const { return (size()*sizeof(T) / refcnt()) + sizeof(*this); }
  protected:
    /* read-write access (ensures the refcount is 1) */
    pointer base() {
      allocator().duplicate_if_aliased(id);
      return static_cast<pointer>(allocator().obj_data(id)); 
    }
    /* read-only access */
    const_pointer const_base() const { 
      SVEC_ASSERT(id == 0 || refcnt()); return static_cast<pointer>(allocator().obj_data(id)); 
    }
    block_allocator& allocator() const { return *palloc; }
    bool allocator_destroyed() const { return palloc == 0; }
    node_id allocate(size_type n) {
      return node_id(allocator().allocate(gmm::uint32_type(n*sizeof(value_type)))); SVEC_ASSERT(refcnt() == 1);
    }
  };

  template<class T> inline bool small_vector<T>::operator<(const small_vector<T>& other) const {
    return std::lexicographical_compare(begin(), end(), other.begin(), other.end());
  }

  template<class T> inline small_vector<T>& small_vector<T>::addmul(T v, const small_vector<T>& other) {
    const_iterator b = other.begin(); iterator it = begin();
    for (size_type i=0; i < size(); ++i) *it++ += v * *b++; 
    return *this;
  }

  template<class T> std::ostream& operator<<(std::ostream& os, const small_vector<T>& v) {
    os << "["; for (size_type i=0; i < v.size(); ++i) { if (i) os << ", "; os << v[i]; }
    os << "]"; return os;
  }

  template<class T> inline small_vector<T> operator *(T x, const small_vector<T>& m)
  { return m*x; }
}

#endif
