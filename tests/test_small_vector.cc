/*===========================================================================
 
 Copyright (C) 2007-2012 Yves Renard, Julien Pommier.
 
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
 
===========================================================================*/


#ifndef DEBUG_SMALL_VECTOR
# define DEBUG_SMALL_VECTOR
#endif

#ifdef GETFEM_HAVE_SYS_TIMES
#  include <sys/times.h>
#endif
#include <valarray>
#include <unistd.h>
#include "getfem/bgeot_small_vector.h"
#include "getfem/getfem_mesh.h"

bool quick = false;

#ifdef GETFEM_HAVE_SYS_TIMES
struct chrono {
  struct ::tms t;
  ::clock_t t_elapsed;
  float cpu_, elapsed_, system_;
  float nbclocktk;
public:
  chrono() { nbclocktk = ::sysconf(_SC_CLK_TCK); init(); }
  chrono& init() { elapsed_=0; cpu_=0; system_ =0; return *this; }
  void tic() { t_elapsed = ::times(&t); }
  chrono& toc() { 
    struct tms t2; ::clock_t t2_elapsed = ::times(&t2); 
    elapsed_ += (t2_elapsed - t_elapsed) / nbclocktk;
    cpu_     += (t2.tms_utime - t.tms_utime) / nbclocktk;
    system_  += (t2.tms_stime - t.tms_stime) / nbclocktk;
    memcpy(&t, &t2, sizeof(struct tms));
    return *this;
  }
  float cpu() const { return cpu_; }
  float elapsed() const { return elapsed_; }
  float system() const { return system_; }
};
#else
struct chrono {
  float t,cpu_;
public:
  chrono() { }
  chrono& init() { cpu_=0; return *this; }
  void tic() { t = float(::clock())/float(CLOCKS_PER_SEC); }
  chrono& toc() {
    float t2 = float(::clock())/float(CLOCKS_PER_SEC);
    cpu_ += t2 - t; t = t2; return *this;
  }
  float cpu() const { return cpu_; }
  float elapsed() const { return cpu_; }
  float system() const { return 0.; }
};
#endif
#define REFCNT

namespace test {
  using bgeot::dim_type;
  using bgeot::size_type;  

  class block_allocator {
  public:
    typedef gmm::uint16_type uint16_type;
    /* number of objects stored in a same block, power of 2 */
    enum { p2_BLOCKSZ = 8, BLOCKSZ = 1<<p2_BLOCKSZ }; 
    struct node_id {
      gmm::uint32_type p;
      bool aliased() const { return (p != (p*2)/2); }// & 0x80000000; }
      void set_aliased() { p += 0x80000000; }
      void unset_aliased() { p = (p*2)/2; }//&= ~0x80000000; }
      gmm::uint32_type id() const { return (p*2)/2; } // & 0x7fffffff; }
      gmm::uint32_type bid() const { return id()/BLOCKSZ; }
      gmm::uint32_type chunkid() const { return id()%BLOCKSZ; }
      explicit node_id(gmm::uint32_type p_) : p(p_) {}
      node_id(gmm::uint32_type b, gmm::uint32_type cid, bool al_) : p((b*BLOCKSZ+cid)+(al_?0x80000000:0)) {}
      bool null() const { return p == 0; }
      void nullify() { p = 0; }
    };


    typedef gmm::uint32_type size_type;
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
      return blocks[id.bid()].data + BLOCKSZ + (id.chunkid())*blocks[id.bid()].objsz;
    }
    dim_type obj_sz(node_id id) {
      return dim_type(blocks[id.bid()].objsz);
    }
    /* reference counting */
    unsigned char& refcnt(node_id id) {
      return blocks[id.bid()].refcnt(id.chunkid());
    }
    node_id inc_ref(node_id& id) {
      if (!id.null()) {
	if (refcnt(id)++==1) id.set_aliased(); 
	else if (refcnt(id) == 0) {
	  --refcnt(id);
	  return duplicate(id);
	}
      }
      return id;
    }
    void dec_ref(node_id& id) {
      SVEC_ASSERT(id.null() || refcnt(id));
      if (!id.null()) {
	--refcnt(id);
	if (refcnt(id) == 0) {
	  ++refcnt(id);
	  deallocate(id);
	} else if (refcnt(id) == 1) {
	  SVEC_ASSERT(id.aliased());
	  id.unset_aliased();
	}
      }
    }
    void duplicate_if_aliased(node_id& id) {
      if (id.aliased()) {
	SVEC_ASSERT(refcnt(id)>=1);
	if (refcnt(id)>1) {
	  --refcnt(id);
	  id = duplicate(id); SVEC_ASSERT(id.null() || (refcnt(id)==1 && !id.aliased()));
	} else id.unset_aliased();
      }
    }
    /* allocation of a chunk */
    node_id allocate(size_type n);
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
    void insert_block_into_unfilled(size_type bid);
    void remove_block_from_unfilled(size_type bid);
  };
  
  /* common class for all mini_vec, provides access to the common static allocator */
  struct static_block_allocator {
    static block_allocator alloc;
  };
  
  /**
    small_vector class: container for small vectors of POD types. Should be as fast as 
    std::vector<T> while beeing smaller and uses copy-on-write. The gain is especially
    valuable on 64 architectures.
  */
  template<typename T> class small_vector : public static_block_allocator {
    typedef block_allocator::node_id node_id;
    mutable node_id id;
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

    void out_of_range_error(void) const { GMM_ASSERT1(false, "out of range"); }
    reference operator[](size_type l) { if (l >= size()) out_of_range_error(); return base()[l]; }
    value_type operator[](size_type l) const { if (l >= size()) out_of_range_error(); return const_base()[l]; }
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
	SVEC_ASSERT(id.null() || refcnt()); 
	swap(other);
	SVEC_ASSERT(refcnt()); SVEC_ASSERT(other.id.null() || other.refcnt());
      } else { allocator().dec_ref(id); id.nullify(); }
    }
    const small_vector<T>& operator=(const small_vector<T>& other) { 
      /* order very important when &other == this */
      node_id id2 = allocator().inc_ref(other.id); 
      allocator().dec_ref(id); id = id2;
      SVEC_ASSERT(id.null() || refcnt()); SVEC_ASSERT(other.id.null() || other.refcnt());
      return *this;
    }
    void swap(small_vector<T> &v) { std::swap(id,v.id); }
    small_vector() : id(0) {}
    small_vector(size_type n) : id(allocate(n)) {}
    small_vector(const small_vector<T>& v) : id(allocator().inc_ref(v.id)) {}
    ~small_vector() { allocator().dec_ref(id); }

    small_vector(T v1, T v2) : id(allocate(2)) 
    { begin()[0] = v1; begin()[1] = v2; }
    small_vector(T v1, T v2, T v3) : id(allocate(3)) 
    { begin()[0] = v1; begin()[1] = v2; begin()[2] = v3; }
    template<class UNOP> small_vector(const small_vector<T>& a, UNOP op) 
      : id(allocate(a.size())) { std::transform(a.begin(), a.end(), begin(), op); }
    template<class BINOP> small_vector(const small_vector<T>& a, const small_vector<T>& b, BINOP op) 
      : id(allocate(a.size())) { std::transform(a.begin(), a.end(), b.begin(), begin(), op); }
    bool empty() const { return this->id==node_id(0); }
    unsigned char refcnt() const { return allocator().refcnt(id); }
    dim_type size() const { return allocator().obj_sz(id)/sizeof(value_type); }
    small_vector<T> operator+(const small_vector<T>& other) const 
    { return small_vector<T>(*this,other,std::plus<T>()); }
    small_vector<T> operator-(const small_vector<T>& other) const 
    { return small_vector<T>(*this,other,std::minus<T>()); }
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
      SVEC_ASSERT(id.null() || refcnt()); return static_cast<const_pointer>(allocator().obj_data(id)); 
    }
    block_allocator& allocator() const { return alloc; }
    node_id allocate(size_type n) {
      return (allocator().allocate(n*sizeof(value_type))); SVEC_ASSERT(refcnt() == 1);
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

  block_allocator static_block_allocator::alloc;

  block_allocator::block_allocator() {
    for (size_type i=0; i < OBJ_SIZE_LIMIT; ++i) 
      first_unfilled[i] = i ? size_type(-1) : 0; 
    /* bloc 0 is reserved for objects of size 0 -- it won't grow */
    blocks.push_back(block(0)); blocks.front().init();
  }
  block_allocator::~block_allocator() {
    for (size_type i=0; i < blocks.size(); ++i) 
      if (!blocks[i].empty()) blocks[i].clear();
  }
  block_allocator::node_id block_allocator::allocate(block_allocator::size_type n) {
    if (n == 0) return node_id(0);
    if (n >= OBJ_SIZE_LIMIT) 
      GMM_ASSERT1(false, 
		  "attempt to allocate a supposedly \"small\" object of " 
		  << n << " bytes\n");
    //cout << "dim = " << n << " ";
    if (first_unfilled[n] == size_type(-1)) {
      blocks.push_back(block(n)); blocks.back().init();
      insert_block_into_unfilled(gmm::uint32_type(blocks.size()-1));
      if (first_unfilled[n] >= (1<<(31-p2_BLOCKSZ))) {//(node_id(1)<<(sizeof(node_id)*CHAR_BIT - p2_BLOCKSZ))) {
	GMM_ASSERT1(false,
		    "allocation slots exhausted for objects of size " << 
		    n << " (" << first_unfilled[n] << " allocated!),\n" << 
		    "either increase the limit, or check for a leak in your code.");
	//cout << "created new block " << first_unfilled[n] << "\n";
      }
    }
    block &b = blocks[first_unfilled[n]]; SVEC_ASSERT(b.objsz == n);
    if (b.empty()) b.init(); /* realloc memory if needed */
    size_type vid = b.first_unused_chunk; SVEC_ASSERT(vid < BLOCKSZ); 
    node_id id(first_unfilled[n], vid, false);
    SVEC_ASSERT(b.refcnt(b.first_unused_chunk)==0);
    b.refcnt(vid) = 1; b.count_unused_chunk--;
    if (b.count_unused_chunk) {
      do b.first_unused_chunk++; while (b.refcnt(b.first_unused_chunk));
    } else {
      b.first_unused_chunk = BLOCKSZ;
      remove_block_from_unfilled(first_unfilled[n]);
    }
    //cout << "allocated " << first_unfilled[n] << ", " << vid << " @" << obj_data(id) << "\n";
    SVEC_ASSERT(obj_data(id));
    //SVEC_ASSERT(refcnt(id) == 0);
    return id;
  }
  void block_allocator::deallocate(block_allocator::node_id nid) {
    if (nid.null()) return;
    size_type bid = nid.bid();
    size_type vid = nid.chunkid();
    block &b = blocks[bid];      
    //cout << "deallocate " << bid << "/dim=" << b.dim << ", " << vid << ", unused=" << b.unused << "\n";
    SVEC_ASSERT(b.refcnt(vid) == 1);
    b.refcnt(vid) = 0;
    if (b.count_unused_chunk++ == 0) {
      insert_block_into_unfilled(bid); 
      b.first_unused_chunk = gmm::uint16_type(vid);
    } else {
      b.first_unused_chunk = gmm::uint16_type(std::min<size_type>(b.first_unused_chunk,vid));
      if (b.count_unused_chunk == BLOCKSZ) b.clear();
    }
  }
  void block_allocator::memstats()  {
    cout << "block_allocator memory statistics:\ntotal number of blocks: " 
	 << blocks.size() << ", each blocks stores " << BLOCKSZ 
	 << " chuncks; size of a block header is " << sizeof(block) << " bytes\n";
    for (size_type d = 0; d < OBJ_SIZE_LIMIT; ++d) {
      size_type total_cnt=0, used_cnt=0, mem_total = 0, bcnt = 0;
      for (size_type i=0; i < blocks.size(); ++i) {
	if (blocks[i].objsz != d) continue; else bcnt++;
	if (!blocks[i].empty()) {
	  total_cnt += BLOCKSZ;
	  used_cnt += BLOCKSZ - blocks[i].count_unused_chunk;
	  mem_total += (BLOCKSZ+1)*blocks[i].objsz;
	}
	mem_total = gmm::uint32_type(mem_total + sizeof(block));
      }
      if (mem_total)
	cout << " sz " << d << ", memory used = " << mem_total << " bytes for " 
	     << total_cnt << " nodes, unused space = " 
	     << (total_cnt == 0 ? 100. : 100. - 100.* used_cnt / total_cnt) 
	     << "%, bcnt=" << bcnt << "\n";
    }
  }
  void block_allocator::insert_block_into_unfilled(block_allocator::size_type bid) {
    SVEC_ASSERT(bid < blocks.size());
    dim_type dim = dim_type(blocks[bid].objsz);
    SVEC_ASSERT(bid != first_unfilled[dim]);
    SVEC_ASSERT(blocks[bid].prev_unfilled+1 == 0);
    SVEC_ASSERT(blocks[bid].next_unfilled+1 == 0);
    blocks[bid].prev_unfilled = size_type(-1);
    blocks[bid].next_unfilled = first_unfilled[dim];
    if (first_unfilled[dim] != size_type(-1)) {
      SVEC_ASSERT(blocks[first_unfilled[dim]].prev_unfilled+1 == 0);
      blocks[first_unfilled[dim]].prev_unfilled = bid;
    }
    first_unfilled[dim] = bid;
    //cout << "** bloc " << bid << " has been INSERTED in unfilled list, which is now"; show_unfilled(bid);
  }
  void block_allocator::remove_block_from_unfilled(block_allocator::size_type bid) {
    SVEC_ASSERT(bid < blocks.size());
    dim_type dim = dim_type(blocks[bid].objsz);
    //cout << "** bloc " << bid << " is going to be REMOVE unfilled list, which is now"; show_unfilled(bid);
    size_type p = blocks[bid].prev_unfilled; blocks[bid].prev_unfilled = size_type(-1);
    size_type n = blocks[bid].next_unfilled; blocks[bid].next_unfilled = size_type(-1);
    if (p != size_type(-1)) { blocks[p].next_unfilled = n; }
    if (n != size_type(-1)) { blocks[n].prev_unfilled = p; }
    if (first_unfilled[dim] == bid) { SVEC_ASSERT(p+1==0); first_unfilled[dim] = n; }
    //cout << "** bloc " << bid << " has been REMOVED in unfilled list, which is now"; show_unfilled(bid);
  }

} /* namespace test */

namespace getfem {
  std::allocator<double> al;
  struct micro_vec {
    typedef double value_type;
    typedef double *  pointer;
    typedef double &  reference;
    typedef double *iterator;
    typedef const double * const_iterator;
    pointer p;
    pointer base() const { return (double*)((unsigned long)p&(~7UL)); }
    reference operator[](size_type i) { return base()[i]; }
    const double& operator[](size_type i) const { return base()[i]; }
    micro_vec() : p(0) {}
    static pointer alloc(size_type n) { 
      if (n==0) return 0;
      double *p = al.allocate(n); //new double[n]; 
      return (double*) ((char*)p + n); 
    }
    micro_vec(size_type n) : p(alloc(n)) {}
    micro_vec(const micro_vec& v) : p(alloc(v.size())) {
      memcpy(base(), v.base(), size()*sizeof(double));
    }
    iterator begin() { return base(); }
    const_iterator begin() const { return base(); }
    iterator end() { return base()+size(); }
    const_iterator end() const { return base()+size(); }

    void resize(size_type) { /* hum */ }
    void deallocate() { 
      if (p) { al.deallocate(base(),size()); } //delete[] base();
    }
    micro_vec& operator=(const micro_vec& other) {
      if (&other != this) {
	if (other.size() != size()) {
	  deallocate();      
	  p = alloc(other.size());
	}
	memcpy(base(), other.base(), size()*sizeof(double));      
      }
      return *this;
    }
    ~micro_vec() { deallocate(); }
    size_type size() const { return (unsigned long)p & 7UL; }
    micro_vec(const micro_vec& va, const micro_vec& vb) : p(alloc(va.size())) { 
      /*double * restrict__ pb = base(), *pe = base()+size(), *pva = va.base(), *pvb = vb.base();
	for (; pb < pe;) *pb++=*pva+++*pvb++;*/
      for (size_type i=0; i < va.size(); ++i) (*this)[i] = va[i]+vb[i];
    }
    micro_vec operator+(const micro_vec& v) const {
      return micro_vec(*this, v);
    }
    
    micro_vec& operator+=(const micro_vec& v) {
      for (unsigned i=0; i < size(); ++i) {
	(*this)[i] += v[i];
      }
      return *this;
    }
    micro_vec& operator-=(const micro_vec& v) {
      for (unsigned i=0; i < v.size(); ++i) {
	(*this)[i] -= v[i];
      }
      return *this;
    }
  };





  template<class V> void init(std::vector<V>& vv) {
    for (size_type i=0; i < vv.size(); ++i) {
      size_type sz = (rand()%5)  +1;
      vv[i] = V(sz);
    }
    
    for (size_type i=0; i < vv.size(); ++i) {
      vv[i].resize(2);
      vv[i] = V(2);
      vv[i][0] = double(i); vv[i][1] = double(i)/double(23); //vv[i][2] = i;
    }
  }

  template <class V> void rrun(std::vector<V>& vv) {
    chrono c;
    c.init().tic();
    cout << GMM_PRETTY_FUNCTION << "\n";
    size_type N = quick ? 10 : 25;
    for (size_type k=0; k < N; ++k) {
      init(vv);
    }    
    cout << " creation : " << c.toc().cpu() << " sec\n, size=" 
	 << vv.capacity()*((char*)&vv[1] - (char*)&vv[0]) << " + "
         << vv.size() * vv[0].size()*((char*)&vv[1][0] - (char*)&vv[0][0]) 
	 << " bytes\n";

    c.init().tic();
    
    for (size_type i=0; i < N*5; ++i) {
      for (size_type j=0; j < vv.size(); ++j) {
	for (size_type k=0; k < vv[j].size(); ++k) {
	  vv[j][k] = double(j+k);
	  assert(vv[j][k] == j+k);
	}
	vv[j] = vv[j];
	for (size_type k=0; k < vv[j].size(); ++k) assert(vv[j][k] == j+k);
      }
    }
    cout << "remplissage : " << c.toc().cpu() << " sec\n";

    c.init().tic();
    for (size_type i=0; i < N*5; ++i) {
      for (size_type j=0; j < vv.size(); ++j) {
	for (size_type k=0; k < vv[j].size(); ++k) {
	  if (gmm::abs(vv[j][k] - double(j+k))>1e-14)
	    { cout << "vv[j][k]=" << vv[j][k] << "!=" <<j+k<< "\n"; abort(); }
	}
      }
    }
    cout << "lecture : " << c.toc().cpu() << " sec\n";

    c.init().tic();
    for (size_type i=0; i < N; ++i) {
      std::vector<V> w, z;
      w = vv;
      std::random_shuffle(vv.begin(), vv.end());
      z = w;
    }
    cout << " copies : " << c.toc().cpu() << " sec\n";

    c.init().tic();
    for (size_type i=0; i < N; ++i) {
      init(vv);
      for (size_type j=1; j < vv.size()-1; ++j) {
	vv[j] = vv[j-1] + vv[j] + vv[j] + vv[j+1];
      }
    }
    cout << " operator + : " << c.toc().cpu() << " sec\n";
    c.init().tic();
    for (size_type i=0; i < N; ++i) {
      init(vv);
      for (size_type j=vv.size()-1; j > 0; --j) {
	vv[j] -= vv[j-1];
      }
    }
    cout << " operator -=: " << c.toc().cpu() << " sec\n";
//     for (size_type i=0; i < (vv.size()<4 ? vv.size() : 4); ++i) {
//       printf("V[%d]@%p, V[%d][0]@%p\n", int(i), &vv[i], int(i), &vv[i][0]);
//     }
    gmm::lexicographical_less<V, gmm::approx_less<typename V::value_type> > comp;
    c.init().tic();
    init(vv); std::random_shuffle(vv.begin(), vv.end());
    for (size_type i=0; i < N*4; ++i) {
      size_type cnt = 0;
      for (size_type j=vv.size()-1; j > 0; --j) {
	cnt += comp(vv[j],vv[j-1]) ? 1:-1;
      }
    }
    cout << " comparaison: " << c.toc().cpu() << " sec\n";
    
    getfem::mesh m;

    cout << "mesh<base_node> : empty size = " << m.memsize() << "\n";
    init(vv);
    c.init().tic();
    std::random_shuffle(vv.begin(), vv.end());
    for (size_type i=0; i < N; ++i) {
      m.clear();
      for (size_type j=0; j < vv.size(); ++j) {
	m.add_point(vv[j]);
      }
    }
    cout << "mesh<base_node> : size for " << m.nb_points() 
	 << " points: " << m.memsize() << " bytes, cpu="
	 << c.toc().cpu() << " sec\n";
    
  }

  template<class MICRO_VEC> void hop() {
    /*std::vector<MICRO_VEC::node_id> w(300);
    for (size_type j=0; j < w.size(); ++j) {
      w[j] = MICRO_VEC::alloc.allocate(3);
    }
    std::random_shuffle(w.begin(), w.end());
    for (size_type j=0; j < w.size(); ++j) {
      MICRO_VEC::alloc.deallocate(w[j]);
    }
    MICRO_VEC::alloc.memstats();*/
  }


  template<class V>void refccheck() {
    { 
      std::vector<V> v;
      v.clear();
      v.resize(259, V(3));
    /*for (size_type i=0; i < v.size(); ++i) { 
      cout << "v[" << i << "]=" << *(unsigned*)(&v[i]) << ", refcnt=" << (int)v[i].refcnt() << "\n"; 
      }*/
    }
    {
      std::vector<V> v(342,V(2)); v[0][0] = 3; v[0][1] = 6;
      for (size_type i=1; i < v.size(); ++i) {
	v[i] = v[0];
	/*	cout << "v[" << 0 << "]=" << *(unsigned*)(&v[0]) << ", refcnt=" << (int)v[0].refcnt() << ", "; 
	cout << "v[" << i << "]=" << *(unsigned*)(&v[i]) << ", refcnt=" << (int)v[i].refcnt() 
	     << "val={" << v[i].const_begin()[0] << "," << v[i].const_begin()[1] << "}\n";
	*/
      }
      for (size_type i=1; i < v.size(); ++i) {
	for (size_type j=0; j < std::max(v[0].size(),v[i].size()); ++j) {
	  typename V::value_type a = v[0][j];
	  typename V::value_type b = v[i][j];
	  assert(a == b);
	}
	scalar_type d = gmm::vect_dist2_sqr(v[i],v[0]);
	assert(d ==0.);
      }
      for (size_type i=1; i < v.size(); ++i) v[i] = v[i-1];
    }
    {
      std::vector<V> w(342);
      std::vector<V> v(w);
      for (size_type i=1; i < w.size(); ++i) v[i] = w[i];
      for (size_type i=1; i < w.size(); ++i) assert(gmm::vect_dist2_sqr(v[i],w[i])==0.);
      for (size_type i=1; i < w.size(); ++i) v[i] = w[0]+w[1];
      for (size_type i=1; i < w.size(); ++i) v[i] = w[0];
      for (size_type i=0; i < 8; ++i) {
	v[0].resize(i);
	v[i] = v[0];
      }
      v[0].clear();
      for (size_type i=0; i < 8; ++i) {
	v[0].push_back(i);
      }
    }
  }

#ifdef REFCNT
#  define IFREFCNT(x) (x)
#else
#  define IFREFCNT(x) 1
#endif


  template<class MICRO_VEC> void hop2() {
    MICRO_VEC v(1); v[0] = 100;
    assert(v.refcnt() == 1);
    MICRO_VEC w(v);
    assert(v.refcnt() == IFREFCNT(2));
    assert(w.refcnt() == IFREFCNT(2));
    w = v;
    assert(v.refcnt() == IFREFCNT(2));
    assert(w.refcnt() == IFREFCNT(2));
    MICRO_VEC z(3);
    w = z;
    assert(v.refcnt() == 1);
    assert(w.refcnt() == IFREFCNT(2));
    assert(z.refcnt() == IFREFCNT(2));
    z[0] = 1;
    assert(v.refcnt() == 1);
    assert(w.refcnt() == 1);
    assert(z.refcnt() == 1);
    MICRO_VEC *pz;
    for (size_type i=0; i < 67; ++i) {
      pz = new MICRO_VEC(z);
      //cout << "i=" << i << ", pz=" << (void*)pz->const_begin() << ", pz->cnt=" << (int)pz->refcnt() << ", z->cnt=" << (int)z.refcnt() << "\n";
      assert(pz->refcnt() == IFREFCNT((i < test::block_allocator::MAXREF-2) ? i+2 : 1));
      assert(z.refcnt() == IFREFCNT(std::min<unsigned>((i+2),test::block_allocator::MAXREF-1)));
    }
    std::vector<MICRO_VEC*> pp;
    for (size_type i=0; i < 67; ++i) {
      pp.push_back(new MICRO_VEC(z));
      for (size_type j=0; j <= i; ++j) {
	pz = pp[j]; assert(pz->size() == 3);
	//cout << (int)pz->refcnt() << "\n";
	(*pz)[0] = 3*(i-j)+1;(*pz)[1] = 3*(i-j)+2;(*pz)[2] = 3*(i-j)+3;
	//cout << (int)pz->refcnt() << "\n";
	assert(pz->refcnt() == 1);
      }
      assert(pz->refcnt() == 1);
      //cout << "i=" << i << ", pp[i]=" << (void*)pz->const_begin() << ", pp[i]->cnt=" << (int)pz->refcnt() << "\n";
    }
    w = MICRO_VEC(2); w[0] = 2; w[1] = 4; //typename MICRO_VEC::set(), 2, 4);
    w *= 2; w /= 2; w += w; w -= w; w = w+w; w = w-w; w = (w*3) - w/2;
    refccheck<MICRO_VEC>();
  }

  /*void runhop() {
    cout << "coucou\n";
    for (size_type i=0; i < 5; ++i) {
      hop2<bgeot::small_vector<double> >();
      //hop2<bgeot::small_vector<unsigned> >();
      hop2<bgeot::small_vector<float> >();
      hop2<bgeot::small_vector<short> >();
      hop2<bgeot::small_vector<char> >();
      hop<bgeot::small_vector<double> >();       
    }
    bgeot::static_block_allocator::palloc->memstats();
    for (size_type i=0; i < 5; ++i) {
      hop2<test::small_vector<double> >();
      //hop2<test::small_vector<unsigned> >();
      hop2<test::small_vector<float> >();
      hop2<test::small_vector<short> >();
      hop2<test::small_vector<char> >();
      hop<test::small_vector<double> >(); 
    }
    test::static_block_allocator::alloc.memstats();
  }
  */

  void run() {
    //runhop();
    size_type N=quick ? 2311 : 20000;
    //std::vector<base_vector> bv(N);
    std::vector<base_node> vv(N);
    //std::vector<micro_vec> sv(N);
    std::vector<std::valarray<double> > av(N);
    //std::vector<test::small_vector<double> > mv(N);
    std::vector<bgeot::small_vector<double> > Sv(N);
    
    //cerr << av[0].size() << "\n";
    //base_node n(3); cerr << n[3] << "\n";
    //rrun(bv);
    rrun(vv);
    //rrun(sv);
    //rrun(mv);
    rrun(Sv);
    //rrun(av);
    bgeot::static_block_allocator::palloc->memstats();
    cout << "sizeof(size_type)=" << sizeof(size_type) 
	 << ", sizeof(base_node)=" << sizeof(base_node) 
	 << ", sizeof(base_small_vector)=" << sizeof(base_small_vector) << "\n";
  }
}

int main(int argc, char **argv) {

  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  if (argc == 2 && strcmp(argv[1],"-quick")==0) quick = true;
  getfem::run();
}
