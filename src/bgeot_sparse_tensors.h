// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2000-2006 Julien Pommier
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

/**@file bgeot_sparse_tensors.h
   @author  Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date January 2003.
   @brief Sparse tensors, used during the assembly.

   "sparse" tensors: these are not handled like sparse matrices
   
   As an example, let say that we have a tensor t(i,j,k,l) of
   dimensions 4x2x3x3, with t(i,j,k,l!=k) == 0. 
   
   Then the tensor shape will be represented by a set of 3 objects of type
     'tensor_mask':
   mask1: {i}, "1111"
   mask2: {j}, "11"
   mask3: {k,l}, "100"
                 "010"
                 "001"
   They contain a binary tensor indicating the non-null elements.
   
   The set of these three masks define the shape of the tensor
   (class tensor_shape)

   If we add information about the location of the non-null elements
   (by mean of strides), then we have an object of type 'tensor_ref'
   
   Iteration on the data of one or more tensor should be done via the
   'multi_tensor_iterator', which can iterate over common non-null
   elements of a set of tensors.


   maximum (virtual) number of elements in a tensor : 2^31
   maximum number of dimensions : 254

   "ought to be enough for anybody"
*/
#ifndef BGEOT_SPARSE_TENSORS
#define BGEOT_SPARSE_TENSORS

#include <bgeot_config.h>
#include <dal_except.h>
#include <dal_bit_vector.h>
// #include <gmm_kernel.h> // for i/o on vectors it is convenient
#include <iostream>
#include <bitset>

namespace bgeot {
  typedef dal::uint32_type index_type;
  typedef dal::int32_type stride_type; /* signé! */

  //  typedef std::vector<index_type> tensor_ranges;
  class tensor_ranges : public std::vector<index_type> {
  public:
    tensor_ranges() : std::vector<index_type>() {}
    tensor_ranges(size_type n) : std::vector<index_type>(n) {}
    tensor_ranges(size_type n, index_type V) : std::vector<index_type>(n,V) {}
  };
  typedef std::vector<stride_type> tensor_strides;
  typedef std::vector<dim_type> index_set;

  typedef scalar_type * TDIter;

  std::ostream& operator<<(std::ostream& o, const tensor_ranges& r); 
  
  /* stupid && inefficient loop structure */
  struct tensor_ranges_loop {
    tensor_ranges sz;
    tensor_ranges cnt;
    bool finished_;
  public:
    tensor_ranges_loop(const tensor_ranges& t) : sz(t), cnt(t.size()), finished_(t.size() == 0) { 
      std::fill(cnt.begin(), cnt.end(), 0); 
    }
    index_type index(dim_type i) { return cnt[i]; }
    bool finished() const { return finished_; }
    bool next() { 
      index_type i = 0;
      while (++cnt[i] >= sz[i]) {
	cnt[i] = 0; i++; if (i >= sz.size()) { finished_ = true; break; }
      }
      return finished_;
    }
  };

  /* handle a binary mask over a given number of indices */
  class tensor_mask {
    tensor_ranges r;
    index_set idxs;
    std::vector<bool> m;
    tensor_strides s; /* strides in m */
    mutable index_type card_; /* finally i should have kept m as a dal::bit_vector ... */
    mutable bool card_uptodate;
  public:
    tensor_mask() { set_card(0); }
    explicit tensor_mask(const tensor_ranges& r_, const index_set& idxs_) {
      assign(r_,idxs_);
    }    
    /* constructeur par fusion */
    explicit tensor_mask(const tensor_mask& tm1, const tensor_mask& tm2, bool and_op);
    explicit tensor_mask(const std::vector<const tensor_mask*>& tm);    
    explicit tensor_mask(const std::vector<const tensor_mask*> tm1, 
			 const std::vector<const tensor_mask*> tm2, bool and_op);
    void swap(tensor_mask &tm) {
      r.swap(tm.r); idxs.swap(tm.idxs);
      m.swap(tm.m); s.swap(tm.s); 
      std::swap(card_, tm.card_);
      std::swap(card_uptodate, tm.card_uptodate);
    }
    void assign(const tensor_ranges& r_, const index_set& idxs_) {
      r = r_; idxs = idxs_; eval_strides(); m.assign(size(),false);
      set_card(0);
    }
    void assign(const tensor_mask& tm) { 
      r = tm.r; 
      idxs = tm.idxs; 
      m = tm.m; 
      s = tm.s;
      card_ = tm.card_; card_uptodate = tm.card_uptodate;
    }
    void assign(const std::vector<const tensor_mask* >& tm);
    void assign(const tensor_mask& tm1, const tensor_mask& tm2, bool and_op);
    
    void clear() { r.resize(0); idxs.resize(0); m.clear(); s.resize(0); set_card(0); }
    const tensor_ranges& ranges() const { return r; }
    const index_set& indexes() const { return idxs; }
    const tensor_strides& strides() const { return s; }
    index_set& indexes() { return idxs; }
    void eval_strides() {
      s.resize(r.size()+1); s[0]=1;
      for (index_type i=0; i < r.size(); ++i) {
	s[i+1]=s[i]*r[i];
      }
    }
    index_type ndim() const { return r.size(); }
    index_type size() const { return s[r.size()]; }
    void set_card(index_type c) const { card_ = c; card_uptodate = true; }
    void unset_card() const { card_uptodate = false; }
    index_type card(bool just_look=false) const {       
      if (!card_uptodate || just_look) {
	index_type c = std::count_if(m.begin(), m.end(), 
				     std::bind2nd(std::equal_to<bool>(),true));
	if (just_look) return c;
	card_ = c;
      }
      return card_;
    }
    index_type pos(tensor_ranges& global_r) const {
      index_type p = 0;
      for (index_type i=0; i < r.size(); ++i) 
	p+= s[i]*global_r[idxs[i]];
      return p;
    }
    index_type lpos(tensor_ranges& local_r) const {
      index_type p = 0;
      for (index_type i=0; i < r.size(); ++i) 
	p+= s[i]*local_r[i];
      return p;
    }
    bool operator()(tensor_ranges& global_r) const {
      return m[pos(global_r)];
    }
    bool operator()(stride_type p) const { return m[p]; }
    void set_mask_val(stride_type p, bool v) { m[p]=v; card_uptodate = false; }
    struct Slice {
      dim_type dim;
      index_type i0;
      Slice(dim_type d, index_type i0_) : dim(d), i0(i0_) {}
    };

    /* cree un masque de tranche */
    void set_slice(index_type dim, index_type range, index_type islice) {
      r.resize(1); r[0] = range;
      idxs.resize(1); idxs[0] = dim;
      m.clear(); m.assign(range,false); m[islice] = 1; set_card(1);
      eval_strides();
    }
    explicit tensor_mask(index_type range, Slice slice) {
      set_slice(slice.dim, range, slice.i0); 
    }

    struct Diagonal {
      dim_type i0, i1;
      Diagonal(dim_type i0_, dim_type i1_) : i0(i0_), i1(i1_) {}
    };

    /* cree un masque diagonal */
    void set_diagonal(index_type n, index_type i0, index_type i1) {
      assert(n);
      r.resize(2); r[0] = r[1] = n;
      idxs.resize(2); idxs[0] = i0; idxs[1] = i1;
      m.assign(n*n, false); 
      for (index_type i=0; i < n; ++i) m[n*i+i]=true;
      set_card(n);
      eval_strides();
    }
    explicit tensor_mask(index_type n, Diagonal diag) {
      set_diagonal(n, diag.i0, diag.i1);
    }
    void set_triangular(index_type n, index_type i0, index_type i1) {
      assert(n);
      r.resize(2); r[0] = r[1] = n;
      idxs.resize(2); idxs[0] = i0; idxs[1] = i1;
      m.assign(n*n,false); unset_card();
      for (index_type i=0; i < n; ++i)
	for (index_type j=i; j < n; ++j) m[i*n+j]=true;
      eval_strides();
    }
    void set_full(index_type dim, index_type range) {
      assert(range);
      r.resize(1); r[0] = range;
      idxs.resize(1); idxs[0] = dim;
      m.assign(range, true); set_card(range);
      eval_strides();
    }
    void set_empty(index_type dim, index_type range) {
      assert(range);
      r.resize(1); r[0] = range;
      idxs.resize(1); idxs[0] = dim;
      m.assign(range,false); set_card(0);
      eval_strides();
    }
    explicit tensor_mask(index_type dim, index_type range) {
      set_full(dim, range);
    }
    void set_zero() {
      m.assign(size(),false); set_card(0);
    }
    void shift_dim_num_ge(dim_type dim, int shift) {
      for (dim_type i=0; i < idxs.size(); ++i) {
	if (idxs[i] >= dim) idxs[i]+=shift;
      }
      check_assertions();
    }
    void gen_mask_pos(tensor_strides& p) const {
      check_assertions();
      p.resize(card());
      index_type i = 0;
      for (tensor_ranges_loop l(r); !l.finished(); l.next()) {
	if (m[lpos(l.cnt)]) p[i++] = lpos(l.cnt);
      }
      assert(i==card());
    }
    void unpack_strides(const tensor_strides& packed, tensor_strides& unpacked) const;

    /* c'est mieux que celle ci renvoie un int ..
       ou alors un unsigned mais dim_type c'est dangereux */
    int max_dim() const {
      index_set::const_iterator it = std::max_element(idxs.begin(),idxs.end());
      return (it == idxs.end() ? -1 : *it);
    }
    void check_assertions() const;
    void print(std::ostream &o) const;
    void print_() const { print(cerr); }
  };



  typedef std::vector<tensor_mask> tensor_mask_container;

  struct tensor_index_to_mask {
    short_type mask_num;
    short_type mask_dim;
    tensor_index_to_mask() : mask_num(short_type(-1)), 
			     mask_dim(short_type(-1)) {}
    bool is_valid() { return mask_num != short_type(-1) && 
	mask_dim != short_type(-1); }
  };


  /* 
     defini une "forme" de tenseur creux 
     la fonction merge permet de faire des unions / intersections entre ces formes
  */
  class tensor_shape {
    mutable std::vector<tensor_index_to_mask> idx2mask;
    tensor_mask_container masks_;

    /* verifie si un masque est completement vide,
       si c'est le cas alors tous les autres masques sont vidés
       (le tenseur est identiquement nul) */
    void check_empty_mask() {
      if (card() == 0) {
	for (dim_type i=0; i < masks_.size(); ++i) {
	  masks_[i].set_zero();
	}
      }
    }

    static void find_linked_masks(dim_type mnum, const tensor_shape &ts1, const tensor_shape &ts2, 
				dal::bit_vector& treated1, dal::bit_vector& treated2, 
				std::vector<const tensor_mask*>& lstA,
				std::vector<const tensor_mask*>& lstB) {
      // gare aux boucles infinies si aucun des indices n'est valide
      assert(mnum < ts1.masks().size());
      assert(!treated1[mnum]);
      treated1.add(mnum);
      lstA.push_back(&ts1.mask(mnum));
      for (dim_type i=0; i < ts1.mask(mnum).indexes().size(); ++i) {
	dim_type ii = ts1.mask(mnum).indexes()[i];
	if (ts2.index_is_valid(ii) && !treated2[ts2.index_to_mask_num(ii)])
	  find_linked_masks(ts2.index_to_mask_num(ii),ts2,ts1,treated2,treated1,lstB,lstA);
      }
    }

  protected:
    dim_type index_to_mask_num(dim_type ii) const { 
      if (index_is_valid(ii)) return idx2mask[ii].mask_num; else return dim_type(-1); 
    }
  public:
    void clear() { masks_.resize(0); idx2mask.resize(0); }
    void swap(tensor_shape& ts) {
      idx2mask.swap(ts.idx2mask);
      masks_.swap(ts.masks_);
    }
    dim_type ndim() const { return idx2mask.size(); }
    bool index_is_valid(dim_type ii) const {  
      assert(ii < idx2mask.size()); return idx2mask[ii].is_valid(); 
    }
    const tensor_mask& index_to_mask(dim_type ii) const { 
      assert(index_is_valid(ii)); return masks_[idx2mask[ii].mask_num]; 
    }
    dim_type index_to_mask_dim(dim_type ii) const { 
      assert(index_is_valid(ii)); return idx2mask[ii].mask_dim; 
    }
    index_type dim(dim_type ii) const 
    { assert(index_is_valid(ii)); return index_to_mask(ii).ranges()[index_to_mask_dim(ii)]; 
    }
    tensor_mask_container& masks() { return masks_; }
    const tensor_mask_container& masks() const { return masks_; }
    const tensor_mask& mask(dim_type i) const { assert(i<masks_.size()); return masks_[i]; }
    stride_type card(bool just_look=false) const { 
      stride_type n = 1; 
      for (dim_type i=0; i < masks().size(); ++i) 
	n *= masks()[i].card(just_look); 
      return n; 
    }    
    void push_mask(const tensor_mask& m) { masks_.push_back(m); update_idx2mask(); }
    void remove_mask(dim_type mdim) { /* be careful with this function.. remove
					 only useless mask ! */
      masks_.erase(masks_.begin()+mdim);
      update_idx2mask();
    }
    void remove_unused_dimensions() {
      dim_type nd = 0;
      for (dim_type i=0; i < ndim(); ++i) {
	if (index_is_valid(i)) {
	  masks_[idx2mask[i].mask_num].indexes()[idx2mask[i].mask_dim] = nd++;
	}
      }
      set_ndim_noclean(nd);
      update_idx2mask();
    }

    void update_idx2mask() const {
      /*
	dim_type N=0;
	for (dim_type i=0; i < masks_.size(); ++i) {
	N = std::max(N, std::max_element(masks_.indexes().begin(), masks_.indexes.end()));
	}
	idx2mask.resize(N); 
      */

      std::fill(idx2mask.begin(), idx2mask.end(), tensor_index_to_mask());
      for (dim_type i=0; i < masks_.size(); ++i) {
	for (dim_type j=0; j < masks_[i].indexes().size(); ++j) {
	  dim_type k = masks_[i].indexes()[j];
	  if (k >= idx2mask.size() || idx2mask[k].is_valid()) DAL_INTERNAL_ERROR("");
	  idx2mask[k].mask_num = i; idx2mask[k].mask_dim = j;
	}
      }
    }
    void assign_shape(const tensor_shape& other) { 
      masks_ = other.masks_;
      idx2mask = other.idx2mask;
      //      update_idx2mask(); 
    }
    void set_ndim(dim_type n) {
      clear();
      idx2mask.resize(n); update_idx2mask();
    }
    void set_ndim_noclean(dim_type n) {idx2mask.resize(n);}

    tensor_shape() {}
    
    /* create an "empty" shape of dimensions nd */
    explicit tensor_shape(dim_type nd) : idx2mask(nd,tensor_index_to_mask()) {
      masks_.reserve(16);
    }
    explicit tensor_shape(const tensor_ranges& r) {
      masks_.reserve(16);
      set_full(r);
    }
    void set_full(const tensor_ranges& r) {
      idx2mask.resize(r.size());
      masks_.resize(r.size());
      for (dim_type i=0; i < r.size(); ++i) masks_[i].set_full(i,r[i]);
      update_idx2mask();
    }

    void set_empty(const tensor_ranges& r) { 
      idx2mask.resize(r.size());
      masks_.resize(r.size());
      for (dim_type i=0; i < r.size(); ++i) masks_[i].set_empty(i,r[i]);
      update_idx2mask();
    }


    /* fusion d'une autre forme */
    void merge(const tensor_shape &ts2, bool and_op = true) {
      /* quelques verifs de base */
      if (ts2.ndim() != ndim()) DAL_INTERNAL_ERROR("");
      if (ts2.ndim()==0) return; /* c'est un scalaire */
      for (index_type i = 0; i < ndim(); ++i) 
	if (index_is_valid(i) && ts2.index_is_valid(i)) 
	  if (ts2.dim(i) != dim(i)) DAL_INTERNAL_ERROR("");

      tensor_mask_container new_mask;
      dal::bit_vector mask_treated1; mask_treated1.sup(0,masks().size());
      dal::bit_vector mask_treated2; mask_treated2.sup(0,ts2.masks().size());
      std::vector<const tensor_mask*> lstA, lstB; lstA.reserve(10); lstB.reserve(10);
      for (index_type i = 0; i < ndim(); ++i) {
	dim_type i1 = index_to_mask_num(i);
	dim_type i2 = ts2.index_to_mask_num(i);
	lstA.clear(); lstB.clear();
	if (index_is_valid(i) && !mask_treated1[i1])
	  find_linked_masks(i1, *this, ts2, mask_treated1, mask_treated2,
			    lstA, lstB);
	else if (ts2.index_is_valid(i) && !mask_treated2[i2])
	  find_linked_masks(i2, ts2, *this, mask_treated2, mask_treated1,
			    lstB, lstA);
	else continue;
	if (!(lstA.size() || lstB.size())) DAL_INTERNAL_ERROR("");
	new_mask.push_back(tensor_mask(lstA,lstB,and_op));
      }
      masks_ = new_mask;
      update_idx2mask();
      check_empty_mask();
    }

    void shift_dim_num_ge(dim_type dim_num, int shift) {
      for (dim_type m = 0; m < masks().size(); ++m) {
	masks()[m].shift_dim_num_ge(dim_num,shift);
      }
    }
    /* the permutation vector might be greater than the current ndim,
       in which case some indexes will be unused (when p[i] == dim_type(-1))
    */
    void permute(const std::vector<dim_type> p, bool revert=false) {
      std::vector<dim_type> invp(ndim()); std::fill(invp.begin(), invp.end(), dim_type(-1));

      /* build the inverse permutation and check that this IS really a permuation */
      for (dim_type i=0; i < p.size(); ++i) {
	if (p[i] != dim_type(-1)) { assert(invp[p[i]] == dim_type(-1)); invp[p[i]] = i; }
      }
      for (dim_type i=0; i < invp.size(); ++i) assert(invp[i] != dim_type(-1));
      
      /* do the permutation (quite simple!) */
      for (dim_type m=0; m < masks().size(); ++m) {
	for (dim_type i=0; i < masks()[m].indexes().size(); ++i) {
	  if (!revert) {
	    masks()[m].indexes()[i] = invp[masks()[m].indexes()[i]];
	  } else {
	    masks()[m].indexes()[i] = p[masks()[m].indexes()[i]];
	  }
	}
      }
      set_ndim_noclean(p.size());
      update_idx2mask();
    }

    /* forme d'une tranche (c'est la forme qu'on applique à un tenseur pour
       en extraire la tranche) */
    tensor_shape slice_shape(tensor_mask::Slice slice) const {
      assert(slice.dim < ndim() && slice.i0 < dim(slice.dim));
      tensor_shape ts(ndim());
      ts.push_mask(tensor_mask(dim(slice.dim), slice));
      ts.merge(*this); /* le masque peut se retrouver brutalement vidé si on a tranché au mauvais endroit! */
      return ts;
    }

    tensor_shape diag_shape(tensor_mask::Diagonal diag) const {
      assert(diag.i1 != diag.i0 && diag.i0 < ndim() && diag.i1 < ndim());
      assert(dim(diag.i0) == dim(diag.i1));
      tensor_shape ts(ndim());
      ts.push_mask(tensor_mask(dim(diag.i0), diag));
      ts.merge(*this);
      return ts;
    }

    /*
      void diag(index_type i0, index_type i1) {
      assert(i0 < idx.size() && i1 < idx.size());
      assert(idx[i0].n == idx[i1].n);
      tensor_shape ts2 = *this;
      ts2.masks.resize(1);
      ts2.masks[0].set_diagonal(idx[i0].n, i0, i1);
      ts2.idx[i0].mask_num = ts2.idx[i1].mask_num = 0;
      ts2.idx[i0].mask_dim = 0; ts2.idx[i1].mask_dim = 1;      
      }
    */
    void print(std::ostream& o) const;
    void print_() const { print(cerr); }
  };


  /* reference to a tensor: 
     - a shape
     - a data pointer
     - a set of strides
  */
  class tensor_ref : public tensor_shape {
    std::vector< tensor_strides > strides_;
    TDIter *pbase_; /* pointeur sur un pointeur qui designe les données
		       ça permet de changer la base pour toute une serie
		       de tensor_ref en un coup */

    stride_type base_shift_;

    void remove_mask(dim_type mdim) {
      tensor_shape::remove_mask(mdim);
      assert(strides_[mdim].size() == 0 ||
	     (strides_[mdim].size() == 1 && strides_[mdim][0] == 0)); /* sanity check.. */
      strides_.erase(strides_.begin()+mdim);
    }
  public:
    void swap(tensor_ref& tr) {
      tensor_shape::swap(tr);
      strides_.swap(tr.strides_);
      std::swap(pbase_, tr.pbase_);
      std::swap(base_shift_, tr.base_shift_);
    }
    const std::vector< tensor_strides >& strides() const { return strides_; }
    std::vector< tensor_strides >& strides() { return strides_; }
    TDIter base() const { return (pbase_ ? (*pbase_) : 0); }
    TDIter *pbase() const { return pbase_; }
    stride_type base_shift() const { return base_shift_; }
    void set_base(TDIter &new_base) { pbase_ = &new_base; base_shift_ = 0; }

    void clear() { strides_.resize(0); pbase_ = 0; base_shift_ = 0; tensor_shape::clear(); }

    

    /* s'assure que le stride du premier indice est toujours bien égal à zéro */
    void  ensure_0_stride() {
      for (index_type i=0; i < strides_.size(); ++i) {
	if (strides_[i].size() >= 1 && strides_[i][0] != 0) {
	  stride_type s = strides_[i][0];
	  base_shift_ += s;
	  for (index_type j=0; j < strides_[i].size(); ++j) strides_[i][j] -= s;
	}
      }
    }

    /* constructeur à partir d'une forme : ATTENTION ce constructeur n'alloue pas la
       mémoire nécessaire pour les données !! */
    explicit tensor_ref(const tensor_shape& ts) : tensor_shape(ts), pbase_(0), base_shift_(0) {
      strides_.reserve(16);
      init_strides();
    }
    explicit tensor_ref(const tensor_ranges& r, TDIter *pbase__=0) 
      : tensor_shape(r), pbase_(pbase__), base_shift_(0) {
      strides_.reserve(16);
      init_strides();
    }
    void init_strides() {
      strides_.resize(masks().size());
      stride_type s = 1;
      for (dim_type i = 0; i < strides_.size(); ++i) {
	index_type n = mask(i).card();
	strides_[i].resize(n);
	for (index_type j=0;j<n;++j) strides_[i][j] = j*s;
	s *= n;
      }
    }
    tensor_ref() : pbase_(0), base_shift_(0) { strides_.reserve(16); }

    void set_sub_tensor(const tensor_ref& tr, const tensor_shape& sub);

    /* constructeur à partir d'un sous-tenseur à partir d'un tenseur et d'une forme 
       hypothese: la forme 'sub' doit être un sous-ensemble de la forme du tenseur
    */
    explicit tensor_ref(const tensor_ref& tr, const tensor_shape& sub) {
      set_sub_tensor(tr,sub);
    }

    /* slices a tensor_ref, at dimension 'dim', position 'islice'
       ... not straightforward for sparse tensors !
    */
    explicit tensor_ref(const tensor_ref& tr, tensor_mask::Slice slice);

    /* create a diagonal of another tensor */
    explicit tensor_ref(const tensor_ref& tr, tensor_mask::Diagonal diag) {
      set_sub_tensor(tr, tr.diag_shape(diag));
      ensure_0_stride();
    }

    void print(std::ostream& o) const;

    void print_() const { print(cerr); }
  };
    
    std::ostream& operator<<(std::ostream& o, const tensor_mask& m);
    std::ostream& operator<<(std::ostream& o, const tensor_shape& ts);
    std::ostream& operator<<(std::ostream& o, const tensor_ref& tr);

  /* minimalistic data for iterations */
  struct packed_range {
    const stride_type *pinc;
    const stride_type *begin, *end;
    index_type n;
    /*    index_type cnt;*/
  };
  /* additionnal data */
  struct packed_range_info {
    index_type range;
    dim_type original_masknum;
    dim_type n;
    std::vector<stride_type> mask_pos; /* pour l'iteration avec maj de la valeur des indices */
    bool operator<(const packed_range_info& pi) const {
      if (n < pi.n) return true;
      else return false;
    }
    stride_type mean_increm; /* valeur moyenne de l'increment (utilisé pour le tri) */
    tensor_strides inc; /* not strides but increments to the next index value,
				     with inc[range-1] == -sum(inc[0..range-2]) (automatic rewinding!) 
				     of course, stride_type MUST be signed
				  */
    std::bitset<32> have_regular_strides;
  };

  /* the big one */
  class multi_tensor_iterator {
    index_type N; /* number of simultaneous tensors */
    std::vector<packed_range> pr;
    std::vector<packed_range_info> pri;

    std::vector<index_type> bloc_rank;
    std::vector<index_type> bloc_nelt;

    std::vector<TDIter> it;
    std::vector<TDIter*> pit0;
    tensor_strides itbase;
    struct  index_value_data {
      dim_type cnt_num;
      const stride_type **ppinc; /* pointe vers pr[cnt_num].pinc, initialisé par rewind()
				  et pas avant (à cause de pbs lors de la copie de multi_tensor_iterator sinon) 
				  permet de déduire la valeur du compteur: (*ppinc - pincbase) (à diviser par nn=(pri[cnt_num].n-N))
			       */
      const stride_type *pincbase;
      const stride_type *pposbase; /* pointe dans pri[cnt_num].mask_pos, retrouve la position dans le masque en fonction
				  du compteur déduit ci-dessus et des champs div et mod ci-dessous */
      index_type div, mod, nn;
      stride_type pos_; /* stores the position when the indexe is not part of the pri array
			  (hence the index only has 1 value, and ppos == &pos_, and pcnt = &zero */
    };
    std::vector<index_value_data> idxval;
    std::vector<stride_type> vectorized_strides_; /* if the tensor have regular strides, the mti might be vectorizable */
    index_type vectorized_size_;                 /* the size of each vectorizable chunk */
    index_type vectorized_pr_dim;                /* all pr[i], i >= vectorized_pr_dim, can be accessed via vectorized_strides */
  public:
    void clear() { 
      N = 0; pr.clear(); pri.clear(); bloc_rank.clear(); bloc_nelt.clear(); 
      it.clear(); pit0.clear(); itbase.clear(); idxval.clear(); 
    }
    void swap(multi_tensor_iterator& m) {
      std::swap(N,m.N);  pr.swap(m.pr);  pri.swap(m.pri);
      bloc_rank.swap(m.bloc_rank); bloc_nelt.swap(m.bloc_nelt);
      it.swap(m.it); pit0.swap(m.pit0); itbase.swap(m.itbase);
      idxval.swap(m.idxval);
    }
    void rewind() { 
      for (dim_type i=0; i < pr.size(); ++i) { 
	pr[i].pinc = pr[i].begin = &pri[i].inc[0]; pr[i].end = pr[i].begin+pri[i].inc.size(); 
      }
      for (dim_type n=0; n < N; ++n) it[n] = *(pit0[n]) + itbase[n];
      for (dim_type i=0; i < idxval.size(); ++i) {
	if (idxval[i].cnt_num != dim_type(-1)) {
	  idxval[i].ppinc = &pr[idxval[i].cnt_num].pinc;
	  idxval[i].pincbase = &pri[idxval[i].cnt_num].inc[0];
	  idxval[i].pposbase = &pri[idxval[i].cnt_num].mask_pos[0];
	  idxval[i].nn = (N-pri[idxval[i].cnt_num].n);
	} else {
	  static const stride_type *null=0;
	  idxval[i].ppinc = &null;
	  idxval[i].pincbase = 0;
	  idxval[i].pposbase = &idxval[i].pos_;
	  idxval[i].nn = 1;
	}
      }
    }
    dim_type ndim() const { return idxval.size(); }
    /* get back the value of an index from then current iterator position */
    index_type index(dim_type ii) {
      index_value_data& iv = idxval[ii];
      index_type cnt = (*iv.ppinc - iv.pincbase)/iv.nn;
      return ((iv.pposbase[cnt]) % iv.mod)/ iv.div;
    }
    index_type vectorized_size() const { return vectorized_size_; }
    const std::vector<stride_type>& vectorized_strides() const { return vectorized_strides_; }
    bool next(unsigned i_stop = unsigned(-1), unsigned i0_ = unsigned(-2)) {//=pr.size()-1) {
      unsigned i0 = (i0_ == unsigned(-2) ? pr.size()-1 : i0_);
      while (i0 != i_stop) {
	for (unsigned n = pr[i0].n; n < N; ++n) {
	  //	  index_type pos = pr[i0].cnt * (N-pri[i0].n) + (n - pri[i0].n);
	  it[n] += *pr[i0].pinc; pr[i0].pinc++; 
	}
	if (pr[i0].pinc != pr[i0].end) {
	  return true;
	} else {
	  pr[i0].pinc = pr[i0].begin; i0--;
	}
      }
      return false;
    }
    bool vnext() { return next(unsigned(-1), vectorized_pr_dim); }
    bool bnext(dim_type b) { return next(bloc_rank[b]-1, bloc_rank[b+1]-1); }
    bool bnext_useful(dim_type b) { return bloc_rank[b] != bloc_rank[b+1]; }
    /* version speciale pour itérer sur des tenseurs de même dimensions
       (doit être un poil plus rapide) */    
    bool qnext1() {
      if (pr.size() == 0) return false;
      std::vector<packed_range>::reverse_iterator p_ = pr.rbegin();
     while (p_!=pr.rend()) {
	it[0] += *(p_->pinc++);
	if (p_->pinc != p_->end) {
	  return true;
	} else {
	  p_->pinc = p_->begin; p_++;
	}
      }
      return false;
    }

    bool qnext2() { 
      if (pr.size() == 0) return false;
      std::vector<packed_range>::reverse_iterator p_ = pr.rbegin();
      while (p_!=pr.rend()) {
	it[0] += *(p_->pinc++);
	it[1] += *(p_->pinc++);
	if (p_->pinc != p_->end) {
	  return true;
	} else {
	  p_->pinc = p_->begin; p_++;
	}
      }
      return false;
    }

    scalar_type& p(dim_type n) { return *it[n]; }

    multi_tensor_iterator() {}
    multi_tensor_iterator(std::vector<tensor_ref> trtab, bool with_index_values) {
      init(trtab, with_index_values);
    }
    void assign(std::vector<tensor_ref> trtab, bool with_index_values) {
      multi_tensor_iterator m(trtab, with_index_values);
      swap(m);
    }
    multi_tensor_iterator(const tensor_ref& tr0, bool with_index_values) {
      std::vector<tensor_ref> trtab(1); trtab[0] = tr0;
      init(trtab, with_index_values);
    }
    void assign(const tensor_ref& tr0, bool with_index_values) {
      multi_tensor_iterator m(tr0, with_index_values);
      swap(m);
    }
    multi_tensor_iterator(const tensor_ref& tr0, 
			  const tensor_ref& tr1,  bool with_index_values) {
      std::vector<tensor_ref> trtab(2); trtab[0] = tr0; trtab[1] = tr1;
      init(trtab, with_index_values);
    }
    void assign(const tensor_ref& tr0, const tensor_ref& tr1,  bool with_index_values) {
      multi_tensor_iterator m(tr0, tr1, with_index_values);
      swap(m);
    }
    multi_tensor_iterator(const tensor_ref& tr0, 
			  const tensor_ref& tr1, 
			  const tensor_ref& tr2, bool with_index_values) {
      std::vector<tensor_ref> trtab(3); trtab[0] = tr0; trtab[1] = tr1; trtab[2] = tr2;
      init(trtab, with_index_values);
    }
    void assign(const tensor_ref& tr0, const tensor_ref& tr1, const tensor_ref& tr2,  bool with_index_values) {
      multi_tensor_iterator m(tr0, tr1, tr2, with_index_values);
      swap(m);
    }
    void init(std::vector<tensor_ref> trtab, bool with_index_values);
    void print() const;
  };


  /* handles a tree of reductions
     The tree is used if more than two tensors are reduced, i.e.
       z(:,:)=t(:,i).u(i,j).v(j,:) 
     in that case, the reduction against j can be performed on u(:,j).v(j,:) = w(:,:)
     and then, z(:,:) = t(:,i).w(i,:) 
  */
  struct tensor_reduction {
    struct tref_or_reduction {
      tensor_ref tr_;
      tensor_reduction *reduction;
      tensor_ref &tr() { return tr_; }
      const tensor_ref &tr() const { return tr_; }
      explicit tref_or_reduction(const tensor_ref &tr__, const std::string& s) 
	: tr_(tr__), reduction(0), ridx(s) {}
      explicit tref_or_reduction(tensor_reduction *p, const std::string& s) 
	: reduction(p), ridx(s) {
	reduction->result(tr_);
      }
      bool is_reduction() const { return reduction != 0; }
      void swap(tref_or_reduction &other) { tr_.swap(other.tr_); std::swap(reduction, other.reduction); }
      std::string ridx;      /* reduction indexes, no index can appear
			      twice in the same tensor */
      std::vector<dim_type> gdim; /* mapping to the global virtual
				     tensor whose range is the
				     union of the ranges of each
				     reduced tensor */
      std::vector<dim_type> rdim; /* mapping to the dimensions of the
				     reduced tensor ( = dim_type(-1) for
				     dimensions i s.t. ridx[i] != ' ' ) */
				     
    };
    tensor_ranges reduced_range;
    std::string reduction_chars; /* list of all indexes used for reduction */
    tensor_ref trres;
    typedef std::vector<tref_or_reduction>::iterator trtab_iterator;
    std::vector<tref_or_reduction> trtab;
    multi_tensor_iterator mti;
    std::vector<scalar_type> out_data; /* optional storage of output */
    TDIter pout_data;
  public:
    tensor_reduction() { clear(); }
    virtual ~tensor_reduction() { clear(); }
    void clear();

    /* renvoie les formes diagonalisées 
       pour bien faire, il faudrait que cette fonction prenne en argument
       le required_shape de l'objet ATN_reducted_tensor, et fasse le merge
       avec ce qu'elle renvoie... non trivial
    */
    static void diag_shape(tensor_shape& ts, const std::string& s) {
      for (index_type i=0; i < s.length(); ++i) {
	size_type pos = s.find(s[i]);
	if (s[i] != ' ' && pos != i) { // ce n'est pas de l'indice => reduction sur la diagonale
	  ts = ts.diag_shape(tensor_mask::Diagonal(pos,i));
	}
      }
    }

    void insert(const tensor_ref& tr_, const std::string& s);
    void prepare(const tensor_ref* tr_out = NULL);
    void do_reduction();
    void result(tensor_ref& res) const {
      res=trres;
      res.remove_unused_dimensions();
    }
  private:
    void insert(const tref_or_reduction& tr_, const std::string& s);
    void update_reduction_chars();
    void pre_prepare();
    void make_sub_reductions();
    size_type find_best_sub_reduction(dal::bit_vector &best_lst, std::string &best_idxset);
  };

} /* namespace bgeot */

#endif
