//========================================================================
//
// Library : Basic GEOmetric Tool (bgeot)
// File    : bgeot_sparse_tensors.cc : tensors with a certain pattern.
//           
// Date    : January 2003.
// Author  : Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2000-2005 Julien Pommier
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

#include <bitset>
#include "gmm_blas_interface.h"
#include "bgeot_sparse_tensors.h"


namespace bgeot {
  std::ostream& operator<<(std::ostream& o, const tensor_ranges& r) {
    for (size_type i=0; i < r.size(); ++i) {
      if (i) o << "x";
      o << "[0.." << r[i] << "]"; 
    }
    return o;
  }

  /*
    strides:
    0 1 2 3 4 5 6 7
    ---------------
    x x   x
    x       x x
    x   x   x x x x
    x   x   x x   x

    => n={0,4,6,8,8}

    basic iterator on a set of full tensors, just used for iterating over binary masks
  */
  template<typename IT> class basic_multi_iterator {
    unsigned N;
    index_set idxnums;
    tensor_ranges ranges;
    tensor_strides strides;
    tensor_ranges cnt;
    index_set ilst2idxnums;
    std::vector<const tensor_strides*> slst;
    std::vector<IT> iter;
    std::vector<int> n;
  public:
    const tensor_ranges& getcnt() const { return cnt; }
    basic_multi_iterator() { 
      N = 0; idxnums.reserve(16); strides.reserve(64); 
      slst.reserve(16); 
      ilst2idxnums.reserve(64); iter.reserve(4);
    }
    const tensor_ranges& all_ranges() const { return ranges; }
    const index_set& all_indexes() const { return idxnums; }
    /* /!!!!!\ attention les strides ont 1 elt de plus que r et idxs (pour les tensor_masks)
       (ça intervient aussi dans prepare()) */
    void insert(const index_set& idxs, const tensor_ranges& r, const tensor_strides& s, IT it_) {
      assert(idxs.size() == r.size()); assert(s.size() == r.size()+1);
      slst.push_back(&s);
      for (unsigned int i=0; i < idxs.size(); ++i) {
	index_set::const_iterator f=std::find(idxnums.begin(), idxnums.end(), idxs[i]);	
	if (f == idxnums.end()) {
	  ilst2idxnums.push_back(idxnums.size());
	  idxnums.push_back(idxs[i]);
	  ranges.push_back(r[i]);
	} else {
	  ilst2idxnums.push_back(f-idxnums.begin());
	  assert(ranges[ilst2idxnums.back()] == r[i]);
	}
      }
      iter.push_back(it_);
      N++;
    }
    void prepare() {
      unsigned c = 0;
      strides.assign(N*idxnums.size(),0);
      for (unsigned int i=0; i < N; ++i) {
	for (unsigned int j=0; j < slst[i]->size()-1; ++j) {
	  stride_type s = (*slst[i])[j];
	  strides[ilst2idxnums[c]*N + i] = s;
	  c++;
	}
      }
      n.assign(N+1,-1);
      for (unsigned int i=0; i < idxnums.size(); ++i) {
	for (unsigned int j=0; j < N; ++j) {
	  if (strides[i*N+j]) n[j+1] = i;
	}
      }
      cnt.assign(idxnums.size(),0);
    }
    /* unfortunatly these two function don't inline
       and are not optimized by the compiler 
       (when b and N are known at compile-time) (at least for gcc-3.2) */
    bool next(unsigned int b) { return next(b,N); }
    bool next(unsigned int b, unsigned int NN) {
      int i0 = n[b+1];
      while (i0 > n[b]) {
	if (++cnt[i0] < ranges[i0]) {
	  for (unsigned int i=b; i < NN; ++i) {
	    iter[i] += strides[i0*NN+i];
	  }
	  return true;
	} else {
 	  for (unsigned int i=b; i < NN; ++i) {
	    iter[i] -= strides[i0*NN+i]*(ranges[i0]-1);
	  }
	  cnt[i0]=0; i0--;
	}
      }
      return false;
    }
    template<unsigned b, unsigned NN> bool qnext() { return next(b,NN); }
    IT it(unsigned int b) const { return iter[b]; }
  };


  /* 
     constructeur par fusion de deux masques binaires
     (avec potentiellement une intersection non vide)
  */
  tensor_mask::tensor_mask(const tensor_mask& tm1, const tensor_mask& tm2, bool and_op) {
    assign(tm1,tm2,and_op); unset_card();
  }
#if 0
  void tensor_mask::assign(const tensor_mask& tm1, const tensor_mask& tm2, bool and_op) {
    clear(); unset_card();
    if (tm1.ndim()==0) { assign(tm2); return; }
    if (tm2.ndim()==0) { assign(tm1); return; }
    r = tm1.ranges();
    idxs = tm1.indexes();
    
    /* ce tableau va permettre de faire une boucle commune sur les
       deux masques, les dimension inutilisées étant mises à 1 pour
       ne pas géner */
    tensor_ranges global_r(std::max(tm1.max_dim(),tm2.max_dim())+1, index_type(1));
 
    for (index_type i = 0; i < tm1.indexes().size(); ++i) 
      global_r[tm1.indexes()[i]] = tm1.ranges()[i];

    /* detect new indices */
    for (index_type i = 0; i < tm2.indexes().size(); ++i) {
      index_set::const_iterator f=std::find(tm1.indexes().begin(), tm1.indexes().end(), tm2.indexes()[i]);
      if (f == tm1.indexes().end()) {
	assert(global_r[tm2.indexes()[i]]==1);
	global_r[tm2.indexes()[i]] = tm2.ranges()[i];
	r.push_back(tm2.ranges()[i]);
	idxs.push_back(tm2.indexes()[i]);	
      }
      else assert(global_r[tm2.indexes()[i]] = tm2.ranges()[i]);
    }
    eval_strides();
    assert(size());
    m.assign(size(),false);
    /* sans doute pas très optimal, mais la fonction n'est pas critique .. */    
    for (tensor_ranges_loop l(global_r); !l.finished(); l.next()) {
      if (and_op) {
	if (tm1(l.cnt) && tm2(l.cnt)) m.add(pos(l.cnt));
      } else {
	if (tm1(l.cnt) || tm2(l.cnt)) m.add(pos(l.cnt));
      }
    }
    unset_card();
  }
#endif


  void tensor_mask::assign(const tensor_mask& tm1, const tensor_mask& tm2, bool and_op) {
    clear(); unset_card();
    /* check simple cases first, since this function can be called often and
       is quite expensive */
    if (tm1.ndim()==0) { assign(tm2); return; }
    if (tm2.ndim()==0) { assign(tm1); return; }

    /* same masks ? */
    if (tm1.indexes() == tm2.indexes() &&
	tm1.ranges() == tm2.ranges() &&
	tm1.strides() == tm2.strides()) {
      r = tm1.ranges(); idxs=tm1.indexes(); s=tm1.strides();
      assert(tm1.m.size() == tm2.m.size()); 
      m = tm1.m;
      if (and_op) {
	for (index_type i=0; i < tm2.m.size(); ++i)
	  if (!tm2.m[i]) m[i] = false;
      } else {
	for (index_type i=0; i < tm2.m.size(); ++i)
	  if (tm2.m[i]) m[i] = true;
      }
      return;
    }

    basic_multi_iterator<unsigned> bmit;
    bmit.insert(tm1.indexes(), tm1.ranges(), tm1.strides(), 0); 
    bmit.insert(tm2.indexes(), tm2.ranges(), tm2.strides(), 0); 
    r = bmit.all_ranges(); idxs = bmit.all_indexes(); eval_strides();
    assert(size());
    m.assign(size(),false);    
    bmit.insert(indexes(), ranges(), strides(), 0);
    bmit.prepare(); 
    //cout << "tm1=" << tm1 << "\ntm2=" << tm2 << endl;
    if (and_op) {
      do {
	if (tm1.m[bmit.it(0)]) {
	  do {
	    if (tm2.m[bmit.it(1)]) {
	      m[bmit.it(2)] = 1;
	    }
	    //	    cout << "at cnt=" << bmit.getcnt() << ", it0=" << bmit.it(0) << "=" << tm1.m[bmit.it(0)]
	    //		 << ", it1=" << bmit.it(1) << "=" << tm2.m[bmit.it(1)] << ", res=" << bmit.it(2) << "=" << m[bmit.it(2)] << endl;
	  } while (bmit.qnext<1,3>());
	}
      } while (bmit.qnext<0,3>());
    } else {
      do {
	bool v1 = tm1.m[bmit.it(0)];
	do {
	  if (v1 || (tm2.m[bmit.it(1)]))
	    m[bmit.it(2)] = 1;
	} while (bmit.qnext<1,3>());
      } while (bmit.qnext<0,3>());
    }
    //    cout << "output: " << *this << endl;
  }


  /* construit un masque formé du produit cartesien de plusieurs masques (disjoints)
     /!\ aucune verif sur la validite des arguments */
  tensor_mask::tensor_mask(const std::vector<const tensor_mask*>& tm) { 
    assign(tm); 
  }
#if 0
  // REMPLACE PAR LA FONCTION SUIVANTE
  void tensor_mask::assign(const std::vector<const tensor_mask*> & tm) {
    index_set mask_start; unset_card();
    clear();
    r.reserve(255); idxs.reserve(255); mask_start.reserve(255);
    for (dim_type i=0; i < tm.size(); ++i) {
      mask_start.push_back(r.size());
      for (dim_type j=0; j < tm[i]->ndim(); ++j) { 
	r.push_back(tm[i]->ranges()[j]);
	idxs.push_back(tm[i]->indexes()[j]);
      }
    }
    eval_strides(); assert(size());
    m.assign(size(),false);
    for (tensor_ranges_loop l(r); !l.finished(); l.next()) {
      bool is_in = true;
      for (dim_type i=0; is_in && i < tm.size(); ++i) {
	index_type s_ = 0; 
	for (dim_type j=0; j < tm[i]->ndim(); ++j) 
	  s_ += l.cnt[j+mask_start[i]]*tm[i]->strides()[j];
	if (!tm[i]->m[s_]) is_in = false;
      }
      if (is_in) m.add(lpos(l.cnt));
    }
  }
#endif
  // PRODUIT CARTESIEN DE MASQUES BINAIRES DISJOINTS
  void tensor_mask::assign(const std::vector<const tensor_mask*> & tm) {
    unset_card();
    if (tm.size() == 0) { clear(); return; }
    if (tm.size() == 1) { assign(*tm[0]); return; }
    clear();
    basic_multi_iterator<unsigned> bmit;
    for (dim_type i=0; i < tm.size(); ++i)
      bmit.insert(tm[i]->indexes(), tm[i]->ranges(), tm[i]->strides(), 0);
    r = bmit.all_ranges(); idxs = bmit.all_indexes(); eval_strides();
    assert(size());
    m.assign(size(), false);
    bmit.insert(indexes(), ranges(), strides(), 0);
    bmit.prepare(); 
    unsigned int N = tm.size();
    bool finished = false;
    while (!finished) {
      bool is_in = true;
      unsigned int b;
      for (b=0; b < N; ++b) {
	if (!tm[b]->m[bmit.it(b)]) {
	  is_in = false; break;
	}
      }
      if (is_in) { m[bmit.it(N)] = 1; b = N-1; }
      while (!finished && !bmit.next(b)) { if (b == 0) finished = true; b--; }
    }
  }

  void tensor_mask::unpack_strides(const tensor_strides& packed, tensor_strides& unpacked) const {
    if (packed.size() != card()) 
      assert(packed.size()==card());
    unpacked.assign(size(),INT_MIN);
    index_type i=0;
    for (tensor_ranges_loop l(r); !l.finished(); l.next()) {
      if (m[lpos(l.cnt)]) unpacked[lpos(l.cnt)] = packed[i++];
    }
  }

  void tensor_mask::check_assertions() const {
    if (r.size() != idxs.size()) DAL_INTERNAL_ERROR("");
    if (s.size() != idxs.size()+1) DAL_INTERNAL_ERROR("");
    if (m.size() != size()) DAL_INTERNAL_ERROR("");
    dal::bit_vector bv;
    for (dim_type i=0; i < idxs.size(); ++i) {
      if (bv.is_in(i)) DAL_INTERNAL_ERROR(""); 
      bv.add(i);
    }
  }

  tensor_mask::tensor_mask(const std::vector<const tensor_mask*> tm1, const std::vector<const tensor_mask*> tm2, bool and_op) {
    assign(tensor_mask(tm1), tensor_mask(tm2), and_op);
  }

  void tensor_ref::set_sub_tensor(const tensor_ref& tr, const tensor_shape& sub) {
    assign_shape(sub);
    /* fusionne sub et tr.shape */
    merge(tr);
    
    /* reserve l'espace pour les strides */
    strides_.resize(masks().size());
    for (dim_type i = 0; i < strides_.size(); ++i)
      strides_[i].resize(mask(i).card());
    
    pbase_ = tr.pbase_; base_shift_ = tr.base_shift();
    
    /*    
    cout << "\n  -> entrée dans set_sub_tensor: " << endl 
	 << "tr.shape=" << (tensor_shape&)(tr) << endl
	 << "     sub=" << sub << endl;
    */
    /* pre-calcul rapide des strides sur tr, pour chacun de ses masques */
    std::vector<tensor_strides> trstrides_unpacked(tr.masks().size());
    for (dim_type im = 0; im < tr.masks().size(); ++im) {
      tr.mask(im).check_assertions();
      tr.mask(im).unpack_strides(tr.strides()[im], trstrides_unpacked[im]);
    }
    
    
    /* on parcours chaque masque (mergé) */
    for (dim_type im = 0; im < masks().size(); ++im) {
      const tensor_mask& m = masks()[im];
      /* par construction, tous les masques de tr sont inclus (ou egaux)
	 dans un masque de l'objet courant 
	 
	 on construit donc la liste des masques de tr qui sont inclus dans m
      */
      std::vector<dim_type> trmasks; trmasks.reserve(tr.masks().size());
      for (dim_type i=0; i < m.indexes().size(); ++i) {
	if (tr.index_is_valid(m.indexes()[i])) { 	/* l'index n'est pas forcément valide pour tr !*/
	  dim_type im2 = tr.index_to_mask_num(m.indexes()[i]);
	  if (std::find(trmasks.begin(), trmasks.end(), im2)==trmasks.end()) trmasks.push_back(im2);
	}
      }
      
      tensor_ranges gcnt(tr.ndim(),0);
      stride_type stcnt = 0;

      for (tensor_ranges_loop l(m.ranges()); !l.finished(); l.next()) {
	/* recopie les indice de la boucle dans les indices globaux */
	for (dim_type i=0; i < m.ranges().size(); ++i) gcnt[m.indexes()[i]] = l.index(i);
	
	bool in_m = m(gcnt);
	bool in_trm = true;
	stride_type tr_s = 0;
	
	/* verifie si l'element est bien marqué non nul dans les masques de tr
	     et met à jour le stride */
	for (dim_type i=0; i < trmasks.size(); ++i) { 
	  const tensor_mask &mm = tr.mask(trmasks[i]);

	  //cout << "  mm=" << mm << endl << "gcnt=" << gcnt << endl;
	  if (mm(gcnt)) {
	    tr_s += trstrides_unpacked[trmasks[i]][mm.pos(gcnt)];
	    assert(trstrides_unpacked[trmasks[i]][mm.pos(gcnt)]>=0); // --- DEBUG --- 
	  } else { in_trm = false; break; }
	}
	  /* verifie que m est un sous-ensemble des masques de trmasks */
	if (in_m) assert(in_trm);
	if (!in_trm) assert(!in_m);
	/* recopie le stride si l'element est non nul dans m */
	if (in_m) {
	  //	  cout << "ajout du " << stcnt << "eme elt @ stride=" << tr_s << endl;
	    strides_[im][stcnt++] = tr_s;
	}
      }
      
      /* verif que yapa bug */
      assert(stcnt == stride_type(m.card()));
      }
    ensure_0_stride(); /* c'est plus propre comme ça */
  }

  /* slices a tensor_ref, at dimension 'dim', position 'islice'
     ... not straightforward for sparse tensors !
  */
  tensor_ref::tensor_ref(const tensor_ref& tr, tensor_mask::Slice slice) {
    set_sub_tensor(tr, tr.slice_shape(slice));

    /* shift the base according to the old stride */
    ensure_0_stride();
    
    /* create a mask m2 with one less dimension than m1 */
    const tensor_mask& m1(index_to_mask(slice.dim));
    dim_type mdim = index_to_mask_dim(slice.dim);
    if (m1.ndim() > 1) { 
      tensor_ranges r(m1.ranges()); r.erase(r.begin()+mdim);
      index_set idx(m1.indexes()); idx.erase(idx.begin()+mdim);
      tensor_mask m2(r,idx);
      index_type pos1 = 0, pos2 = 0;
      for (tensor_ranges_loop l(m1.ranges()); !l.finished(); l.next()) {
	if (l.index(mdim) == slice.i0) {
	  m2.set_mask_val(pos2++, m1(pos1));
	} else assert(m1(pos1) == 0);
	pos1++;
      }
      
      
      /* replace the old mask by the new one */
      assert(index_to_mask_num(slice.dim) < masks().size());
      masks()[index_to_mask_num(slice.dim)] = m2;
    } else {
      /* simply remove the mask since it only contained the dimension 'dim' */
      remove_mask(index_to_mask_num(slice.dim));
    }
    /* shift all indexes greater than dim */
    shift_dim_num_ge(slice.dim,-1);
    set_ndim_noclean(ndim()-1);
    update_idx2mask();
  }

  
  struct compare_packed_range {
    const std::vector<packed_range_info>& pri;
    std::vector<stride_type> mean_inc;
    compare_packed_range(const std::vector<packed_range_info>& pri_) : pri(pri_) {}
    bool operator()(dim_type a, dim_type b) {
      if (pri[a].n < pri[b].n) return true;
      else if (pri[a].n > pri[b].n) return false;
      else { /* c'est la que ça devient interessant */
	if (pri[a].mean_increm > pri[b].mean_increm) return true;
      }
      return false;      
    }
  };

  void multi_tensor_iterator::init(std::vector<tensor_ref> trtab, bool with_index_values) {
    N = trtab.size();
    index_type N_packed_idx = 0;

    /* non null elements across all tensors */
      
    tensor_shape ts(trtab[0]);
    for (dim_type i=1; i < trtab.size(); ++i)
      ts.merge(trtab[i]);
    
    
    /* apply the mask to all tensors, 
       updating strides according to it */
    for (dim_type i = 0; i < N; ++i) {
      tensor_ref tmp = trtab[i];
      trtab[i].set_sub_tensor(tmp,ts);
    }

    /* find significant indexes (ie remove indexes who only address 1 element) */
    dal::bit_vector packed_idx; packed_idx.sup(0,ts.ndim());
    for (index_type mi=0; mi < ts.masks().size(); ++mi) {
      if (ts.masks()[mi].card() != 1) {
	packed_idx.add(mi);
	N_packed_idx++;
      } else {
	/* sanity check (TODO: s'assurer que sub_tensor_ref appelle bien ensure_0_stride) */
	for (dim_type j=0; j < N; ++j) {
	  if (trtab[j].strides()[mi].size() != 0) {
	    assert(trtab[j].strides()[mi].size() == 1);
	    assert(trtab[j].strides()[mi][0] == 0);
	  }
	}
      }
    }

    pr.resize(N_packed_idx);
    pri.resize(N_packed_idx);

    /* evaluation of "ranks" per indice */

    index_type pmi = 0;
    for (index_type mi = 0; mi < ts.masks().size(); ++mi) {
      if (packed_idx[mi]) {
	index_type n;
	pri[pmi].original_masknum = mi;
	pri[pmi].range = ts.masks()[mi].card();
	for (n = 0; n < N; n++)
	  if (trtab[n].index_is_valid(mi)) break;
	pri[pmi].n = pr[pmi].n = n;
	pmi++;
      }
    }

    /* sort the packed_range_info according to the "ranks" */
    std::sort(pri.begin(), pri.end());


    /* eval bloc ranks */      
    bloc_rank.resize(N+1); bloc_rank[0] = 0;
    bloc_nelt.resize(N+1); bloc_nelt[0] = 0;
    for (index_type i=1; i <= N; ++i) {
      bloc_rank[i] = bloc_rank[i-1];
      while (bloc_rank[i] < pri.size() && pri[bloc_rank[i]].n == i-1) bloc_rank[i]++;
      bloc_nelt[i] = bloc_rank[i] - bloc_rank[i-1];
    }

    /* "package" the strides in structure pri */
    for (pmi = 0; pmi < pri.size(); ++pmi) {
      index_type mi = pri[pmi].original_masknum;
      pri[pmi].mean_increm = 0;
      pri[pmi].inc.assign(pri[pmi].range*(N-pri[pmi].n), 0);
      pri[pmi].have_regular_strides.reset();
      pri[pmi].have_regular_strides = std::bitset<32>((1 << N)-1);
      for (dim_type n=pri[pmi].n; n < N; ++n) {
        index_type pos0 = (n - pri[pmi].n);
	for (index_type i = 0; i < pri[pmi].range; ++i) {
	  index_type pos = i * (N-pri[pmi].n) + pos0;
	  if (i != pri[pmi].range-1) {
	    stride_type increm = trtab[n].strides()[mi][i+1] - trtab[n].strides()[mi][i];
	    pri[pmi].inc[pos]      = increm;
            if (pri[pmi].inc[pos] != pri[pmi].inc[pos0]) 
              pri[pmi].have_regular_strides[n] = false;
	    pri[pmi].mean_increm += increm;
	  } else { pri[pmi].inc[pos] = -trtab[n].strides()[mi][i]; }
	}
      }
      if (pri[pmi].n!=N)
	pri[pmi].mean_increm /= ((N-pri[pmi].n)*(pri[pmi].range-1));
    }

    /* optimize them w/r to mean strides (without modifying the "ranks") */
    index_set pr_reorder(pri.size());
    for (pmi = 0; pmi < pri.size(); ++pmi) {
      pr_reorder[pmi]=pmi;
    }
    std::sort(pr_reorder.begin(), pr_reorder.end(), compare_packed_range(pri));
    {
      std::vector<packed_range> tmppr(pr);
      std::vector<packed_range_info> tmppri(pri);
      for (dim_type i =0; i < pri.size(); ++i) {
	pr[i]  = tmppr [pr_reorder[i]];
	pri[i] = tmppri[pr_reorder[i]];
      }
    }

    /* setup data necessary to get back (quite quickly) the index values while iterating */
    if (with_index_values) {
      idxval.resize(ts.ndim());

      for (dim_type i=0; i < ts.ndim(); ++i) {
	idxval[i].ppinc = NULL; idxval[i].pposbase = NULL; idxval[i].pos_ = 0;
      }

      for (index_type mi = 0; mi < ts.masks().size(); ++mi) {
	tensor_strides v;
	pmi = index_type(-1);
	for (dim_type i=0; i < pr.size(); ++i)
	  if (pri[i].original_masknum == mi) pmi = i;

	if (pmi != index_type(-1)) {
	  ts.masks()[mi].gen_mask_pos(pri[pmi].mask_pos);
	} else { /* not very nice .. */
	  ts.masks()[mi].gen_mask_pos(v); assert(v.size()==1);
	}
	for (dim_type i=0; i < ts.masks()[mi].indexes().size(); ++i) {
	  dim_type ii = ts.masks()[mi].indexes()[i];
	  idxval[ii].cnt_num = pmi; //packed_idx[mi] ? pmi : dim_type(-1);
	  idxval[ii].pos_ = (pmi != index_type(-1)) ? 0 : v[0];
	  idxval[ii].mod = ts.masks()[mi].strides()[i+1];
	  idxval[ii].div = ts.masks()[mi].strides()[i];
	}
      }
    }

    /* check for opportunities to vectorize the loops with the mti 
       (assuming regular strides etc.)
     */
    vectorized_strides_.resize(N); vectorized_size_ = 0;
    std::fill(vectorized_strides_.begin(), vectorized_strides_.end(), 0);
    vectorized_pr_dim = pri.size();
    for (vectorized_pr_dim = pri.size()-1; vectorized_pr_dim != index_type(-1); vectorized_pr_dim--) {
      std::vector<packed_range_info>::const_iterator pp = pri.begin() + vectorized_pr_dim;
      if (vectorized_pr_dim == pri.size()-1) {
        if (pp->have_regular_strides.count() == N) vectorized_size_ = pp->range;
        for (dim_type n=pp->n; n < N; ++n) {
          vectorized_strides_[n] = pp->inc[n];
        }
      } else {
        if (pp->have_regular_strides.count() != N) break;
        bool still_ok = true;
        for (dim_type n=pp->n; n < N; ++n) {
          if (stride_type(vectorized_strides_[n]*vectorized_size_) != pp->inc[n]) still_ok = false;
        }
        if (still_ok) {
          vectorized_size_ *= pp->range;
        } else break;
      }
    }

    it.resize(N); pit0.resize(N); itbase.resize(N);
    for (dim_type n=0; n < N; ++n) {
      pit0[n]=trtab[n].pbase();
      itbase[n]=trtab[n].base_shift();
    }
    rewind();
  }

  void multi_tensor_iterator::print() const {
    cout << "MTI(N=" << N << "): "; 
    for (dim_type i=0; i < pr.size(); ++i) 
      cout << "  pri[" << int(i) << "]: n=" << int(pri[i].n) 
           << ", range=" << pri[i].range << ", mean_increm=" 
           << pri[i].mean_increm << ", regular = " << pri[i].have_regular_strides 
           << ", inc=" << pri[i].inc << "\n";
    cout << "bloc_rank: " << bloc_rank << ", bloc_nelt: " << bloc_nelt << "\n";    
    cout << "vectorized_size : " << vectorized_size_ << ", strides = " << vectorized_strides_ << ", pr_dim=" << vectorized_pr_dim << "\n";
  }

  void tensor_reduction::insert(const tensor_ref& tr_, const std::string& s) {
    tensor_shape ts(tr_); 
    diag_shape(ts,s);
    trtab.push_back(tref_or_reduction(tensor_ref(tr_, ts), s));
    //cout << "reduction.insert('" << s << "')\n";
    //cout << "Reduction: INSERT tr(ndim=" << int(tr_.ndim()) << ", s='" << s << "')\n";
    //cout << "Content: " << tr_ << endl;
    //cout << "shape: " << (tensor_shape&)tr_ << endl;
  }
  
  void tensor_reduction::insert(const tref_or_reduction& t, const std::string& s) {
    if (!t.is_reduction()) {
      insert(t.tr(),s);
    } else {
      trtab.push_back(t); trtab.back().ridx = s;
      //cout << "reduction.insert_reduction('" << s << "')\n";
      //cout << "Reduction: INSERT REDUCTION tr(ndim=" << int(t.tr().ndim()) << ", s='" << s << "')\n";   
    }
  }

  /* ensure that  r(i,i).t(j,:,j,j,k,o) 
     becomes      r(i,A).t(j,:,B,C,k,o)
     and updates reduction_chars accordingly
  */
  void tensor_reduction::update_reduction_chars() {
    reduction_chars.clear();
    for (trtab_iterator it = trtab.begin(); it != trtab.end(); ++it) {
      assert(it->ridx.size() == it->tr().ndim());
      for (unsigned i =0; i < it->ridx.size(); ++i) {
	if (it->ridx[i] != ' ' &&
	    reduction_chars.find(it->ridx[i]) == std::string::npos)
	  reduction_chars.push_back(it->ridx[i]);
      }
    }
    /* for each tensor, if a diagonal reduction inside the tensor is used,
       the mask of the tensor has been 'and'-ed with a diagonal mask
       and a second 'virtual' reduction index is used */    
    for (trtab_iterator it = trtab.begin(); it != trtab.end(); ++it) {
      it->gdim.resize(it->ridx.size());
      for (unsigned i =0; i < it->ridx.size(); ++i) {
        char c = it->ridx[i];
	if (c != ' ' && it->ridx.find(c) != i) { 
          for (c = 'A'; c <= 'Z'; ++c)
            if (reduction_chars.find(c) == std::string::npos) break;
          it->ridx[i] = c;
          reduction_chars.push_back(it->ridx[i]);
        }
      }
    }
  }

  /* 
     initialize 'reduced_range' and it->rdim 
  */
  void tensor_reduction::pre_prepare() {
    for (trtab_iterator it = trtab.begin(); it != trtab.end(); ++it) {
      assert(it->ridx.size() == it->tr().ndim());
      it->rdim.resize(it->ridx.size());
      //cout << " rix = '" << it->ridx << "'\n";
      for (unsigned i =0; i < it->ridx.size(); ++i) {
	if (it->ridx[i] == ' ') {
	  reduced_range.push_back(it->tr().dim(i));
	  it->rdim[i] = reduced_range.size()-1;
	} else it->rdim[i] = dim_type(-1);
      }
    }
  }

  /* look for a possible sub-reduction on a subset of the tensors.
     returns the subset, and the list of concerned reduction indexes */
  size_type tensor_reduction::find_best_sub_reduction(dal::bit_vector &best_lst, std::string &best_idxset) {
    dal::bit_vector lst;
    std::string idxset;
    best_lst.clear(); best_idxset.clear();

    update_reduction_chars();

    //cout << "find_best_reduction: reduction_chars='" << reduction_chars << "'\n";
    if (trtab.size() > 32) 
      DAL_INTERNAL_ERROR("wow it was assumed that nobody would ever need a reduction on more than 32 tensors..");

    std::vector<std::bitset<32> > idx_occurences(reduction_chars.size());

    for (unsigned ir=0; ir < reduction_chars.size(); ++ir) {
      char c = reduction_chars[ir];
      for (unsigned tnum=0; tnum < trtab.size(); ++tnum)
	idx_occurences[ir][tnum] = (trtab[tnum].ridx.find(c) != std::string::npos); 
      //cout << "find_best_reduction: idx_occurences[" << ir << "] = " << idx_occurences[ir] << "\n";
    }
    size_type best_redsz = 100000000;
    for (unsigned ir=0; ir < reduction_chars.size(); ++ir) {
      lst.clear(); lst.add(ir);
      idxset.resize(0); idxset.push_back(reduction_chars[ir]);
      /* add other possible reductions */
      for (unsigned ir2=0; ir2 < reduction_chars.size(); ++ir2) {
	if (ir2 != ir) {
	  if ((idx_occurences[ir2] | idx_occurences[ir]) == idx_occurences[ir]) {
	    lst.add(ir2);
	    idxset.push_back(reduction_chars[ir2]);
	  }
	}
      }
      /* evaluate the cost */
      size_type redsz = 1;
      for (unsigned tnum=0; tnum < trtab.size(); ++tnum) {
	if (!idx_occurences[ir][tnum]) continue;
        std::bitset<32> once(reduction_chars.size());
	for (size_type i=0; i < trtab[tnum].tr().ndim(); ++i) {
	  bool ignore = false;
	  for (dal::bv_visitor j(lst); !j.finished(); ++j) {
	    if (trtab[tnum].ridx[i] == reduction_chars[j]) 
              if (once[j]) ignore = true; else once[j] = true;
          }
	  if (!ignore)
	    redsz *= trtab[tnum].tr().dim(i);
	}
      }
      //cout << "   test " << reduction_chars[ir] << ": lst=" << lst << ", redsz=" << redsz << "\n";
      if (redsz < best_redsz) {
	best_redsz = redsz;
        best_lst.clear();
        for (unsigned i=0; i < trtab.size(); ++i)
          if (idx_occurences[ir][i]) best_lst.add(i);
	best_idxset = idxset;
      }
    }
    /*cout << "find_best_reduction: lst = " << best_lst << " [nt=" 
	 << trtab.size() << "], idx_set='" << best_idxset 
	 << "', redsz=" << best_redsz << "\n";
    */
    return best_redsz;
  }

  void tensor_reduction::make_sub_reductions() {
    dal::bit_vector bv; std::string red;
    int iter = 1;
    do {
      find_best_sub_reduction(bv,red);
      if (bv.card() < trtab.size() && red.size()) {
        //cout << "making sub reduction\n";
        tensor_reduction *sub = new tensor_reduction();
        std::vector<dim_type> new_rdim; new_rdim.reserve(16);
	std::string new_reduction;
        for (dal::bv_visitor tnum(bv); !tnum.finished(); ++tnum) {
	  tref_or_reduction &t = trtab[tnum];
          std::string re = t.ridx; t.ridx.clear();
          for (unsigned i = 0; i < re.size(); ++i) {
            bool reduced = false;
	    char c = re[i];
            if (c != ' ') {
              if (red.find(re[i]) == std::string::npos)  c = ' ';
              else reduced = true; 
            }
            if (!reduced) { 
              t.ridx.push_back(re[i]);
	      new_rdim.push_back(t.rdim[i]);
	      new_reduction.push_back(re[i]);
	    }
	    re[i] = c;
          }
          //cout << "  sub-";
          sub->insert(trtab[tnum], re);
        }
	//cout << "  new_reduction = '" << new_reduction << "'\n";
        sub->prepare();
	/*cout << "  " << new_reduction.size() << " == " << int(sub->trres.ndim()) << "?\n";
	  assert(new_reduction.size() == sub->trres.ndim());*/
        trtab[bv.first_true()] = tref_or_reduction(sub, new_reduction); 
        trtab[bv.first_true()].rdim = new_rdim;
        std::vector<tref_or_reduction> trtab2; trtab2.reserve(trtab.size() - bv.card() + 1);
        for (size_type i=0; i < trtab.size(); ++i)
          if (!bv.is_in(i) || i == bv.first_true())
            trtab2.push_back(trtab[i]);
        trtab.swap(trtab2);
	//cout << "make_sub_reductions[" << iter << "] : still " << trtab.size() << " tensors\n";
	/*for (size_type i=0; i < trtab.size(); ++i)
	  cout << "    dim = " << trtab[i].tr().ndim() << " : '" << trtab[i].ridx << "'\n";
        */
	++iter;
      } else {
	//cout << "Ok, no more reductions (bv.card() == " << bv.card() << ")\n\n";
	break;
      }
    } while (1);
  }

  void tensor_reduction::prepare(const tensor_ref* tr_out) {
    pre_prepare();
    make_sub_reductions();

    /* create the output tensor */
    if (tr_out == NULL) {
      trres = tensor_ref(reduced_range);
      out_data.resize(trres.card());
      pout_data = &out_data[0];
      trres.set_base(pout_data);
    } else {
      if (tr_out->ndim() != reduced_range.size()) 
	DAL_INTERNAL_ERROR("");
      for (dim_type i=0; i < tr_out->ndim(); ++i) if (tr_out->dim(i) != reduced_range[i]) 
	DAL_INTERNAL_ERROR("");
      trres = *tr_out;
    }

    /* prepare the mapping from each dimension of each tensor to the global range 
       (the first dimensions are reserved for non-reduced dimensions, i.e. those
       of 'reduced_range'
    */
    tensor_ranges global_range; /* global range across all tensors of the 
				   reduction */
    std::string   global_chars; /* name of indexes (or ' ') for each dimension
				   of global_range */
    global_range.reserve(16); 
    global_range.assign(reduced_range.begin(), reduced_range.end());
    global_chars.insert(size_type(0), reduced_range.size(), ' ');    
    for (trtab_iterator it = trtab.begin(); it != trtab.end(); ++it) {
      assert(it->rdim.size() == it->tr().ndim());
      it->gdim = it->rdim;
      for (unsigned i=0; i < it->ridx.size(); ++i) {
        if (it->rdim[i] == dim_type(-1)) {
          assert(it->ridx[i] != ' '); 
          std::string::size_type p = global_chars.find(it->ridx[i]);
          if (p == std::string::npos) {
            global_chars.push_back(it->ridx[i]);
            global_range.push_back(it->tr().dim(i));
            it->gdim[i] = global_range.size() - 1;
          } else {
            if (it->tr().dim(i) != global_range[p]) 
              DAL_THROW(std::invalid_argument, 
                        "inconsistent dimensions for reduction index " 
                        << it->ridx[i] << "(" << int(it->tr().dim(i)) 
                        << " != " << int(global_range[p]) << ")");
            it->gdim[i] = p;
          }
        }
      }
      //cout << " rdim = " << it->rdim << ", gdim = " << it->gdim << "\n";
    }
    //cout << "global_range = " << global_range << ", global_chars = '" << global_chars << "'\n";
    
    std::vector<dim_type> reorder(global_chars.size(), dim_type(-1));
    /* increase the dimension of the tensor holding the result */
    for (dim_type i=0; i < reduced_range.size(); ++i) reorder[i] = i;
    //cout << "reorder = '" << reorder << "'\n";
    trres.permute(reorder);
    std::vector<tensor_ref> tt; tt.reserve(trtab.size()+1);
    tt.push_back(trres);

    /* permute all tensors (and increase their number of dimensions) */
    for (trtab_iterator it = trtab.begin(); it != trtab.end(); ++it) {
      std::fill(reorder.begin(), reorder.end(), dim_type(-1));
      for (dim_type i=0; i < it->gdim.size(); ++i) {
        reorder[it->gdim[i]] = i;
      }
      //cout << "reorder = '" << reorder << "'\n";
      it->tr().permute(reorder);
      tt.push_back(it->tr());
      //cout << "MTI[" << it-trtab.begin() << "/" << trtab.size() << "] : " << (tensor_shape&)it->tr();
    }
    
    /* now, the iterator can be built */
    mti.init(tt,false);
  }
  
  static void do_reduction1(bgeot::multi_tensor_iterator &mti) {
    do {
      scalar_type s1 = 0;
      do { 
        s1 += mti.p(1); 
      } while (mti.bnext(1));
      mti.p(0) += s1;
    } while (mti.bnext(0));
  }

  static bool do_reduction2v(bgeot::multi_tensor_iterator &mti) {
    int n = mti.vectorized_size();
    const std::vector<stride_type> &s = mti.vectorized_strides();
    if (n && s[0] && s[1] && s[2] == 0) {
      int incx = s[1], incy = s[0];
      /*mti.print();
        scalar_type *b[3]; 
        for (int i=0; i < 3; ++i)       b[i] = &mti.p(i);*/
      do {
        /*cout << "vectorized_ reduction2a : n=" << n << ", s = " << s << " mti.p=" << &mti.p(0)-b[0] << "," 
          << &mti.p(1)-b[1] << "," << &mti.p(2)-b[2] << "\n";*/
	gmm::daxpy_(&n, &mti.p(2), &mti.p(1), &incx, &mti.p(0), &incy);
      } while (mti.vnext());
      return true;
    } else return false;
  }
  static void do_reduction2a(bgeot::multi_tensor_iterator &mti) {
    if (!do_reduction2v(mti)) {
      do {
        mti.p(0) += mti.p(1)*mti.p(2);
      } while (mti.bnext(0));
    }
  }

  static void do_reduction2b(bgeot::multi_tensor_iterator &mti) {
    do {
      scalar_type s1 = 0;
      do { 
        scalar_type s2 = 0;
        do {
          s2 += mti.p(2);
        } while (mti.bnext(2));
        s1 += mti.p(1)*s2; 
      } while (mti.bnext(1));
      mti.p(0) += s1;
    } while (mti.bnext(0));    
  }

  static bool do_reduction3v(bgeot::multi_tensor_iterator &mti) {
    int n = mti.vectorized_size();
    const std::vector<stride_type> &s = mti.vectorized_strides();
    if (n && s[0] && s[1] && s[2] == 0 && s[3] == 0) {
      int incx = s[1], incy = s[0];
      do {
        double v = mti.p(2)*mti.p(3);
	gmm::daxpy_(&n, &v, &mti.p(1), &incx, &mti.p(0), &incy);
      } while (mti.vnext());
      return true;
    } else return false;
  }

  static void do_reduction3a(bgeot::multi_tensor_iterator &mti) {
    if (!do_reduction3v(mti)) {
      do {
        mti.p(0) += mti.p(1)*mti.p(2)*mti.p(3);
      } while (mti.bnext(0));
    }
  }

  static void do_reduction3b(bgeot::multi_tensor_iterator &mti) {
    do {
      scalar_type s1 = 0;
      do { 
        scalar_type s2 = 0;
        do {
          scalar_type s3 = 0;
          do { 
            s3 += mti.p(3);
          } while (mti.bnext(3));
          s2 += mti.p(2)*s3;
        } while (mti.bnext(2));
        s1 += mti.p(1)*s2; 
      } while (mti.bnext(1));
      mti.p(0) += s1;
    } while (mti.bnext(0));	
  }

  void tensor_reduction::do_reduction() {
    /* on s'assure que le tenseur destination a bien été remis à zero
       avant le calcul (c'est obligatoire malheureusement, conséquence
       de l'utilisation de masque qui ne s'arrêtent pas forcement sur les 
       'frontieres' entre les differents tenseurs reduits) */
    //std::fill(out_data.begin(), out_data.end(), 0.);
    memset(&out_data[0], 0, out_data.size()*sizeof(out_data[0]));
    for (unsigned i=0; i < trtab.size(); ++i) {
      if (trtab[i].is_reduction()) { 
        trtab[i].reduction->do_reduction(); 
	trtab[i].reduction->result(trtab[i].tr());
	//cout << "resultat intermediaire: " << trtab[i].tr() << endl;
      }
    }
    mti.rewind();
    dim_type N = trtab.size();
    if (N == 1) {
      do_reduction1(mti);
    } else if (N == 2) {
      if (!mti.bnext_useful(2) && !mti.bnext_useful(1)) {
        do_reduction2a(mti);
      } else {
        do_reduction2b(mti);
      }
    } else if (N == 3) {
      if (!mti.bnext_useful(1) && (!mti.bnext_useful(2)) && !mti.bnext_useful(3)) {
        do_reduction3a(mti);
      } else {
        do_reduction3b(mti);
      }
    } else if (N == 4) {
      do {
	scalar_type s1 = 0;
	do { 
	  scalar_type s2 = 0;
	  do {
	    scalar_type s3 = 0;
	    do { 
              scalar_type s4 = 0;
              do {
                s4 += mti.p(4);
              } while (mti.bnext(4));
              s3 += mti.p(3)*s4;
	    } while (mti.bnext(3));
	    s2 += mti.p(2)*s3;
	  } while (mti.bnext(2));
	  s1 += mti.p(1)*s2; 
	} while (mti.bnext(1));
	mti.p(0) += s1;
      } while (mti.bnext(0));	
    } else {
      DAL_THROW(std::invalid_argument, "unhandled reduction case ! (N=" << int(N) << ")");
    }
  }

  void tensor_reduction::clear() {
    for (unsigned i=0; i < trtab.size(); ++i) {
      if (trtab[i].is_reduction()) delete trtab[i].reduction;
      trtab[i].reduction = 0;
    }
    trtab.clear(); 
    trres.clear(); 
    reduced_range.resize(0);
    reduction_chars.clear();

    out_data.resize(0);
    pout_data = 0; trtab.reserve(10);
    mti.clear();
  }


  void tensor_mask::print(std::ostream &o) const {
    index_type c=card(true);
    check_assertions();
    o << "   mask : card=" << c << "(card_=" << card_ << ", uptodate=" << card_uptodate << "), indexes=";
    for (dim_type i=0; i < idxs.size(); ++i) 
      o << (i==0?"":", ") << int(idxs[i]) << ":" << int(r[i]);
    o << "   ";
    if (c == size()) o << " FULL" << endl;
    else {
      o << "m={";
      if (idxs.size() == 1) {
	for (index_type i=0; i < m.size(); ++i) o << (m[i] ? 1 : 0);
      } else {
	for (tensor_ranges_loop l(r); !l.finished(); l.next()) {
	  if (l.cnt[0] == 0 && l.cnt[1] == 0 && r.size()>2) {
	    o << "\n   -> (:,:";
	    for (dim_type i=2; i < r.size(); ++i) o << "," << l.cnt[i]; 
	    o << ")={";
	  }
	  o << (m[lpos(l.cnt)] ? 1 : 0);
	  if (l.cnt[0] == r[0]-1) 
	    if (l.cnt[1] != r[1]-1) o << ","; 
	    else if (idxs.size() > 2) o << "}";
	}
      }
      o << "}" << endl;
    }
  }



  void tensor_shape::print(std::ostream& o) const {
    o << "  tensor_shape: n=" << idx2mask.size() << ", idx2mask=";
    for (dim_type i=0; i < idx2mask.size(); ++i) {
      if (i) o << ",";
      if (idx2mask[i].is_valid()) {
	o << "r" << dim(i) << ":m" << int(idx2mask[i].mask_num) << "/" << int(idx2mask[i].mask_dim);
      } else o << " (na) ";
    }
    o << endl;
    //      o << "     masks[1.."<< masks_.size() << "]={" << endl;
    for (dim_type i=0; i < masks_.size(); ++i) o << masks_[i];
    o << "  ^-- end tensor_shape" << endl;
  }

  void tensor_ref::print(std::ostream& o) const {
    o << "tensor_ref, n=" << int(ndim()) << ", card=" << card(true) << ", base=" << base() << endl;
    for (dim_type i=0; i < strides().size(); ++i) {
      o << " * strides["<<i<<"]={";
      for (size_type j=0; j < strides()[i].size(); ++j) o << (j>0?",":"") << strides()[i][j];
      o << "}" << endl;
    }
    multi_tensor_iterator mti(*this, true);
    do {
      for (dim_type i = 0; i < mti.ndim(); ++i) {
	o << (i==0?"(":",");
	if (index_is_valid(i))
	  o << mti.index(i);
	else o << "*";
      }
      o << ")";
      if (base()) {
	o << " = " << mti.p(0) << "\t@base+" << &mti.p(0) - base();
      } else o << "\t@" << size_t(&mti.p(0))/sizeof(scalar_type);
      o << endl;
    } while (mti.qnext1());

    o << "^---- end tensor_ref" << endl;
  }

  std::ostream& operator<<(std::ostream& o, const tensor_mask& m) { 
    m.print(o); return o; 
  }
  std::ostream& operator<<(std::ostream& o, const tensor_shape& ts) { 
    ts.print(o); return o; 
  }
  std::ostream& operator<<(std::ostream& o, const tensor_ref& tr) {
    tr.print(o); return o;
  }
  void print_tm(const tensor_mask& tm) { tm.print_(); }
  void print_ts(const tensor_shape& ts) { ts.print_(); }
  void print_tr(const tensor_ref& tr) { tr.print_(); }
}
