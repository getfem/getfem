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
    void insert(const index_set& idxs, const tensor_ranges& r, const tensor_strides& s, IT _it) {
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
      iter.push_back(_it);
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
    //cerr << "tm1=" << tm1 << "\ntm2=" << tm2 << endl;
    if (and_op) {
      do {
	if (tm1.m[bmit.it(0)]) {
	  do {
	    if (tm2.m[bmit.it(1)]) {
	      m[bmit.it(2)] = 1;
	    }
	    //	    cerr << "at cnt=" << bmit.getcnt() << ", it0=" << bmit.it(0) << "=" << tm1.m[bmit.it(0)]
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
    //    cerr << "output: " << *this << endl;
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
	index_type _s = 0; 
	for (dim_type j=0; j < tm[i]->ndim(); ++j) 
	  _s += l.cnt[j+mask_start[i]]*tm[i]->strides()[j];
	if (!tm[i]->m[_s]) is_in = false;
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


  tensor_mask::tensor_mask(const std::vector<const tensor_mask*> tm1, const std::vector<const tensor_mask*> tm2, bool and_op) {
    assign(tensor_mask(tm1), tensor_mask(tm2), and_op);
  }

  void tensor_ref::set_sub_tensor(const tensor_ref& tr, const tensor_shape& sub) {
    assign_shape(sub);
    /* fusionne sub et tr.shape */
    merge(tr);
    
    /* reserve l'espace pour les strides */
    _strides.resize(masks().size());
    for (dim_type i = 0; i < _strides.size(); ++i)
      _strides[i].resize(mask(i).card());
    
    _pbase = tr._pbase; _base_shift = tr.base_shift();
    
    /*    
    cerr << "\n  -> entrée dans set_sub_tensor: " << endl 
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

	  //cerr << "  mm=" << mm << endl << "gcnt=" << gcnt << endl;
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
	  //	  cerr << "ajout du " << stcnt << "eme elt @ stride=" << tr_s << endl;
	    _strides[im][stcnt++] = tr_s;
	}
      }
      
      /* verif que yapa bug */
      assert(stcnt == stride_type(m.card()));
      }
    ensure_0_stride(); /* c'est plus propre comme ça */
  }
  
  struct compare_packed_range {
    const std::vector<packed_range_info>& pri;
    std::vector<stride_type> mean_inc;
    compare_packed_range(const std::vector<packed_range_info>& _pri) : pri(_pri) {}
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
	pri[pmi].n = n;
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

      for (dim_type n=pri[pmi].n; n < N; ++n) {
	for (index_type i = 0; i < pri[pmi].range; ++i) {
	  index_type pos = i * (N-pri[pmi].n) + (n - pri[pmi].n);
	  if (i != pri[pmi].range-1) {
	    stride_type increm = trtab[n].strides()[mi][i+1] - trtab[n].strides()[mi][i];
	    pri[pmi].inc[pos]      = increm;
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
	idxval[i].ppinc = NULL; idxval[i].pposbase = NULL; idxval[i]._pos = 0;
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
	  idxval[ii]._pos = (pmi != index_type(-1)) ? 0 : v[0];
	  idxval[ii].mod = ts.masks()[mi].strides()[i+1];
	  idxval[ii].div = ts.masks()[mi].strides()[i];
	}
      }
    }
    it.resize(N); pit0.resize(N); itbase.resize(N);
    for (dim_type n=0; n < N; ++n) {
      pit0[n]=trtab[n].pbase();
      itbase[n]=trtab[n].base_shift();
    }
    rewind();
  }

  void tensor_reduction::insert(const tensor_ref& _tr, const std::string& s) {
    tensor_shape ts(_tr); 
    diag_shape(ts,s);
    trtab.push_back(tensor_ref(_tr, ts));
    tensor_ref& tr = trtab.back();
    /*
    cerr << "Reduction: INSERT tr(ndim=" << int(_tr.ndim()) << ", s='" << s << "')\n";
    cerr << "shape: " << (tensor_shape&)_tr << endl;
    */
    if (s.length() != tr.ndim()) DAL_INTERNAL_ERROR("");
    tr2r_dim.push_back(std::vector<dim_type>(s.length()));
    for (index_type i=0; i < s.length(); ++i) {
      size_type pos = (s[i] == ' ') ? std::string::npos : redidx.find(s[i]);

      /* si on detecte un réduction sur la diagonale, alors il faut compter
	 la deuxième (troisième..) occurence de l'indice comme un nouvel indice de reduc 
	 (dans l'ideal il aurait fallu lui donner un nouveau nom
	 mais sans doit marcher sans rien toucher
      */
      if (s.find(s[i]) != i) pos = std::string::npos; 
      if (pos == std::string::npos) {
	redidx += s[i];
	r.push_back(tr.dim(i));
	tr2r_dim.back()[i] = r.size()-1;
	if (s[i] == ' ') {
	  tr2r_dim[0].push_back(r.size()-1);
	}
      } else {
	if (tr.dim(i) != r[pos]) 
	  DAL_THROW(std::invalid_argument, 
		    "inconsistent dimensions for reduction index " << int(s[i]));
	tr2r_dim.back()[i] = pos;
      }
    }
  }

  void tensor_reduction::prepare(const tensor_ref* tr_out) {
    std::vector<dim_type> reorder(r.size());

    /* create the output tensor */
    tensor_ranges outr(tr2r_dim[0].size());
    for (dim_type i=0; i < tr2r_dim[0].size(); ++i) outr[i] = r[tr2r_dim[0][i]];
    if (tr_out == NULL) {
      trtab[0] = tensor_ref(outr);
      out_data.resize(trtab[0].card());
      pout_data = &out_data[0];
      trtab[0].set_base(pout_data);
      cerr << "resultat de reduction: nb de dimensions=" << int(trtab[0].ndim()) << endl;
    } else {
      if (tr_out->ndim() != outr.size()) 
	DAL_INTERNAL_ERROR("");
      for (dim_type i=0; i < tr_out->ndim(); ++i) if (tr_out->dim(i) != outr[i]) 
	DAL_INTERNAL_ERROR("");
      trtab[0] = *tr_out;
    }

    /* permute all tensors (and increase their number of dimensions) */
    for (dim_type tnum=0; tnum < tr2r_dim.size(); tnum++) {
      std::fill(reorder.begin(), reorder.end(), dim_type(-1));
      for (dim_type i=0; i < tr2r_dim[tnum].size(); ++i) reorder[tr2r_dim[tnum][i]] = i;
      trtab[tnum].permute(reorder);
    }
    mti.init(trtab,false);
  }

  void tensor_reduction::do_reduction() {
    /* on s'assure que le tenseur destination a bien été remis à zero
       avant le calcul (c'est obligatoire malheureusement, conséquence
       de l'utilisation de masque qui ne s'arrêtent pas forcement sur les 
       'frontieres' entre les differents tenseurs reduits) */
    std::fill(out_data.begin(), out_data.end(), 0.);
    mti.rewind();
    dim_type N = trtab.size();
    if (N == 2) {
      do {
	scalar_type s1 = 0;
	do { 
	  s1 += mti.p(1); 
	} while (mti.bnext(1));
	mti.p(0) += s1;
      } while (mti.bnext(0));
    } else if (N == 3) {
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
    } else if (N == 4) {
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
    } else {
      DAL_THROW(std::invalid_argument, "unhandled reduction case ! (N=" << N << ")");
    }
  }

  void tensor_mask::print(std::ostream &o) const {
    index_type c=card(true);
    check_assertions();
    o << "   mask : card=" << c << "(_card=" << _card << ", uptodate=" << card_uptodate << "), indexes=";
    for (dim_type i=0; i < idxs.size(); ++i) 
      o << (i==0?"":", ") << int(idxs[i]) << ":" << int(r[i]);
    o << "   ";
    if (c == size()) o << " FULL" << endl;
    else {
      o << "m={";
      if (idxs.size() == 1) {
	for (index_type i=0; i < m.size(); ++i) o << m[i] ? 1 : 0;
      } else {
	for (tensor_ranges_loop l(r); !l.finished(); l.next()) {
	  if (l.cnt[0] == 0 && l.cnt[1] == 0 && r.size()>2) {
	    o << "\n   -> (:,:"; for (dim_type i=2; i < r.size(); ++i) o << "," << l.cnt[i]; 
	    o << ")={";
	  }
	  o << m[lpos(l.cnt)] ? 1 : 0;
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
    //      o << "     masks[1.."<< _masks.size() << "]={" << endl;
    for (dim_type i=0; i < _masks.size(); ++i) o << _masks[i];
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
  void print_tm(const tensor_mask& tm) { tm._print(); }
  void print_ts(const tensor_shape& ts) { ts._print(); }
  void print_tr(const tensor_ref& tr) { tr._print(); }
}
