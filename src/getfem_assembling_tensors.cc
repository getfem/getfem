// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_assembling_tensors.cc : 
//           
// Date    : January 2003.
// Author  : Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2003-2005 Julien Pommier
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

#include <getfem_assembling_tensors.h>
#include <getfem_mat_elem.h>

extern "C" void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);

namespace getfem {
  size_type vdim_specif_list::nb_mf() const { 
    return std::count_if(begin(),end(),std::mem_fun_ref(&vdim_specif::is_mf_ref));
  }
  size_type vdim_specif_list::nbelt() const { 
    size_type sz = 1;
    for (const_iterator it = begin(); it != end(); ++it) sz *= (*it).dim;
    return sz;
  }
  void vdim_specif_list::build_strides_for_cv(size_type cv, tensor_ranges& r, 
					      std::vector<tensor_strides >& str) const {
    stride_type s = 1, cnt = 0;
    str.resize(size());
    r.resize(size());      
    for (const_iterator it = begin(); it != end(); ++it, ++cnt) {
      if ((*it).is_mf_ref()) {
	r[cnt] = (*it).pmf->nb_dof_of_element(cv); 
	ref_mesh_dof_ind_ct::const_iterator ii = (*it).pmf->ind_dof_of_element(cv).begin();
	str[cnt].resize(r[cnt]);
	for (index_type j=0; j < r[cnt]; ++j) {
	  str[cnt][j] = (*it).pmf->ind_dof_of_element(cv)[j]*s;
	}
      } else {
	r[cnt] = (*it).dim;
	str[cnt].resize(r[cnt]);
	for (index_type j=0; j < (*it).dim; ++j) {
	  str[cnt][j] = j*s;
	}
      }
      s *= (*it).dim;
    }
  }

  void ATN::update_childs_required_shape() {
    for (dim_type i=0; i < nchilds(); ++i) {
      child(i).merge_required_shape(tensor_shape(child(i).ranges()));
    }
  }
  void ATN::set_number(unsigned &gcnt) { 
    if (number_ == unsigned(-1)) {
      for (unsigned i=0; i < nchilds(); ++i) child(i).set_number(gcnt);
      number_ = ++gcnt;
    }
  }

  /* 
     general class for tensor who store their data
  */
  class ATN_tensor_w_data : public ATN_tensor {
    TDIter data_base;
  protected:
    std::vector<scalar_type> data;
  public:
    void reinit();
    void reinit0() { ATN_tensor_w_data::reinit(); std::fill(data.begin(), data.end(),0); }
  };

  /* note that the data is NOT filled with zeros */
  void ATN_tensor_w_data::reinit() {
    tr.assign_shape(req_shape);
    tr.init_strides();
    if (tr.card() > 10000000) {
      cerr << "warning, a tensor of size " << tr.card()
	   << " will be created, it needs " 
	   << tr.card()*sizeof(scalar_type) << " bytes of memory\n";
    }
    if (tr.card() == 0) {
      cerr << "WARNING: tensor " << name() << " will be created with a size of " 
	   << ranges() << " and 0 non-nul elements!" << endl;
    }
    data.resize(tr.card());
    data_base = &data[0];
    tr.set_base(data_base);
  }


  /* 
     general class for the computation of a reduced product of tensors
     (templated by the number of product tensors)
   
     should be very effective.
  */
  typedef std::vector<std::pair<ATN_tensor*,std::string> > reduced_tensor_arg_type;

  class ATN_reduced_tensor : public ATN_tensor_w_data {
    /* used for specification of tensors and reduction indices , see below */
    reduced_tensor_arg_type red;
    bgeot::tensor_reduction tred;    
  public:
    void check_shape_update(size_type , dim_type) {
      shape_updated_ = false;
      for (dim_type i=0; i < nchilds(); ++i) {
	if (child(i).is_shape_updated()) {
	  shape_updated_ = true;
	}
      }
      if (shape_updated_) {
	r_.resize(0);
	for (dim_type i=0; i < nchilds(); ++i) {
	  std::string s = red_n(i);
          if (s.size() != child(i).ranges().size()) {
            ASM_THROW_TENSOR_ERROR("wrong number of indexes for the " << int(i+1) 
                                   << "th argument of the reduction " << name() 
                                   << " (ranges=" << child(i).ranges() << ")");
          }
	  for (size_type j=0; j < s.length(); ++j) {
	    if (s[j] == ' ') r_.push_back(child(i).ranges()[j]);
	  }
	}
      }
    }
    void update_childs_required_shape() {
      /* pourrait être mieux, cf les commentaires de la fonction 
	 tensor_reduction::required_shape */
      for (dim_type n=0; n < nchilds(); ++n) {
	tensor_shape ts(child(n).ranges());
	tensor_ranges rn(child(n).ranges());
	const std::string &s = red[n].second;
	assert(rn.size() == s.size());
	for (unsigned i=0; i < rn.size(); ++i)
	  if (s[i] != ' ') {
	    size_type p = s.find(s[i]);
	    if (p != size_type(-1) && p < i && rn[p] != rn[i]) 
	      ASM_THROW_TENSOR_ERROR("can't reduce the diagonal of a tensor of size " << rn << " with '" << s << "'");
	  }
	//cerr << "ts = " << child(n).ranges() << ", red=" << red[n].second << "\n";
	bgeot::tensor_reduction::diag_shape(ts, red[n].second);
	/*cerr << "REDUCTION '" << red[n].second << "' -> sending required to child#" << int(n) << " " << child(n).name() << ":" << endl;
          cerr << ts << endl;*/
	child(n).merge_required_shape(ts);
	/*cerr << "------>required shape is now: " << child(n).required_shape() << endl;*/
      }
    }

    /* 
       r is a container of pair<vtensor&,std::string> 
       where the strings specify the reduction indices:
     
       if a_{ik}b_{kjl} is reduced against k and l, then the strings are 
       " k" and "k l"
    */
    ATN_reduced_tensor(reduced_tensor_arg_type& r) : red(r) { 
      for (size_type i=0; i < r.size(); ++i) add_child(*red[i].first);
    }
   
    std::string red_n(size_type n) {
      std::string s = red[n].second;
      if (s.length() == 0)
	s.append(red[n].first->ranges().size(), ' ');
      return s;
    }
    void reinit() {
      tred.clear();
      for (dim_type i=0; i < red.size(); ++i) {
	/*	cerr << "ATN_reduced_tensor::reinit : insertion of r(" << red_n(i) 
		<< "), tr[" << red[i].first->ranges() << "\n" << red[i].first->tensor() << endl;*/
	/*if (red[i].first->ranges().size() != red_n(i).length()) {
	  ASM_THROW_TENSOR_ERROR("wrong number of indexes for the " << int(i+1) 
				 << "th argument of the reduction " << name() 
				 << " (ranges=" << red[i].first->ranges() << ")");
                                 }*/
	tred.insert(red[i].first->tensor(), red_n(i));
      }
      /* reserve the memory for the output 
	 the memory is set to zero since the reduction may only affect a subset 
	 of this tensor hence a part of it would not be initialized
      */
      ATN_tensor_w_data::reinit0();
      /* on fournit notre propre tenseur pour stocker les resultats */
      tred.prepare(&tensor());
    }
  private:
    
    void exec_(size_type , dim_type ) {
      std::fill(data.begin(), data.end(), 0.); /* do_reduction ne peut pas le faire puisque ce n'est pas lui 
						  le proprietaire du tenseur de sortie */
      tred.do_reduction();
    }
  };


  /* slice tensor: 
     no need of a temporary storage for the slice, since direct access
     can be provided via strides.
  */
  class ATN_sliced_tensor : public ATN_tensor {
    dim_type slice_dim;
    size_type slice_idx;
  public:
    ATN_sliced_tensor(ATN_tensor& a, dim_type slice_dim_, size_type slice_idx_) : 
      slice_dim(slice_dim_), slice_idx(slice_idx_)  { add_child(a); }
    void check_shape_update(size_type , dim_type) {
      if ((shape_updated_ = child(0).is_shape_updated())) {
	if (slice_dim >= child(0).ranges().size() ||
	    slice_idx >= child(0).ranges()[slice_dim]) {
	  ASM_THROW_TENSOR_ERROR("can't slice tensor " << child(0).ranges() << 
				 " at index " << int(slice_idx) << " of dimension " << int(slice_dim));
	}
	r_ = child(0).ranges(); r_.erase(r_.begin()+slice_dim);
      }
    }
    void update_childs_required_shape() {
      tensor_shape ts = req_shape; 
      ts.set_ndim_noclean(ts.ndim()+1);
      ts.shift_dim_num_ge(slice_dim,+1);
      ts.push_mask(tensor_mask(child(0).ranges()[slice_dim],tensor_mask::Slice(slice_dim, slice_idx)));
      child(0).merge_required_shape(ts);
    }
    void reinit() {
      tensor() = tensor_ref(child(0).tensor(),  tensor_mask::Slice(slice_dim, slice_idx));
    }
  private:
    void exec_(size_type, dim_type) {}
  };

  /* tensor with reoderer indices:
     t{i,j,k} -> t{j,i,k}
     reorder=   0 1 2       1 0 2
  */
  class ATN_permuted_tensor : public ATN_tensor {
    std::vector<dim_type> reorder;
  public:
    /* attention on ne s'assure pas que reorder est une permutation */
    ATN_permuted_tensor(ATN_tensor& a, const std::vector<dim_type>& reorder_) :
      reorder(reorder_)  { add_child(a); }
  private:
    void check_shape_update(size_type , dim_type) {
      if ((shape_updated_ = child(0).is_shape_updated())) {
	if (reorder.size() != child(0).ranges().size())
	  ASM_THROW_TENSOR_ERROR("can't reorder tensor '" << name() << "' of dimensions " 
				 << child(0).ranges() << 
				 " with this permutation: " << reorder);
	r_.resize(reorder.size());
	std::fill(r_.begin(), r_.end(), dim_type(-1));

	/*
	  --- TODO: A VERIFIER !!!!! ---
	*/
	for (size_type i=0; i < reorder.size(); ++i)
	  r_[i] = child(0).ranges()[reorder[i]];
      }
    }
    void update_childs_required_shape() {
      tensor_shape ts = req_shape;
      ts.permute(reorder, true);
      child(0).merge_required_shape(ts);
    }
    void reinit() {
      tensor() = child(0).tensor();
      tensor().permute(reorder);
    }
    void exec_(size_type, dim_type) {}
  };

  /* diagonal tensor: take the "diagonal" of a tensor
     (ie diag(t(i,j,k), {i,k}) == t(i,j,i)) 

     /!\ the number of dimensions DO NOT change 
  */
  class ATN_diagonal_tensor : public ATN_tensor {
    dim_type i1, i2;
  public:
    ATN_diagonal_tensor(ATN_tensor& a, dim_type i1_, dim_type i2_) : 
      i1(i1_), i2(i2_) { add_child(a); }
  private:
    void check_shape_update(size_type , dim_type) {
      if ((shape_updated_ = child(0).is_shape_updated())) {
	if (i1 >= child(0).ranges().size() || i2 >= child(0).ranges().size() ||
	    i1 == i2 || child(0).ranges()[i1] != child(0).ranges()[i2])
	  ASM_THROW_TENSOR_ERROR("can't take the diagonal of a tensor of sizes " << child(0).ranges() << 
				 " at indexes " << int(i1) << " and " << int(i2));
	r_ = child(0).ranges();
      }
    }
    void update_childs_required_shape() {
      tensor_shape ts = req_shape.diag_shape(tensor_mask::Diagonal(i1,i2));
      child(0).merge_required_shape(ts);
    }
    void reinit() {
      tensor() = tensor_ref(child(0).tensor(), tensor_mask::Diagonal(i1,i2));
    }
    void exec_(size_type, dim_type) {}
  };

  /* called (if possible, i.e. if not an exact integration) for each integration point
     during mat_elem->compute() */    
  struct computed_tensor_integration_callback : public mat_elem_integration_callback {
    bgeot::tensor_reduction red;
    bool was_called;
    std::vector<TDIter> tensor_bases; /* each tref of 'red' has a reference into this vector */
    virtual void exec(bgeot::base_tensor &t, bool first, scalar_type c) {
      if (first) { 
        resize_t(t);
        std::fill(t.begin(), t.end(), 0.);
        was_called = true;
      }      
      assert(t.size());
      for (unsigned k=0; k != eltm.size(); ++k) { /* put in the 'if (first)' ? */
        tensor_bases[k] = const_cast<TDIter>(&(*eltm[k]->begin()));
      }
      red.do_reduction();
      int one = 1, n = red.out_data.size(); assert(n);
      daxpy_(&n, &c, const_cast<double*>(&(red.out_data[0])), &one, (double*)&(t[0]), &one);
    }
    void resize_t(bgeot::base_tensor &t) {
      bgeot::multi_index r; 
      if (red.reduced_range.size())
	r.assign(red.reduced_range.begin(), red.reduced_range.end());
      else { r.resize(1); r[0]=1; }
      t.adjust_sizes(r);      
    }
  };

  /*
    ATN_computed_tensor , the largest of all

    This object has become quite complex. It is the glue with the
    mat_elem_*.  It is able to perform an inline reduction (i.e. a
    reduction applied during the mat_elem->compute()) when it is
    allowed (i.e. no exact integration), or do the same reduction
    after the mat_elem->compute().
    The reduction may also involve other ATN_tensors.
  */

  struct mf_comp_vect;

  struct mf_comp {
    pnonlinear_elem_term nlt;
    const mesh_fem* pmf; /* always defined except when op_type == NORMAL */
    mf_comp_vect *owner;
				     
    ATN_tensor *data;
    std::vector<const mesh_fem*> auxmf; /* used only by nonlinear terms */
    typedef enum { BASE=1, GRAD=2, HESS=3, NORMAL=4, GRADGT=5, GRADGTINV=6, NONLIN=7, DATA=8 } op_type;
    op_type op; /* the numerical values indicates the number 
		   of dimensions in the tensor */
    bool vectorize; /* true if vectorization was required (adds
		       an addiational dimension to the tensor which
		       represents the component number
		    */
    std::string reduction;
    /*
      vectorization of non-vector FEM:

      phi1  0    0
      0   phi1   0
      0     0  phi1
      phi2  0    0
      0   phi2   0
      0     0  phi2
      ...
    */
    mf_comp(mf_comp_vect *ow, const mesh_fem* pmf_, op_type op_, bool vect) :
      nlt(0), pmf(pmf_), owner(ow), data(0), op(op_), vectorize(vect) { }
    mf_comp(mf_comp_vect *ow, const std::vector<const mesh_fem*> vmf, pnonlinear_elem_term nlt_) : 
      nlt(nlt_), pmf(vmf[0]), owner(ow), data(0), auxmf(vmf.begin()+1, vmf.end()), op(NONLIN), vectorize(false) { }
    mf_comp(mf_comp_vect *ow, ATN_tensor *t) :
      nlt(0), pmf(0), owner(ow), data(t), op(DATA), vectorize(0) {}
    void push_back_dimensions(size_type cv, tensor_ranges &rng, bool only_reduced=false) const;
    bool reduced(unsigned i) const
    { if (i >= reduction.size()) return false; else return reduction[i] != ' '; }
  };

  struct mf_comp_vect : public std::vector<mf_comp> {
    const mesh_im *main_im;
  public:
    mf_comp_vect() : std::vector<mf_comp>(), main_im(0) { }
    mf_comp_vect(const mf_comp_vect &other) : std::vector<mf_comp>(other), main_im(other.main_im) {
      for (size_type i=0; i < size(); ++i) (*this)[i].owner = this;
    }
    void set_im(const mesh_im &mim) {
      main_im = &mim;
    }
    const mesh_im& get_im() const { 
      return *main_im; 
    }
  private:
    mf_comp_vect& operator=(const mf_comp_vect &other);
  };

  void mf_comp::push_back_dimensions(size_type cv, tensor_ranges &rng,
				     bool only_reduced) const {
    switch (op) {
      case NONLIN:
	for (unsigned j=0; j < nlt->sizes().size(); ++j)
	  if (!only_reduced || !reduced(j)) 
	    rng.push_back(nlt->sizes()[j]);
	break;
      case DATA:
	for (unsigned i=0; i < data->ranges().size(); ++i) 
	  if (!only_reduced || !reduced(i)) 
	    rng.push_back(data->ranges()[i]);
	break;
      case NORMAL:
	assert(pmf==0);
	assert(&owner->get_im());
	assert(owner->get_im().linked_mesh().dim() != dim_type(-1));
	rng.push_back(owner->get_im().linked_mesh().dim());
	break;
      case GRADGT:
      case GRADGTINV:
	assert(pmf==0);
	assert(&owner->get_im());
	rng.push_back(owner->get_im().linked_mesh().dim()); // N
	rng.push_back(owner->get_im().linked_mesh().structure_of_convex(cv)->dim()); // P
	break;
      default:
	unsigned d = 0;
	if (!only_reduced || !reduced(d)) rng.push_back(pmf->nb_dof_of_element(cv));
	++d; 
	if (vectorize) {
	  if (!only_reduced || !reduced(d)) rng.push_back(pmf->get_qdim());
	  ++d;
	}
        
	if (op == GRAD || op == HESS) {
	  if (!only_reduced || !reduced(d)) rng.push_back(pmf->linked_mesh().dim());
	  ++d;
	}
	if (op == HESS) {
	  if (!only_reduced || !reduced(d)) rng.push_back(pmf->linked_mesh().dim());
	  ++d;
	}
	break;
    }
  }


  class ATN_computed_tensor : public ATN_tensor {
    mf_comp_vect mfcomp;
    pmat_elem_computation pmec;
    pmat_elem_type pme;
    pintegration_method pim;
    bgeot::pgeometric_trans pgt;
    base_tensor t;
    std::vector<scalar_type> data;     
    TDIter data_base;
    stride_type tsize;
    dal::bit_vector req_bv;  /* bit_vector of values the mat_elem has to compute
                               (useful when only a subset is required from the 
                               possibly very large elementary tensor) */
    bool has_inline_reduction; /* true if used with reductions inside the comp, for example:
                                  "comp(Grad(#1)(:,i).Grad(#2)(:,i))" */
    computed_tensor_integration_callback icb; /* callback for inline reductions */

    /* if inline reduction are to be done, but were not possible (i.e. if exact
       integration was used) then a fallback is used: apply the reduction
       afterward, on the large expanded tensor */
    bgeot::tensor_reduction fallback_red; 
    bool fallback_red_uptodate;
    TDIter fallback_base;

    size_type cv_shape_update;
    //mat_elem_inline_reduction inline_red;
  public:
    ATN_computed_tensor(const mf_comp_vect &mfcomp_) : 
      mfcomp(mfcomp_), pmec(0), pme(0), pim(0), pgt(0), data_base(0) { 
      has_inline_reduction = false;
      bool in_data = false;
      for (size_type i=0; i < mfcomp.size(); ++i) {
        if (mfcomp[i].reduction.size() || mfcomp[i].op == mf_comp::DATA) {
          has_inline_reduction = true; 
          if (mfcomp[i].op == mf_comp::DATA) { add_child(*mfcomp[i].data); in_data = true; }
        }
        if (mfcomp[i].op != mf_comp::DATA && in_data) {
          /* constraint of fallback 'do_post_reduction' */
          ASM_THROW_ERROR("data tensors inside comp() cannot be intermixed with Grad() and Base() etc., they must appear LAST");
        }
      }
    }

  private:    
    /* mostly for non-linear terms, such as a 3x3x3x3 tensor which may have
       many symmetries or many null elements..  in that case, it is preferable
       for getfem_mat_elem to handle only a sufficient subset of the tensor,
       and build back the full tensor via adequate strides and masks */
    /*
    stride_type add_packed_dims(const bgeot::tensor<unsigned> &packed_idx, stride_type s) {
      size_type d = r_.size();
      const bgeot::multi_index sz = packed_idx.sizes();
      bgeot::index_set dims(sz.size());
      bgeot::tensor_ranges rngs(sz.size());
      for (dim_type i=0; i < dims.size(); ++i) {
        dims[i] = d+i; // new dimension numbers 
        rngs[i] = sz[i];
      }

      r_.insert(r_.end(), rngs.begin(), rngs.end());
      tensor_strides v;
      tensor_mask m(rngs,dims);
      for (index_type i=0; i < packed_idx.size(); ++i)
        if (packed_idx[i] != unsigned(-1)) m.set_mask_val(i,true);
      v.resize(m.card());
      stride_type max_stride = 0;
      for (index_type i=0, p=0; i < packed_idx.size(); ++i) {
        if (m(i)) {
          v[p++] = s*packed_idx[i];
          max_stride = std::max<stride_type>(max_stride, packed_idx[i]*s);
        }
      }
      cout << "add_packed_dims : sz=" << sz << ", r_=" << r_ << ", s=" << s << ", card=" << m.card() << ", max_stride=" << max_stride << "\n";
      assert(tensor().masks().size() == tensor().strides().size());
      tensor().set_ndim_noclean(tensor().ndim()+sz.size());
      tensor().push_mask(m);
      tensor().strides().push_back(v);
      return max_stride;
    }
    */

    /* append a dimension (full) to tref */
    stride_type add_dim(const tensor_ranges& rng, dim_type d, stride_type s, tensor_ref &tref) {
      assert(d < rng.size());
      tensor_strides v;
      index_type r = rng[d];
      tensor_mask m; m.set_full(d, r);	    
      v.resize(r);
      for (index_type i=0; i < r; ++i) v[i] = s*i;
      assert(tref.masks().size() == tref.strides().size());
      tref.set_ndim_noclean(tref.ndim()+1);
      tref.push_mask(m);
      tref.strides().push_back(v);
      return s*r;
    }
    
    /* append a vectorized dimension to tref */
    stride_type add_vdim(const tensor_ranges& rng, dim_type d, stride_type s, tensor_ref &tref) {
      assert(d < rng.size()-1);
      index_type r = rng[d], q=rng[d+1];
      tensor_strides v;
      tensor_ranges trng(2); trng[0] = q; trng[1] = r;
      index_set ti(2); ti[0] = d+1; ti[1] = d;
      tensor_mask m(trng,ti);
      assert(r%q == 0);
      v.resize(r);
      tensor_ranges cnt(2);
      for (cnt[1]=0; cnt[1] < r; cnt[1]++) {
	cnt[0] = cnt[1] % q;
	m.set_mask_val(m.lpos(cnt), true);
	v[cnt[1]] = s*(cnt[1]/q);
      }
      assert(tref.masks().size() == tref.strides().size());
      tref.set_ndim_noclean(tref.ndim()+2);
      tref.push_mask(m);
      tref.strides().push_back(v);
      return s*(r/q);
    }

    /* called when the FEM has changed */
    void update_pmat_elem(size_type cv) {
      pme = 0;
      for (size_type i=0; i < mfcomp.size(); ++i) {
        if (mfcomp[i].op == mf_comp::DATA) continue;
	pfem fem = (mfcomp[i].pmf ? mfcomp[i].pmf->fem_of_element(cv) : NULL);
	pmat_elem_type pme2 = 0;
	switch (mfcomp[i].op) {
	  case mf_comp::BASE: pme2 = mat_elem_base(fem); break;
	  case mf_comp::GRAD: pme2 = mat_elem_grad(fem); break;
	  case mf_comp::HESS: pme2 = mat_elem_hessian(fem); break;
	  case mf_comp::NORMAL: pme2 = mat_elem_unit_normal(); break;
	  case mf_comp::GRADGT:
	  case mf_comp::GRADGTINV: 
	    pme2 = mat_elem_grad_geotrans(mfcomp[i].op == mf_comp::GRADGTINV); break;
	  case mf_comp::NONLIN: {
	      std::vector<pfem> ftab(1+mfcomp[i].auxmf.size()); 
	      ftab[0] = fem;
	      for (unsigned k=0; k < mfcomp[i].auxmf.size(); ++k) 
		ftab[k+1] = mfcomp[i].auxmf[k]->fem_of_element(cv);
	      pme2 = mat_elem_nonlinear(mfcomp[i].nlt, ftab);
	    } break;
          case mf_comp::DATA: /*ignore*/;
        } 
        if (pme == 0) pme = pme2;
        else pme = mat_elem_product(pme, pme2);
      }
      
      if (pme == 0) pme = mat_elem_empty(); //ASM_THROW_ERROR("no Base() or Grad() or etc!");
    }

   

    size_type 
    push_back_mfcomp_dimensions(size_type cv, const mf_comp& mc, 
				unsigned &d, const bgeot::tensor_ranges &rng,
				bgeot::tensor_ref &tref, size_type tsz=1) {
      if (mc.op == mf_comp::NONLIN) {
	for (size_type j=0; j < mc.nlt->sizes().size(); ++j)
	  tsz = add_dim(rng, d++, tsz, tref);
      } else if (mc.op == mf_comp::DATA) {
        assert(tsz == 1);
	tref = mc.data->tensor();
	tsz *= tref.card();
	d += tref.ndim();
      } else if (mc.op == mf_comp::NORMAL) {
	tsz = add_dim(rng, d++, tsz, tref);
      } else if (mc.op == mf_comp::GRADGT || 
		 mc.op == mf_comp::GRADGTINV) {
	tsz = add_dim(rng, d++, tsz, tref);
	tsz = add_dim(rng, d++, tsz, tref);
      } else {
        size_type target_dim = mc.pmf->fem_of_element(cv)->target_dim();
        size_type qdim = mc.pmf->get_qdim();
	if (mc.vectorize) {
	  if (target_dim == qdim) {
	    tsz = add_dim(rng, d++, tsz, tref);
	    tsz = add_dim(rng, d++, tsz, tref);
	  } else {
	    tsz = add_vdim(rng, d, tsz, tref);
	    d += 2;
	  }
	} else tsz = add_dim(rng, d++, tsz, tref);
	if (mc.op == mf_comp::GRAD || mc.op == mf_comp::HESS) {
	  tsz = add_dim(rng, d++, tsz, tref);
	}
	if (mc.op == mf_comp::HESS) {
	  tsz = add_dim(rng, d++, tsz, tref);
	}
      }
      return tsz;
    }

    void update_shape_with_inline_reduction(size_type cv) {
      fallback_red_uptodate = false;
      icb.tensor_bases.resize(mfcomp.size()); /* todo : resize(nb_mfcomp_not_data) */
      for (size_type i=0; i < mfcomp.size(); ++i) {
	tensor_ref tref;
	tensor_ranges rng;
	unsigned d = 0;
	mfcomp[i].push_back_dimensions(cv,rng);
	push_back_mfcomp_dimensions(cv,mfcomp[i], d, rng, tref); 
	assert(tref.ndim() == rng.size() && d == rng.size());
	if (mfcomp[i].reduction.size() == 0) 
	  mfcomp[i].reduction.insert(size_type(0), tref.ndim(), ' ');
	if (mfcomp[i].op != mf_comp::DATA) /* should already have the correct base */
          tref.set_base(icb.tensor_bases[i]);
        tref.update_idx2mask();
	if (mfcomp[i].reduction.size() != tref.ndim()) {
	  ASM_THROW_TENSOR_ERROR("wrong number of indexes for the " << int(i+1) 
				 << "th argument of the reduction " << name() 
				 << " (expected " << int(tref.ndim()) 
                                 << " indexes, got " << mfcomp[i].reduction.size());
	}
	icb.red.insert(tref, mfcomp[i].reduction);
      }
      //cout << "  build inline reduction : " << name() << "\n";
      icb.red.prepare();
      //cout << "\n done inline reduction : " << name() << "\n";
      icb.red.result(tensor());
      r_.resize(tensor().ndim()); 
      for (unsigned i=0; i < tensor().ndim(); ++i) r_[i] = tensor().dim(i);
      tsize = tensor().card();
      //cerr << "update_shape_with_inline_reduction: tensor=" << tensor() << "\nr_=" << r_ << ", tsize=" << tsize << "\n";
    }
    
    void update_shape_with_expanded_tensor(size_type cv) {
      icb.red.clear();
      unsigned d = 0;
      for (size_type i=0; i < mfcomp.size(); ++i) {
	tsize = push_back_mfcomp_dimensions(cv, mfcomp[i], d, r_, tensor(), tsize);
      }
      assert(d == r_.size());
      tensor().update_idx2mask();
    }

    void check_shape_update(size_type cv, dim_type) {
      const mesh_im& mi = mfcomp.get_im();
      pintegration_method pim2;
      bgeot::pgeometric_trans pgt2;
      bool fem_changed = false;
      pgt2 = mi.linked_mesh().trans_of_convex(cv);
      pim2 = mi.int_method_of_element(cv);
      // cerr << "computed tensor cv=" << cv << " f=" << int(face) << "\n";
      /* shape is considered for update only if the FEM changes,
	 changes of pgt or integration method does not affect shape
	 (only the mat_elem) */
      cv_shape_update = cv;
      shape_updated_ = (pme == 0); //false;
      for (size_type i=0; i < nchilds(); ++i)
        shape_updated_ = shape_updated_ || child(i).is_shape_updated();
      for (size_type i=0; shape_updated_ == false && i < mfcomp.size(); ++i) {
	if (mfcomp[i].pmf == 0) continue;
	if  (current_cv == size_type(-1)) {
	  shape_updated_ = true; fem_changed = true;
	} else {
	  fem_changed = fem_changed || 
	    (mfcomp[i].pmf->fem_of_element(current_cv) != 
	     mfcomp[i].pmf->fem_of_element(cv));
	  /* for FEM with non-constant nb_dof.. */ 
	  shape_updated_ = shape_updated_ ||
	    (mfcomp[i].pmf->nb_dof_of_element(current_cv) != 
	     mfcomp[i].pmf->nb_dof_of_element(cv)); 
	}
      }
      if (shape_updated_) {
        r_.resize(0);
        /* get the new ranges */
        for (unsigned i=0; i < mfcomp.size(); ++i) 
          mfcomp[i].push_back_dimensions(cv, r_, true);
      }
      if (fem_changed || shape_updated_) {
        /* build the new mat_elem structure */      
        update_pmat_elem(cv);
      }
      if (shape_updated_ || fem_changed || pgt != pgt2 || pim != pim2) {
	pgt = pgt2; pim = pim2;
	pmec = mat_elem(pme, pim, pgt, has_inline_reduction);
      }
    }

    void reinit() {
      if (!shape_updated_) return;
      tensor().clear();
      tsize = 1;
      if (has_inline_reduction) {
        update_shape_with_inline_reduction(cv_shape_update);
      } else {
        update_shape_with_expanded_tensor(cv_shape_update);
      }
      data_base = 0;
      tensor().set_base(data_base);
    }
    void update_childs_required_shape() {
    }

    /* fallback when inline reduction is not possible:
       do the reduction after evaluation of the mat_elem */
    void do_post_reduction(size_type cv) {
      if (!fallback_red_uptodate) {
        fallback_red_uptodate = true;
        std::string s;
        size_type sz = 1;
        tensor_ref tref; tensor_ranges r;
        unsigned cnt, d=0;
        /* insert the tensorial product of Grad() etc */
        for (cnt=0; cnt < mfcomp.size() && mfcomp[cnt].op != mf_comp::DATA; ++cnt) {
	  mfcomp[cnt].push_back_dimensions(cv,r);
          sz = push_back_mfcomp_dimensions(cv, mfcomp[cnt], d, r, tref, sz);
          s += mfcomp[cnt].reduction;
        }
        fallback_red.clear();
        tref.set_base(fallback_base);
        fallback_red.insert(tref, s);
        /* insert the optional data tensors */
        for ( ; cnt < mfcomp.size(); ++cnt) {
          assert(mfcomp[cnt].op == mf_comp::DATA);
          fallback_red.insert(mfcomp[cnt].data->tensor(), mfcomp[cnt].reduction);
        }
        fallback_red.prepare();
        fallback_red.result(tensor()); /* this SHOULD NOT, IN ANY CASE change the shape or strides of tensor() .. */
        assert(icb.red.reduced_range == fallback_red.reduced_range);
      }
      icb.resize_t(t);
      fallback_base = &(*t.begin());
      fallback_red.do_reduction();
    }
    
    void exec_(size_type cv, dim_type face) {
      const mesh_im& mim = mfcomp.get_im();
      for (unsigned i=0; i < mfcomp.size(); ++i) {
        if (mfcomp[i].op == mf_comp::DATA) {
          size_type fullsz = 1;
          for (unsigned j=0; j < mfcomp[i].data->ranges().size(); ++j)
	    fullsz *= mfcomp[i].data->ranges()[j];
          if (fullsz != size_type(mfcomp[i].data->tensor().card())) 
            ASM_THROW_TENSOR_ERROR("aaarg inline reduction will explode with non-full tensors. "
                                   "Complain to the author, I was too lazy to do that properly");
        }
      }
      icb.was_called = false;
      if (face == dim_type(-1)) {
	  pmec->gen_compute(t, mim.linked_mesh().points_of_convex(cv), cv, 
			    has_inline_reduction ? &icb : 0);
      } else {
	pmec->gen_compute_on_face(t, mim.linked_mesh().points_of_convex(cv), face, cv, 
				  has_inline_reduction ? &icb : 0);
      }
      

      if (has_inline_reduction && icb.was_called == false) {
        do_post_reduction(cv);
        data_base = &fallback_red.out_data[0];
      } else data_base = &(*t.begin());
      if (t.size() != size_type(tsize)) DAL_INTERNAL_ERROR(t.size() << " != " << tsize);
    }
  };


  /* extract data for each dof of the convex */
  class ATN_tensor_from_dofs_data : public ATN_tensor_w_data {
    const base_asm_data *basm; //scalar_type* global_array;
    vdim_specif_list vdim;
    multi_tensor_iterator mti;
  public:
    ATN_tensor_from_dofs_data(const base_asm_data *basm_, 
			      const vdim_specif_list& d) :
      basm(basm_), vdim(d) {
    }
    void check_shape_update(size_type cv, dim_type) {
      shape_updated_ = false;
      r_.resize(vdim.size());
      for (dim_type i=0; i < vdim.size(); ++i) {
	if (vdim[i].is_mf_ref()) {
	  if (vdim[i].pmf->nb_dof_of_element(cv) != ranges()[i]) {
	    r_[i] = vdim[i].pmf->nb_dof_of_element(cv);
	    shape_updated_ = true;
	  }
	} else if (vdim[i].dim != ranges()[i]) {
	  r_[i] = vdim[i].dim;
	  shape_updated_ = true;
	}
      }
    }
    virtual void init_required_shape() { req_shape = tensor_shape(ranges()); }
    void reinit() {
      ATN_tensor_w_data::reinit();
      mti.assign(tensor(), true);
    }

  private:
    void exec_(size_type cv, dim_type ) {
      tensor_ranges r;
      std::vector< tensor_strides > str;
      vdim.build_strides_for_cv(cv, r, str);
      assert(r == ranges());
      mti.rewind();
      basm->copy_with_mti(str, mti);
    }
  };
  
  /* enforce symmetry of a 2D tensor 
     (requiring only the upper-triangle of its child and
     duplicating it) */
  class ATN_symmetrized_tensor : public ATN_tensor_w_data {
    multi_tensor_iterator mti;
  public:
    ATN_symmetrized_tensor(ATN_tensor& a) { add_child(a); }
    void check_shape_update(size_type , dim_type) {
      if ((shape_updated_ = child(0).is_shape_updated())) {
	if (child(0).ranges().size() != 2 || 
	    child(0).ranges()[0] != child(0).ranges()[1])
	  ASM_THROW_TENSOR_ERROR("can't symmetrize a non-square tensor of sizes " << child(0).ranges());
	r_ = child(0).ranges();
      }
    }
    void update_childs_required_shape() {
      tensor_shape ts = req_shape;
      tensor_shape ts2 = req_shape; 
      index_set perm(2); perm[0] = 1; perm[1] = 0; ts2.permute(perm);
      ts.merge(ts2, false);
      tensor_mask dm; dm.set_triangular(ranges()[0],0,1);
      tensor_shape tsdm(2); tsdm.push_mask(dm);
      ts.merge(tsdm, true);
      child(0).merge_required_shape(ts);
    }

    void reinit() {
      req_shape.set_full(ranges()); // c'est plus simple comme ça
      ATN_tensor_w_data::reinit0();
      mti.assign(child(0).tensor(),true);
    }
  private:
    void exec_(size_type, dim_type) {
      std::fill(data.begin(), data.end(), 0.);
      mti.rewind();
      index_type n = ranges()[0];
      do {
	index_type i=mti.index(0), j=mti.index(1);
	data[i*n+j]=data[j*n+i]=mti.p(0);
      } while (mti.qnext1());
    }
  };


  template<class UnaryOp> class ATN_unary_op_tensor : public ATN_tensor_w_data {
    multi_tensor_iterator mti;
  public:
    ATN_unary_op_tensor(ATN_tensor& a) { add_child(a); }
    void check_shape_update(size_type , dim_type) {
      if ((shape_updated_ = (ranges() != child(0).ranges())))      
	r_ = child(0).ranges();
    }
    void reinit() {
      ATN_tensor_w_data::reinit0();
      mti.assign(tensor(), child(0).tensor(),false);
    }
  private:
    void update_cv_(size_type, dim_type) {
      mti.rewind();
      do {
	mti.p(0) = UnaryOp()(mti.p(1));
      } while (mti.qnext2());
    }
  };

  /* sum AND scalar scaling */
  class ATN_tensors_sum_scaled : public ATN_tensor_w_data {
    std::vector<multi_tensor_iterator> mti;
    std::vector<scalar_type> scales; /* utile pour des somme "scalées" du genre 0.5*t1 + 0.5*t2 */
  public:
    ATN_tensors_sum_scaled(ATN_tensor& t1, scalar_type s1) {
      add_child(t1);
      scales.resize(1); scales[0]=s1;
    }
    void push_scaled_tensor(ATN_tensor& t, scalar_type s) {
      add_child(t); scales.push_back(s);
    }
    void check_shape_update(size_type , dim_type) {
      if ((shape_updated_ = child(0).is_shape_updated()))
	r_ = child(0).ranges();
      for (size_type i=1; i < nchilds(); ++i)
	if (ranges() != child(i).ranges()) 
	  ASM_THROW_TENSOR_ERROR("can't add two tensors of sizes " << 
				 ranges() << " and " << child(i).ranges());
    }
    void reinit() {
      ATN_tensor_w_data::reinit0();
      mti.resize(nchilds());
      for (size_type i=0; i < nchilds(); ++i)
	mti[i].assign(tensor(), child(i).tensor(),false);
    }
    void apply_scale(scalar_type s) { 
      for (size_type i=0; i < scales.size(); ++i) scales[i] *= s; 
    }
    ATN_tensors_sum_scaled* is_tensors_sum_scaled() { return this; }
  private:
    void exec_(size_type, dim_type) {
      //if (cv == 0) {
      //cerr << "ATN_tensors_sum["<< name() << "] req_shape=" << req_shape << endl;
      //}
      std::fill(data.begin(), data.end(), 0.);
      mti[0].rewind();
      do {
	mti[0].p(0) = mti[0].p(1)*scales[0];
      } while (mti[0].qnext2());
      for (size_type i=1; i < nchilds(); ++i) {
	mti[i].rewind();
	do {
	  mti[i].p(0) = mti[i].p(0)+mti[i].p(1)*scales[i];
	} while (mti[i].qnext2());
      }
    }
  };

  class ATN_tensor_scalar_add : public ATN_tensor_w_data {
    scalar_type v;
    multi_tensor_iterator mti;
    int sgn; /* v+t or v-t ? */
  public:
    ATN_tensor_scalar_add(ATN_tensor& a, scalar_type v_, int sgn_) : 
      v(v_), sgn(sgn_) { add_child(a); }
    void check_shape_update(size_type , dim_type) {
      if ((shape_updated_ = child(0).is_shape_updated()))
	r_ = child(0).ranges();
    }
    void reinit() {
      ATN_tensor_w_data::reinit();
      mti.assign(tensor(), child(0).tensor(),false);
    }
  private:
    void exec_(size_type, dim_type) {
      std::fill(data.begin(), data.end(), v);
      mti.rewind();
      do {
	if (sgn > 0)
	  mti.p(0) += mti.p(1);
	else mti.p(0) -= mti.p(1);
      } while (mti.qnext2());
    }
  };

  class ATN_print_tensor : public ATN {
    std::string name;
  public:
    ATN_print_tensor(ATN_tensor& a, std::string n_) 
      : name(n_) { add_child(a); }
    void reinit() {}
  private:
    void exec_(size_type cv, dim_type face) {
      multi_tensor_iterator mti(child(0).tensor(), true);
      cout << "------- > evaluation of " << name << ", at" << endl;
      cout << "convex " << cv; 
      if (face != dim_type(-1)) cout << ", face " << int(face);
      cout << endl;
      cout << "  size   = " << child(0).ranges() << endl;
      mti.rewind();
      do {
	cout << " @[";
	for (size_type i=0; i < child(0).ranges().size(); ++i)
	  cout <<(i>0?",":"") << mti.index(i);
	cout << "] = " << mti.p(0) << endl;
      } while (mti.qnext1());
    }
  };


  /* 
     -------------------
     analysis of the supplied string
     -----------------
  */

  std::string asm_tokenizer::syntax_err_print() {      
    std::string s;
    if (tok_pos - err_msg_mark > 80) err_msg_mark = tok_pos - 40;
    if (str.length() - err_msg_mark < 80) s = tok_substr(err_msg_mark, str.length());
    else { s = tok_substr(err_msg_mark,err_msg_mark+70); s.append(" ... (truncated)"); }
    s += "\n" + std::string(std::max(int(tok_pos - err_msg_mark),0), '-') + "^^";
    return s;
  }
  
  void asm_tokenizer::get_tok() {
    std::string sub = str.substr(tok_pos, 10);
    //cerr << "tok_pos=" << tok_pos << "... '" << sub << "'..." << endl;
    curr_tok_ival = -1;
    while (tok_pos < str.length() && isspace(str[tok_pos])) ++tok_pos;
    if (tok_pos == str.length()) { 
      curr_tok_type = END; tok_len = 0;
    } else if (strchr("{}(),;:=-.*/+", str[tok_pos])) { 
      curr_tok_type = tok_type_enum(str[tok_pos]); tok_len = 1;
    } else if (str[tok_pos] == '$' || str[tok_pos] == '#' || str[tok_pos] == '%') { 
      curr_tok_type = str[tok_pos] == '$' ? ARGNUM_SELECTOR : 
	(str[tok_pos] == '#' ? MFREF : IMREF);
      tok_len = 1; 
      curr_tok_ival = 0;
      while (isdigit(str[tok_pos+tok_len])) { 
	curr_tok_ival*=10; 
	curr_tok_ival += str[tok_pos+tok_len] - '0';
	++tok_len;
      }
      curr_tok_ival--;
    } else if (isalpha(str[tok_pos])) {
      curr_tok_type = IDENT;
      tok_len = 0;
      while (isalnum(str[tok_pos+tok_len]) || str[tok_pos+tok_len] == '_') ++tok_len;
    } else if (isdigit(str[tok_pos])) {
      curr_tok_type = NUMBER;
      char *p;
      curr_tok_dval = strtod(&str[0]+tok_pos, &p);
      tok_len = p - &str[0] - tok_pos;
    }
    curr_tok = str.substr(tok_pos, tok_len);
    /*
      cerr << "TOK[ ";
      switch (curr_tok_type) {
      case OPEN_PAR: case CLOSE_PAR: case EQUAL:
      case COMMA: case SEMICOLON: case COLON:
      case PLUS: case MINUS: case DIVIDE:
      case OPEN_BRACE: case CLOSE_BRACE:
      case PRODUCT: case MULTIPLY: cerr << "'" << char(curr_tok_type) << "'"; break;
      case IDENT: cerr << "IDENT: '" << tok() << "'"; break;
      case NUMBER: cerr << "NUMBER: " << curr_tok_dval; break;
      case MFREF: cerr << "MFREF: " << tok_mfref_num(); break;
      case  ARGNUM_SELECTOR: cerr << "ARGNUM: " << tok_argnum(); break;
      case END: cerr << "END"; break;
      }
      cerr << "]" << endl;*/
  }

  const mesh_fem& generic_assembly::do_mf_arg(std::vector<const mesh_fem*> * multimf) {
    if (!multimf) advance(); // special hack for NonLin$i(#a,#b,..)
    accept(OPEN_PAR,"expecting '('");
    if (tok_type() != MFREF) ASM_THROW_PARSE_ERROR("expecting mesh_fem reference");
    if (tok_mfref_num() >= mftab.size()) 
      ASM_THROW_PARSE_ERROR("reference to a non-existant mesh_fem #" << tok_mfref_num()+1);
    const mesh_fem& mf_ = *mftab[tok_mfref_num()]; advance();
    if (multimf) {
      multimf->resize(1); (*multimf)[0] = &mf_;
      while (advance_if(COMMA)) {
	if (tok_type() != MFREF) ASM_THROW_PARSE_ERROR("expecting mesh_fem reference");
	if (tok_mfref_num() >= mftab.size()) 
	  ASM_THROW_PARSE_ERROR("reference to a non-existant mesh_fem #" << tok_mfref_num()+1);
	multimf->push_back(mftab[tok_mfref_num()]); advance();
      }
    }
    accept(CLOSE_PAR, "expecting ')'");
    return mf_;
  }

  /* "inline" reduction operations inside comp(..) */
  std::string generic_assembly::do_comp_red_ops() {
    std::string s;
    if (advance_if(OPEN_PAR)) {
      size_type j = 0;
      do {
	if (tok_type() == COLON) {
	  s.push_back(' '); advance(); j++;
	} else if (tok_type() == IDENT) {
	  if (tok().length()==1 && isalpha(tok()[0]) || islower(tok()[0])) {
	    s.push_back(tok()[0]); advance(); j++;
	  } else ASM_THROW_PARSE_ERROR("invalid reduction index '" << tok() << 
				       "', only lower case characters allowed");
	}
      } while (advance_if(COMMA));
      accept(CLOSE_PAR, "expecting ')'");
    }
    return s;
  }

  ATN_tensor* generic_assembly::do_comp() {
    accept(OPEN_PAR, "expecting '('");
    mf_comp_vect what;
    bool in_data = false;
    /* the first optinal argument is the "main" mesh_im, i.e. the one
       whose integration methods are used, (and whose linked_mesh is
       used for mf_comp::NORMAL, mf_comp::GRADGT etc computations). If
       not given, then the first mesh_im pushed is used (then expect
       problems when assembling simultaneously on two different
       meshes).
    */
    if (tok_type() == IMREF) {
      if (tok_imref_num() >= imtab.size()) 
	ASM_THROW_PARSE_ERROR("reference to a non-existant mesh_im #" << tok_imref_num()+1);
      what.set_im(*imtab[tok_imref_num()]); advance();
      accept(COMMA, "expecting ','");
    } else {
      what.set_im(*imtab[0]);
    }
    do {
      if (tok_type() == CLOSE_PAR) break;
      if (tok_type() != IDENT) ASM_THROW_PARSE_ERROR("expecting Base or Grad or Hess, Normal, etc..");
      std::string f = tok(); 
      const mesh_fem *pmf = 0;
      if (f.compare("Base")==0 || f.compare("vBase")==0) {
	pmf = &do_mf_arg(); what.push_back(mf_comp(&what, pmf, mf_comp::BASE, f[0] == 'v'));
      } else if (f.compare("Grad")==0 || f.compare("vGrad")==0) {
	pmf = &do_mf_arg(); what.push_back(mf_comp(&what, pmf, mf_comp::GRAD, f[0] == 'v'));
      } else if (f.compare("Hess")==0 || f.compare("vHess")==0) {
	pmf = &do_mf_arg(); what.push_back(mf_comp(&what, pmf, mf_comp::HESS, f[0] == 'v'));
      } else if (f.compare("NonLin")==0) {	
	size_type num = 0; /* default value */
	advance();
	if (tok_type() == ARGNUM_SELECTOR) { num = tok_argnum(); advance(); }
	if (num >= innonlin.size()) ASM_THROW_PARSE_ERROR("NonLin$" << num << " does not exist");
	std::vector<const mesh_fem*> allmf;
	pmf = &do_mf_arg(&allmf); what.push_back(mf_comp(&what, allmf, innonlin[num]));
      } else if (f.compare("Normal") == 0) {
	advance();
	accept(OPEN_PAR,"expecting '('"); accept(CLOSE_PAR,"expecting ')'");
	pmf = 0; what.push_back(mf_comp(&what, pmf, mf_comp::NORMAL, false));
      } else if (f.compare("GradGT") == 0 ||
		 f.compare("GradGTInv") == 0) {
	advance();
	accept(OPEN_PAR,"expecting '('"); accept(CLOSE_PAR,"expecting ')'");
	pmf = 0; 
	what.push_back(mf_comp(&what, pmf, 
			       f.compare("GradGT") == 0 ?
			       mf_comp::GRADGT : 
			       mf_comp::GRADGTINV, false));
      } else {
        if (vars.find(f) != vars.end()) {
          what.push_back(mf_comp(&what, vars[f]));
          in_data = true;
          advance();
        } else {
          ASM_THROW_PARSE_ERROR("expecting Base, Grad, vBase, NonLin ...");
        }
      }
	
      if (!in_data && f[0] != 'v' && pmf && pmf->get_qdim() != 1 && f.compare("NonLin")) {
	ASM_THROW_PARSE_ERROR("Attempt to use a vector mesh_fem as a scalar mesh_fem");
      }
      what.back().reduction = do_comp_red_ops();
    } while (advance_if(PRODUCT));
    accept(CLOSE_PAR, "expecting ')'");

    return record(new ATN_computed_tensor(what));
  }

  void generic_assembly::do_dim_spec(vdim_specif_list& lst) {
    lst.resize(0);
    accept(OPEN_PAR, "expecting '('");    
    while (true) {
      if (tok_type() == IDENT) {
	if (tok().compare("mdim")==0) lst.push_back(vdim_specif(do_mf_arg().linked_mesh().dim()));
	else if (tok().compare("qdim")==0) lst.push_back(vdim_specif(do_mf_arg().get_qdim()));
	else ASM_THROW_PARSE_ERROR("expecting mdim(#mf) or qdim(#mf) or a number or a mesh_fem #id");
      } else if (tok_type() == NUMBER) {
	lst.push_back(vdim_specif(tok_number_ival()+1));
	advance();
      } else if (tok_type() == MFREF) {
	if (tok_mfref_num() >= mftab.size()) 
	  ASM_THROW_PARSE_ERROR("reference to a non-existant mesh_fem #" << tok_mfref_num()+1);
	lst.push_back(vdim_specif(mftab[tok_mfref_num()]));
	advance();
      } else if (tok_type() != CLOSE_PAR) ASM_THROW_PARSE_ERROR("expecting mdim(#mf) or qdim(#mf) or a number or a mesh_fem #id");
      /*      if (mfcnt && !lst.back().is_mf_ref())
	      ASM_THROW_PARSE_ERROR("#mf argument must be given after numeric dimensions");*/
      if (advance_if(CLOSE_PAR)) break;
      accept(COMMA,"expecting ',' or ')'");
    }
  }


  ATN_tensor* generic_assembly::do_data() {
    //    ATN_tensor *t;
    size_type datanum = 0; /* par défaut */
    if (tok_type() != OPEN_PAR) { /* on peut oublier le numero de dataset */
      if (tok_type() != ARGNUM_SELECTOR) ASM_THROW_PARSE_ERROR("expecting dataset number");
      datanum = tok_argnum();
      advance(); 
    }
    if (datanum >= indata.size()) ASM_THROW_PARSE_ERROR("wong dataset number: " << datanum);
    
    vdim_specif_list sz;
    do_dim_spec(sz);
   
    if (sz.nbelt() != indata[datanum]->vect_size()) 
      ASM_THROW_PARSE_ERROR("invalid size for data argument " << datanum+1 << 
			    " real size is " << indata[datanum]->vect_size() << 
			    " expected size is " << sz.nbelt());
    return record(new ATN_tensor_from_dofs_data(indata[datanum], sz));
  }

  std::pair<ATN_tensor*, std::string> generic_assembly::do_red_ops(ATN_tensor* t) {
    std::string s;

    if (advance_if(OPEN_PAR)) {
      size_type j = 0;
      do {
	if (tok_type() == COLON) {
	  s.push_back(' '); advance(); j++;
	} else if (tok_type() == NUMBER) {
	  t = record(new ATN_sliced_tensor(*t, j, tok_number_ival())); advance();
	} else if (tok_type() == IDENT) {
	  if (tok().length()==1 && isalpha(tok()[0]) || islower(tok()[0])) {
	    s.push_back(tok()[0]); advance(); j++;
	  } else ASM_THROW_PARSE_ERROR("invalid reduction index '" << tok() << 
				       "', only lower case characters allowed");
	}
      } while (advance_if(COMMA));
      accept(CLOSE_PAR, "expecting ')'");
    }
    return std::pair<ATN_tensor*,std::string>(t,s);
  }

  /*
    ( expr ) 
    variable
    comp(..)
    data(data)
  */
  tnode generic_assembly::do_tens() {
    tnode t;
    push_mark();
    if (advance_if(OPEN_PAR)) {
      t = do_expr();
      accept(CLOSE_PAR, "expecting ')'");
    } else if (tok_type() == NUMBER) {
      t.assign(tok_number_dval()); advance();
    } else if (tok_type() == IDENT) {
      if (vars.find(tok()) != vars.end()) {
	t.assign(vars[tok()]); advance();
      } else if (tok().compare("comp")==0) {
	advance(); t.assign(do_comp());
      } else if (tok().compare("data")==0) {
	advance(); t.assign(do_data());
      } else if (tok().compare("sym")==0) {
	advance(); 
	tnode t2 = do_expr();
	if (t2.type() != tnode::TNTENSOR) ASM_THROW_PARSE_ERROR("can't symmetrise a scalar!");
	t.assign(record(new ATN_symmetrized_tensor(*t2.tensor())));
      } else ASM_THROW_PARSE_ERROR("unknown identifier: " << tok());
    } else ASM_THROW_PARSE_ERROR("unexpected token: " << tok());
    pop_mark();
    return t;
  }

  /*
    handle tensorial product/reduction

    a(:,i).b(j,i).(c)(1,:,i)
  */
  tnode generic_assembly::do_prod() {
    reduced_tensor_arg_type ttab;
    
    do {
      tnode t = do_tens();
      if (t.type() == tnode::TNCONST) {
	if (ttab.size() == 0) return t;
	else ASM_THROW_PARSE_ERROR("can't mix tensor and scalar into a reduction expression");
      }
      ttab.push_back(do_red_ops(t.tensor()));
    } while (advance_if(PRODUCT));
    if (ttab.size() == 1 && ttab[0].second.length() == 0) {
      return tnode(ttab[0].first);
    } else {
      return tnode(record(new ATN_reduced_tensor(ttab)));
    }
  }

  /* calls do_prod() once,
     and handle successive reordering/diagonals transformations */
  tnode generic_assembly::do_prod_trans() {
    tnode t = do_prod();
    while (advance_if(OPEN_BRACE)) {
      index_set reorder;
      size_type j = 0;
      dal::bit_vector check_permut;
      if (t.type() != tnode::TNTENSOR) ASM_THROW_PARSE_ERROR("can't use reorder braces on a constant!");
      for (;; ++j) {
	size_type i;
	if (tok_type() == COLON) i = j;
	else if (tok_type() == NUMBER) i = tok_number_ival(1000);
	else ASM_THROW_PARSE_ERROR("only numbers or colons allowed here");
	if (check_permut.is_in(i)) { /* on prend la diagonale du tenseur */
	  t = tnode(record(new ATN_diagonal_tensor(*t.tensor(), dim_type(i), dim_type(j))));
	  check_permut.add(j);
	  reorder.push_back(j);
	} else { 
	  check_permut.add(i);
	  reorder.push_back(i);
	}
	advance();
	if (advance_if(CLOSE_BRACE)) break;
	accept(COMMA, "expecting ','");
      }
      if (check_permut.first_false() != reorder.size()) {
	cerr << check_permut << endl;
	cerr << reorder << endl;
	ASM_THROW_PARSE_ERROR("you did not give a real permutation:" << reorder);
      }
      t = tnode(record(new ATN_permuted_tensor(*t.tensor(), reorder)));
    }
    return t;
  }

  /* 
     term := prod_trans*prod_trans/prod_trans ...
  */
  tnode generic_assembly::do_term() {
    push_mark();
    err_set_mark();
    tnode t = do_prod_trans();
    while (true) {
      bool mult;
      if (advance_if(MULTIPLY)) mult = true;
      else if (advance_if(DIVIDE)) mult = false;
      else break;
      tnode t2 = do_prod();
      if (mult == false && t.type() == tnode::TNCONST && t2.type() == tnode::TNTENSOR) 
	ASM_THROW_PARSE_ERROR("can't divide a constant by a tensor");
      if (t.type() == tnode::TNTENSOR && t2.type() == tnode::TNTENSOR) {
	ASM_THROW_PARSE_ERROR("tensor term-by-term productor division not implemented yet! are you sure you need it ?");
      } else if (t.type() == tnode::TNCONST && t2.type() == tnode::TNCONST) {
	if (mult) t.assign(t.xval()*t2.xval());
	else 	{ t2.check0(); t.assign(t.xval()/t2.xval()); }
      } else {
	if (t.type() != tnode::TNTENSOR) std::swap(t,t2);
	scalar_type v = t2.xval();
	if (!mult) {
	  if (v == 0.) { ASM_THROW_PARSE_ERROR("can't divide by zero"); }
	  else v = 1./v;
	}
	if (t.tensor()->is_tensors_sum_scaled() && !t.tensor()->is_frozen()) {
	  t.tensor()->is_tensors_sum_scaled()->apply_scale(v);
	} else {
	  t.assign(record(new ATN_tensors_sum_scaled(*t.tensor(), v)));
	}
      }
    }
    pop_mark();
    return t;
  }

  /* 
     expr := term + term - term + ...
     suboptimal for things like t1+1-2-1 (which gives (((t1+1)-2)-1) )
     ... could be fixed but noone needs that i guess
  */
  tnode generic_assembly::do_expr() {
    bool negt=false;
    push_mark();
    if (advance_if(MINUS)) negt = true;
    tnode t = do_term();
    if (negt) {
      if (t.type() == tnode::TNCONST) t.assign(-t.xval());
      else t.assign(record(new ATN_tensor_scalar_add(*t.tensor(), 0., -1)));
    }
    while (true) {
      int plus;
      if (advance_if(PLUS)) plus = +1;
      else if (advance_if(MINUS)) plus = -1;
      else break;
      tnode t2 = do_term();
      if (t.type() == tnode::TNTENSOR && t2.type() == tnode::TNTENSOR) {
	if (!t.tensor()->is_tensors_sum_scaled() || t.tensor()->is_frozen()) {
	  t.assign(record(new ATN_tensors_sum_scaled(*t.tensor(), +1))); 
	}
	t.tensor()->is_tensors_sum_scaled()->push_scaled_tensor(*t2.tensor(), scalar_type(plus));
      } else if (t.type() == tnode::TNCONST && t2.type() == tnode::TNCONST) {
	t.assign(t.xval()+t2.xval()*plus);
      } else {
	int tsgn = 1;	
	if (t.type() != tnode::TNTENSOR) { std::swap(t,t2); if (plus<0) tsgn = -1; }
	else if (plus<0) t2.assign(-t2.xval());
	t.assign(record(new ATN_tensor_scalar_add(*t.tensor(), t2.xval(), tsgn)));
      } 
    }
    pop_mark();
    return t;
  }

  /* instr := ident '=' expr |
     print expr |
     M(#mf,#mf) '+=' expr |
     V(#mf) '+=' expr */
  void generic_assembly::do_instr() {
    enum { wALIAS, wOUTPUT_ARRAY, wOUTPUT_MATRIX, wPRINT, wERROR } what = wERROR;
    std::string ident;

    /* get the rhs */
    if (tok_type() != IDENT) ASM_THROW_PARSE_ERROR("expecting identifier");
    if (vars.find(tok()) != vars.end()) ASM_THROW_PARSE_ERROR("redefinition of identifier " << tok());

    push_mark();
    ident = tok();
    advance();
     
    size_type print_mark = 0;
    size_type arg_num = size_type(-1);

    vdim_specif_list vds;

    if (ident.compare("print") == 0) {
      print_mark = tok_mark();
      what = wPRINT;
    } else if (tok_type() == ARGNUM_SELECTOR ||
	       tok_type() == OPEN_PAR) {
      if (tok_type() == ARGNUM_SELECTOR) {
	arg_num = tok_argnum();
	advance(); 
      } else { arg_num = 0; }
      
      do_dim_spec(vds);

      /* check the validity of the output statement */
      if (ident.compare("V")==0) {
	what = wOUTPUT_ARRAY;
	if (arg_num >= outvec.size()) { outvec.resize(arg_num+1); outvec[arg_num] = 0; }
	/* if we are allowed to dynamically create vectors */
	if (outvec[arg_num] == 0) {
	  if (vec_fact != 0) {
	    tensor_ranges r(vds.size());
	    for (size_type i=0; i < vds.size(); ++i) r[i] = vds[i].dim;
	    outvec[arg_num] = vec_fact->create_vec(r);
	  }
	  else ASM_THROW_PARSE_ERROR("output vector $" << arg_num+1 << " does not exist");
	}
      } else if (vds.nb_mf()==2 && vds.size() == 2 && ident.compare("M")==0) {
	what = wOUTPUT_MATRIX;
	if (arg_num >= outmat.size()) { outmat.resize(arg_num+1); outmat[arg_num] = 0; }
	/* if we are allowed to dynamically create matrices */
	if (outmat[arg_num] == 0) {
	  if (mat_fact != 0)
	    outmat[arg_num] = mat_fact->create_mat(vds[0].pmf->nb_dof(),
						   vds[1].pmf->nb_dof());
	  else ASM_THROW_PARSE_ERROR("output matrix $" << arg_num+1 << " does not exist");
	}
      } else ASM_THROW_PARSE_ERROR("not a valid output statement");

      accept(PLUS); 
      accept(EQUAL);
    } else if (advance_if(EQUAL)) {
      what = wALIAS;
    } else ASM_THROW_PARSE_ERROR("missing '=' or ':='");
      
    tnode t = do_expr();
    if (t.type() != tnode::TNTENSOR)
      ASM_THROW_PARSE_ERROR("left hand side is a constant, not a tensor!");
      
    switch (what) {
    case wPRINT: {
      record_out(new ATN_print_tensor(*t.tensor(), tok_substr(print_mark, 
							      tok_mark())));
    } break;
    case wOUTPUT_ARRAY: {
      ATN *pout = outvec[arg_num]->build_output_tensor(*t.tensor(), vds);
      record_out(pout);
    } break;
    case wOUTPUT_MATRIX: {
      ATN *pout = outmat[arg_num]->build_output_tensor(*t.tensor(), *vds[0].pmf, *vds[1].pmf);
      record_out(pout);
    } break;
    case wALIAS: {
      vars[ident] = t.tensor(); t.tensor()->freeze();
    } break;
    default: DAL_INTERNAL_ERROR(""); break;
    }
    pop_mark();
  }

  struct atn_number_compare {
    bool operator()(ATN *a, ATN *b) {
      assert(a && b);
      return (a->number() < b->number());
    }
  };

  void generic_assembly::parse() {
    if (parse_done) return;
    do {
      if (tok_type() == END) break;
      do_instr();
    } while (advance_if(SEMICOLON));
    if (tok_type() != END) ASM_THROW_PARSE_ERROR("unexpected token: '" << tok() << "'");
    if (outvars.size() == 0) cerr << "warning: assembly without output\n";

    /* reordering of atn_tensors and outvars */
    unsigned gcnt = 0;
    for (size_type i=0; i < outvars.size(); ++i) 
      outvars[i]->set_number(gcnt);

    std::sort(atn_tensors.begin(), atn_tensors.end(), atn_number_compare());
    std::sort(outvars.begin(), outvars.end(), atn_number_compare());

    /* remove non-numbered (ie unused) atn_tensors */
    while (atn_tensors.size() && atn_tensors.back()->number() == unsigned(-1)) {
      cerr << "warning: the expression " << atn_tensors.back()->name() << " won't be evaluated since it is not used!\n";
      delete atn_tensors.back(); atn_tensors.pop_back();
    }
    parse_done = true;
  }

  /* caution: the order of the loops is really important */
  void generic_assembly::exec(size_type cv, dim_type face) {
    bool update_shapes = false;
    for (size_type i=0; i < atn_tensors.size(); ++i) {
      atn_tensors[i]->check_shape_update(cv,face);
      update_shapes =  (update_shapes || atn_tensors[i]->is_shape_updated());
      /*if (atn_tensors[i]->is_shape_updated()) {
	cerr << "[cv=" << cv << ",f=" << int(face) << "], shape_updated: " //<< typeid(*atn_tensors[i]).name() 
             << " [" << atn_tensors[i]->name() << "]\n  -> r=" << atn_tensors[i]->ranges() << "\n    ";
	     }*/
    }

    if (update_shapes) {

      /*cerr << "updated shapes: cv=" << cv << " face=" << int(face) << ": ";
      for (size_type k=0; k < mftab.size(); ++k)
	cerr << mftab[k]->nb_dof_of_element(cv) << " "; cerr << "\n";
      */

      for (size_type i=0; i < atn_tensors.size(); ++i) {
	atn_tensors[i]->init_required_shape();
      }
      for (size_type i=0; i < outvars.size(); ++i) {
	outvars[i]->update_childs_required_shape();
      }
      for (size_type i=atn_tensors.size()-1; i!=size_type(-1); --i) {
	atn_tensors[i]->update_childs_required_shape();
      }
      for (size_type i=0; i < atn_tensors.size(); ++i) {
	atn_tensors[i]->reinit();
      }
      for (size_type i=0; i < outvars.size(); ++i) {
	outvars[i]->reinit();
      }
    }
    for (size_type i=0; i < atn_tensors.size(); ++i)
      atn_tensors[i]->exec(cv,face);
    for (size_type i=0; i < outvars.size(); ++i)
      outvars[i]->exec(cv, face);
  }

  struct cv_fem_compare {
    const std::deque<const mesh_fem *>& mf;
    cv_fem_compare(const std::deque<const mesh_fem *>& mf_) : mf(mf_) {}
    bool operator()(size_type a, size_type b) {
      for (size_type i=0; i < mf.size(); ++i) {
	/*
	  if (mf[i]->fem_of_element(a) < mf[i]->fem_of_element(b))
	  return true;
	if (mf[i]->fem_of_element(a) == mf[i]->fem_of_element(b) &&
	    mf[i]->nb_dof_of_element(a) < mf[i]->nb_dof_of_element(b))
	  return true;
	*/

	/* sort by nb_dof and then by fem */
	if (mf[i]->nb_dof_of_element(a) < mf[i]->nb_dof_of_element(b)) {
	  return true;
	}
	if (mf[i]->nb_dof_of_element(a) == mf[i]->nb_dof_of_element(b) && 
	    mf[i]->fem_of_element(a) < mf[i]->fem_of_element(b)) {
	  return true;
	}
      }
      return false;
    }
  };
  /* reorder the convexes in order to minimize the number of
     shape modifications during the assembly (since this can be
     very expensive) */
  static void get_convex_order(const dal::bit_vector& cvlst,
			       const std::deque<const mesh_im *>& imtab, 
			       const std::deque<const mesh_fem *>& mftab, 
			       const dal::bit_vector& candidates, 
			       std::vector<size_type>& cvorder) {
    cvorder.reserve(candidates.card()); cvorder.resize(0);
    for (dal::bv_visitor cv(candidates); !cv.finished(); ++cv) {
      if (cvlst.is_in(cv) && imtab[0]->int_method_of_element(cv) != im_none()) {
	bool ok = true;
        for (size_type i=0; i < mftab.size(); ++i)
          if (!mftab[i]->convex_index().is_in(cv)) {
	    ok = false;
            //ASM_THROW_ERROR("the convex " << cv << " has no FEM for the #" << i+1 << " mesh_fem");	  
	  }
        if (ok) cvorder.push_back(cv);
      } else if (!imtab[0]->linked_mesh().convex_index().is_in(cv)) {
        ASM_THROW_ERROR("the convex " << cv << " is not part of the mesh");
      } else {
	/* skip convexes without integration method */
      }
    }
    std::sort(cvorder.begin(), cvorder.end(), cv_fem_compare(mftab));
  }
  
  void generic_assembly::consistency_check() {
    //if (mftab.size() == 0) ASM_THROW_ERROR("no mesh_fem for assembly!");
    if (imtab.size() == 0) ASM_THROW_ERROR("no mesh_im (integration methods) given for assembly!");
    const getfem_mesh& m = imtab[0]->linked_mesh();
    for (unsigned i=0; i < mftab.size(); ++i) {
      if (&mftab[i]->linked_mesh() != &m)
	ASM_THROW_ERROR("the mesh_fem/mesh_im live on different meshes!");
    }
    for (unsigned i=0; i < imtab.size(); ++i) {
      if (&imtab[i]->linked_mesh() != &m)
	ASM_THROW_ERROR("the mesh_fem/mesh_im live on different meshes!");
    }
  }

  void generic_assembly::assembly(const mesh_region &r) {
    std::vector<size_type> cv;
    r.from_mesh(imtab.at(0)->linked_mesh());
    r.error_if_not_homogeneous();

    consistency_check();
    get_convex_order(imtab.at(0)->convex_index(), imtab, mftab, r.index(), cv);
    parse();
    for (size_type i=0; i < cv.size(); ++i) {
      mesh_region::face_bitset nf = r[cv[i]];
      dim_type f = dim_type(-1);
      while (nf.any()) {
	cerr << "generic_assembly::exec(" << cv[i] << ")\n";
	if (nf[0]) exec(cv[i],f);
	nf >>= 1; f+=1;
      }
    }
  }
} /* end of namespace */
