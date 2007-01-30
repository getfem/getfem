// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2003-2007 Julien Pommier
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

/**@file getfem_assembling_tensors.h
   @author  Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date January 2003.
   @brief Generic assembly implementation.
*/
#ifndef GETFEM_ASSEMBLING_TENSORS_H__
#define GETFEM_ASSEMBLING_TENSORS_H__

#include "gmm/gmm_kernel.h"
#include "getfem_mesh_fem.h"
#include "getfem_mesh_im.h"
#include "bgeot_sparse_tensors.h"
#include "getfem_mat_elem_type.h"
#include "getfem_mat_elem.h"
#include <map>

#define ASM_THROW_PARSE_ERROR(x) DAL_THROW(std::invalid_argument, "parse error: " << x << endl << "found here:\n " << syntax_err_print());
#define ASM_THROW_TENSOR_ERROR(x) DAL_THROW(std::invalid_argument, "tensor error: " << x);
#define ASM_THROW_ERROR(x) DAL_THROW(std::invalid_argument, "error: " << x);

namespace getfem {
  using bgeot::stride_type;
  using bgeot::index_type;
  using bgeot::index_set;
  using bgeot::tensor_ranges;
  using bgeot::tensor_strides;
  using bgeot::tensor_mask;
  using bgeot::tensor_shape;
  using bgeot::tensor_ref;
  using bgeot::multi_tensor_iterator;
  using bgeot::TDIter;

  class ATN_tensor;

  /*
     base class for the tree built from the expression of the tensor assembly
     (ATN == Assembly Tree Node)
  */
  class ATN {
    std::deque< ATN_tensor* > childs_;
    std::string name_;   /* the name is a part of the parsed string */
    unsigned number_;    /* a unique number, which is used for the ordering of the tree */
  protected:
    size_type current_cv;
    dim_type current_face;
  public:
    ATN(const std::string& n=std::string("unnamed")) : 
      name_(n), number_(unsigned(-1)), current_cv(size_type(-1)), current_face(dim_type(-1)) {}
    virtual ~ATN() {}

    void add_child(ATN_tensor& a) { childs_.push_back(&a); }
    ATN_tensor& child(size_type n) { return *childs_[n]; }
    size_type nchilds() { return childs_.size(); }
    /* reinit is called each time the object need to reset itself
       (when the shape of one of its childs has changed) */
    virtual void reinit() = 0;   
    /* do the computations for a given convex */
    void exec(size_type cv, dim_type face) {      
      if (cv != current_cv || face != current_face) {
	exec_(cv,face);
	current_cv = cv;
	current_face = face;
      }
    }
    const std::string& name() { return name_; }
    void set_name(const std::string& n) { name_ = n; }
    /* the "root" nodes expect to get all tensor values 
       others nodes have a more specific behavior
     */
    virtual void update_childs_required_shape();

    /* numbering og tensors, such that if i < j then tensor(j)
       cannot be in the sub-tree of tensor(i) */
    void set_number(unsigned &gcnt);
    unsigned number() const { return number_; }
  private:
    virtual void exec_(size_type , dim_type ) {}
  };
  
  class ATN_tensors_sum_scaled;

  /* Base class for every node except the "final" ones */
  class ATN_tensor : public ATN {
  protected:
    tensor_ranges r_;
    bool shape_updated_;
    tensor_ref   tr;
    tensor_shape req_shape;
    bool frozen_; /* used to recognize intermediate results of
		     computations stored in a temporary variable: they
		     cannot be modified a posteriori (like it could
		     happen with an ATN_tensors_sum_scaled) */
  public:
    ATN_tensor() { shape_updated_ = false; frozen_ = false; }
    bool is_shape_updated() const { return shape_updated_; }
    void freeze() { frozen_ = true; }
    bool is_frozen() const { return frozen_; }
    const tensor_ranges& ranges() const { return r_; }
    const tensor_shape& required_shape() const { return req_shape; }
    /* check_shape_update is called for each node of the tree
       if the shape of the tensor has been modified, the flag
       shape_updated_ should be set, and r_ should contain the
       new dimensions. This function is called in such an order
       that the shape updates are automatically propagated in the tree */
    virtual void check_shape_update(size_type , dim_type) {}
    /* if the shape was updated, the node should initialise its req_shape */
    virtual void init_required_shape() { req_shape.set_empty(r_); }
    /* then each node update the req_shape of its childs.
     */
    virtual void update_childs_required_shape() {
      for (dim_type i=0; i < nchilds(); ++i) {
	child(i).merge_required_shape(req_shape);
      }
    }
    /* ... then reserve some memory if necessary for tensor storage 
       in 'reinit' (inherited here from ATN)
     */
    tensor_ref& tensor() { 
      return tr; 
    }

    void merge_required_shape(const tensor_shape& shape_from_parent) {
      req_shape.merge(shape_from_parent, false);
    }
    /* recognize sums of scaled tensors
       (in order to stack those sums on the same object)
       (dynamic_cast prohibited for the moment: crashes matlab) */
    virtual ATN_tensors_sum_scaled* is_tensors_sum_scaled() { return 0; }
  };


  /* simple list of "virtual" dimensions, i.e. which may be constant
     or be given by a mesh_fem */
  struct vdim_specif {
    size_type dim;
    const mesh_fem *pmf;
    bool is_mf_ref() const { return (pmf != 0); }
    vdim_specif() { dim = size_type(-1); pmf = 0; }
    vdim_specif(size_type i) { dim = i; pmf = 0; }
    vdim_specif(const mesh_fem *pmf_) { dim = pmf_->nb_dof(); pmf = pmf_; }
  };
  class vdim_specif_list : public std::vector< vdim_specif > {
  public:
    vdim_specif_list() { reserve(8); }
    size_type nb_mf() const;
    size_type nbelt() const;
    void build_strides_for_cv(size_type cv, tensor_ranges& r, 
			      std::vector<tensor_strides >& str) const;
  };

  /* final node for array output: array means full array of 0,1,2,3 or more dimensions,
     stored in a vector VEC in fortran order
  */     
  template< typename VEC > class ATN_array_output : public ATN {
    VEC& v;
    vdim_specif_list vdim;
    multi_tensor_iterator mti;
    tensor_strides strides;
  public:
    ATN_array_output(ATN_tensor& a, VEC& v_, vdim_specif_list &d) : v(v_), vdim(d) {
      strides.resize(vdim.size()+1);
      add_child(a);
      strides[0] = 1;
      for (size_type i=0; i < vdim.size(); ++i) strides[i+1] = strides[i]*vdim[i].dim;
      if (gmm::vect_size(v) != size_type(strides[vdim.size()])) 
	ASM_THROW_TENSOR_ERROR("wrong size for output vector: supplied vector size is " << 
			       gmm::vect_size(v) << " while it should be " << strides[vdim.size()]);      
    }
  private:
    void reinit() {
      mti = multi_tensor_iterator(child(0).tensor(),true);
    }
    void exec_(size_type cv, dim_type) {
      tensor_ranges r;
      std::vector< tensor_strides > str;
      vdim.build_strides_for_cv(cv, r, str);
      if (child(0).ranges() != r) {
	ASM_THROW_TENSOR_ERROR("can't output a tensor of dimensions " << child(0).ranges() << 
			       " into an output array of size " << r);
      }
      mti.rewind();
      do {
	typename gmm::linalg_traits<VEC>::iterator it = gmm::vect_begin(v);
	for (dim_type i = 0; i < mti.ndim(); ++i) it+=str[i][mti.index(i)];
	*it += mti.p(0);
      } while (mti.qnext1());
    }
  };

  /* final node for sparse matrix output */
  template< typename MAT > class ATN_smatrix_output : public ATN {
    const mesh_fem &mf_r, &mf_c;
    MAT& m;
    multi_tensor_iterator mti;
    struct ijv { // just a fast cache for the mti output (yes it makes a small difference)
      scalar_type *p;
      unsigned i,j;
    };
    std::vector<ijv> it;
  public:
    ATN_smatrix_output(ATN_tensor& a, const mesh_fem& mf_r_, 
		       const mesh_fem& mf_c_, MAT& m_) 
      : mf_r(mf_r_), mf_c(mf_c_), m(m_) { 
      add_child(a); 
      it.reserve(100);
    }
  private:
    void reinit() {
      mti = multi_tensor_iterator(child(0).tensor(),true);
      it.resize(0);
    }
    void exec_(size_type cv, dim_type) {
      size_type nb_r = mf_r.nb_dof_of_element(cv);
      size_type nb_c = mf_c.nb_dof_of_element(cv);
      if (child(0).tensor().ndim() != 2)
	ASM_THROW_TENSOR_ERROR("cannot write a " << 
			       int(child(0).tensor().ndim()) << 
			       "D-tensor into a matrix!");
      if (child(0).tensor().dim(0) != nb_r ||
	  child(0).tensor().dim(1) != nb_c) {
	ASM_THROW_TENSOR_ERROR("size mismatch for sparse matrix output:"
			       " tensor dimension is " << 
			       child(0).ranges()
			       << ", while the elementary matrix for convex "
			       << cv << " should have " << nb_r << "x"
			       << nb_c << " elements");
      }
      std::vector<size_type> cvdof_r(mf_r.ind_dof_of_element(cv).begin(), mf_r.ind_dof_of_element(cv).end());
      std::vector<size_type> cvdof_c(mf_c.ind_dof_of_element(cv).begin(), mf_c.ind_dof_of_element(cv).end());
      /*mti.rewind();
      do {
	if (mti.p(0)) {
	  size_type dof_i = cvdof_r[mti.index(0)];
	  size_type dof_j = cvdof_c[mti.index(1)];
	  m(dof_i, dof_j) += mti.p(0);
	}
      } while (mti.qnext1());
      */
      if (it.size() == 0) {
	mti.rewind();
	do {
	  ijv v;
	  v.p = &mti.p(0);
	  v.i = mti.index(0);
	  v.j = mti.index(1);
	  it.push_back(v);
	} while (mti.qnext1());
      }
      for (unsigned i=0; i < it.size(); ++i)
	if (*it[i].p) 
	  m(cvdof_r[it[i].i], cvdof_c[it[i].j]) += *it[i].p;
    }
  };



  /* some wrappers : their aim is to provide a better genericity,
     and to avoid the whole templatization of the 'generic_assembly' class,
     which is quite(!) big
  */
  class base_asm_data {
  public:
    virtual size_type vect_size() const = 0;
    virtual void copy_with_mti(const std::vector<tensor_strides> &, multi_tensor_iterator &) const = 0;
    virtual ~base_asm_data() {}
  };

  template< typename VEC > class asm_data : public base_asm_data {
    const VEC &v;
  public:
    asm_data(const VEC *v_) : v(*v_) {}
    size_type vect_size() const {
      return gmm::vect_size(v); 
    }
    /* used to transfert the data for the current convex to the mti of ATN_tensor_from_dofs_data */
    void copy_with_mti(const std::vector<tensor_strides> &str, multi_tensor_iterator &mti) const {
      size_type ppos;
      do {
	ppos = 0;
	for (unsigned i = 0; i < mti.ndim(); ++i) ppos+=str[i][mti.index(i)];
	mti.p(0) = v[ppos];
      } while (mti.qnext1());
    }
  };

  class base_asm_vec {
  public:
    virtual ATN* build_output_tensor(ATN_tensor &a, 
					    vdim_specif_list& vdim)=0;
    virtual ~base_asm_vec() {}
  };

  template< typename VEC > class asm_vec : public base_asm_vec {
    VEC *v;
  public:
    asm_vec(VEC *v_) : v(v_) {}
    virtual ATN* build_output_tensor(ATN_tensor &a, 
				     vdim_specif_list& vdim) {
      ATN *t = new ATN_array_output<VEC>(a, *v, vdim); return t;
    }
    VEC *vec() { return v; }
    ~asm_vec() {}
  };

  /* the "factory" is only useful for the matlab interface,
     since the number of output arrays and sparse matrices is unknown
     for user-supplied assemblies. Hence they are created "on-the-fly" */
  class base_vec_factory {
  public:
    virtual base_asm_vec* create_vec(const tensor_ranges& r) = 0;
    virtual ~base_vec_factory() {}
  };

  template< typename VEC > class vec_factory : public base_vec_factory, private std::deque<asm_vec<VEC> > {
  public:
    base_asm_vec* create_vec(const tensor_ranges& r) {
      size_type sz = 1; for (size_type i=0; i < r.size(); ++i) sz *= r[i];
      if (sz == 0) ASM_THROW_TENSOR_ERROR("can't create a vector of size " << r);
      asm_vec<VEC> v(new VEC(sz));
      push_back(v); return &this->back();
    }
    ~vec_factory() { 
      for (size_type i=0; i < this->size(); ++i) {
	delete (*this)[i].vec();
      }
    }
  };


  /* matrix wrappers */
  class base_asm_mat {
  public:
    virtual ATN*
    build_output_tensor(ATN_tensor& a, const mesh_fem& mf1, const mesh_fem& mf2) = 0;
    virtual ~base_asm_mat() {}
  };

  template< typename MAT > class asm_mat : public base_asm_mat {
    MAT *m;
  public:
    asm_mat(MAT* m_) : m(m_) {}
    ATN*
    build_output_tensor(ATN_tensor& a, const mesh_fem& mf1, const mesh_fem& mf2) {
      return new ATN_smatrix_output<MAT>(a, mf1, mf2, *m);
    }
    MAT *mat() { return m; }
    ~asm_mat() {}
  };

  class base_mat_factory {
  public:
    virtual base_asm_mat* create_mat(size_type m, size_type n) = 0;
    virtual ~base_mat_factory() {};
  };

  template< typename MAT > class mat_factory : public base_mat_factory, private std::deque<asm_mat<MAT> > {
  public:
    base_asm_mat* create_mat(size_type m, size_type n) { 
      push_back(asm_mat<MAT>(new MAT(m, n))); return &this->back();
    }
    ~mat_factory() { 
      for (size_type i=0; i < this->size(); ++i) {
	delete ((*this)[i]).mat(); 
      }
    }
  };


  class tnode {
  public:
    typedef enum { TNCONST, TNTENSOR, TNNONE } node_type;
  private:
    node_type type_;
    scalar_type x;
    ATN_tensor *t;
  public:
    tnode() : type_(TNNONE), x(1e300), t(NULL) {}
    tnode(scalar_type x_) { assign(x_); }
    tnode(ATN_tensor *t_) { assign(t_); }
    void assign(scalar_type x_) { type_ = TNCONST; t = NULL; x = x_; }
    void assign(ATN_tensor *t_) { type_ = TNTENSOR; t = t_; x = 1e300; }
    ATN_tensor* tensor() { assert(type_ == TNTENSOR); return t; }
    scalar_type xval() { assert(type_ == TNCONST); return x; }
    node_type type() { return type_; }
    void check0() { if (xval() == 0) ASM_THROW_ERROR("division by zero"); }
  };

  class asm_tokenizer {
  public:
    typedef enum { OPEN_PAR='(', CLOSE_PAR=')', COMMA=',', 
		   SEMICOLON=';', COLON=':', EQUAL='=', MFREF='#', IMREF='%', 
		   PLUS='+',MINUS='-', PRODUCT='.',MULTIPLY='*', 
		   DIVIDE='/', ARGNUM_SELECTOR='$', 
		   OPEN_BRACE='{', CLOSE_BRACE='}', 
		   END=0, IDENT=1, NUMBER=2 } tok_type_enum;
  private:
    std::string str;
    size_type tok_pos, tok_len;
    tok_type_enum curr_tok_type;
    std::string curr_tok;
    int curr_tok_ival;
    double curr_tok_dval;
    size_type err_msg_mark;
    std::deque<size_type> marks;
  public:
    asm_tokenizer() {}
    void set_str(const std::string& s_) {
      str = s_; tok_pos = 0; tok_len = size_type(-1); curr_tok_type = END;
      err_msg_mark = 0; get_tok(); 
    }
    std::string tok() const { return curr_tok; }
    tok_type_enum tok_type() const { return curr_tok_type; }
    size_type tok_mark() { return tok_pos; }
    std::string tok_substr(size_type i1, size_type i2) { return str.substr(i1, i2-i1); }
    void err_set_mark() {
      err_msg_mark = tok_pos;
    }
    void push_mark() { marks.push_back(tok_pos); }
    void pop_mark() { assert(marks.size()); marks.pop_back(); }
    std::string mark_txt() { assert(marks.size()); return tok_substr(marks.back(),tok_pos); }

    /* returns a friendly message indicated the location of the syntax error */
    std::string syntax_err_print();
    void accept(tok_type_enum t, const char *msg_="syntax error") { 
      if (tok_type() != t) ASM_THROW_PARSE_ERROR(msg_); advance();
    }
    void accept(tok_type_enum t, tok_type_enum t2, const char *msg_="syntax error") {
      if (tok_type() != t && tok_type() != t2) ASM_THROW_PARSE_ERROR(msg_); advance();
    }
    bool advance_if(tok_type_enum t) { 
      if (tok_type() == t) { advance(); return true; } else return false; 
    }
    void advance() {
      tok_pos += tok_len;
      get_tok();
    }
    void get_tok();
    double tok_number_dval() { assert(tok_type()==NUMBER); return curr_tok_dval; }
    int tok_number_ival(int maxval=10000000) {
      int n=int(tok_number_dval());
      if (n != curr_tok_dval) ASM_THROW_PARSE_ERROR("non an integer"); 
      if (n > maxval) ASM_THROW_PARSE_ERROR("out of bound integer");
      return n-1; /* -1 pour un indicage qui commence à 1! */ 
    }
    size_type tok_mfref_num() { assert(tok_type()==MFREF); return curr_tok_ival; }
    size_type tok_imref_num() { assert(tok_type()==IMREF); return curr_tok_ival; }
    size_type tok_argnum() { assert(tok_type()==ARGNUM_SELECTOR); return curr_tok_ival; }
  };




  /** Generic assembly of vectors, matrices. 
      
      Many examples of use available @link asm here@endlink.
   */
  class generic_assembly : public asm_tokenizer {
    std::vector<const mesh_fem *> mftab;/* list of the mesh_fem used in the computation */
    std::vector<const mesh_im *> imtab;/* list of the mesh_im used in the computation */
    std::vector<pnonlinear_elem_term> innonlin;  /* alternatives to base, grad, hess in comp() for non-linear computations) */
    std::vector<base_asm_data*> indata;              /* data sources */
    std::vector<base_asm_vec*> outvec;               /* vectors in which is done the assembly */
    std::vector<base_asm_mat*> outmat;               /* matrices in which is done the assembly */

    base_vec_factory *vec_fact; /* if non null, used to fill the outvec list with a given vector class */
    base_mat_factory *mat_fact; /* if non null, used to fill the outmat list with a given matrix class */

    std::vector<ATN*> outvars;        /* the list of "final tensors" which produce some
						   output in outvec and outmat */
    
    std::map<std::string, ATN_tensor *> vars; /* the list of user variables */
    std::vector<ATN_tensor*> atn_tensors;     /* keep track of all tensors objects (except the ones listed in 'outvars')
						for deallocation when all is done 
						note that they are not stored in a random order, but are reordered such that
						the childs of the i-th ATN_tensor are all stored at indices j < i
						this assumption is largely used for calls to shape updates and exec(cv,f)
					     */
    bool parse_done;

  public:
    generic_assembly() : vec_fact(0), mat_fact(0), parse_done(false) {}
    generic_assembly(const std::string& s_) :
      vec_fact(0), mat_fact(0), parse_done(false)
    { set_str(s_); }
    generic_assembly(const std::string& s_,
	       std::vector<const mesh_fem*>& mftab_, 
	       std::vector<const mesh_im*>& imtab_, 
	       std::vector<base_asm_data*> indata_,
	       std::vector<base_asm_mat*> outmat_,
	       std::vector<base_asm_vec*> outvec_) : 
      mftab(mftab_), imtab(imtab_),
      indata(indata_), outvec(outvec_), outmat(outmat_),
      vec_fact(0), mat_fact(0), parse_done(false)
    { set_str(s_); }    
    ~generic_assembly() {
      for (size_type i = 0; i < atn_tensors.size(); ++i) delete atn_tensors[i];
      for (size_type i = 0; i < outvars.size(); ++i) delete outvars[i];
      for (size_type i = 0; i < indata.size(); ++i) delete indata[i];
      /* the destruction of outvec and outmat is assured, if necessary
	 by the vec_fact and asm_fact (since they derive from deque<asm_mat>) */
      if (vec_fact==0) for (size_type i = 0; i < outvec.size(); ++i) delete outvec[i];
      if (mat_fact==0) for (size_type i = 0; i < outmat.size(); ++i) delete outmat[i];
    }

    void set(const std::string& s_) { set_str(s_); }
    const std::vector<const mesh_fem*>& mf() const { return mftab; }
    const std::vector<const mesh_im*>& im() const { return imtab; }
    const std::vector<pnonlinear_elem_term> nonlin() const { return innonlin; }
    const std::vector<base_asm_data*>& data() const { return indata; }
    const std::vector<base_asm_vec*>& vec() const { return outvec; }
    const std::vector<base_asm_mat*>& mat() const { return outmat; }
    /// Add a new mesh_fem
    void push_mf(const mesh_fem& mf_) { mftab.push_back(&mf_); }
    /// Add a new mesh_im
    void push_mi(const mesh_im& im_) { imtab.push_back(&im_); }
    /// Add a new non-linear term
    void push_nonlinear_term(pnonlinear_elem_term net) {
      innonlin.push_back(net);
    }
    /// Add a new data (dense array)
    template< typename VEC > void push_data(const VEC& d) { 
      indata.push_back(new asm_data<VEC>(&d)); 
    }
    /// Add a new output vector
    template< typename VEC > void push_vec(VEC& v) { 
      asm_vec<VEC> *pv = new asm_vec<VEC>(&(gmm::linalg_cast(v)));
      outvec.push_back(pv);
    }
    /// Add a new output vector (fake const version..)
    template< typename VEC > void push_vec(const VEC& v) { 
      asm_vec<VEC> *pv = new asm_vec<VEC>(&(gmm::linalg_cast(v)));
      outvec.push_back(pv);
    }
    /// Add a new output matrix (fake const version..)
    template< typename MAT > void push_mat(const MAT& m) { 
      outmat.push_back(new asm_mat<MAT>(&(gmm::linalg_cast(m)))); 
    }
    /// Add a new output matrix
    template< typename MAT > void push_mat(MAT& m) { 
      outmat.push_back(new asm_mat<MAT>(&(gmm::linalg_cast(m)))); 
    }

    template <typename T> void push_mat_or_vec(T &v) {
      push_mat_or_vec(v, typename gmm::linalg_traits<T>::linalg_type());
    }

    /// used by the getfem_interface..
    void set_vec_factory(base_vec_factory *fact) { vec_fact = fact; }
    void set_mat_factory(base_mat_factory *fact) { mat_fact = fact; }

  private:
    ATN_tensor* record(ATN_tensor *t) {
      t->set_name(mark_txt());
      atn_tensors.push_back(t); return t;
    }
    ATN* record_out(ATN *t) {
      t->set_name(mark_txt());
      outvars.push_back(t); return t;
    }
    const mesh_fem& do_mf_arg_basic();
    const mesh_fem& do_mf_arg(std::vector<const mesh_fem*> *multimf = 0);
    void do_dim_spec(vdim_specif_list& lst);
    std::string do_comp_red_ops();
    ATN_tensor* do_comp();
    ATN_tensor* do_data();
    std::pair<ATN_tensor*, std::string> do_red_ops(ATN_tensor* t);
    tnode do_tens();
    tnode do_prod();
    tnode do_prod_trans();
    tnode do_term();
    tnode do_expr();
    void do_instr();
    void exec(size_type cv, dim_type face);
    void consistency_check();
    template <typename T> void push_mat_or_vec(T &v, gmm::abstract_vector) {
      push_vec(v);
    }
    template <typename T> void push_mat_or_vec(T &v, gmm::abstract_matrix) {
      push_mat(v);
    }
  public:
    /* parse the string 'str' and build the tree of vtensors */
    void parse();
    /* do the assembly on the whole mesh */
    //void volumic_assembly();
    /* do the assembly on the specified convexes */
    //void volumic_assembly(const dal::bit_vector& cvlst);
    /* do the assembly on the specified boundary */
    //void boundary_assembly(size_type boundary_number);

    /** do the assembly on the specified region (boundary or set of convexes) */
    void assembly(const mesh_region &region = 
		  mesh_region::all_convexes());
  };
}  /* end of namespace getfem.                                             */


#endif /* GETFEM_ASSEMBLING_TENSORS_H__  */
