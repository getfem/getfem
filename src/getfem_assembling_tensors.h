/* -*- c++ -*- (enables emacs c++ mode)                                    */

#ifndef __GETFEM_ASSEMBLING_TENSORS_H
#define __GETFEM_ASSEMBLING_TENSORS_H

#include <gmm.h>
#include <getfem_mesh_fem.h>
#include <bgeot_sparse_tensors.h>
#include <numeric>
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

  /* 
     base class for the tree built from the expression of the tensor assembly
     (ATN == Assembly Tree Node)
  */
  class ATN_tensor;

  class ATN {
    std::deque< ATN_tensor* > _childs;
    std::string _name;   /* the name is a part of the parsed string */
    unsigned _number;    /* a unique number, which is used for the ordering of the tree */
  protected:
    size_type current_cv;
    dim_type current_face;
  public:
    ATN(const std::string& n=std::string("unnamed")) : 
      _name(n), _number(unsigned(-1)), current_cv(size_type(-1)), current_face(dim_type(-1)) {}
    virtual ~ATN() {}

    void add_child(ATN_tensor& a) { _childs.push_back(&a); }
    ATN_tensor& child(size_type n) { return *_childs[n]; }
    size_type nchilds() { return _childs.size(); }
    /* reinit is called each time the object need to reset itself
       (when the shape of one of its childs has changed) */
    virtual void reinit() = 0;   
    /* do the computations for a given convex */
    void exec(size_type cv, dim_type face) {      
      if (cv != current_cv || face != current_face) {
	_exec(cv,face);
	current_cv = cv;
	current_face = face;
      }
    }
    const std::string& name() { return _name; }
    void set_name(const std::string& n) { _name = n; }
    /* the "root" nodes expect to get all tensor values 
       others nodes have a more specific behavior
     */
    virtual void update_childs_required_shape();

    /* numérotation des tenseurs, telle que si i < j alors le tenseur(j) 
       ne peut pas etre dans le sous-arbre du tenseur(i) */
    void set_number(unsigned &gcnt);
    unsigned number() const { return _number; }
  private:
    virtual void _exec(size_type , dim_type ) {}
  };
  
  class ATN_tensors_sum_scaled;

  /* Base class for every node except the "final" ones */
  class ATN_tensor : public ATN {
  protected:
    tensor_ranges _r;
    bool _shape_updated;
    tensor_ref   tr;
    tensor_shape req_shape;
    bool _frozen; /* pour reconnaitre les resultats intermédiaires de calculs
		     stockés dans une variable temporaire: ils ne peuvent pas être
		     modifiés à posteriori (comme ça pourrait arriver avec un 
		     ATN_tensors_sum_scaled) */
  public:
    ATN_tensor() { _shape_updated = false; _frozen = false; }
    bool is_shape_updated() const { return _shape_updated; }
    void freeze() { _frozen = true; }
    bool is_frozen() const { return _frozen; }
    const tensor_ranges& ranges() const { return _r; }
    const tensor_shape& required_shape() const { return req_shape; }
    /* check_shape_update is called for each node of the tree
       if the shape of the tensor has been modified, the flag
       _shape_updated should be set, and _r should contain the
       new dimensions. This function is called in such an order
       that the shape updates are automatically propagated in the tree */
    virtual void check_shape_update(size_type , dim_type) {}
    /* if the shape was updated, the node should initialise its req_shape */
    virtual void init_required_shape() { req_shape.set_empty(_r); }
    /* then each node update the req_shape of its childs.
     */
    virtual void update_childs_required_shape() {
      for (dim_type i=0; i < nchilds(); ++i) {
	child(i).merge_required_shape(req_shape);
      }
    }
    /* ... puis reserve de la memoire si necessaire pour le stockage du tenseur 
       dans 'reinit' (heritee ici de la class ATN)
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
    vdim_specif(const mesh_fem *_pmf) { dim = _pmf->nb_dof(); pmf = _pmf; }
  };
  class vdim_specif_list : public std::deque< vdim_specif > {
  public:
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
    ATN_array_output(ATN_tensor& a, VEC& _v, vdim_specif_list &d) : v(_v), vdim(d) {
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
    void _exec(size_type cv, dim_type) {
      tensor_ranges r;
      std::vector< tensor_strides > str;
      vdim.build_strides_for_cv(cv, r, str);
      if (child(0).ranges() != r) {
	ASM_THROW_TENSOR_ERROR("can't output a tensor of dimensions " << child(0).ranges() << 
			       " into an output array of size " << r);
      }
      mti.rewind();
      do {
	typename VEC::iterator it = v.begin();
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
  public:
    ATN_smatrix_output(ATN_tensor& a, const mesh_fem& _mf_r, 
		       const mesh_fem& _mf_c, MAT& _m) 
      : mf_r(_mf_r), mf_c(_mf_c), m(_m) { 
      add_child(a); 
    }
  private:
    void reinit() {
      mti = multi_tensor_iterator(child(0).tensor(),true);
    }
    void _exec(size_type cv, dim_type) {
      size_type nb_r = mf_r.nb_dof_of_element(cv);
      size_type nb_c = mf_c.nb_dof_of_element(cv);
      if (child(0).tensor().dim(0) != nb_r ||
	  child(0).tensor().dim(1) != nb_c) {
	ASM_THROW_TENSOR_ERROR("size mismatch for sparse matrix output: tensor dimension is " << 
			       child(0).ranges() << ", while the elementary matrix for convex " << cv << 
			       " should have " << nb_r << "x" << nb_c << " elements");
      }
      std::vector<size_type> cvdof_r(nb_r);
      std::vector<size_type> cvdof_c(nb_c);
      std::copy(mf_r.ind_dof_of_element(cv).begin(), mf_r.ind_dof_of_element(cv).end(), cvdof_r.begin());
      std::copy(mf_c.ind_dof_of_element(cv).begin(), mf_c.ind_dof_of_element(cv).end(), cvdof_c.begin());
      mti.rewind();
      do {
	if (mti.p(0)) {
	  size_type dof_i = cvdof_r[mti.index(0)];
	  size_type dof_j = cvdof_c[mti.index(1)];
	  m(dof_i, dof_j) += mti.p(0);
	}
      } while (mti.qnext1());
    }
  };



  /* some wrappers : their aim is to provide a better genericity,
     and to avoid the whole templatization of the 'generic_assembly' class,
     which is quite(!) big
  */
  class base_asm_data {
  public:
    virtual size_type vect_size() const = 0;
    virtual const scalar_type *get_data_ptr() const = 0;
    virtual ~base_asm_data() {}
  };

  template< typename VEC > class asm_data : public base_asm_data {
    VEC *v;
  public:
    asm_data(VEC *_v) : v(_v) {}
    size_type vect_size() const {
      return gmm::vect_size(*v); 
    }
    const scalar_type *get_data_ptr() const { return &((*v)[0]); }
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
    asm_vec(VEC *_v) : v(_v) {}
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
    asm_mat(MAT* _m) : m(_m) {}
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
    typedef enum { CONST, TENSOR, NONE } node_type;
  private:
    node_type _type;
    scalar_type x;
    ATN_tensor *t;
  public:
    tnode() : _type(NONE), x(1e300), t(NULL) {}
    tnode(scalar_type _x) { assign(_x); }
    tnode(ATN_tensor *_t) { assign(_t); }
    void assign(scalar_type _x) { _type = CONST; t = NULL; x = _x; }
    void assign(ATN_tensor *_t) { _type = TENSOR; t = _t; x = 1e300; }
    ATN_tensor* tensor() { assert(_type == TENSOR); return t; }
    scalar_type xval() { assert(_type == CONST); return x; }
    node_type type() { return _type; }
    void check0() { if (xval() == 0) ASM_THROW_ERROR("division by zero"); }
  };

  class asm_tokenizer {
  public:
    typedef enum { OPEN_PAR='(', CLOSE_PAR=')', COMMA=',', 
		   SEMICOLON=';', COLON=':', EQUAL='=', MFREF='#', 
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
    void set_str(const std::string& _s) {
      str = _s; tok_pos = 0; tok_len = size_type(-1); curr_tok_type = END;
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
    void accept(tok_type_enum t, const char *_msg="syntax error") { 
      if (tok_type() != t) ASM_THROW_PARSE_ERROR(_msg); advance();
    }
    void accept(tok_type_enum t, tok_type_enum t2, const char *_msg="syntax error") {
      if (tok_type() != t && tok_type() != t2) ASM_THROW_PARSE_ERROR(_msg); advance();
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
    size_type tok_argnum() { assert(tok_type()==ARGNUM_SELECTOR); return curr_tok_ival; }
  };




  /* main class for generic assembly */
  class generic_assembly : public asm_tokenizer {
    std::deque<const mesh_fem *> mftab;/* list of the mesh_fem used in the computation */
  
    std::deque<base_asm_data*> indata;              /* data sources */
    std::deque<base_asm_vec*> outvec;               /* vectors in which is done the assembly */
    std::deque<base_asm_mat*> outmat;               /* matrices in which is done the assembly */

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
    generic_assembly(const std::string& _s) :
      vec_fact(0), mat_fact(0), parse_done(false)
    { set_str(_s); }
    generic_assembly(const std::string& _s,
	       std::deque<const mesh_fem*>& _mftab, 
	       std::deque<base_asm_data*> _indata,
	       std::deque<base_asm_mat*> _outmat,
	       std::deque<base_asm_vec*> _outvec) : 
      mftab(_mftab), 
      indata(_indata), outvec(_outvec), outmat(_outmat),
      vec_fact(0), mat_fact(0), parse_done(false)
    { set_str(_s); }    
    ~generic_assembly() {
      for (size_type i = 0; i < atn_tensors.size(); ++i) delete atn_tensors[i];
      for (size_type i = 0; i < outvars.size(); ++i) delete outvars[i];
      for (size_type i = 0; i < indata.size(); ++i) delete indata[i];
      /* the destruction of outvec and outmat is assured, if necessary
	 by the vec_fact and asm_fact (since they derive from deque<asm_mat>) */
      if (vec_fact==0) for (size_type i = 0; i < outvec.size(); ++i) delete outvec[i];
      if (mat_fact==0) for (size_type i = 0; i < outmat.size(); ++i) delete outmat[i];
    }

    void set(const std::string& _s) { set_str(_s); }
    const std::deque<const mesh_fem*>& mf() const { return mftab; }
    const std::deque<base_asm_data*>& data() const { return indata; }
    const std::deque<base_asm_vec*>& vec() const { return outvec; }
    const std::deque<base_asm_mat*>& mat() const { return outmat; }
    void push_mf(const mesh_fem& _mf) { mftab.push_back(&_mf); }
    template< typename VEC > void push_data(VEC& d) { 
      indata.push_back(new asm_data<VEC>(&d)); 
    }
    template< typename VEC > void push_vec(VEC& v) { 
      asm_vec<VEC> *pv = new asm_vec<VEC>(&v);
      outvec.push_back(pv);
    }
    template< typename MAT > void push_mat(MAT& m) { 
      outmat.push_back(new asm_mat<MAT>(&m)); 
    }

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
    const mesh_fem& do_mf_arg();
    void do_dim_spec(vdim_specif_list& lst);
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
  public:
    /* parse the string 'str' and build the tree of vtensors */
    void parse();
    /* do the assembly on the whole mesh */
    void volumic_assembly();
    /* do the assembly on the specified convexes */
    void volumic_assembly(const dal::bit_vector& cvlst);
    /* do the assembly on the specified boundary */
    void boundary_assembly(size_type boundary_number);
  };
}  /* end of namespace getfem.                                             */


#endif /* __GETFEM_ASSEMBLING_TENSORS_H  */
