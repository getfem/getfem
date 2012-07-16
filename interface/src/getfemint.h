/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2001-2012 Y. Renard, J. Pommier.
 
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

/**\file getfemint.h
   \brief main getfem-interface header.
*/

#ifndef GETFEMINT_H__
#define GETFEMINT_H__

#include <getfemint_std.h>
#include <set>
#include <getfem/getfem_mesh_fem.h>
#include <gfi_array.h>
#include <getfem/dal_shared_ptr.h>

namespace getfem {
  class stored_mesh_slice;
  class mesh;
  class mesh_fem;
  class mat_elem_type;
  typedef boost::intrusive_ptr<const mat_elem_type> pmat_elem_type;
  class level_set;
  class mesh_level_set;
  class abstract_xy_function;
  class mesher_signed_distance;
  class cont_struct_getfem_model;
}

namespace getfemint
{
  /* exception-throwing version of the allocation functions of gfi_array.h */
  gfi_array* checked_gfi_array_create(int ndim, const int *dims, gfi_type_id type, gfi_complex_flag is_complex = GFI_REAL);
  gfi_array* checked_gfi_array_create_0(gfi_type_id type, gfi_complex_flag is_complex = GFI_REAL);
  gfi_array* checked_gfi_array_create_1(int M, gfi_type_id type, gfi_complex_flag is_complex = GFI_REAL);
  gfi_array* checked_gfi_array_create_2(int M, int N, gfi_type_id type, gfi_complex_flag is_complex = GFI_REAL);
  gfi_array* checked_gfi_array_from_string(const char*s);
  gfi_array* checked_gfi_create_sparse(int m, int n, int nzmax, gfi_complex_flag is_complex);

  typedef bgeot::dim_type dim_type;
  typedef bgeot::scalar_type scalar_type;
  typedef std::complex<double> complex_type;

  inline int *gfi_get_data(gfi_array *g, int) { return gfi_int32_get_data(g); }
  inline unsigned *gfi_get_data(gfi_array *g, unsigned) { return gfi_uint32_get_data(g); }
  inline char *gfi_get_data(gfi_array *g, char) { return gfi_char_get_data(g); }
  inline double *gfi_get_data(gfi_array *g, double) { return gfi_double_get_data(g); }
  inline complex_type *gfi_get_data(gfi_array *g, complex_type) { return (complex_type*)gfi_double_get_data(g); }


  typedef gmm::row_matrix<gmm::wsvector<scalar_type> >  gf_real_sparse_by_row;
  typedef gmm::col_matrix<gmm::wsvector<scalar_type> >  gf_real_sparse_by_col;
  typedef gmm::csc_matrix_ref<const double *, const unsigned int *, const unsigned int *>
  gf_real_sparse_csc_const_ref;
  typedef gmm::row_matrix<gmm::wsvector<complex_type> >  gf_cplx_sparse_by_row;
  typedef gmm::col_matrix<gmm::wsvector<complex_type> >  gf_cplx_sparse_by_col;
  typedef gmm::csc_matrix_ref<const complex_type *, const unsigned int *, const unsigned int *>
  gf_cplx_sparse_csc_const_ref;

  class getfem_object;
  class getfemint_mesh;
  class getfemint_mesh_fem;
  class getfemint_mesh_im;
  class getfemint_mdbrick;
  class getfemint_mdstate;
  class getfemint_model;
  class getfemint_mesh_slice;
  class getfemint_precond;
  class getfemint_gsparse;
  class getfemint_pfem;
  class getfemint_levelset;
  class getfemint_mesh_levelset;
  class getfemint_global_function;
  class getfemint_mesher_object;
  class getfemint_cont_struct;
  class gsparse;

  class sub_index : public gmm::unsorted_sub_index{
  public:
    sub_index() {}
    template <class IT> sub_index(IT b, IT e) : gmm::unsorted_sub_index(b,e) {}
    template <class CONT> sub_index(CONT c) : gmm::unsorted_sub_index(c.begin(),c.end()) {}
    const sub_index &check_range(size_type n) const;
  };

  /* fake full vector/matrix class for matlab/python/etc arrays (necessary for template functions) */
  class array_dimensions {
    unsigned sz;
    unsigned ndim_;
#define ARRAY_DIMENSIONS_MAXDIM 5
    unsigned sizes_[ARRAY_DIMENSIONS_MAXDIM];
  public:
    array_dimensions(): sz(0), ndim_(0) {}
    void push_back(unsigned d) {
      GMM_ASSERT1(ndim_ != ARRAY_DIMENSIONS_MAXDIM-1,
                  " max. nb of dimensions for an output argument exceeded!");
      if (ndim_ == 0) sz = 1;
      sizes_[ndim_++] = d; sz *= d;
    }
    size_type push_back(const array_dimensions& other,
                        unsigned d0, unsigned n,
                        bool matlab_row_matrix_is_a_vector);
    void assign_dimensions(const gfi_array *mx);
    array_dimensions(const gfi_array *mx) {
      assign_dimensions(mx);
    }
    explicit array_dimensions(unsigned sz_) { sz = sz_; ndim_=1; sizes_[0]=sz_; }
    template <typename IVECT> void assign(const IVECT &v) {
      for (unsigned i=0; i < v.size(); ++i) push_back(v[i]);
    }
    unsigned size() const { return sz; }
    unsigned ndim() const { return ndim_; }
    /* dim(0) = first dimension, dim(-1) = last dimension, ..*/
    unsigned dim(int d) const
    { if (d < 0) d = ndim_ + d; return d >= 0 && d < int(ndim_) ? sizes_[d] : 1; }
    unsigned getm() const { return dim(0); }
    unsigned getn() const { return dim(1); }
    unsigned getp() const { return dim(2); }
    unsigned getq() const { return dim(3); }
    const unsigned *sizes() const { return sizes_; }
    void reshape(unsigned n_, unsigned m_, unsigned p_=1);
    void opt_transform_col_vect_into_row_vect();
    friend std::ostream& operator << (std::ostream &os,
                                      const array_dimensions &d) {
      os << d.dim(0);
      for (unsigned i=1; i < d.ndim(); ++i) os << "x" << d.dim(i);
      return os;
    }
  };


  /* fake full vector/matrix class for matlab/python/etc arrays (necessary for template functions) */
  template<typename T> class garray : public array_dimensions {
  public:
    typedef T value_type;
    typedef value_type* iterator;
    typedef const value_type* const_iterator;
  protected:
    dal::shared_array<T> data;
  public:
    value_type& operator[](size_type i)
    { if (i >= size()) THROW_INTERNAL_ERROR; return data[unsigned(i)]; }
    const value_type& operator[](size_type i)
      const { if (i >= size()) THROW_INTERNAL_ERROR; return data[unsigned(i)]; }
    value_type& operator()(size_type i, size_type j, size_type k=0) {
      if (i+j*getm()+k*getm()*getn() >= size()) THROW_INTERNAL_ERROR;
      return data[unsigned(i+j*getm()+k*getm()*getn())];
    }
    const value_type& operator()(size_type i, size_type j, size_type k=0) const {
      if (i+j*getm()+k*getm()*getn() >= size()) THROW_INTERNAL_ERROR;
      return data[unsigned(i+j*getm()+k*getm()*getn())];
    }
    iterator begin() { return data.get(); }
    iterator end() { return data.get()+size(); }
    const_iterator begin() const { return data.get(); }
    const_iterator end() const { return data.get()+size(); }

    /* with this contructor, the data is refcounted */
    garray(value_type *v, int sz_) : array_dimensions(sz_), data(v,true) {}

    /* copies the array vector into a new VECT object */
    template<class VECT> VECT to_vector() const {
      VECT v(begin(), end()); /*v(size());
                                std::copy(begin(), end(), v.begin());*/
      return v;
    }
    garray() {}
    garray(const gfi_array *mx) : array_dimensions(mx) {}
    bool in_range(value_type vmin, value_type vmax) {
      for (size_type i = 0; i < size(); i++)
        if (data[i] < vmin || data[i] > vmax) return false;
      return true;
    }
  };

  class darray : public garray<double> {
  public:
    darray() {}
    /* set data point to the given gfi_array. Freeing data is still the responsability of the gfi_array owner */
    void assign(const gfi_array *mx) {
      if (gfi_array_get_class(mx) == GFI_DOUBLE) {
        /* creation from an array of doubles : just store a ref to the array */
        assign_dimensions(mx);
        data.reset(gfi_double_get_data(mx), false);
      } else if (gfi_array_get_class(mx) == GFI_UINT32 || gfi_array_get_class(mx) == GFI_INT32) {
        /* creation from an array of int : allocation of new storage, and copy of the content */
        assign_dimensions(mx);
        data.reset(new double[size()], true);
        if (gfi_array_get_class(mx) == GFI_INT32)
          std::copy(gfi_int32_get_data(mx), gfi_int32_get_data(mx)+size(), data.get());
        else
          std::copy(gfi_uint32_get_data(mx), gfi_uint32_get_data(mx)+size(), data.get());
      } else THROW_INTERNAL_ERROR;
    }
    darray(const gfi_array *mx) { assign(mx); }
    /* constructor for refcounted arrays */
    darray(value_type *v, int sz_) : garray<double>(v,sz_) {}

    getfem::base_node col_to_bn(unsigned j, unsigned k=0) const {
      getfem::base_node P(getm());
      for (unsigned i=0; i < getm(); i++) P[i] = operator()(i,j,k);
      return P;
    }
    getfem::base_matrix row_col_to_bm(unsigned k=0) const {
      getfem::base_matrix M(getm(),getn());
      for (unsigned i=0; i < getm(); i++)
        for (unsigned j=0; j < getn(); j++) M(i,j) = operator()(i,j,k);
      return M;
    }
  };

  class carray : public garray<complex_type> {
  public:
    carray() {}
    /* set data point to the given gfi_array. Freeing data is still the responsability of the gfi_array owner */
    void assign(const gfi_array *mx) {
      if (gfi_array_get_class(mx) == GFI_DOUBLE && gfi_array_is_complex(mx)) {
        /* creation from an array of complexes : just store a ref to the array */
        assign_dimensions(mx);
        data.reset(reinterpret_cast<complex_type*>(gfi_double_get_data(mx)), false);
      } else if (gfi_array_get_class(mx) == GFI_DOUBLE || gfi_array_get_class(mx) == GFI_UINT32 || gfi_array_get_class(mx) == GFI_INT32) {
        /* creation from an array of int or doubles : allocation of new storage, and copy of the content */
        assign_dimensions(mx);
        data.reset(new complex_type[size()], true);
        if (gfi_array_get_class(mx) == GFI_DOUBLE)
          std::copy(gfi_double_get_data(mx), gfi_double_get_data(mx)+size(), data.get());
        else if (gfi_array_get_class(mx) == GFI_INT32)
          std::copy(gfi_int32_get_data(mx), gfi_int32_get_data(mx)+size(), data.get());
        else if (gfi_array_get_class(mx) == GFI_UINT32)
          std::copy(gfi_uint32_get_data(mx), gfi_uint32_get_data(mx)+size(), data.get());
      } else THROW_INTERNAL_ERROR;
    }
    carray(const gfi_array *mx) { assign(mx); }
    /* constructor for refcounted arrays */
    carray(value_type *v, int sz_) : garray<complex_type>(v,sz_) {}
  };

  /* reference to a real or complex gfi_array */
  class rcarray {
  public:
    typedef enum { REAL, COMPLEX } value_type;
    darray &real() { if (v != REAL) THROW_INTERNAL_ERROR; return *d; }
    carray &cplx() { if (v != COMPLEX) THROW_INTERNAL_ERROR; return *c; }

    rcarray() : mx(0), d(0), c(0), v(REAL) {}
    rcarray(const gfi_array *mx_, int v_ = -1) : d(0), c(0) { assign(mx_,v_); }
    ~rcarray() { clear(); }
    void assign(const gfi_array *mx_, int v_ = -1) {
      mx = mx_;
      v = v_; if (v == -1) v = gfi_array_is_complex(mx) ? COMPLEX : REAL;
      clear();
      if (v == REAL)
        d.reset(new darray(mx));
      else c.reset(new carray(mx));
    }
    bool is_complex() const { return v == COMPLEX; }
    carray &to_complex() {
      if (v == REAL)
        { c.reset(new carray(mx)); d.reset(0); v = COMPLEX; }
      return cplx();
    }
    darray &to_real() {
      if (v == COMPLEX)
        THROW_BADARG("expected a real (not a complex) array in this context");
      return real();
    }
    void clear() { d.reset(0); c.reset(0); }
    array_dimensions& sizes() { if (d.get()) return *d; else return *c; }
    const array_dimensions& sizes() const { if (d.get()) return *d; else return *c; }
  private:
    const gfi_array *mx;
    dal::shared_ptr<darray> d;
    dal::shared_ptr<carray> c;
    int v;
  };


  /* fake full vector/matrix class for matlab/python/etc arrays (necessary for template functions) */
  class iarray : public garray<int> {
  public:
    iarray() {}
    void assign(const gfi_array *mx) {
      if (gfi_array_get_class(mx) == GFI_INT32)
        data.reset(gfi_int32_get_data(mx), false);
      else if (gfi_array_get_class(mx) == GFI_UINT32)
        data.reset((int*)gfi_uint32_get_data(mx), false);
      else THROW_INTERNAL_ERROR;
      assign_dimensions(mx);
    }
    iarray(const gfi_array *mx) { assign(mx); }
    /* constructor for refcounted arrays */
    iarray(value_type *v, int sz_) : garray<int>(v,sz_) {}
  };
}

/* gmm interface for darrays and carrays */
namespace gmm {
  template<typename T> struct linalg_traits<getfemint::garray<T> > {
    typedef getfemint::garray<T> this_type;
    typedef linalg_modifiable is_reference; /* is a writable reference */
    typedef abstract_vector linalg_type;
    typedef T value_type;
    typedef value_type origin_type;
    typedef value_type& reference;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef gmm::abstract_dense storage_type;
    typedef linalg_true index_sorted;
    static size_type size(const this_type &v) { return v.size(); }
    static iterator begin(this_type &v) { return v.begin(); }
    static const_iterator begin(const this_type &v) { return v.begin(); }
    static iterator end(this_type &v) { return v.end(); }
    static const_iterator end(const this_type &v) { return v.end(); }
    static const origin_type* origin(const this_type &v) { return v.begin(); }
    static origin_type* origin(this_type &v) { return v.begin(); }
    static void clear(origin_type* o, const iterator &it, const iterator &ite)
    { std::fill(it, ite, value_type(0)); }
    static void do_clear(this_type &v) { std::fill(v.begin(), v.end(), 0.); }
    static value_type access(const origin_type *, const const_iterator &it,
                             const const_iterator &, size_type i)
    { return it[i]; }
    static reference access(origin_type *, const iterator &it,
                            const iterator &, size_type i)
    { return it[i]; }
  };

  template<> struct linalg_traits<getfemint::darray> :
    public linalg_traits<getfemint::garray<double> > {};
  template<> struct linalg_traits<getfemint::carray> :
    public linalg_traits<getfemint::garray<getfemint::complex_type> > {};

  template <typename T> struct temporary_dense_vector<getfemint::garray<T> > {
    typedef std::vector<T> vector_type;
  };
} /* fin de l'intermède */


namespace getfemint {
  //bool is_static_object(id_type id, id_type cid);
  template<typename T> struct select_int_double_or_complex {
    typedef int value_type;
  };
  template<> struct select_int_double_or_complex<double> {
    typedef double value_type;
  };
  template<> struct select_int_double_or_complex<float> {
    typedef double value_type;
  };
  template<typename V> struct select_int_double_or_complex<std::complex<V> > {
    typedef std::complex<typename select_int_double_or_complex<V>::value_type> value_type;
  };

  /* input argument of the gf_* mex-files */
  class mexarg_in {
    //void check_int_values(int vmin=INT_MIN, int vmax=INT_MAX);
    double to_scalar_(bool isint);
    void error_if_nonwritable(getfem_object *o, bool want_writeable);
  public:
    const gfi_array *arg;
    int argnum;

    mexarg_in(const gfi_array *arg_, int num_) { arg = arg_; argnum = num_; }
    bool                                 is_string() { return (gfi_array_get_class(arg) == GFI_CHAR); }
    bool                                 is_cell() { return (gfi_array_get_class(arg) == GFI_CELL); }
    bool                                 is_object_id(id_type *pid=0, id_type *pcid=0);
    bool                                 is_mesh();
    bool                                 is_mesh_fem();
    bool                                 is_mesh_im();
    bool                                 is_mdbrick();
    bool                                 is_mdstate();
    bool                                 is_model();
    bool                                 is_mesh_slice();
    bool                                 is_levelset();
    bool                                 is_mesh_levelset();
    bool                                 is_global_function();
    bool                                 is_mesher_object();
    bool                                 is_cont_struct();
    bool                                 is_sparse() { return (gfi_array_get_class(arg) == GFI_SPARSE || is_gsparse()); };
    bool                                 is_gsparse();
    bool                                 is_complex(); /* true for complex garrays AND complex sparse matrices (native or gsparse) */
    bool                                 is_integer();
    bool                                 is_bool();
    bool                                 to_bool();
    int                                  to_integer(int min_val=INT_MIN, int max_val=INT_MAX);
    size_type                            to_convex_number(const getfem::mesh &m);
    size_type                            to_face_number(size_type nbf);
    double                               to_scalar(double min_val=-1e300, double max_val=1e300);
    complex_type                         to_scalar(complex_type);
    std::string                          to_string();
    id_type                              to_object_id(id_type *pid=0, id_type *pcid=0);
    bgeot::base_poly *                   to_poly();
    const getfem::mesh_fem *             to_const_mesh_fem();
    getfem::mesh_fem *                   to_mesh_fem();
    const getfem::mesh_im *              to_const_mesh_im();
    getfem::mesh_im *                    to_mesh_im();
    const getfem::mesh *                 to_const_mesh();
    const getfem::mesh *                 to_const_mesh(id_type& id);
    getfemint_mesh *                     to_getfemint_mesh(bool writeable=false);
    getfemint_mesh_fem *                 to_getfemint_mesh_fem(bool writeable=false);
    getfemint_mesh_im *                  to_getfemint_mesh_im(bool writeable=false);
    getfemint_mdbrick *                  to_getfemint_mdbrick(bool writeable=false);
    getfemint_mdstate *                  to_getfemint_mdstate(bool writeable=false);
    getfemint_model *                    to_getfemint_model(bool writeable=false);
    getfem::mesh *                       to_mesh();
    getfemint_mesh_slice *               to_getfemint_mesh_slice(bool writeable=false);
    getfemint_levelset *                 to_getfemint_levelset(bool writeable=false);
    // const getfem::level_set *            to_const_levelset();
    getfem::level_set *                  to_levelset();
    getfemint_mesh_levelset *            to_getfemint_mesh_levelset(bool writeable=false);
    getfem::mesh_level_set *             to_mesh_levelset();
    const getfem::abstract_xy_function * to_const_global_function();
    const getfem::mesher_signed_distance * to_const_mesher_object();
    getfem::abstract_xy_function *       to_global_function();
    getfem::mesher_signed_distance *     to_mesher_object();
    const getfem::cont_struct_getfem_model * to_const_cont_struct();
    getfem::cont_struct_getfem_model *   to_cont_struct();
    getfemint_global_function *          to_getfemint_global_function(bool writeable=false);
    getfemint_mesher_object *            to_getfemint_mesher_object(bool writeable=false);
    getfemint_cont_struct *              to_getfemint_cont_struct(bool writable=false);
    getfem::pintegration_method          to_integration_method();
    getfemint_pfem*                      to_getfemint_pfem();
    getfem::pfem                         to_fem();
    getfem::pmat_elem_type               to_mat_elem_type();
    bgeot::pgeometric_trans              to_pgt();
    bgeot::pconvex_structure             to_convex_structure();
    getfemint_precond *                  to_precond();
    getfem::mesh_region                  to_mesh_region();

    carray to_carray();
    carray to_carray(int expected_dim);
    carray to_carray(int expected_n, int expected_m,
                     int expected_k=1, int expected_q=1);

    /* do not perform any check on the number of dimensions of the array */
    darray to_darray();

    /* expect the argument to be a row or column vector of given dimension */
    darray to_darray(int expected_dim);

    /* expect the argument to be a matrix (or possibly a 3D array)
       if any of the arguments has a value of -1, the corresponding dimension
       is not checked
     */
    darray to_darray(int expected_m, int expected_n,
                                          int expected_k=1, int expected_q=1);

    /* convertion to a real or complex array */
    rcarray to_rcarray();
    rcarray to_rcarray(int expected_dim);
    rcarray to_rcarray(int expected_m, int expected_n,
                                           int expected_k=1, int expected_q=1);

    iarray to_iarray();
    iarray to_iarray(int expected_dim);
    iarray to_iarray(int expected_m, int expected_n,
                                          int expected_k=1, int expected_q=1);

    /* template friendly version */
    garray<double>            to_garray(double) { return to_darray(); }
    garray<double>            to_garray(int expected_dim, double) { return to_darray(expected_dim); }
    garray<double>            to_garray(int expected_m, int expected_n, double) { return to_darray(expected_m, expected_n); }
    garray<complex_type>      to_garray(complex_type) { return to_carray(); }
    garray<complex_type>      to_garray(int expected_dim, complex_type) { return to_carray(expected_dim); }
    garray<complex_type>      to_garray(int expected_m, int expected_n, complex_type) { return to_carray(expected_m, expected_n); }

    dal::bit_vector           to_bit_vector(const dal::bit_vector *subsetof = NULL, int shiftvals=-config::base_index());
    sub_index                 to_sub_index();
    getfem::base_node         to_base_node() { return to_base_node(-1); }
    getfem::base_node         to_base_node(int expected_dim);
    void                      to_sparse(gf_real_sparse_csc_const_ref& M);
    void                      to_sparse(gf_cplx_sparse_csc_const_ref& M);
    dal::shared_ptr<gsparse>  to_sparse();
    getfemint_gsparse *       to_getfemint_gsparse();

    mexarg_in &check_trailing_dimension(int expected_dim);
    void check_dimensions(array_dimensions &v, int expected_m, int expected_n, int expected_p=-1, int expected_q=-1);
    void check_dimensions(const array_dimensions &v, int expected_dim);
private:
    friend class mexargs_in;
    mexarg_in() {}
  };

  /* output argument of the gf_* mex-files */
  class mexargs_out;
  class mexarg_out {
  public:
    gfi_array *& arg;
    int argnum;
    mexarg_out(gfi_array * &arg_, int num_) : arg(arg_), argnum(num_) {}

    void from_object_id(id_type id, id_type class_id);
    void from_object_id(std::vector<id_type> ids, id_type class_id);
    void from_integer(int i);
    void from_scalar(double v);
    void from_string(const char *s);
    template<class STR_CONT> void from_string_container(const STR_CONT& s);
    void from_bit_vector(const dal::bit_vector& bv, int shift=config::base_index());
    void from_mesh_region(const getfem::mesh_region &region);
    typedef enum { USE_NATIVE_SPARSE, USE_GSPARSE, USE_DEFAULT_SPARSE } output_sparse_fmt;


    /**
       BIG CAUTION:
       the from_sparse function MAY erase the sparse matrix M (in order to avoid unnecessary copies of sparse matrices)
    */
    void from_sparse(gf_real_sparse_by_col& M, output_sparse_fmt fmt = USE_DEFAULT_SPARSE);
    void from_sparse(gf_cplx_sparse_by_col& M, output_sparse_fmt fmt = USE_DEFAULT_SPARSE);
    void from_sparse(gsparse& M, output_sparse_fmt fmt = USE_DEFAULT_SPARSE);

    void from_tensor(const getfem::base_tensor& t);
    carray create_carray_v(unsigned dim);
    carray create_carray_h(unsigned dim);
    carray create_carray(unsigned n, unsigned m);
    carray create_carray(unsigned n, unsigned m, unsigned p);
    darray create_darray_v(unsigned dim);
    darray create_darray_h(unsigned dim);
    darray create_darray(unsigned n, unsigned m);
    darray create_darray(unsigned n, unsigned m, unsigned p);
    iarray create_iarray_v(unsigned dim);
    iarray create_iarray_h(unsigned dim);
    iarray create_iarray(unsigned n, unsigned m);
    iarray create_iarray(unsigned n, unsigned m, unsigned p);
    /* overloaded functions, useful for templates */
    darray create_array_v(unsigned n, double) { return create_darray_v(n); }
    darray create_array_h(unsigned n, double) { return create_darray_h(n); }
    darray create_array(unsigned n, unsigned m, double) { return create_darray(n,m); }
    darray create_array(unsigned n, unsigned m, unsigned p, double) { return create_darray(n,m,p); }
    carray create_array_v(unsigned n, complex_type) { return create_carray_v(n); }
    carray create_array_h(unsigned n, complex_type) { return create_carray_h(n); }
    carray create_array(unsigned n, unsigned m, complex_type) { return create_carray(n,m); }
    carray create_array(unsigned n, unsigned m, unsigned p, complex_type) { return create_carray(n,m,p); }
    darray create_array(const array_dimensions &dims, double);
    carray create_array(const array_dimensions &dims, complex_type);
    /* allocates a gsparse object in the workspace */
    gsparse &create_gsparse();

    template<class VECT> void from_dcvector(VECT& v) {
      typedef typename VECT::value_type T;
      create_array_h(unsigned(v.size()), T());
      std::copy(v.begin(), v.end(), gfi_get_data(arg, T()));
    }
    template<class VECT> void from_ivector(VECT& v) {
      create_iarray_h(unsigned(v.size()));
      std::copy(v.begin(), v.end(), gfi_int32_get_data(arg));
    }
    template<class VEC_CONT> void from_vector_container(const VEC_CONT& vv);
  };

  template<class STR_CONT> void
  mexarg_out::from_string_container(const STR_CONT& s)
  {
    arg = checked_gfi_array_create_2(int(s.size()), 1, GFI_CELL);
    gfi_array **c = gfi_cell_get_data(arg);
    typename STR_CONT::const_iterator it;
    size_type cnt = 0;
    for (it = s.begin(); it != s.end(); ++it, ++cnt) {
      c[cnt] = checked_gfi_array_from_string(it->c_str());
    }
  }

  /* output a matrix from a container of vectors (of same dimensions) */
  template<class VEC_CONT> void
  mexarg_out::from_vector_container(const VEC_CONT& vv)
  {
    size_type n = vv.size();
    size_type m = (n == 0) ? 0 : vv[0].size();
    darray w  = create_darray(unsigned(m),unsigned(n));
    for (size_type j=0; j < n; ++j)
      std::copy(vv[j].begin(), vv[j].end(), &w(0,j));
  }

  gfi_array *create_object_id(int nid, id_type *ids, id_type cid, bool not_as_a_vector=false);
  inline gfi_array *create_object_id(id_type id, id_type cid) {
    return create_object_id(1, &id, cid, true);
  }

  /* handles the list of input arguments */
  class mexargs_in {
    const gfi_array **in;
    dal::bit_vector idx;
    int nb_arg;
    bool use_cell;
    mexarg_in last;
    /* copy forbidden */
    mexargs_in(const mexargs_in& );
    mexargs_in& operator=(const mexargs_in& );
  public:
    mexargs_in(int n, const gfi_array *p[], bool use_cell);
    ~mexargs_in();
    void check() const { if (idx.card() == 0) THROW_INTERNAL_ERROR; }
    const gfi_array *pop_gfi_array(size_type decal=0, int *out_idx = NULL) {
      size_type i = idx.first_true();
      check();
      if (decal >= idx.card()) THROW_INTERNAL_ERROR;
      while (decal>0) { i++; check(); if (idx.is_in(i)) decal--; }
      idx.sup(i);
      if (out_idx) *out_idx = int(i);
      return in[i];
    }
    mexarg_in &pop(size_type decal=0) {
      int i;
      const gfi_array *m = pop_gfi_array(decal, &i);
      last = mexarg_in(m,i+1);
      return last;
    }
    mexarg_in &last_popped() { return last; }
    void restore(size_type i) { idx.add(i); }
    mexarg_in front() const {
      check();
      return mexarg_in(in[idx.first_true()],
                       int(idx.first_true()));
    }
    int narg() const { return nb_arg; }
    int remaining() const { return int(idx.card()); }
    void get_array(const gfi_array **& m) {
      if (remaining()) {
        m = new const gfi_array *[remaining()];
        for (size_type i=0; remaining(); i++) {
          m[i] = pop_gfi_array();
        }
      } else m = NULL;
    }
    void pop_all() { idx.clear(); }
  };

  /* handles the list of output arguments */
  class mexargs_out {
    mutable std::deque<gfi_array *> out; /* deque because mexarg_out hold a reference to this array content */
    int nb_arg; /* if equal to -1, the number of output arguments is unknown */
    int idx;
    int okay; /* if 0, the destructor will destroy the allacted arrays in 'out'
                 and will call workspace().destroy_newly_created_objects */
    bool scilab_flag;
    /* copy forbidden */
    mexargs_out(const mexargs_out& );
    mexargs_out& operator=(const mexargs_out& );
  public:
    mexargs_out(int nb_arg_);
    ~mexargs_out();
    void check() const;
    mexarg_out pop();
    mexarg_out front() const { check(); return mexarg_out(out[idx], idx+1); }
    bool narg_in_range(int min, int max) const {
      if ((scilab_flag) && (max==0) && (min==0)) return true;
      return (nb_arg == -1 || (nb_arg >= min && (nb_arg <= max || max == -1)));
    }
    bool narg_known() const { return nb_arg != -1; }
    bool remaining() const { return !narg_known() || (std::max(nb_arg,1) - idx); }
    void return_packed_obj_ids(const std::vector<id_type>& ids, id_type class_id);
    std::deque<gfi_array *>& args() { return out; }
    void set_okay(bool ok) { okay = ok; }
    void set_scilab(bool _scilab_flag) {scilab_flag = _scilab_flag;}
    bool get_scilab() const {return scilab_flag;}
  };

  std::string cmd_normalize(const std::string& a);
  bool cmd_strmatch(const std::string& a, const char *s);
  bool cmd_strmatchn(const std::string& a, const char *s, unsigned n);
  bool check_cmd(const std::string& cmdname, const char *s,
                 const mexargs_in& in,
                 int min_argin=0, int max_argin=-1);
  bool check_cmd(const std::string& cmdname, const char *s,
                 const mexargs_out& out,
                 int min_argout=0, int max_argout=-1);
  bool check_cmd(const std::string& cmdname, const char *s,
                 const mexargs_in& in, const mexargs_out& out,
                 int min_argin=0, int max_argin=-1,
                 int min_argout=0, int max_argout=-1);

  inline void bad_cmd(std::string& cmd) {
    THROW_BADARG("Bad command name: " << cmd); }

  void check_cv_fem(const getfem::mesh_fem& mf, size_type cv);
  void check_cv_im(const getfem::mesh_im& mim, size_type cv);

  const double& get_NaN();
  bool is_NaN(const double&);
}  /* end of namespace getfemint.                                          */

#endif /* GETFEMINT_H__                                                    */
