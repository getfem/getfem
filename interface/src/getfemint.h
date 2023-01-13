/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2001-2020 Y. Renard, J. Pommier.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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
#include <getfem/dal_static_stored_objects.h>
#include <getfem/dal_bit_vector.h>
#include <getfem/getfem_config.h>
#include <getfem/getfem_mesh.h>
#include <gfi_array.h>

namespace getfemint {

  /* This is very important that the following classes respects the
     alphabetic order. The order has to be the same as in getfem.py !!!
     Do not forget to modify also 'name_of_getfemint_class_id' just after
  */
  typedef enum { CONT_STRUCT_CLASS_ID,
                 CVSTRUCT_CLASS_ID,
                 ELTM_CLASS_ID,
                 FEM_CLASS_ID,
                 GEOTRANS_CLASS_ID,
                 GLOBAL_FUNCTION_CLASS_ID,
                 INTEG_CLASS_ID,
                 LEVELSET_CLASS_ID,
                 MESH_CLASS_ID,
                 MESHFEM_CLASS_ID,
                 MESHIM_CLASS_ID,
                 MESHIMDATA_CLASS_ID,
                 MESH_LEVELSET_CLASS_ID,
                 MESHER_OBJECT_CLASS_ID,
                 MODEL_CLASS_ID,
                 PRECOND_CLASS_ID,
                 SLICE_CLASS_ID,
                 SPMAT_CLASS_ID,
                 POLY_CLASS_ID,    /* Not fully interfaced. Remain at the end */
                 GETFEMINT_NB_CLASS } getfemint_class_id;

  
  /* Associate the class ID found in the matlab structures referencing
     getfem object to a class name which coincides with the class name
     given by matlab to the structure.
     
     IMPORTANT: Should correspond to the getfemint_class_id
                In particular, it should be in alphabetic order.
  */
  const char *name_of_getfemint_class_id(id_type cid);

  /* exception-throwing version of the allocation functions of gfi_array.h */
  gfi_array* checked_gfi_array_create(int ndim, const int *dims,
                                      gfi_type_id type,
                                      gfi_complex_flag is_complex = GFI_REAL);
  gfi_array* checked_gfi_array_create_0(gfi_type_id type,
                                        gfi_complex_flag is_complex = GFI_REAL);
  gfi_array* checked_gfi_array_create_1(int M, gfi_type_id type,
                                        gfi_complex_flag is_complex = GFI_REAL);
  gfi_array* checked_gfi_array_create_2(int M, int N, gfi_type_id type,
                                        gfi_complex_flag is_complex = GFI_REAL);
  gfi_array* checked_gfi_array_from_string(const char*s);
  gfi_array* checked_gfi_create_sparse(int m, int n, int nzmax,
                                       gfi_complex_flag is_complex);

  typedef bgeot::dim_type dim_type;
  typedef bgeot::scalar_type scalar_type;
  typedef std::complex<double> complex_type;

  inline int *gfi_get_data(gfi_array *g, int)
  { return gfi_int32_get_data(g); }
  inline unsigned *gfi_get_data(gfi_array *g, unsigned)
  { return gfi_uint32_get_data(g); }
  inline char *gfi_get_data(gfi_array *g, char)
  { return gfi_char_get_data(g); }
  inline double *gfi_get_data(gfi_array *g, double)
  { return gfi_double_get_data(g); }
  inline complex_type *gfi_get_data(gfi_array *g, complex_type)
  { return (complex_type*)gfi_double_get_data(g); }


  typedef gmm::row_matrix<gmm::wsvector<scalar_type> >  gf_real_sparse_by_row;
  typedef gmm::col_matrix<gmm::wsvector<scalar_type> >  gf_real_sparse_by_col;
  typedef gmm::csc_matrix_ref<const double *, const unsigned int *,
                              const unsigned int *>
  gf_real_sparse_csc_const_ref;
  typedef gmm::row_matrix<gmm::wsvector<complex_type> >  gf_cplx_sparse_by_row;
  typedef gmm::col_matrix<gmm::wsvector<complex_type> >  gf_cplx_sparse_by_col;
  typedef gmm::csc_matrix_ref<const complex_type *, const unsigned int *,
                              const unsigned int *>
  gf_cplx_sparse_csc_const_ref;
  class gsparse;

  class sub_index : public gmm::unsorted_sub_index {
  public:
    sub_index() {}
    template <class IT> sub_index(IT b, IT e)
      : gmm::unsorted_sub_index(b,e) {}
    template <class CONT> sub_index(CONT c)
      : gmm::unsorted_sub_index(c.begin(),c.end()) {}
    const sub_index &check_range(size_type n) const;
  };

  /* fake full vector/matrix class for matlab/python/etc arrays
     (necessary for template functions) */
# define ARRAY_DIMENSIONS_MAXDIM 5

  class array_dimensions {
    unsigned sz;
    unsigned ndim_;
    unsigned sizes_[ARRAY_DIMENSIONS_MAXDIM];
  public:
    array_dimensions(): sz(0), ndim_(0) { sizes_[0]=0; sizes_[1]=0; }
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
    explicit array_dimensions(unsigned sz_)
    { sz = sz_; ndim_=1; sizes_[0]=sz_; }
    template <typename IVECT> void assign(const IVECT &v) {
      for (unsigned i=0; i < v.size(); ++i) push_back(unsigned(v[i]));
    }
    unsigned size() const { return sz; }
    unsigned ndim() const { return ndim_; }
    /* dim(0) = first dimension, dim(-1) = last dimension, ..*/
    unsigned dim(int d) const {
      if (d < 0) d = ndim_ + d;
      return d >= 0 && d < int(ndim_) ? sizes_[d] : 1;
    }
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


  /* fake full vector/matrix class for matlab/python/etc arrays
     (necessary for template functions) */
  template<typename T> class garray : public array_dimensions {
  public:
    typedef T value_type;
    typedef value_type* iterator;
    typedef const value_type* const_iterator;
  protected:
    std::shared_array_ptr<T> data; // Should be replaced by std::shared_ptr<T[]>
                                   // when it will be suported.
  public:
    value_type& operator[](size_type i)
    { if (i >= size()) THROW_INTERNAL_ERROR; return (data.get())[unsigned(i)]; }
    const value_type& operator[](size_type i) const
    { if (i >= size()) THROW_INTERNAL_ERROR; return (data.get())[unsigned(i)]; }
    value_type& operator()(size_type i, size_type j, size_type k=0) {
      if (i+j*getm()+k*getm()*getn() >= size()) THROW_INTERNAL_ERROR;
      return (data.get())[unsigned(i+j*getm()+k*getm()*getn())];
    }
    const value_type& operator()(size_type i, size_type j,
                                 size_type k=0) const {
      if (i+j*getm()+k*getm()*getn() >= size()) THROW_INTERNAL_ERROR;
      return (data.get())[unsigned(i+j*getm()+k*getm()*getn())];
    }
    iterator begin() { return data.get(); }
    iterator end() { return data.get()+size(); }
    const_iterator begin() const { return data.get(); }
    const_iterator end() const { return data.get()+size(); }

    /* with this contructor, the data is refcounted */
    garray(value_type *v, int sz_) : array_dimensions(sz_), data(v) {}

    /* copies the array vector into a new VECT object */
    template<class VECT> VECT to_vector()
      const { VECT v(begin(), end()); return v; }
    garray() {}
    garray(const gfi_array *mx) : array_dimensions(mx) {}
    bool in_range(value_type vmin, value_type vmax) {
      for (size_type i = 0; i < size(); i++)
        if ((data.get())[i] < vmin || (data.get())[i] > vmax) return false;
      return true;
    }
  };

  class darray : public garray<double> {
  public:
    darray() {}
    /* set data point to the given gfi_array.
       Freeing data is still the responsability of the gfi_array owner */
    void assign(const gfi_array *mx) {
      if (gfi_array_get_class(mx) == GFI_DOUBLE) {
        /* creation from an array of doubles : just store a ref to the array */
        assign_dimensions(mx);
        data = std::shared_array_ptr<double>(std::shared_ptr<double>(),
                                             gfi_double_get_data(mx));
      } else if (gfi_array_get_class(mx) == GFI_UINT32 ||
                 gfi_array_get_class(mx) == GFI_INT32) {
        /* creation from an array of int : allocation of new storage,
           and copy of the content */
        assign_dimensions(mx);
        data = std::make_shared_array<double>(size());
        if (gfi_array_get_class(mx) == GFI_INT32)
          std::copy(gfi_int32_get_data(mx), gfi_int32_get_data(mx)+size(),
                    data.get());
        else
          std::copy(gfi_uint32_get_data(mx), gfi_uint32_get_data(mx)+size(),
                    data.get());
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
    /* set data point to the given gfi_array.
       Freeing data is still the responsability of the gfi_array owner */
    void assign(const gfi_array *mx) {
      if (gfi_array_get_class(mx) == GFI_DOUBLE && gfi_array_is_complex(mx)) {
        /* creation from an array of complexes: just store a ref to the array */
        assign_dimensions(mx);
        data = std::shared_array_ptr<complex_type>
               (std::shared_ptr<complex_type>(),
                reinterpret_cast<complex_type*>(gfi_double_get_data(mx)));
      } else if (gfi_array_get_class(mx) == GFI_DOUBLE ||
                 gfi_array_get_class(mx) == GFI_UINT32 ||
                 gfi_array_get_class(mx) == GFI_INT32) {
        /* creation from an array of int or doubles: allocation of new storage,
           and copy of the content */
        assign_dimensions(mx);
        data = std::make_shared_array<complex_type>(size());
        if (gfi_array_get_class(mx) == GFI_DOUBLE)
          std::copy(gfi_double_get_data(mx), gfi_double_get_data(mx)+size(),
                    data.get());
        else if (gfi_array_get_class(mx) == GFI_INT32)
          std::copy(gfi_int32_get_data(mx), gfi_int32_get_data(mx)+size(),
                    data.get());
        else if (gfi_array_get_class(mx) == GFI_UINT32)
          std::copy(gfi_uint32_get_data(mx), gfi_uint32_get_data(mx)+size(),
                    data.get());
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
      if (v == REAL) d = std::make_shared<darray>(mx);
      else c = std::make_shared<carray>(mx);
    }
    bool is_complex() const { return v == COMPLEX; }
    carray &to_complex() {
      if (v == REAL) {
        c = std::make_shared<carray>(mx);
        d = std::shared_ptr<darray>();
        v = COMPLEX;
      }
      return cplx();
    }
    darray &to_real() {
      if (v == COMPLEX)
        THROW_BADARG("expected a real (not a complex) array in this context");
      return real();
    }
    void clear() { c=std::shared_ptr<carray>(); d=std::shared_ptr<darray>(); }
    array_dimensions& sizes() { if (d.get()) return *d; else return *c; }
    const array_dimensions& sizes() const
    { if (d.get()) return *d; else return *c; }
  private:
    const gfi_array *mx;
    std::shared_ptr<darray> d;
    std::shared_ptr<carray> c;
    int v;
  };


  /* fake full vector/matrix class for matlab/python/etc arrays
     (necessary for template functions) */
  class iarray : public garray<int> {
  public:
    iarray() {}
    void assign(const gfi_array *mx) {
      if (gfi_array_get_class(mx) == GFI_INT32)
        data = std::shared_array_ptr<int>(std::shared_ptr<int>(),
                                          gfi_int32_get_data(mx));
      else if (gfi_array_get_class(mx) == GFI_UINT32)
        data = std::shared_array_ptr<int>(std::shared_ptr<int>(),
                                          (int*)gfi_uint32_get_data(mx));
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
}


namespace getfemint {
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
    typedef std::complex<typename select_int_double_or_complex<V>::value_type>
    value_type;
  };

  /* input argument of the gf_* files */
  class mexarg_in {
    double to_scalar_(bool isint);
  public:
    const gfi_array *arg;
    int argnum;

    mexarg_in(const gfi_array *arg_, int num_) { arg = arg_; argnum = num_; }
    bool is_string() { return (gfi_array_get_class(arg) == GFI_CHAR); }
    bool is_cell() { return (gfi_array_get_class(arg) == GFI_CELL); }
    bool is_object_id(id_type *pid=0, id_type *pcid=0) const;
    bool is_sparse();
    bool is_complex(); /* true for complex garrays AND complex sparse matrices
                          (native or gsparse) */
    bool is_integer();
    bool is_bool();


    bool                 to_bool();
    int                  to_integer(int min_val=INT_MIN, int max_val=INT_MAX);
    size_type            to_convex_number(const getfem::mesh &m);
    short_type           to_face_number(short_type nbf);
    double               to_scalar(double min_val=-1e300, double max_val=1e300);
    complex_type         to_scalar(complex_type);
    std::string          to_string();
    id_type              to_object_id(id_type *pid=0, id_type *pcid=0);
    getfem::mesh_region  to_mesh_region();
    carray               to_carray();
    carray               to_carray(int expected_dim);
    carray               to_carray(int expected_n, int expected_m,
                                   int expected_k=1, int expected_q=1);
    darray               to_darray(); /* do not perform any check on the
                                         number of dimensions of the array */
    /* expect the argument to be a row or column vector of given dimension */
    darray               to_darray(int expected_dim);
    /* expect the argument to be a matrix (or possibly a 3D array)
       if any of the arguments has a value of -1, the corresponding dimension
       is not checked
     */
    darray               to_darray(int expected_m, int expected_n,
                                   int expected_k=1, int expected_q=1);
    /* convertion to a real or complex array */
    rcarray              to_rcarray();
    rcarray              to_rcarray(int expected_dim);
    rcarray              to_rcarray(int expected_m, int expected_n,
                                    int expected_k=1, int expected_q=1);
    iarray               to_iarray();
    iarray               to_iarray(int expected_dim);
    iarray               to_iarray(int expected_m, int expected_n,
                                   int expected_k=1, int expected_q=1);

    /* template friendly version */
    garray<double>            to_garray(double) { return to_darray(); }
    garray<double>            to_garray(int expected_dim, double)
    { return to_darray(expected_dim); }
    garray<double>            to_garray(int expected_m, int expected_n, double)
    { return to_darray(expected_m, expected_n); }
    garray<complex_type>      to_garray(complex_type) { return to_carray(); }
    garray<complex_type>      to_garray(int expected_dim, complex_type)
    { return to_carray(expected_dim); }
    garray<complex_type>      to_garray(int expected_m, int expected_n,
                                        complex_type)
    { return to_carray(expected_m, expected_n); }
    dal::bit_vector           to_bit_vector(const dal::bit_vector *subsetof = 0,
                                            int shiftvs=-config::base_index());
    sub_index                 to_sub_index();
    getfem::base_node         to_base_node() { return to_base_node(-1); }
    getfem::base_node         to_base_node(int expected_dim);
    void                      to_sparse(gf_real_sparse_csc_const_ref& M);
    void                      to_sparse(gf_cplx_sparse_csc_const_ref& M);
    std::shared_ptr<gsparse>  to_sparse();
    
    mexarg_in &check_trailing_dimension(int expected_dim);
    void check_dimensions(array_dimensions &v, int expected_m,
                          int expected_n, int expected_p=-1, int expected_q=-1);
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
    void from_bit_vector(const dal::bit_vector& bv,
                         int shift=config::base_index());
    void from_mesh_region(const getfem::mesh_region &region);
    typedef enum { USE_NATIVE_SPARSE, USE_GSPARSE,
                   USE_DEFAULT_SPARSE } output_sparse_fmt;

    /**
       BIG CAUTION:
       the from_sparse function MAY erase the sparse matrix M
       (in order to avoid unnecessary copies of sparse matrices)
    */
    void from_sparse(gf_real_sparse_by_col& M,
                     output_sparse_fmt fmt = USE_DEFAULT_SPARSE);
    void from_sparse(gf_cplx_sparse_by_col& M,
                     output_sparse_fmt fmt = USE_DEFAULT_SPARSE);
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
    darray create_array(unsigned n, unsigned m, double)
    { return create_darray(n,m); }
    darray create_array(unsigned n, unsigned m, unsigned p, double)
    { return create_darray(n,m,p); }
    carray create_array_v(unsigned n, complex_type)
    { return create_carray_v(n); }
    carray create_array_h(unsigned n, complex_type)
    { return create_carray_h(n); }
    carray create_array(unsigned n, unsigned m, complex_type)
    { return create_carray(n,m); }
    carray create_array(unsigned n, unsigned m, unsigned p, complex_type)
    { return create_carray(n,m,p); }
    darray create_array(const array_dimensions &dims, double);
    carray create_array(const array_dimensions &dims, complex_type);

    template<class VECT> void from_dcvector(VECT& v) {
      typedef typename VECT::value_type T;
      create_array_h(unsigned(v.size()), T());
      std::copy(v.begin(), v.end(), gfi_get_data(arg, T()));
    }
    template<class VECT> void from_dlvector(VECT& v) {
      typedef typename VECT::value_type T;
      create_array_v(unsigned(v.size()), T());
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

  gfi_array *create_object_id(int nid, id_type *ids, id_type cid,
                              bool not_as_a_vector=false);
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
    void pop_all() { idx.clear(); }
  };

  /* handles the list of output arguments */
  class mexargs_out {
    mutable std::deque<gfi_array *> out; /* deque because mexarg_out hold a
                                            reference to this array content */
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
    bool remaining() const
    { return !narg_known() || (std::max(nb_arg,1) - idx); }
    void return_packed_obj_ids(const std::vector<id_type>& ids,
                               id_type class_id);
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


  // Gives the class id of an object
  // To be completed when an object class is added.
  id_type class_id_of_object(const dal::pstatic_stored_object &p,
                             const void **q = 0);

# define getfemint_declare_getfem_class(CLASS)                          \
  } namespace getfem { class CLASS; } namespace getfemint {
# define getfemint_delare_bgeot_class(CLASS)                            \
  } namespace bgeot { class CLASS; } namespace getfemint {
  
  // Functions for CONT_STRUCT_CLASS_ID
  getfemint_declare_getfem_class(cont_struct_getfem_model)
  bool is_cont_struct_object(const mexarg_in &p);
  id_type store_cont_struct_object
  (const std::shared_ptr<getfem::cont_struct_getfem_model> &shp);
  getfem::cont_struct_getfem_model *to_cont_struct_object(const mexarg_in &p);

  // Functions for CVSTRUCT_CLASS_ID
  getfemint_delare_bgeot_class(convex_structure)
  typedef std::shared_ptr<const bgeot::convex_structure> pconvex_structure;
  bool is_cvstruct_object(const mexarg_in &p);
  id_type store_cvstruct_object(const pconvex_structure &shp);
  pconvex_structure to_cvstruct_object(const mexarg_in &p);

  // Functions for ELTM_CLASS_ID
  getfemint_declare_getfem_class(mat_elem_type)
  typedef std::shared_ptr<const getfem::mat_elem_type> pmat_elem_type;
  bool is_eltm_object(const mexarg_in &p);
  id_type store_eltm_object(const pmat_elem_type &shp);
  pmat_elem_type to_eltm_object(const mexarg_in &p);

  // Functions for FEM_CLASS_ID
  getfemint_declare_getfem_class(virtual_fem)
  typedef std::shared_ptr<const getfem::virtual_fem> pfem;
  bool is_fem_object(const mexarg_in &p);
  id_type store_fem_object(const pfem &shp);
  pfem to_fem_object(const mexarg_in &p);

  // Functions for GEOTRANS_CLASS_ID
  getfemint_delare_bgeot_class(geometric_trans)
  typedef std::shared_ptr<const bgeot::geometric_trans> pgeometric_trans;
  bool is_geotrans_object(const mexarg_in &p);
  id_type store_geotrans_object(const pgeometric_trans &shp);
  pgeometric_trans to_geotrans_object(const mexarg_in &p);

  // Functions for GLOBAL_FUNCTION_CLASS_ID
  getfemint_declare_getfem_class(abstract_xy_function)
  typedef std::shared_ptr<const getfem::abstract_xy_function> pxy_function;
  bool is_global_function_object(const mexarg_in &p);
  id_type store_global_function_object(const pxy_function &shp);
  pxy_function to_global_function_object(const mexarg_in &p);

  // Functions for INTEG_CLASS_ID
  getfemint_declare_getfem_class(integration_method)
  typedef std::shared_ptr<const getfem::integration_method> pintegration_method;
  bool is_integ_object(const mexarg_in &p);
  id_type store_integ_object(const pintegration_method &shp);
  pintegration_method to_integ_object(const mexarg_in &p);

  // Functions for LEVELSET_CLASS_ID
  getfemint_declare_getfem_class(level_set)
  bool is_levelset_object(const mexarg_in &p);
  id_type store_levelset_object(const std::shared_ptr<getfem::level_set> &shp);
  getfem::level_set *to_levelset_object(const mexarg_in &p);

  // Functions for MESH_CLASS_ID
  getfemint_declare_getfem_class(mesh)
  bool is_mesh_object(const mexarg_in &p);
  id_type store_mesh_object(const std::shared_ptr<getfem::mesh> &shp);
  getfem::mesh *to_mesh_object(const mexarg_in &p);
  bool has_mesh_object(const mexarg_in &p); // mesh or mesh_fem or mesh_im ...
  getfem::mesh *extract_mesh_object(const mexarg_in &p);

  // Functions for MESHFEM_CLASS_ID
  getfemint_declare_getfem_class(mesh_fem)
  bool is_meshfem_object(const mexarg_in &p);
  id_type store_meshfem_object(const std::shared_ptr<getfem::mesh_fem> &shp);
  getfem::mesh_fem *to_meshfem_object(const mexarg_in &p);

  // Functions for MESHIM_CLASS_ID
  getfemint_declare_getfem_class(mesh_im)
  bool is_meshim_object(const mexarg_in &p);
  id_type store_meshim_object(const std::shared_ptr<getfem::mesh_im> &shp);
  getfem::mesh_im *to_meshim_object(const mexarg_in &p);

  // Functions for MESHIMDATA_CLASS_ID
  getfemint_declare_getfem_class(im_data)
  bool is_meshimdata_object(const mexarg_in &p);
  id_type store_meshimdata_object(const std::shared_ptr<getfem::im_data> &shp);
  getfem::im_data *to_meshimdata_object(const mexarg_in &p);

  // Functions for MESH_LEVELSET_CLASS_ID
  getfemint_declare_getfem_class(mesh_level_set)
  bool is_mesh_levelset_object(const mexarg_in &p);
  id_type store_mesh_levelset_object
  (const std::shared_ptr<getfem::mesh_level_set> &shp);
  getfem::mesh_level_set *to_mesh_levelset_object(const mexarg_in &p);
   
  // Functions for MESHER_OBJECT_CLASS_ID
  getfemint_declare_getfem_class(mesher_signed_distance)
  typedef std::shared_ptr<const getfem::mesher_signed_distance>
    pmesher_signed_distance;
  bool is_mesher_object(const mexarg_in &p);
  id_type store_mesher_object(const pmesher_signed_distance &shp);
  pmesher_signed_distance to_mesher_object(const mexarg_in &p);

  // Functions for MODEL_CLASS_ID
  getfemint_declare_getfem_class(model)
  bool is_model_object(const mexarg_in &p);
  id_type store_model_object(const std::shared_ptr<getfem::model> &shp);
  getfem::model *to_model_object(const mexarg_in &p);

  // Functions for PRECOND_CLASS_ID
  class gprecond_base;
  bool is_precond_object(const mexarg_in &p);
  id_type store_precond_object(const std::shared_ptr<gprecond_base> &shp);
  gprecond_base *to_precond_object(const mexarg_in &p);

  // Functions for SLICE_CLASS_ID
  getfemint_declare_getfem_class(stored_mesh_slice)
  bool is_slice_object(const mexarg_in &p);
  id_type store_slice_object
  (const std::shared_ptr<getfem::stored_mesh_slice> &shp);
  getfem::stored_mesh_slice *to_slice_object(const mexarg_in &p);

  // Functions for SPMAT_CLASS_ID
  class gsparse;
  bool is_spmat_object(const mexarg_in &p);
  id_type store_spmat_object
  (const std::shared_ptr<gsparse> &shp);
  gsparse *to_spmat_object(const mexarg_in &p);

  // Functions for POLY_CLASS_ID
  struct getfemint_poly : virtual public dal::static_stored_object
  { bgeot::base_poly p; };
  bool is_poly_object(const mexarg_in &p);
  id_type store_poly_object(const std::shared_ptr<getfemint_poly> &shp);
  getfemint_poly *to_poly_object(const mexarg_in &p);




}  /* end of namespace getfemint.                                          */

#endif /* GETFEMINT_H__                                                    */
