/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_vector.h : plain vectors.                              */
/*     									   */
/*                                                                         */
/* Date : June 01, 1995.                                                   */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; version 2 of the License.                 */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software Foundation, */
/* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.         */
/*                                                                         */
/* *********************************************************************** */


#ifndef __BGEOT_VECTOR_H
#define __BGEOT_VECTOR_H

#include <bgeot_config.h>
#include <vector>


namespace bgeot
{

  /* ******************************************************************** */
  /*	Basic Operations for fixed size vectors (inspired from M.T.L.).   */
  /* ******************************************************************** */

  template <int NN> struct bgeot_count { enum { N = NN }; };
  template <> struct bgeot_count<0> { enum { N = 0 }; };

  template <int N, class InIter, class OutIter> inline
    void copy(InIter first, bgeot_count<N>, OutIter result)
  { *result = *first; copy(++first, bgeot_count<N-1>(), ++result); }

  template <class InIter, class OutIter> inline
    void copy(InIter first, bgeot_count<0>, OutIter result)
  { }
  
  template <int N, class OutIter, class T> inline
    void fill(OutIter first, bgeot_count<N>, const T& a)
  { *first = a; fill(++first, bgeot_count<N-1>(), a); }

  template <class OutIter, class T> inline
    void fill(OutIter first, bgeot_count<0>, const T& a)
  { }

  template <int N, class InIter, class OutIter, class T> inline
    void addmul(InIter first, bgeot_count<N>, OutIter result, const T& a)
  { *result += *first * a; addmul(++first, bgeot_count<N-1>(), ++result, a); }

  template <class InIter, class OutIter, class T> inline 
    void addmul(InIter first, bgeot_count<0>, OutIter result, const T& a)
  { }

  template <int N, class OutIter, class T> inline
    void timeegal(OutIter first, bgeot_count<N>, const T& a)
  { *first *= a; timeegal(++first, bgeot_count<N-1>(), a); }

  template <class OutIter, class T> inline
    void timeegal(OutIter first, bgeot_count<0>, const T& a)
  { }

  template <int N, class OutIter, class T> inline
    void overegal(OutIter first, bgeot_count<N>, const T& a)
  { *first /= a; overegal(++first, bgeot_count<N-1>(), a); }

  template <class OutIter, class T> inline
    void overegal(OutIter first, bgeot_count<0>, const T& a)
  { }

  template <int N, class InIter, class OutIter> inline
    void plusegal(InIter first, bgeot_count<N>, OutIter result)
  { *result += *first; plusegal(++first, bgeot_count<N-1>(), ++result); }

  template <class InIter, class OutIter> inline 
    void plusegal(InIter first, bgeot_count<0>, OutIter result)
  { }

  template <int N, class InIter, class OutIter> inline
    void minusegal(InIter first, bgeot_count<N>, OutIter result)
  { *result -= *first; minusegal(++first, bgeot_count<N-1>(), ++result); }

  template <class InIter, class OutIter> inline 
    void minusegal(InIter first, bgeot_count<0>, OutIter result)
  { }

  template <int N, class InIter, class OutIter> inline
    bool egalegal(InIter first, bgeot_count<N>, OutIter result)
  { 
    return (*result == *first) ?
      egalegal(++first, bgeot_count<N-1>(), ++result) : false;
  }

  template <class InIter, class OutIter> inline 
    bool egalegal(InIter first, bgeot_count<0>, OutIter result)
  { return true; }

  template <int N, class InIter, class OutIter> inline
    void vectsp(InIter first, bgeot_count<N>, OutIter result, double& res)
  { 
    res += (*result) * (*first);
    vectsp(++first, bgeot_count<N-1>(), ++result, res);
  }

  template <class InIter, class OutIter> inline 
    void vectsp(InIter first, bgeot_count<0>, OutIter result, double& res)
  { }


  /* ******************************************************************** */
  /*		Class fsvector<T,n>: fixed size vector.	        	  */
  /* ******************************************************************** */

  /** Plain vectors of arbitrary base type
   * ($<$T$>$). It as been build especially for very small vectors.
   * The copy operator is a true copy (i.e. no "sharing" copy).
   * The classes are very similar to std::vector, with additional
   * classical linear algebraic operations.
   */
  template<class T, int N> class fsvector
  {
    public:
    
      typedef T value_type;
      typedef value_type* pointer;
      typedef const value_type* const_pointer;
      typedef value_type* iterator;
      typedef const value_type* const_iterator;
      typedef value_type& reference;
      typedef const value_type& const_reference;
      typedef size_t size_type;
      typedef ptrdiff_t difference_type;

      typedef dal::reverse_iter<const_iterator> const_reverse_iterator;
      typedef dal::reverse_iter<iterator> reverse_iterator;

    protected:

      T array[N];

    public:

      /// Return the iterator on the first component.
      inline iterator begin() { return &(array[0]); }
      /// Return the  constant iterator on the first component.
      inline const_iterator begin() const { return &(array[0]); }
      ///  Return the iterator over the last component.   
      inline iterator end() { return &(array[N]); }
      ///  Return the constant iterator over the last component.
      inline const_iterator end() const { return &(array[N]); }
      inline reverse_iterator rbegin() { return reverse_iterator(end()); }
      inline const_reverse_iterator rbegin() const
      { return const_reverse_iterator(end()); }
      inline reverse_iterator rend() { return reverse_iterator(begin()); }
      inline const_reverse_iterator rend() const
      { return const_reverse_iterator(begin()); }
      /// Return the size of the vector (i.e. n).       
      inline size_type size() const { return N; }
      inline size_type max_size() const { return N; }
      inline size_type capacity() const { return N; }
      /// Return true if the vector has no component (n = 0)
      inline bool empty() const { return (N == 0); }
      void out_of_range_error(void) const;

      #ifdef __GETFEM_VERIFY
        reference operator[](size_type n) { 
	  if (n >= N) out_of_range_error();
	  return array[n];
	}
        const_reference operator[](size_type n) const {
	  if (n >= N) out_of_range_error();
	  return array[n];
	}
      #else
	/// Return a reference on the component n of the vector
        inline reference operator[](size_type n) { return array[n]; }
	/// Return a constant reference on the component n of the vector
        inline const_reference operator[](size_type n) const {return array[n];}
      #endif

      fsvector(void) {}
      /// Constructor. For 2 components initialized with a0 and a1.
      fsvector(T a0, T a1) {
	if (N != 2) DAL_THROW(std::invalid_argument, "bad argument");
	iterator p = begin(); *p++ = a0; *p++ = a1;
      }
      /// Constructor. For 3 components initialized with a0, a1 and a2.
      fsvector(T a0, T a1, T a2) {
	if (N != 3) DAL_THROW(std::invalid_argument, "bad argument");
	iterator p = begin(); *p++ = a0; *p++ = a1; *p++ = a2;
      }
      /// Constructor. For 4 components initialized with a0, a1, a2 and a3.
      fsvector(T a0, T a1, T a2, T a3) {
	if (N != 4) DAL_THROW(std::invalid_argument, "bad argument");
	iterator p = begin(); *p++=a0;*p++=a1;*p++=a2;*p++=a3;
      }
      fsvector(size_type l) { 
	if (N != l) DAL_THROW(std::invalid_argument, "bad argument");
      }
      fsvector(const fsvector<T,N> &v)
      { copy(v.begin(), bgeot_count<N>(), begin()); }
      fsvector &operator =(const fsvector<T, N> &v)
      { copy(v.begin(), bgeot_count<N>(), begin()); return *this; }

      void fill(const T &a) { bgeot::fill(begin(), bgeot_count<N>(), a); }
      void addmul(const T &a, const fsvector<T,N>& v)
      { bgeot::addmul(v.begin(), bgeot_count<N>(), begin(), a); }

      ///  Multiply the vector with the scalar a.
      fsvector &operator *=(const T &a)
      { timeegal(begin(), bgeot_count<N>(), a); return *this; }
      /// Divide the vector with the scalar a.
      fsvector &operator /=(const T &a)
      { overegal(begin(), bgeot_count<N>(), a); return *this; }
      /// Add vector v to current vector.
      fsvector &operator +=(const fsvector<T,N>& v)
      { plusegal(v.begin(), bgeot_count<N>(), begin()); return *this; }
      /// Substract vector v to current vector V.
      fsvector &operator -=(const fsvector<T,N>& v)
      { minusegal(v.begin(), bgeot_count<N>(), begin()); return *this; }
      bool operator ==(const fsvector<T,N> &v) const
      { return egalegal(v.begin(), bgeot_count<N>(), begin()); }
      bool operator !=(const fsvector<T,N> &v) const
      { return !(egalegal(v.begin(), bgeot_count<N>(), begin())); }
      
    size_type memsize() const { return sizeof(fsvector<T,N>); }
  };

  template<class T, int N>  void fsvector<T,N>::out_of_range_error(void) const
  { DAL_THROW(std::out_of_range, "out of range"); }

  /** @name functions on fsvector (fixed size vector)
   */
  //@{

  /// Multiply m by the scalar x.
  template<class T, int N>
    inline fsvector<T,N> operator *(const fsvector<T,N>& m, const T &x)
  { fsvector<T,N> p = m; p *= x; return p; }

  /// Multiply m by the scalar x.
  template<class T, int N>
    inline fsvector<T,N> operator *(T x, const fsvector<T,N>& m)
  { fsvector<T,N> p = m; p *= x; return p; }
  
  /// Divide m by the scalar x.
  template<class T, int N>
    inline fsvector<T,N> operator /(const fsvector<T,N>& m, const T &x)
  { fsvector<T,N> p = m; p /= x; return p; }
  
  /// Add m and n.
  template<class T, int N>
    inline fsvector<T,N> operator +(const fsvector<T,N>& m, const fsvector<T,N>& n)
  { fsvector<T,N> p = m; p += n; return p; }
  
  /// Substract n to m.
  template<class T, int N>
    inline fsvector<T,N> operator -(const fsvector<T,N>& m, const fsvector<T,N>& n)
  { fsvector<T,N> p = m; p -= n; return p; }
  
  /// Gives the opposite of m.
  template<class T, int N> inline fsvector<T,N> operator -(const fsvector<T,N>& m)
  { fsvector<T,N> p = m; p *= -1; return p; }
  
  /// Gives  m.
  template<class T, int N> inline fsvector<T,N> operator +(const fsvector<T,N>& p)
  { return p; }
  
  /// Scalar product between v and w.
  template<class T, int N>
    inline double vect_sp(const fsvector<T,N>&v, const fsvector<T,N>& w)
  {
    double res = 0.0;
    vectsp(v.begin(),bgeot_count<N>(), w.begin(), res);
    return res;
  }

  //@}

  /* ******************************************************************** */
  /*		Class vsvector<T>: variable size vector.     		  */
  /* ******************************************************************** */

  /** Plain vectors of arbitrary base type
   * ($<$T$>$) and with arbitrary number of components.
   *  It as been build especially for small vectors.
   * The copy operator is a true copy (i.e. no "sharing" copy).
   * The classes is directly derived from std::vector, with additional
   * classical linear algebraic operations.
   */
  template<class T> class vsvector : public std::vector<T>
  {
   
    public:

      typedef typename std::vector<T>::size_type size_type;
      typedef typename std::vector<T>::iterator iterator;
      typedef typename std::vector<T>::const_iterator const_iterator;

      void out_of_range_error(void) const;

      #ifdef __GETFEM_VERIFY
      inline const T& operator [](size_type l) const {
	if (l >= size()) out_of_range_error();
	return *(begin()+l);
      }
      inline T& operator [](size_type l) { 
	if (l >= size()) out_of_range_error();
	return *(begin()+l);
      }
      #endif

      void fill(const T &);
      void addmul(const T &, const vsvector<T>&);

      /// Add vector w to current vector.
      vsvector<T>& operator +=(const vsvector<T>& w);
      /// Substract vector w to current vector v.
      vsvector<T>& operator -=(const vsvector<T>& w);
      /// Multiply the current vector with the scalar x.
      vsvector<T>& operator *=(const T &x);
      /// Divide the current vector with the scalar x.
      vsvector<T>& operator /=(const T &x);
    
      /// Constructor. For 2 components initialized with a0 and a1.
      vsvector(T a0, T a1) : std::vector<T>(size_type(2))
      { iterator p = begin(); *p++ = a0; *p++ = a1; }
      /// Constructor. For 3 components initialized with a0, a1 and a2.
      vsvector(T a0, T a1, T a2) : std::vector<T>(size_type(3)) 
      { iterator p = begin(); *p++ = a0; *p++ = a1; *p++ = a2;}
      /// Constructor. For 4 components initialized with a0, a1, a2 and a3.
      vsvector(T a0, T a1, T a2, T a3) : std::vector<T>(size_type(4))
      { iterator p = begin(); *p++ = a0; *p++ = a1; *p++ = a2; *p++ = a3;}
      /// Constructor. A vector with l components.
      vsvector(size_type l) : std::vector<T>(l) {}
      /// Constructor.
      vsvector(void) : std::vector<T>() {}

    size_type memsize() const { return std::vector<T>::capacity()*sizeof(T) + sizeof(vsvector<T>); }
  };

  template<class T>  void vsvector<T>::out_of_range_error(void) const
  { DAL_THROW(std::out_of_range, "out of range"); }


  template<class T>  void vsvector<T>::addmul(const T &a, const vsvector<T> &v)
  {                             
    typename vsvector<T>::iterator d1 = begin(), e = end();
    const_iterator d2 = v.begin();
    if ( v.size() != this->size())
      DAL_THROW(dimension_error, "dimensions mismatch");
    while (d1 != e) *d1++ += (*d2++) * a;
  }

  template<class T>  void vsvector<T>::fill(const T &x)
  {
    typename vsvector<T>::iterator d = begin(), e = end();
    while (d != e) *d++ = x;
  }
 
  template<class T>  vsvector<T>& vsvector<T>::operator *=(const T &x)
  {
    typename vsvector<T>::iterator d = begin(), e = end();
    while (d != e) *d++ *= x;
    return *this;
  }

  template<class T>  vsvector<T>& vsvector<T>::operator /=(const T &x)
  {
    typename vsvector<T>::iterator d = begin(), e = end();
    while (d != e) *d++ /= x;
    return *this;
  }

  template<class T> vsvector<T>& vsvector<T>::operator +=(const vsvector<T>& w)
  {                             
    typename vsvector<T>::iterator d1 = begin(), e = end();
    typename vsvector<T>::const_iterator d2 = w.begin();
    if (size() != w.size()) DAL_THROW(dimension_error, "dimensions mismatch");
    while (d1 != e) *d1++ += (*d2++);
    return *this;
  }

  template<class T> vsvector<T>& vsvector<T>::operator -=(const vsvector<T>& w)
  {    
    typename vsvector<T>::iterator d1 = begin(), e = end();
    typename vsvector<T>::const_iterator d2 = w.begin();
    if (size() != w.size()) DAL_THROW(dimension_error, "dimensions mismatch");

    while (d1 != e) *d1++ -= (*d2++);
    return *this;
  }

  /** @name functions on vsvector (variable size vector)
   */
  //@{

  /// Multiply m by the scalar x.
  template<class T> inline vsvector<T> operator *(const vsvector<T>& m, const T &x)
  { vsvector<T> p = m; p *= x; return p; }

  /// Multiply m by the scalar x.
  template<class T> inline vsvector<T> operator *(T x, const vsvector<T>& m)
  { vsvector<T> p = m; p *= x; return p; }

  /// Divide m by the scalar x.
  template<class T> inline vsvector<T> operator /(const vsvector<T>& m, const T &x)
  { vsvector<T> p = m; p /= x; return p; }
  
  /// Add m and n.
  template<class T>
    inline vsvector<T> operator +(const vsvector<T>& m, const vsvector<T>& n)
  { vsvector<T> p = m; p += n; return p; }
  
  /// Substract n to m.
  template<class T>
    inline vsvector<T> operator -(const vsvector<T>& m, const vsvector<T>& n)
  { vsvector<T> p = m; p -= n; return p; }


  // Opposite of m.
  template<class T> inline vsvector<T> operator -(const vsvector<T>& m)
  { vsvector<T> p = m; p *= T(-1); return p; }

  /// Gives  p.
  template<class T> inline vsvector<T> operator +(const vsvector<T>& p)
  { return p; }

  /// Scalar product between v and w.
  template<class T> double vect_sp(const vsvector<T>&v, const vsvector<T>&w)
  { /* pas bon sur les complexes.                                           */
    typename vsvector<T>::const_iterator d1 = v.begin(), e = v.end();
    typename vsvector<T>::const_iterator d2 = w.begin();
    double res = 0.0;
    if (v.size() != w.size()) DAL_THROW(dimension_error,"dimensions mismatch");
   
    while (d1 != e) res += (*d1++) * (*d2++);
    return res;
  }

  //@}

  /* ******************************************************************** */
  /*		Generic functions on vectors.            		  */
  /* ******************************************************************** */
  
  /** @name Generic functions on vectors
   */
  //@{
  
  template<class VECT> void vect_out(std::ostream &o, const VECT &v)
  {
    o << '[';
    for (typename VECT::size_type l=1; l < v.size(); l++ ) o << v[l-1] << ' ';
    o << v[v.size()-1] << ']';
  }

  /// Print v on the output stream o. For instance, cout $<<$ v;
  template<class T>
    inline std::ostream &operator <<(std::ostream &o, const vsvector<T>& v)
  { vect_out(o,v); return o; }

  /// Print v on the output stream o. For instance, cout $<<$ v;
  template<class T, int N>
    inline std::ostream &operator <<(std::ostream &o, const fsvector<T,N>& v)
  { vect_out(o,v); return o; }

  /// Gives $\displaystyle \sum_{i=0..(n-1)} |v_i|$.
  template<class VEC> double vect_norm1(const VEC &v)
  {
    typename VEC::const_iterator d1 = v.begin(), e = v.end();
    double res = 0.0;
    while (d1 != e) res += (double)dal::abs(*d1++);
    return v;
  }

  /// Gives $\displaystyle (\sum_{i=0..(n-1)} (v_i)^2)^{1/2}$.
  template<class VEC> double vect_norm2(const VEC &v)
  {
    typename VEC::const_iterator d1 = v.begin(), e = v.end();
    double res = 0.0;
    while (d1 != e) res += dal::sqr((double)dal::abs(*d1++));
    return sqrt(res);
  }

  /// Gives $\displaystyle \sup_{i=0..(n-1)} |v_i|$.
  template<class VEC> double vect_norminf(const VEC &v)
  {
    typename VEC::const_iterator d1 = v.begin(), e = v.end();
    double res = 0.0;
    while (d1 != e) res = std::max(res, (double)dal::abs(*d1++));
    return res;
  }

  /// Gives $\displaystyle \sum_{i=0..(n-1)} |v_i - w_i|$.
  template<class VEC> double vect_dist1(const VEC &v, const VEC &w)
  { 
    typename VEC::const_iterator d1 = v.begin(), e = v.end();
    typename VEC::const_iterator d2 = w.begin();
    double res = 0;
    if (v.size() != w.size()) DAL_THROW(dimension_error,"dimensions mismatch");
    while (d1 != e) res += (double)dal::abs(*d1++ - *d2++);
    return res;
  }

  /// Gives $\displaystyle (\sum_{i=0..(n-1)} |v_i - w_i|^2)^{1/2}$.
  template<class VEC> double vect_dist2(const VEC &v, const VEC &w)
  {
    typename VEC::const_iterator d1 = v.begin(), e = v.end();
    typename VEC::const_iterator d2 = w.begin();
    double res = 0;
    if (v.size() != w.size()) DAL_THROW(dimension_error,"dimensions mismatch");
    while (d1 != e) res += dal::sqr((double)dal::abs(*d1++ - *d2++));
    return sqrt(res);
  }

  /// Gives $\displaystyle \sup_{i=0..(n-1)} |v_i - w_i|$.
  template<class VEC> double vect_distinf(const VEC &v, const VEC &w)
  {
    typename VEC::const_iterator d1 = v.begin(), e = v.end();
    typename VEC::const_iterator d2 = w.begin();
    double res = 0;
    if (v.size() != w.size()) DAL_THROW(dimension_error,"dimensions mismatch");

    while (d1 != e) res += std::max(res, (double)dal::abs(*d1++ - *d2++));
    return res;
  }

  //@}


  /* ******************************************************************** */
  /*		Points.                                     		  */
  /* ******************************************************************** */


  /** For any class of vector VECT, the class PT$<$VECT$>$ 
   *  represent the corresponding class of point.
   */
  template<class VECT> class PT : public VECT
  {
    public :
      typedef VECT vector_type;
      typedef typename VECT::value_type value_type;
      /// Constructor.
      PT(const VECT &v) : VECT(v) { }
      /// Constructor.
      PT(const PT &pt, const VECT &v) : VECT(pt)  { *this += v; }
      /// Constructor.
      PT(void) : VECT() {}
      /// Constructor.
      PT(int n) : VECT(n) {}

      PT &operator *=(const value_type &a)
      { *((VECT *)this) *= a; return *this; }
      PT &operator /=(const value_type &a)
      { *((VECT *)this) /= a; return *this; }
      PT &operator +=(const PT &v) { *((VECT *)this) += v; return *this; }
      PT &operator -=(const PT &v) { *((VECT *)this) -= v; return *this; }
  };

  /** @name generic functions on points
   */
  //@{

  /// Return the vector AB, where A and B are some points.
  template<class PT> typename PT::vector_type vector_from(const PT &A,
							  const PT &B)
  { typename PT::vector_type v = B; v -= A; return v; }

  /// Compute in v the vector AB, where A and B are some points.
  template<class PT> void vector_from(const PT &A,
				      const PT &B, typename PT::vector_type &v)
  { v = B; v -= A; }

  
  //@}

  typedef vsvector<scalar_type> base_vector;
  typedef PT<base_vector> base_node;


}  /* end of namespace bgeot.                                           */


#endif  /* __BGEOT_VECTOR_H */
