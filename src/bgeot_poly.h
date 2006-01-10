// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_poly.h : plain polynomials with several variables.
//           
// Date    : December 01, 2000.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2000-2006 Yves Renard
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


#ifndef BGEOT_POLY_H__
#define BGEOT_POLY_H__

/** @file bgeot_poly.h
    @brief Multivariate polynomials.
*/

#include <bgeot_config.h>
#include <vector>

namespace bgeot
{
  /// used as the common size type in the library
  typedef size_t size_type;
  ///
  /// used as the common short type integer in the library
  typedef dal::uint16_type short_type;
  ///
  
  /** Return the value of @f$ \frac{(n+p)!}{n!p!} @f$ which
   * is the number of monomials of a polynomial of @f$n@f$
   * variables and degree @f$d@f$.
   */
  size_type alpha(short_type n, short_type d);
  
  /** Vector of integer (16 bits type) which represent the powers
   *  of a monomial
   */
  class power_index {
    std::vector<short_type> v;
    mutable short_type degree_;
    mutable size_type global_index_;
    void dirty() const { degree_ = short_type(-1); global_index_ = size_type(-1); }
  public :
    typedef std::vector<short_type>::iterator iterator;
    typedef std::vector<short_type>::const_iterator const_iterator;
    short_type operator[](size_type idx) const { return v[idx]; }
    short_type& operator[](size_type idx) { dirty(); return v[idx]; }
    iterator begin() { dirty(); return v.begin(); }
    const_iterator begin() const { return v.begin(); }
    iterator end() { dirty(); return v.end(); }
    const_iterator end() const { return v.end(); }
    size_type size() const { return v.size(); }
    /// Gives the next power index 
    const power_index &operator ++();
    /// Gives the next power index 
    const power_index operator ++(int)
      { power_index res = *this; ++(*this); return res; } 
    /// Gives the next previous index 
    const power_index &operator --();
    /// Gives the next previous index 
    const power_index operator --(int)
      { power_index res = *this; --(*this); return res; }
    /**  Gives the global number of the index (i.e. the position of
     *   the corresponding monomial
     */
    size_type global_index(void) const;
    /// Gives the degree.
    short_type degree(void) const; 
    /// Constructor
    power_index(short_type nn);
    /// Constructor
    power_index(void) { dirty(); }
  };
  
  /**
   * This class deals with plain polynomials with
   * several variables. 
   *
   * A polynomial of @f$n@f$ variables and degree @f$d@f$ is stored in a vector
   * of @f$\alpha_d^n@f$ components.
   *
   * <h3>Example of code</h3>
   *
   *   the following code is valid :
   *   @code
   *   #include<bgeot_poly.h>
   *   bgeot::polynomial<double> P, Q;
   *   P = bgeot::polynomial<double>(2,2,1); // P = x
   *   Q = bgeot::polynomial<double>(2,2,2); // Q = y
   *   P += Q; // P is equal to x+y. 
   *   P *= Q; // P is equal to xy + y^2
   *   bgeot::power_index pi(P.dim()); 
   *   bgeot::polynomial<double>::const_iterator ite = Q.end();
   *   bgeot::polynomial<double>::const_iterator itb = Q.begin();
   *   for ( ; itb != ite; ++itb, ++pi)
   *     if (*itq != double(0))
   *       cout "there is x to the power " << pi[0]
   *             << " and y to the power "
   *             << pi[1] << " and a coefficient " << *itq << endl;
   *  @endcode
   *
   *  <h3>Monomials ordering.</h3>
   *
   *       The constant coefficient is placed first with the index 0.
   *       Two monomials of different degrees are ordered following
   *       there respective degree.
   *
   *       If two monomials have the same degree, they are ordered with the
   *       degree of the mononomials without the n firsts variables which
   *       have the same degree. The index of the monomial
   *       @f$ x_0^{i_0}x_1^{i_1} ... x_{n-1}^{i_{n-1}} @f$
   *       is then
   *       @f$ \alpha_{d-1}^{n} + \alpha_{d-i_0-1}^{n-1} 
   *          + \alpha_{d-i_0-i_1-1}^{n-2} + ... + \alpha_{i_{n-1}-1}^{1}, @f$
   *       where @f$d = \sum_{l=0}^{n-1} i_l@f$ is the degree of the monomial.
   *       (by convention @f$\alpha_{-1}^{n} = 0@f$).
   *
   *  <h3>Dealing with the vector of power.</h3>
   *
   *        The answer to the question : what is the next and previous
   *        monomial of @f$x_0^{i_0}x_1^{i_1} ... x_{n-1}^{i_{n-1}}@f$ in the
   *        vector is the following :
   *
   *        To take the next coefficient, let @f$l@f$ be the last index between 0
   *        and @f$n-2@f$ such that @f$i_l \ne 0@f$ (@f$l = -1@f$ if there is not), then
   *        make the operations @f$a = i_{n-1}; i_{n-1} = 0; i_{l+1} = a+1;
   *        \mbox{ if } l \ge 0 \mbox{ then } i_l = i_l - 1@f$.
   *
   *        To take the previous coefficient, let @f$l@f$ be the last index 
   *        between 0 and @f$n-1@f$ such that @f$i_l \ne 0@f$ (if there is not, there
   *        is no previous monomial) then make the operations @f$a = i_l;
   *        i_l = 0; i_{n-1} = a - 1; \mbox{ if } l \ge 1 \mbox{ then } 
   *        i_{l-1} = i_{l-1} + 1@f$.
   *
   *  <h3>Direct product multiplication.</h3>
   *
   *        This direct product multiplication of P and Q is the
   *        multiplication considering that the variables of Q follow the
   *        variables of P. The result is a polynomial with the number of
   *        variables of P plus the number of variables of Q.
   *        The resulting polynomials have a smaller degree.
   *
   *  @todo <h3>Horner scheme to evaluate polynomials.</h3>
   *
   */
  template<typename T> class polynomial : public std::vector<T> {
  protected :
    
    short_type n, d;
    
  public :
    
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    
    /// Gives the degree of the polynomial
    short_type degree(void) const { return d; }
    /**  gives the degree of the polynomial, considering only non-zero
     * coefficients
     */
    short_type real_degree(void) const;
    ///     Gives the dimension (number of variables)
    short_type dim(void) const { return n; }
    /// Change the degree of the polynomial to d.
    void change_degree(short_type dd);
    /** Add to the polynomial a monomial of coefficient a and
     * correpsonding to the power index pi.
     */
    void add_monomial(const T &coeff, const power_index &power);
    ///  Add Q to P. P contains the result.
    polynomial &operator +=(const polynomial &Q);
    /// Substract Q to P. P contains the result.
    polynomial &operator -=(const polynomial &Q);
    /// Add Q to P.
    polynomial operator +(const polynomial &Q) const
      { polynomial R = *this; R += Q; return R; }
    /// Substract Q to P.
    polynomial operator -(const polynomial &Q) const
      { polynomial R = *this; R -= Q; return R; }
    polynomial operator -(void) const;
    /// Multiply P with Q. P contains the result.
    polynomial &operator *=(const polynomial &Q);
    /// Multiply P with Q. 
    polynomial operator *(const polynomial &Q) const;
    /** Product of P and Q considering that variables of Q come after
     * variables of P. P contains the result
     */
    void direct_product(const polynomial &Q);
    /// Multiply P with the scalar a. P contains the result.
    polynomial &operator *=(const T &e);
    /// Multiply P with the scalar a.
    polynomial operator *(const T &e) const;
    /// Divide P with the scalar a. P contains the result.
    polynomial &operator /=(const T &e);
    /// Divide P with the scalar a.
    polynomial operator /(const T &e) const
      { polynomial res = *this; res /= e; return res; }   
    /// operator ==.
    bool operator ==(const polynomial &Q) const; 
    /// operator !=.
    bool operator !=(const polynomial &Q) const
    { return !(operator ==(*this,Q)); }   
    /// Derivative of P with respect to the variable k. P contains the result.
    void derivative(short_type k);
    /// Makes P = 1.
    void one(void) { change_degree(0); (*this)[0] = T(1); }
    void clear(void) { change_degree(0); (*this)[0] = T(0); }
    template <typename ITER> T horner(power_index &mi, short_type k,
				   short_type de, const ITER &it) const;
    /** Evaluate the polynomial. "it" is an iterator pointing to the list
     * of variables. A Horner scheme is used.
     */
    template <typename ITER> T eval(const ITER &it) const;

    /// Constructor.
    polynomial(void) : std::vector<T>(1)
      { n = 0; d = 0; (*this)[0] = 0.0; }
    /// Constructor.
    polynomial(short_type dim_, short_type degree_);
    /// Constructor for the polynomial 'x' (k=0), 'y' (k=1), 'z' (k=2) etc.
    polynomial(short_type dim_, short_type degree_, short_type k);
  };


  template<typename T> polynomial<T>::polynomial(short_type nn, short_type dd)
    : std::vector<T>(alpha(nn,dd))
  { n = nn; d = dd; std::fill(this->begin(), this->end(), T(0)); }

  template<typename T> polynomial<T>::polynomial(short_type nn,
						 short_type dd, short_type k)
    : std::vector<T>(alpha(nn,dd)) {
    n = nn; d = std::max(short_type(1), dd);
    std::fill(this->begin(), this->end(), T(0));
    (*this)[k+1] = T(1);
  }

  template<typename T>
  polynomial<T> polynomial<T>::operator *(const polynomial &Q) const
  { polynomial res = *this; res *= Q; return res; }

  template<typename T>
  bool polynomial<T>::operator ==(const polynomial &Q) const { 
    if (dim() != Q.dim()) return false;
    const_iterator it1 = this->begin(), ite1 = this->end();
    const_iterator it2 = Q.begin(), ite2 = Q.end();
    for ( ; it1 != ite1 && it2 != ite2; ++it1, ++it2)
      if (*it1 != *it2) return false;
    for ( ; it1 != ite1; ++it1) if (*it1 != T(0)) return false;
    for ( ; it2 != ite2; ++it2) if (*it2 != T(0)) return false;
    return true;
  }
  
  template<typename T>
  polynomial<T> polynomial<T>::operator *(const T &e) const
  { polynomial res = *this; res *= e; return res; }
  
  template<typename T> short_type polynomial<T>::real_degree(void) const {
    const_iterator it = this->end() - 1, ite = this->begin() - 1;
    size_type l = this->size();
    for ( ; it != ite; --it, --l) { if (*it != T(0)) break; }
    short_type dd = degree();
    while (dd > 0 && alpha(n, dd-1) > l) --dd;
    return dd;
  }
  
  template<typename T> void polynomial<T>::change_degree(short_type dd) {
    this->resize(alpha(n,dd));
    if (dd > d) std::fill(this->begin() + alpha(n,d), this->end(), T(0));
    d = dd;
  }
  
  template<typename T>
  void polynomial<T>::add_monomial(const T &coeff, const power_index &power) {
    size_type i = power.global_index();
    if (n != power.size()) DAL_THROW(dimension_error, "dimensions mismatch");
    if (i >= this->size()) { change_degree(power.degree()); }
    ((*this)[i]) += coeff;
  }
  
  template<typename T> 
  polynomial<T> &polynomial<T>::operator +=(const polynomial &Q) {
    if (Q.dim() != dim())
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    if (Q.degree() > degree()) change_degree(Q.degree());
    iterator it = this->begin();
    const_iterator itq = Q.begin(), ite = Q.end();
    for ( ; itq != ite; ++itq, ++it) *it += *itq;
    return *this;
  }
  
  template<typename T> 
  polynomial<T> &polynomial<T>::operator -=(const polynomial &Q) {
    if (Q.dim() != dim() || dim() == 0)
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    if (Q.degree() > degree()) change_degree(Q.degree());
    iterator it = this->begin();
    const_iterator itq = Q.begin(), ite = Q.end();
    for ( ; itq != ite; ++itq, ++it) *it -= *itq;
    return *this;
  }
  
  template<typename T> 
  polynomial<T> polynomial<T>::operator -(void) const {
    polynomial<T> Q = *this;
    iterator itq = Q.begin(), ite = Q.end();
    for ( ; itq != ite; ++itq) *itq = -(*itq);
    return Q;
  }
  
  template<typename T>
  polynomial<T> &polynomial<T>::operator *=(const polynomial &Q) {
    if (Q.dim() != dim())
      DAL_THROW(dimension_error, "dimensions mismatch");
    
    polynomial aux = *this;
    change_degree(0); (*this)[0] = T(0);
    
    power_index miq(Q.dim()), mia(dim()), mitot(dim());
    if (dim() > 0) miq[dim()-1] = Q.degree();
    const_iterator itq = Q.end() - 1, ite = Q.begin() - 1;
    for ( ; itq != ite; --itq, --miq)
      if (*itq != T(0)) {
	iterator ita = aux.end() - 1, itae = aux.begin() - 1;
	std::fill(mia.begin(), mia.end(), 0);
	if (dim() > 0) mia[dim()-1] = aux.degree();
	for ( ; ita != itae; --ita, --mia)
	  if (*ita != T(0)) {
	    power_index::iterator mita = mia.begin(), mitq = miq.begin();
	    power_index::iterator mit = mitot.begin(), mite = mia.end();
	    for ( ; mita != mite; ++mita, ++mitq, ++mit)
	      *mit = (*mita) + (*mitq); /* on pourrait calculer
					   directement l'index global. */
	    //	     cerr << "*= : " << *this << ", itq*ita=" << (*itq) * (*ita) << endl;
	    //	     cerr << " itq = " << *itq << endl;
	    //	     cerr << " ita = " << *ita << endl;
	    add_monomial((*itq) * (*ita), mitot);
	   
	  }
      }
    return *this;
  }

  template<typename T>
    void polynomial<T>::direct_product(const polynomial &Q) { 
    polynomial aux = *this;

    change_degree(0); n += Q.dim(); (*this)[0] = T(0);
    
    power_index miq(Q.dim()), mia(aux.dim()), mitot(dim());
    if (Q.dim() > 0) miq[Q.dim()-1] = Q.degree();
    const_iterator itq = Q.end() - 1, ite = Q.begin() - 1;
    for ( ; itq != ite; --itq, --miq)
      if (*itq != T(0)) {
	iterator ita = aux.end() - 1, itae = aux.begin() - 1;
	std::fill(mia.begin(), mia.end(), 0); 
	if (aux.dim() > 0) mia[aux.dim()-1] = aux.degree();
	for ( ; ita != itae; --ita, --mia)
	  if (*ita != T(0)) {
	    std::copy(mia.begin(), mia.end(), mitot.begin());
	    std::copy(miq.begin(), miq.end(), mitot.begin() + aux.dim());
	    add_monomial((*itq) * (*ita), mitot); /* on pourrait calculer
					   directement l'index global. */
	  }
      }
  }

  template<typename T>
    polynomial<T> &polynomial<T>::operator *=(const T &e) {
    iterator it = this->begin(), ite = this->end();
    for ( ; it != ite; ++it) (*it) *= e;
    return *this;
  }

  template<typename T>
    polynomial<T> &polynomial<T>::operator /=(const T &e) {
    iterator it = this->begin(), ite = this->end();
    for ( ; it != ite; ++it) (*it) /= e;
    return *this;
  }

  template<typename T>
  inline void polynomial<T>::derivative(short_type k) {
    if (k >= n)
      DAL_THROW(std::out_of_range, "index out of range");
    
     iterator it = this->begin(), ite = this->end();
     power_index mi(dim()); 
     for ( ; it != ite; ++it, ++mi) {
       if ((*it) != T(0) && mi[k] > 0)
         { mi[k]--; (*this)[mi.global_index()] = (*it) * T(mi[k] + 1); mi[k]++; }
       *it = T(0);
     }
     if (d > 0) change_degree(d-1);
  }

   template<typename T> template<typename ITER>
  inline T polynomial<T>::horner(power_index &mi, short_type k, 
				 short_type de, const ITER &it) const {
    if (k == 0)
      return (*this)[mi.global_index()];
    else {
      T v = (*(it+k-1)), res = T(0);
      for (mi[k-1] = degree() - de; mi[k-1] != short_type(-1); (mi[k-1])--)
	res = horner(mi, k-1, de + mi[k-1], it) + v * res;
      mi[k-1] = 0;
      return res;
    }
  }


  template<typename T> template<typename ITER>
  T polynomial<T>::eval(const ITER &it) const {
    /* direct evaluation for common low degree polynomials */
    unsigned deg = degree();
    const_iterator P = this->begin();
    if (deg == 0) return P[0];
    else if (deg == 1) {
      T s = P[0];
      for (size_type i=0; i < dim(); ++i) s += it[i]*P[i+1];
      return s;
    }
 
    switch (dim()) {
      case 1: {
	T x = it[0];
	if (deg == 2)     return P[0] + x*(P[1] + x*(P[2]));
	if (deg == 3)     return P[0] + x*(P[1] + x*(P[2] + x*(P[3])));
	if (deg == 4)     return P[0] + x*(P[1] + x*(P[2] + x*(P[3] + x*(P[4]))));
	if (deg == 5)     return P[0] + x*(P[1] + x*(P[2] + x*(P[3] + x*(P[4] + x*(P[5])))));
	if (deg == 6)     return P[0] + x*(P[1] + x*(P[2] + x*(P[3] + x*(P[4] + x*(P[5] + x*(P[6]))))));
      } break;
      case 2: {
	T x = it[0];
	T y = it[1];
	if (deg == 2)     return P[0] + x*(P[1] + x*(P[3])) + y*(P[2] + x*(P[4]) + y*(P[5]));
	if (deg == 3)     return P[0] + x*(P[1] + x*(P[3] + x*(P[6]))) + y*(P[2] + x*(P[4] + x*(P[7])) + y*(P[5] + x*(P[8]) + y*(P[9])));
	if (deg == 4)     return P[0] + x*(P[1] + x*(P[3] + x*(P[6] + x*(P[10])))) + y*(P[2] + x*(P[4] + x*(P[7] + x*(P[11]))) + y*(P[5] + x*(P[8] + x*(P[12])) + y*(P[9] + x*(P[13]) + y*(P[14]))));
	if (deg == 5)     return P[0] + x*(P[1] + x*(P[3] + x*(P[6] + x*(P[10] + x*(P[15]))))) + y*(P[2] + x*(P[4] + x*(P[7] + x*(P[11] + x*(P[16])))) + y*(P[5] + x*(P[8] + x*(P[12] + x*(P[17]))) + y*(P[9] + x*(P[13] + x*(P[18])) + y*(P[14] + x*(P[19]) + y*(P[20])))));
	if (deg == 6)     return P[0] + x*(P[1] + x*(P[3] + x*(P[6] + x*(P[10] + x*(P[15] + x*(P[21])))))) + y*(P[2] + x*(P[4] + x*(P[7] + x*(P[11] + x*(P[16] + x*(P[22]))))) + y*(P[5] + x*(P[8] + x*(P[12] + x*(P[17] + x*(P[23])))) + y*(P[9] + x*(P[13] + x*(P[18] + x*(P[24]))) + y*(P[14] + x*(P[19] + x*(P[25])) + y*(P[20] + x*(P[26]) + y*(P[27]))))));
      } break;
      case 3: {
	T x = it[0];
	T y = it[1];
	T z = it[2];
	if (deg == 2)     return P[0] + x*(P[1] + x*(P[4])) + y*(P[2] + x*(P[5]) + y*(P[7])) + z*(P[3] + x*(P[6]) + y*(P[8]) + z*(P[9]));
	if (deg == 3)     return P[0] + x*(P[1] + x*(P[4] + x*(P[10]))) + y*(P[2] + x*(P[5] + x*(P[11])) + y*(P[7] + x*(P[13]) + y*(P[16]))) + z*(P[3] + x*(P[6] + x*(P[12])) + y*(P[8] + x*(P[14]) + y*(P[17])) + z*(P[9] + x*(P[15]) + y*(P[18]) + z*(P[19])));
	if (deg == 4)     return P[0] + x*(P[1] + x*(P[4] + x*(P[10] + x*(P[20])))) + y*(P[2] + x*(P[5] + x*(P[11] + x*(P[21]))) + y*(P[7] + x*(P[13] + x*(P[23])) + y*(P[16] + x*(P[26]) + y*(P[30])))) + z*(P[3] + x*(P[6] + x*(P[12] + x*(P[22]))) + y*(P[8] + x*(P[14] + x*(P[24])) + y*(P[17] + x*(P[27]) + y*(P[31]))) + z*(P[9] + x*(P[15] + x*(P[25])) + y*(P[18] + x*(P[28]) + y*(P[32])) + z*(P[19] + x*(P[29]) + y*(P[33]) + z*(P[34]))));
	if (deg == 5)     return P[0] + x*(P[1] + x*(P[4] + x*(P[10] + x*(P[20] + x*(P[35]))))) + y*(P[2] + x*(P[5] + x*(P[11] + x*(P[21] + x*(P[36])))) + y*(P[7] + x*(P[13] + x*(P[23] + x*(P[38]))) + y*(P[16] + x*(P[26] + x*(P[41])) + y*(P[30] + x*(P[45]) + y*(P[50]))))) + z*(P[3] + x*(P[6] + x*(P[12] + x*(P[22] + x*(P[37])))) + y*(P[8] + x*(P[14] + x*(P[24] + x*(P[39]))) + y*(P[17] + x*(P[27] + x*(P[42])) + y*(P[31] + x*(P[46]) + y*(P[51])))) + z*(P[9] + x*(P[15] + x*(P[25] + x*(P[40]))) + y*(P[18] + x*(P[28] + x*(P[43])) + y*(P[32] + x*(P[47]) + y*(P[52]))) + z*(P[19] + x*(P[29] + x*(P[44])) + y*(P[33] + x*(P[48]) + y*(P[53])) + z*(P[34] + x*(P[49]) + y*(P[54]) + z*(P[55])))));
	if (deg == 6)     return P[0] + x*(P[1] + x*(P[4] + x*(P[10] + x*(P[20] + x*(P[35] + x*(P[56])))))) + y*(P[2] + x*(P[5] + x*(P[11] + x*(P[21] + x*(P[36] + x*(P[57]))))) + y*(P[7] + x*(P[13] + x*(P[23] + x*(P[38] + x*(P[59])))) + y*(P[16] + x*(P[26] + x*(P[41] + x*(P[62]))) + y*(P[30] + x*(P[45] + x*(P[66])) + y*(P[50] + x*(P[71]) + y*(P[77])))))) + z*(P[3] + x*(P[6] + x*(P[12] + x*(P[22] + x*(P[37] + x*(P[58]))))) + y*(P[8] + x*(P[14] + x*(P[24] + x*(P[39] + x*(P[60])))) + y*(P[17] + x*(P[27] + x*(P[42] + x*(P[63]))) + y*(P[31] + x*(P[46] + x*(P[67])) + y*(P[51] + x*(P[72]) + y*(P[78]))))) + z*(P[9] + x*(P[15] + x*(P[25] + x*(P[40] + x*(P[61])))) + y*(P[18] + x*(P[28] + x*(P[43] + x*(P[64]))) + y*(P[32] + x*(P[47] + x*(P[68])) + y*(P[52] + x*(P[73]) + y*(P[79])))) + z*(P[19] + x*(P[29] + x*(P[44] + x*(P[65]))) + y*(P[33] + x*(P[48] + x*(P[69])) + y*(P[53] + x*(P[74]) + y*(P[80]))) + z*(P[34] + x*(P[49] + x*(P[70])) + y*(P[54] + x*(P[75]) + y*(P[81])) + z*(P[55] + x*(P[76]) + y*(P[82]) + z*(P[83]))))));
      } break;
    }

    /*
    switch (deg) {
      case 0: return (*this)[0];
      case 1: { 
	T s = (*this)[0];
	for (size_type i=0; i < dim(); ++i) s += it[i]*(*this)[i+1];
	return s; 
      }
      case 2:
      case 3: { 
	if (dim() == 1) {
	  const T &x = it[0]; 
	  if      (deg == 2) return p[0] + x*(p[1] + x*p[2]);
	  else if (deg == 3) return p[0] + x*(p[1] + x*(p[2]+x*p[3]));
	} else if (dim() == 2) {
	  const T &x = it[0]; 
	  const T &y = it[1]; 
	  if      (deg == 2) 
	    return p[0] + p[1]*x + p[2]*y + p[3]*x*x + p[4]*x*y + p[5]*y*y;
	  else if (deg == 3)
	    return p[0] + p[1]*x + p[2]*y + p[3]*x*x + p[4]*x*y + p[5]*y*y + 
	      p[6]*x*x*x + p[7]*x*x*y + p[8]*x*y*y + p[9]*y*y*y;
	} else if (dim() == 3) {
	  const T &x = it[0]; 
	  const T &y = it[1]; 
	  const T &z = it[2]; 
	  if (deg == 2)
	    return p[0] + p[1]*x + p[2]*y + p[3]*z + p[4]*x*x + p[5]*x*y + p[6]*x*z +
	      p[7]*y*y + p[8]*y*z + p[9]*z*z;
	  else if (deg == 3)
	    return p[0] + p[1]*x + p[2]*y + p[3]*z + p[4]*x*x + p[5]*x*y + p[6]*x*z +
	      p[7]*y*y + p[8]*y*z + p[9]*z*z + 
	      p[10]*x*x*x + p[11]*x*x*y + p[12]*x*x*z + p[13]*x*y*y + p[14]*x*y*z + p[15]*x*z*z +
	      p[16]*y*y*y + p[17]*y*y*z + p[18]*y*z*z + 
	      p[19]*z*z*z;
	}
      }
      }*/
    /* for other polynomials, Horner evaluation (quite slow..) */
    power_index mi(dim());
    return horner(mi, dim(), 0, it);
  }

  template<typename ITER>
    typename std::iterator_traits<ITER>::value_type
        eval_monomial(const power_index &mi, ITER it) {
    typename std::iterator_traits<ITER>::value_type res
      = typename std::iterator_traits<ITER>::value_type(1);
    power_index::const_iterator mit = mi.begin(), mite = mi.end();
    for ( ; mit != mite; ++mit, ++it)
      for (short_type l = 0; l < *mit; ++l)
	res *= *it;
    return res;
  }


  /// Print P to the output stream o. for instance cout << P;
  template<typename T>  std::ostream &operator <<(std::ostream &o,
						     const polynomial<T>& P) { 
    bool first = true; size_type n = 0;
    typename polynomial<T>::const_iterator it = P.begin(), ite = P.end();
    power_index mi(P.dim());
    if (it != ite && *it != T(0))
      { o << *it; first = false; ++it; ++n; ++mi; }
    for ( ; it != ite ; ++it, ++mi ) {
      if (*it != T(0)) {
	bool first_var = true;
	if (!first) { if (*it < T(0)) o << " - "; else o << " + "; }
	else if (*it < T(0)) o << "-";
	if (gmm::abs(*it)!=T(1)) { o << gmm::abs(*it); first_var = false; }
	for (short_type j = 0; j < P.dim(); ++j)
	  if (mi[j] != 0) {
	    if (!first_var) o << "*"; first_var = false;
            if (P.dim() <= 7) o << "xyzwvut"[j];
            else o << "x_" << j; 
	    if (mi[j] > 1) o << "^" << mi[j];
	  }
	first = false; ++n;
      }
    }
    if (n == 0) o << "0";
    return o;
  }

  /**
     polynomial variable substitution
     @param P the original polynomial
     @param S the substitution poly (not a multivariate one)
     @param subs_dim : which variable is substituted
     example: poly_subs(x+y*x^2, x+1, 0) = x+1 + y*(x+1)^2
  */
  template<typename T>    
  polynomial<T> poly_substitute_var(const polynomial<T>& P,
				    const polynomial<T>& S,
				    size_type subs_dim) {
    if (S.dim()!=1 || subs_dim >= P.dim())
      DAL_THROW(failure_error, "wrong arguments for polynomial substitution");
    polynomial<T> res(P.dim(),0);
    bgeot::power_index pi(P.dim());
    std::vector< polynomial<T> > Spow(1);
    Spow[0] = polynomial<T>(1, 0); Spow[0].one(); // Spow stores powers of S
    for (size_type k=0; k < P.size(); ++k, ++pi) {
      if (P[k] == T(0)) continue;
      while (pi[subs_dim] >= Spow.size()) 
	Spow.push_back(S*Spow.back());
      const polynomial<T>& p = Spow[pi[subs_dim]];
      bgeot::power_index pi2(pi);
      for (size_type i=0; i < p.size(); ++i) {
	pi2[subs_dim] = i;
	res.add_monomial(p[i]*P[k],pi2);
      }
    }
    return res;
  }
  
  template<typename U, typename T>
  polynomial<T> operator *(U c, const polynomial<T> &p)
  { polynomial<T> q = p; q *= T(c); return q; }

  typedef polynomial<opt_long_scalar_type> base_poly;

  /* usual constant polynomials  */
  
  inline base_poly null_poly(short_type n) { return base_poly(n, 0); }
  inline base_poly one_poly(short_type n)
  { base_poly res=base_poly(n, 0); res.one(); return res;  }
  inline base_poly one_var_poly(short_type n, short_type k)
  { return base_poly(n, 1, k); }

  /** read a base_poly on the stream ist. */
  base_poly read_base_poly(short_type n, std::istream &f);

  /** read a base_poly on the string s. */
  base_poly read_base_poly(short_type n, const std::string &s);

}  /* end of namespace bgeot.                                           */


#endif  /* BGEOT_POLY_H__ */
