/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_sub_index.h : sub indexes                                */
/*     									   */
/* Date : October 13, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GMM++                                            */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */

#ifndef __GMM_SUB_INDEX_H
#define __GMM_SUB_INDEX_H

namespace gmm {

  /* ******************************************************************** */
  /*		sub indexes                               		  */
  /* ******************************************************************** */

  struct basic_index : public std::vector<size_t> {
    
    mutable size_type nb_ref;
    // size_type key1; faire la somme des composantes
    // const basic_index *rind; rindex s'il existe
    

    size_t operator[](size_type i) const {
      return (i < size()) ? std::vector<size_t>::operator[](i) : size_type(-1);
    }
    
    basic_index() : nb_ref(1) {}
    basic_index(size_type j) : std::vector<size_t>(j), nb_ref(1) {}
    template <typename IT> basic_index(IT b, IT e)
      : std::vector<size_t>(e-b), nb_ref(1) { std::copy(b, e, begin()); }
    basic_index(const basic_index *pbi) : nb_ref(1) {
      const_iterator it = pbi->begin(), ite = pbi->end();
      size_type i = 0;
      for ( ; it != ite; ++it) i = std::max(i, *it);
      resize(i+1); std::fill(begin(), end(), size_type(-1));
      for (it = pbi->begin(), i = 0; it != ite; ++it, ++i)
	std::vector<size_t>::operator[](*it) = i;
    }

  };

  typedef const basic_index *pbasic_index;

  struct index_generator {

    template <typename IT> static pbasic_index create_index(IT begin, IT end)
    { return new basic_index(begin, end); }
    static pbasic_index create_rindex(pbasic_index pbi)
    { return new basic_index(pbi); }
    static void attach(pbasic_index pbi) { if (pbi) pbi->nb_ref++; }
    static void unattach(pbasic_index pbi)
      { if (pbi && --(pbi->nb_ref) == 0) delete pbi; }

  };


  struct sub_index {

    typedef basic_index base_type;
    typedef base_type::const_iterator const_iterator;

    pbasic_index ind;
    mutable pbasic_index rind;

    inline void test_rind(void) const
      { if (!rind) rind = index_generator::create_rindex(ind); }
    size_type size(void) const { return ind->size(); }
    size_type index(size_type i) const { return (*ind)[i]; }
    size_type rindex(size_type i) const {
      test_rind();
      if (i < rind->size()) return (*rind)[i]; else return size_type(-1);
    }
   
    const_iterator  begin(void) const { return  ind->begin(); }
    const_iterator    end(void) const { return  ind->end();   }
    const_iterator rbegin(void) const { test_rind(); return rind->begin(); }
    const_iterator   rend(void) const { test_rind(); return rind->end();   }

    sub_index() : ind(0), rind(0) {}
    template <typename IT> sub_index(IT it, IT ite)
      : ind(index_generator::create_index(it, ite)), rind(0) {}
    template <typename CONT> sub_index(const CONT &c)
      : ind(index_generator::create_index(c.begin(), c.end())), rind(0) {}
    ~sub_index()
      { index_generator::unattach(rind); index_generator::unattach(ind); }
    sub_index(const sub_index &si) : ind(si.ind), rind(si.rind)
      { index_generator::attach(rind); index_generator::attach(ind); }
    sub_index &operator =(const sub_index &si) {
      index_generator::unattach(rind); index_generator::unattach(ind);
      ind = si.ind; rind = si.rind; index_generator::attach(rind);
      index_generator::attach(ind);
      return *this;
    }
  };

  inline std::ostream &operator << (std::ostream &o, const sub_index &si) { 
    o << "sub_index(";
    if (si.size()) o << si.index(0);
    for (size_type i = 0; i < si.size(); ++i) o << ", " << si.index(i);
    o << ")";
    return o;
  }

  struct sub_interval {
    size_type min, max; 

    size_type size(void) const { return max - min; }
    size_type index(size_type i) const { return min + i; }
    size_type step(void) const { return 1; }
    size_type rindex(size_type i) const
    { if (i >= min && i < max) return i - min; return size_type(-1); }
    sub_interval(size_type mi, size_type l) : min(mi), max(mi+l) {}
    sub_interval() {}
  };

  inline std::ostream &operator << (std::ostream &o, const sub_interval &si)
  { o << "sub_interval(" << si.min << ", " << si.size() << ")"; return o; }

  struct sub_slice {
    size_type min, max, N; 

    size_type size(void) const { return (max - min) / N; }
    size_type step(void) const { return N; }
    size_type index(size_type i) const { return min + N * i; }
    size_type rindex(size_type i) const { 
      if (i >= min && i < max)
	{ size_type j = (i - min); if (j % N == 0) return j / N; }
      return size_type(-1);
    }
    sub_slice(size_type mi, size_type l, size_type n)
      : min(mi), max(mi+l*n), N(n) {}
    sub_slice(void) {}
  };

  inline std::ostream &operator << (std::ostream &o, const sub_slice &si) {
    o << "sub_slice(" << si.min << ", " << si.size() << ", " << si.step() 
      << ")"; return o;
  }


}

#endif //  __GMM_SUB_INDEX_H
