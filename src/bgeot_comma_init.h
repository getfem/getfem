/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_comma_init.h : basic configuration.                    */
/*     									   */
/* Date : March 27, 2003.                                                  */
/* Author : Julien Pommier, pommier@gmm.insa-tlse.fr                       */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2003  Julien Pommier.                                     */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
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

#ifndef COMMA_INIT
#define COMMA_INIT

/* 
   highly inspired by the boost init.hpp (C) Thorsten Ottosen
   (http://www.cs.auc.dk/~nesotto/init/)    
*/

/**
 *  Template class which forwards insertions to the
 *  container class.
 */ 
namespace bgeot {

  template<typename Container> class Comma_initializer {
    typedef typename Container::value_type       value_type;
    Container& c_;
  public: 
    explicit Comma_initializer( Container& c ) : c_( c ) {}
    
    Comma_initializer& operator,(const value_type& v) {
      c_.push_back(v);
      return *this;
    }
    
    /**
     *  Should only be used with first value. The operator
     *  gives a nice syntax for initializing the container.
     */
    Comma_initializer& operator=(const value_type v) {
      c_.clear();
      c_.push_back(v);
      return *this;
    }
    /**
     *  Should only be used with first value. The operator
     *  gives a nice syntax for appending to the container.
     */
    Comma_initializer& operator+=(value_type v) {
      c_.push_back(v);
      return *this;
    }
  };

  template<typename T> Comma_initializer<T> sc( T& c )
    { return Comma_initializer<T>(c); }

}

#endif
