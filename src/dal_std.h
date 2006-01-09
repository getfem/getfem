// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Dynamic Array Library (dal)
// File    : dal_std.h : Compatibility Header.
//           
// Date    : June 01, 1995.
// Authors : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1995-2006 Yves Renard
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

/**@file dal_std.h
   @brief basic setup for dal (includes, typedefs etc.)
*/
#ifndef DAL_STD_H__
#define DAL_STD_H__

#ifndef NOGETFEM_VERIFY
# define GETFEM_VERIFY
#endif

#ifndef __USE_STD_IOSTREAM
# define __USE_STD_IOSTREAM
#endif

#ifndef __USE_BSD
# define __USE_BSD
#endif

#ifndef __USE_ISOC99
# define __USE_ISOC99
#endif

#if !defined(GMM_USES_MPI) && GETFEM_PARA_LEVEL > 0
# define GMM_USES_MPI
#endif

/* ********************************************************************** */
/*	Compilers detection.						  */
/* ********************************************************************** */

/* for sun CC 5.0 ...
#if defined(__SUNPRO_CC) && __SUNPRO_CC >= 0x500
# include <stdcomp.h>
# undef _RWSTD_NO_CLASS_PARTIAL_SPEC
# undef _RWSTD_NO_NAMESPACE
#endif 
*/
/* for VISUAL C++ ...
   #if defined(_MSC_VER) //  && !defined(__MWERKS__)
   #define _GETFEM_MSVCPP_ _MSC_VER
   #endif
*/

#if defined(__GNUC__)
#  if (__GNUC__ < 3)
#    error : PLEASE UPDATE g++ TO AT LEAST 3.0 VERSION
#  endif
#endif

/* ********************************************************************** */
/*	C++ Standard Headers.						  */
/* ********************************************************************** */
#include <cstdlib>
#include <cstddef>
#include <cmath>
#include <cstring>
#include <cctype>
#include <cassert>
#include <iostream>
//#include <ios> 
#include <fstream>
#include <ctime>
#include <exception>
#include <typeinfo>
#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <vector>
#include <deque>
#include <string>
#include <complex>
#include <limits>
#include <sstream>
#include <numeric>


using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

namespace dal {

  /* ******************************************************************* */
  /*       Clock functions.                                              */
  /* ******************************************************************* */
  
# if  defined(HAVE_SYS_TIMES)
  inline double uclock_sec(void) {
    static double ttclk = 0.;
    if (ttclk == 0.) ttclk = sysconf(_SC_CLK_TCK);
    tms t; times(&t); return double(t.tms_utime) / ttclk;
  }
# else
  inline double uclock_sec(void)
  { return double(clock())/double(CLOCKS_PER_SEC); }
# endif
  
  /* ******************************************************************** */
  /*	Fixed size integer types.                     			  */
  /* ******************************************************************** */
  // Remark : the test program dynamic_array tests the lenght of
  //          resulting integers

  template <int s> struct fixed_size_integer_generator {
    typedef void int_base_type;
    typedef void uint_base_type;  
  };

  template <> struct fixed_size_integer_generator<sizeof(char)> {
    typedef signed char int_base_type;
    typedef unsigned char uint_base_type;
  };

  template <> struct fixed_size_integer_generator<sizeof(short int)
    - ((sizeof(short int) == sizeof(char)) ? 78 : 0)> {
    typedef signed short int int_base_type;
    typedef unsigned short int uint_base_type;
  };

  template <> struct fixed_size_integer_generator<sizeof(int)
    - ((sizeof(int) == sizeof(short int)) ? 59 : 0)> {
    typedef signed int int_base_type;
    typedef unsigned int uint_base_type;
  };
 
  template <> struct fixed_size_integer_generator<sizeof(long)
    - ((sizeof(int) == sizeof(long)) ? 93 : 0)> {
    typedef signed long int_base_type;
    typedef unsigned long uint_base_type;
  };

  template <> struct fixed_size_integer_generator<sizeof(long long)
    - ((sizeof(long long) == sizeof(long)) ? 99 : 0)> {
    typedef signed long long int_base_type;
    typedef unsigned long long uint_base_type;
  };
 
  typedef fixed_size_integer_generator<1>::int_base_type int8_type;
  typedef fixed_size_integer_generator<1>::uint_base_type uint8_type;
  typedef fixed_size_integer_generator<2>::int_base_type int16_type;
  typedef fixed_size_integer_generator<2>::uint_base_type uint16_type;
  typedef fixed_size_integer_generator<4>::int_base_type int32_type;
  typedef fixed_size_integer_generator<4>::uint_base_type uint32_type;
  typedef fixed_size_integer_generator<8>::int_base_type int64_type;
  typedef fixed_size_integer_generator<8>::uint_base_type uint64_type;

// #if INT_MAX == 32767
//   typedef signed int    int16_type;
//   typedef unsigned int uint16_type;
// #elif  SHRT_MAX == 32767
//   typedef signed short int    int16_type;
//   typedef unsigned short int uint16_type;
// #else
// # error "impossible to build a 16 bits integer"
// #endif

// #if INT_MAX == 2147483647
//   typedef signed int    int32_type;
//   typedef unsigned int uint32_type;
// #elif  SHRT_MAX == 2147483647
//   typedef signed short int    int32_type;
//   typedef unsigned short int uint32_type;
// #elif LONG_MAX == 2147483647
//   typedef signed long int    int32_type;
//   typedef unsigned long int uint32_type;
// #else
// # error "impossible to build a 32 bits integer"
// #endif

// #if INT_MAX == 9223372036854775807L || INT_MAX == 9223372036854775807
//   typedef signed int    int64_type;
//   typedef unsigned int uint64_type;
// #elif LONG_MAX == 9223372036854775807L || LONG_MAX == 9223372036854775807
//   typedef signed long int    int64_type;
//   typedef unsigned long int uint64_type;
// #elif LLONG_MAX == 9223372036854775807LL || LLONG_MAX == 9223372036854775807L || LLONG_MAX == 9223372036854775807
//   typedef signed long long int int64_type;
//   typedef unsigned long long int uint64_type;
// #else
// # error "impossible to build a 64 bits integer"
// #endif

#ifdef __GNUC__
/* 
   g++ can issue a warning at each usage of a function declared with this special attribute 
   (also works with typedefs and variable declarations)
*/
# define IS_DEPRECATED __attribute__ ((__deprecated__))
/*
   the specified function is inlined at any optimization level 
*/
# define ALWAYS_INLINE __attribute__((always_inline))
#else
# define IS_DEPRECATED
# define ALWAYS_INLINE
#endif

  /* Pour forcer l'instanciation dans libgetfem de tous les
     type qu'on est suceptible d' afficher. Si il en manque un (c'est
     le cas pour unsigned et float) dans libgetfem.so qui est utilise
     dans libgetfemint.so alors on retrouve un warning
     comme quoi le symbole T_Q13std8ios_base__ et T_Q13std11logic_error__
     sont dupliqués
     visiblement il reinstancie un peu trop de choses .. 

     pour s'assurer que le probleme n'est pas de retour, faire
     nm libgetfem.so | grep 'td::basic_ostream<char,
                              std::char_traits<char> >::operator <<'
     et faire de même dans tous les .o de matlabint, et s'assurer que
     tous ceux de matlabint sont bien inclus dans libgetfem.so

     enfin c'est pas tout a fait ça puisque tous ces trucs sont deja
     instancies dans libcxxstd.a, et c'est le noeud du probleme j'imagine:
     nm /usr/lib/cmplrs/cxx/V6.3-008/libcxxstd.a
       | grep 'td::basic_ostream<char, std::char_traits<char> >::operator <<' 
     ---> std::basic_ostream<char, std::char_traits<char> >::operator
            <<(long long) | 0000000000006368 | T | 0000000000000008
     std::basic_ostream<char, std::char_traits<char> >::operator
         <<(const void*) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(std::basic_ostream<char, std::char_traits<char> >&
	 (*)(std::basic_ostream<char, std::char_traits<char> >&))
	 | 0000000000000000 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(unsigned long long) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(unsigned int) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(unsigned long) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(unsigned short) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator 
         <<(double) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(float) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(int) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(bool) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(long double128) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(long) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(short) | 0000000000006368 | T | 0000000000000008
    std::basic_ostream<char, std::char_traits<char> >::operator
         <<(long double64) | 0000000000006368 | T | 0000000000000008
  */
  struct just_for_the_fine_cxx {
    static void f() {
      long double z(1.0);
      std::stringstream s;
      s << int(1) << double(2.0) 
	<< "hello" << std::string("hello") << unsigned(1) << float(2) 
	<< char('a') << (unsigned char)('b') << short(1)
	<< (unsigned short)(2)  << long(1) << (unsigned long)(2)
	<< (const void*)NULL << bool(1) << z;
    }
  };
  
} /* end of namespace dal.                                                */

#endif /* DAL_STD_H__ */

