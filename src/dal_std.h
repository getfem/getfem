/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Dynamic Array Library (dal)                                  */
/* File    :  dal_std.h : Compatibility Header.                            */
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


#ifndef __DAL_STD_H
#define __DAL_STD_H

#define __GETFEM_VERSION 1
#define __GETFEM_REVISION 0

#define __DAL_VERSION 1
#define __DAL_REVISION 0

#define __GETFEM_VERIFY

#ifndef __USE_STD_IOSTREAM
  #define __USE_STD_IOSTREAM
#endif

/* ********************************************************************** */
/*	C++ Standard Headers.						  */
/* ********************************************************************** */

#include <stddef.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <ctype.h>
#include <assert.h>
#include <limits.h>

/* ********************************************************************** */
/*	S.T.L. Headers.						          */
/* ********************************************************************** */

#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <vector>
#include <string>
#include <strstream>

using std::endl;
using std::cout;
using std::cerr;
using std::ends;
using std::cin;

/* ********************************************************************** */
/*	S.T.L. Reverse iterator definition.				  */
/* ********************************************************************** */

#if defined(_MSC_VER) && !defined(__MWERKS__)
#define _GETFEM_MSVCPP_ _MSC_VER
#endif
#if !defined ( _GETFEM_MSVCPP_ )
#define GETFEM_REVERSE_ITER 1
#else
#define GETFEM_REVERSE_ITER 0
#endif

/** Dynamic Array Library. \\
 *   The Dynamic Array Library (dal) is a library of containers
 *   and algorithms on containers. Thoses containers are initialy adapted
 *   to store data for meshes, but they are also of larger
 *   interest. The library is build on a container which is \\ \\
 *   dal::dynamic\_array$<$T$>$. \\ \\
 *   This container is very similar to std::vector$<$T$>$ of the Standard
 *   Template Library. The major difference
 *   is that memory is allocated by block and
 *   the allocation is automatic when an acces to a non existing element 
 *   is called. \\ \\
 *   Others containers like dal::dynamic\_tas$<$T$>$ or
 *   dal::dynamic\_tree\_sorted$<$T$>$ allow to add or delete elements
 *   in an array. Particularily, dal::dynamic\_tree\_sorted$<$T$>$ combines
 *   the random access of an array and the logarithmic search and insertion
 *   in a balanced sorted tree. \\ \\
 *   The file dal_std.h loads the very standard c++ header files and solve
 *    as much
 *     as possible the incompatibility between differents
 *     configurations. \\ \\
 *     - It assures the existence of the functions \\ \\
 *       template$<$class T$>$ std::abs(T) \\
 *       template$<$class T$>$ std::sqr(T) \\ \\
 *     - It defines fixed size type of integers : \\ \\
 *       int8\_type; uint8\_type; \\
 *       int16\_type; uint16\_type; \\
 *       int32\_type; uint32\_type; \\ \\
 *     - The macro
 *       {\tt __GETFEM_VERIFY -CODE }
 *       allows to switch on or off verifications on libraries,
 *       such as verification of range on arrays, on vectors ...
 */
namespace dal
{

#if GETFEM_REVERSE_ITER
  template <class Iter>
    class reverse_iter : public std::reverse_iterator<Iter>
  {
    typedef std::reverse_iterator<Iter> super; 
#else
  template <class Iter>
    class reverse_iter : public std::reverse_iterator<Iter,
                                           typename Iter::value_type,
                           typename Iter::reference, typename Iter::pointer>
  {
    typedef std::reverse_iterator<Iter, typename Iter::value_type,
      typename Iter::reference, typename Iter::pointer> super; 
#endif
  public:
    typedef typename super::value_type value_type;

#if defined(_GETFEM_MSVCPP_)
    typedef typename super::distance_type difference_type;
    typedef difference_type distance_type;
    typedef typename super::reference_type reference;
#else
    typedef typename super::difference_type difference_type;
    typedef typename super::reference reference;
#endif

    typedef typename super::iterator_category iterator_category;

    inline reverse_iter() {}

    inline reverse_iter(const reverse_iter& x) : super(x) { }    

    inline explicit reverse_iter(Iter x) : super(x) {}

#if GETFEM_REVERSE_ITER
  };
#else
  };
#endif
}

/* ********************************************************************** */
/*	Math functions.                     			          */
/* ********************************************************************** */

namespace dal
{
  template <class T> inline T sqr(T a) { return a * a; }
  template <class T> inline T abs(T a) { return (a < 0) ? -a : a; }
  template <class T> inline T pos(T a) { return (a < 0) ? 0 : a; }
  template <class T> inline T neg(T a) { return (a < 0) ? -a : 0; }
}


#ifndef M_PI   
# define	M_E		2.7182818284590452354       /* e          */
# define	M_LOG2E		1.4426950408889634074       /* 1/ln(2)    */
# define	M_LOG10E	0.43429448190325182765      /* 1/ln(10)   */
# define	M_LN2		0.69314718055994530942      /* ln(2)      */
# define	M_LN10		2.30258509299404568402      /* ln(10)     */
# define	M_PI		3.14159265358979323846      /* pi         */
# define	M_PI_2		1.57079632679489661923      /* pi/2       */
# define	M_PI_4		0.78539816339744830962      /* pi/4       */
# define	M_1_PI		0.31830988618379067154      /* 1/pi       */
# define	M_2_PI		0.63661977236758134308      /* 2/pi       */
# define	M_2_SQRTPI	1.12837916709551257390      /* 2/sqrt(pi) */
# define	M_SQRT2		1.41421356237309504880      /* sqrt(2)    */
# define	M_SQRT1_2	0.70710678118654752440      /* sqrt(2)/2  */

# define M_PIl       3.1415926535897932384626433832795029L  /* pi         */
# define M_PI_2l     1.5707963267948966192313216916397514L  /* pi/2       */
# define M_PI_4l     0.7853981633974483096156608458198757L  /* pi/4       */
# define M_1_PIl     0.3183098861837906715377675267450287L  /* 1/pi       */
# define M_2_PIl     0.6366197723675813430755350534900574L  /* 2/pi       */
# define M_2_SQRTPIl 1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
#endif

/* ********************************************************************** */
/*	Fixed size integer types.                     			  */
/* ********************************************************************** */

namespace dal
{

typedef signed char    int8_type;
typedef unsigned char uint8_type;

#if INT_MAX == 32767
  typedef signed int    int16_type;
  typedef unsigned int uint16_type;
#elif  SHRT_MAX == 32767
  typedef signed short int    int16_type;
  typedef unsigned short int uint16_type;
#else
# error "impossible to build a 16bits integer"
#endif

#if INT_MAX == 2147483647
  typedef signed int    int32_type;
  typedef unsigned int uint32_type;
#elif  SHRT_MAX == 2147483647
  typedef signed short int    int32_type;
  typedef unsigned short int uint32_type;
#elif LONG_MAX == 2147483647
  typedef signed long int    int32_type;
  typedef unsigned long int uint32_type;
#else
# error "impossible to build a 32bits integer"
#endif

#if INT_MAX == 9223372036854775807L || INT_MAX == 9223372036854775807
  typedef signed int    int64_type;
  typedef unsigned int uint64_type;
#elif LONG_MAX == 9223372036854775807L || INT_MAX == 9223372036854775807
  typedef signed long int    int64_type;
  typedef unsigned long int uint64_type;
#else

  // struct int64_type
  // {
  //   uint16_type p[4];

  //  int64_type &operator +=(const int64_type &i)
  //  { 
  //    uint32_type a = p[0]; a += i.p[0]; p[0] = a; a >>= 16;
  //    a += p[1]; a += i.p[1]; p[1] = a; a >>= 16;
  //    a += p[2]; a += i.p[2]; p[2] = a; a >>= 16;
  //    a += p[3]; a += i.p[3]; p[3] = a; a >>= 16;
  //  }

    /* a continuer, passage au complement a 2 pour la soustraction        */
    /* multiplication, transtypage, ...                                   */
    /* a faire si besoin.                                                 */

  // };

  // struct uint64_type : int64_type
  // {
    /* regle differente d'ecriture et de transformation en int ...        */


  // };

#endif

/* ********************************************************************** */
/*	Getfem++ generic errors.                     			  */
/* ********************************************************************** */

  /* errors definitions  */

  class dimension_error : public std::logic_error {
  public:
    dimension_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  class internal_error : public std::logic_error {
  public:
    internal_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  class failure_error : public std::logic_error {
  public:
    failure_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  class not_linear_error : public std::logic_error {
  public:
    not_linear_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  class to_be_done_error : public std::logic_error {
  public:
    to_be_done_error(const std::string& what_arg): std::logic_error (what_arg)
      { }
  };

  #define DAL_STANDARD_CATCH_ERROR   catch(std::logic_error e) \
  { \
    cerr << "============================================\n";\
    cerr << "|      An error has been detected !!!      |\n";\
    cerr << "============================================\n";\
    cerr << e.what() << endl << endl;\
    exit(1);\
  }\
  catch(std::runtime_error e)\
  {\
    cerr << "============================================\n";\
    cerr << "|      An error has been detected !!!      |\n";\
    cerr << "============================================\n";\
    cerr << e.what() << endl << endl;\
    exit(1);\
  }\
  catch(std::bad_alloc) {\
    cerr << "============================================\n";\
    cerr << "|  A bad allocation has been detected !!!  |\n";\
    cerr << "============================================\n";\
  }\
  catch(...) {\
    cerr << "============================================\n";\
    cerr << "|  An unknown error has been detected !!!  |\n";\
    cerr << "============================================\n";\
  }

  #define DAL_THROW(type, thestr) { \
    std::strstream msg; \
    msg << "in "__FILE__ << ", line " << __LINE__ << ": \n" << thestr << ends;\
    throw (type)(msg.str()); \
  }

} /* end of namespace dal.                                                */



#endif /* __DAL_STD_H */
