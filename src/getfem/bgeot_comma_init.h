/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2003-2012 Julien Pommier
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
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

/**@file bgeot_comma_init.h
   @author  Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date March 27, 2003.
   @brief convenient initialization of vectors via overload of "operator,".
   @code
   std::vector<double> foo; bgeot::sc(foo) += 1,3,4;
   @endcode

   highly inspired by the boost init.hpp (C) Thorsten Ottosen
   (http://www.cs.auc.dk/~nesotto/init/)    
*/
#ifndef COMMA_INIT
#define COMMA_INIT

namespace bgeot {

  /**
   *  Template class which forwards insertions to the
   *  container class.
   */ 
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
