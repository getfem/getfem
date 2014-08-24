//  Boost noncopyable.hpp header file  --------------------------------------//

//
//  Copyright (c) 1999-2003 Beman Dawes
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
//  See http://www.boost.org/libs/utility for documentation.
//

#if !defined(BOOST_NONCOPYABLE_HPP_INCLUDED) && !defined(BOOST_NONCOPYABLE_HPP) && !defined(BOOST_CORE_NONCOPYABLE_HPP)

#define BOOST_NONCOPYABLE_HPP_INCLUDED
#define BOOST_NONCOPYABLE_HPP
#define BOOST_CORE_NONCOPYABLE_HPP



namespace boost {

//  Private copy constructor and copy assignment ensure classes derived from
//  class noncopyable cannot be copied.

//  Contributed by Dave Abrahams

namespace noncopyable_  // protection from unintended ADL
{
  class noncopyable
  {
   protected:
      noncopyable() {}
      ~noncopyable() {}
   private:  // emphasize the following members are private
      noncopyable( const noncopyable& );
      const noncopyable& operator=( const noncopyable& );
  };
}

typedef noncopyable_::noncopyable noncopyable;

} // namespace boost

#endif
