#ifndef COMMA_INIT
#define COMMA_INIT

/* 
   highly inspired by the boost init.hpp (C) Thorsten Ottosen (http://www.cs.auc.dk/~nesotto/init/)    
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
  template<typename T> Comma_initializer<T> sc( T& c ) { return Comma_initializer<T>(c); }
}

#endif
