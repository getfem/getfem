#include <dal_basic.h>
#include <deque>

struct titi {
  virtual void toto(void) { cout << "bibi\n"; }
};

struct tutu : public titi {
  virtual void toto(void) { cout << "bubu\n"; }
};


int main(void) {
  try {


    cout << "size of int           : " << sizeof(int)         << endl;
    cout << "size of size_t        : " << sizeof(size_t)      << endl;
    cout << "size of (int *)       : " << sizeof(int *)       << endl;
    cout << "size of short int     : " << sizeof(short int)   << endl;
    cout << "size of long int      : " << sizeof(long int)    << endl;
    cout << "size of long long int : " << sizeof(long long int)    << endl;
    cout << "size of char          : " << sizeof(char)        << endl;
    cout << "size of float         : " << sizeof(float)       << endl;
    cout << "size of double        : " << sizeof(double)      << endl;
    cout << "size of long double   : " << sizeof(long double) << endl;
    
    dal::dynamic_array<int, 4> t;

    try {
      t[(unsigned)(-5)] = 8;
      DAL_THROW(dal::internal_error,
		"negative index does not produce an error");
    }
    catch(std::out_of_range e) {
      cout << "Out of range error successfully catched, ok\n";
    }

    t[64] = 13;
    // cout << "capacity : (should be 80) " << t.capacity() << endl;
    if (t.capacity() != 80)
      DAL_THROW(dal::internal_error, " bad capacity");
    
    dal::dynamic_array<int, 4>::iterator itb = t.begin(), ite = t.end();
    dal::dynamic_array<int, 4>::iterator ita;
    ita = itb++;
    // cout << "range : " << (ita - t.begin()) << endl;
    
    while (itb != ite)  *itb++ = int(3);
    
    
    // std::fill(t.begin(), t.end(), int(3));
    
    // cout << "capacity : (should be 80) " << t.capacity() << endl;
    if (t.capacity() != 80)
      DAL_THROW(dal::internal_error, " bad capacity");
    // cout << "t[64] = (should be 3) " << t[64] << endl;
    if (t[64] != 3) DAL_THROW(dal::internal_error, "iterators don't work");
    
    t.clear();
    // cout << "capacity : (should be 0) " << t.capacity() << endl;
    if (t.capacity() != 0) 
      DAL_THROW(dal::internal_error, " clear does not work");
   
    std::fill(t.begin(), t.end(), int(3));
    // cout << "capacity : (should be 0) " << t.capacity() << endl;
    if (t.capacity() != 0)   
      DAL_THROW(dal::internal_error, " clear does not work");
    t[64] = 6;
    
    dal::dynamic_array<int, 4> t2, t3;
    
    t2[64] = 12;
    t3 = t2 = t;
    
    // cout << "capacity : (should be 80) " << t.capacity() << endl;
    if (t.capacity() != 80)
      DAL_THROW(dal::internal_error, " bad capacity");

    {
      dal::dynamic_array<int, 4>::const_iterator
	it1 = ((const dal::dynamic_array<int, 4> *)(&t))->begin(),
	it2 = ((const dal::dynamic_array<int, 4> *)(&t2))->begin(),
	it3 = ((const dal::dynamic_array<int, 4> *)(&t3))->begin(),
	ite = ((const dal::dynamic_array<int, 4> *)(&t))->end(),
	itb = ((const dal::dynamic_array<int, 4> *)(&t))->begin();
      
      for ( ; it1 != ite; it1++, it2++, it3++)
      {
	size_t ind = it1 - itb;
	if 
	( ( (&(*it1)) != &(t[ind]) ) ||
	  ( (&(*it2)) != &(t2[ind]) ) ||
	  ( (&(*it3)) != &(t3[ind]) ) ||
	  ( (&(*it1)) == (&(*it2)) ) ||
	  ( (&(*it2)) == (&(*it3)) ) ||
	  ( (&(*it1)) == (&(*it3)) ) )
	  DAL_THROW(dal::internal_error, " copy does not work");
      }
    }
    
    {
      dal::dynamic_array<int, 4>::iterator
	it1 = t.begin(),
	it2 = t2.begin(),
	it3 = t3.begin(),
	ite = t.end(),
	itb = t.begin();
      
      for ( ; it1 != ite; it1++, it2++, it3++)
	{
	  size_t ind = it1 - itb;
	  
	  if 
	  ( ( (&(*it1)) != &(t[ind]) ) ||
	    ( (&(*it2)) != &(t2[ind]) ) ||
	    ( (&(*it3)) != &(t3[ind]) ) ||
	    ( (&(*it1)) == (&(*it2)) ) ||
	    ( (&(*it2)) == (&(*it3)) ) ||
	    ( (&(*it1)) == (&(*it3)) ) )
	    DAL_THROW(dal::internal_error, " copy does not work");
	}
    }
    
    
    t[64] = 11;
    
    // cout << "t2[64] = (should be 6 6) " << t2[64] << " " << t3[64]<< endl;
    if (t2[64] != 6)
       DAL_THROW(dal::internal_error, " copy does not work");
    // cout << "capacity : (should be 80) " << t3.capacity() << endl;
    if (t.capacity() != 80)
      DAL_THROW(dal::internal_error, " bad capacity");

  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;
}
