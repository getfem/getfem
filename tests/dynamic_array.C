#include <dal_basic.h>
#include <deque>


int main(void)
{
  
  cout << "size of int : " << sizeof(int) << endl;
  cout << "size of size_t : " << sizeof(size_t) << endl;
  cout << "size of (int *) : " << sizeof(int *) << endl;
  cout << "size of short int : " << sizeof(short int) << endl;
  cout << "size of long int : " << sizeof(long int) << endl;
  cout << "size of char : " << sizeof(char) << endl;
  cout << "size of float : " << sizeof(float) << endl;
  cout << "size of double : " << sizeof(double) << endl;
  cout << "size of long double : " << sizeof(long double) << endl;

  dal::dynamic_array<int, 4> t;

  t[64] = 13;
  cout << "capacity : (should be 80) " << t.capacity() << endl;
  assert(t.capacity() == 80);

  dal::dynamic_array<int, 4>::iterator itb = t.begin(), ite = t.end();
  dal::dynamic_array<int, 4>::iterator ita;
  ita = itb++;
  cout << "range : " << (ita - t.begin()) << endl;

  while (itb != ite)
  {
    // cout << (itb - t.begin()) << " : " << endl;
    *itb++ = int(3);
    //   itb++;
    // cout << (itb - t.begin()) << endl; 
  }

  // std::fill(t.begin(), t.end(), int(3));

  cout << "capacity : (should be 80) " << t.capacity() << endl;
  assert(t.capacity() == 80);
  cout << "t[64] = (should be 3) " << t[64] << endl;
  assert(t[64] == 3);


  t.clear();
  cout << "capacity : (should be 0) " << t.capacity() << endl;
  assert(t.capacity() == 0);
  std::fill(t.begin(), t.end(), int(3));

  cout << "capacity : (should be 0) " << t.capacity() << endl;
  assert(t.capacity() == 0);

  t[64] = 6;
 
  dal::dynamic_array<int, 4> t2, t3;
    
  t2[64] = 12;
  t3 = t2 = t;

  cout << "capacity : (should be 80) " << t.capacity() << endl;
  assert(t.capacity() == 80);

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
//       cout << "ind : " << ind
// 	   << " pointers : " << &(*it1) << " : " << &(*it2)
// 	   << " : " << &(*it3) << endl;
      
      assert( (&(*it1)) == &(t[ind]) );
      assert( (&(*it2)) == &(t2[ind]) );
      assert( (&(*it3)) == &(t3[ind]) );
      assert ( (&(*it1)) != (&(*it2)) );
      assert ( (&(*it2)) != (&(*it3)) );
      assert ( (&(*it1)) != (&(*it3)) );
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
   //    cout << "ind : " << ind
// 	   << " pointers : " << &(*it1) << " : " << &(*it2)
// 	   << " : " << &(*it3) << endl;
      
      assert( (&(*it1)) == &(t[ind]) );
      assert( (&(*it2)) == &(t2[ind]) );
      assert( (&(*it3)) == &(t3[ind]) );
      assert ( (&(*it1)) != (&(*it2)) );
      assert ( (&(*it2)) != (&(*it3)) );
      assert ( (&(*it1)) != (&(*it3)) );
    }
  }

  


  t[64] = 11;
  
  cout << "t2[64] = (should be 6 6) " << t2[64] << " " << t3[64] << endl;
  assert(t2[64] == 6); assert(t3[64] == 6);
  cout << "capacity : (should be 80) " << t3.capacity() << endl;
  assert(t.capacity() == 80);
  
  

  
  return 0;

}
