#include <dal_bit_vector.h>


int main(void)
{
  try {
  dal::bit_vector nn;

  cout << "cardinal : " << nn.card() << endl; assert(nn.card() == 0);
  cout << "first true : " << nn.first_true() << endl;
  assert(nn.first_true() == 1);
  cout << "last  true : " << nn.last_true()  << endl;
  assert(nn.last_true() == 0);
  cout << "first false : " << nn.first_false() << endl;
  assert(nn.first_false() == 0);
  cout << "last  false : " << nn.last_false()  << endl;
  assert(nn.last_false() == 0);
  cout << "first true : " << nn.first_true() << endl;
  assert(nn.first_true() == 1);
  cout << "last  true : " << nn.last_true()  << endl;
  assert(nn.last_true() == 0);
  cout << "first false : " << nn.first_false() << endl;
  assert(nn.first_false() == 0);
  cout << "last  false : " << nn.last_false()  << endl;
  assert(nn.last_false() == 0);


  nn.add(5);
  cout << "cardinal : " << nn.card() << endl; assert(nn.card() == 1);
  cout << nn << endl;

  cout << "first true : " << nn.first_true() << endl;
  assert(nn.first_true() == 5);
  cout << "last  true : " << nn.last_true()  << endl;
  assert(nn.last_true() == 5);
  cout << "first false : " << nn.first_false() << endl;
  assert(nn.first_false() == 0);
  cout << "last  false : " << nn.last_false()  << endl;
  assert(nn.last_false() == 4);
  cout << "first true : " << nn.first_true() << endl;
  assert(nn.first_true() == 5);
  cout << "last  true : " << nn.last_true()  << endl;
  assert(nn.last_true() == 5);
  cout << "first false : " << nn.first_false() << endl;
  assert(nn.first_false() == 0);
  cout << "last  false : " << nn.last_false()  << endl;
  assert(nn.last_false() == 4);


  cout << "nn[5] = " << nn[5] << endl; assert(nn[5] == true);

  nn.add(1); nn.add(10); nn.add(5); nn.add(11);
  cout << "cardinal : " << nn.card() << endl; assert(nn.card() == 4);

  cout << nn << endl;
  cout << "first element : " << nn.first() << endl; assert(nn.first() == 1);
  cout << "last  element : " << nn.last()  << endl; assert(nn.last() == 11);

  nn.sup(8,11); nn.add(3);
  cout << nn << endl;
  cout << "first element : " << nn.first() << endl; assert(nn.first() == 1);
  cout << "last  element : " << nn.last()  << endl; assert(nn.last() == 5);

  dal::bit_vector mm;

  nn.add(31); nn.add(32); nn.add(33);
  mm.add(31); mm.add(35); mm.add(36);
  mm.add(2); mm.add(6);

  nn &= mm;

  cout << nn << endl;
  cout << "first element : " << nn.first() << endl; assert(nn.first() == 31);
  cout << "last  element : " << nn.last()  << endl; assert(nn.last() == 31);

  nn = mm;

  mm.add(28);

  cout << nn << endl;
  cout << "first element : " << nn.first() << endl; assert(nn.first() == 2);
  cout << "last  element : " << nn.last()  << endl; assert(nn.last() == 36);
  cout << "card = " << nn.card() << endl; assert(nn.card() == 5);

  nn.add(64);
  cout << nn << endl;

  nn.add(256);
  cout << nn << endl;

  nn.add(512);
  cout << nn << endl;

  nn.add(1024);
  cout << nn << endl;

  cout << "card = " << nn.card() << endl; assert(nn.card() == 9);

  nn &= mm;
  cout << nn << endl;
  cout << "card = " << nn.card() << endl;

  nn.swap(1024, 36);
  cout << nn << endl;
  cout << "card = " << nn.card() << endl;
  cout << "first element : " << nn.first() << endl; assert(nn.first() == 2);
  cout << "last  element : " << nn.last()  << endl; assert(nn.last() == 1024);

  nn.swap(1024, 36);
  cout << nn << endl;
  cout << "card = " << nn.card() << endl;
  cout << "first element : " << nn.first() << endl; assert(nn.first() == 2);
  cout << "last  element : " << nn.last()  << endl; assert(nn.last() == 36);

  nn.swap(4096, 36);
  cout << nn << endl;
  cout << "card = " << nn.card() << endl;
  cout << "first element : " << nn.first() << endl; assert(nn.first() == 2);
  cout << "last  element : " << nn.last()  << endl; assert(nn.last() == 4096);

  mm |= nn;
  nn = mm;
  cout << nn << endl;

  nn.clear();
  cout << nn << endl;
  nn.add(2048);
  cout << nn << endl;

  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;

}
