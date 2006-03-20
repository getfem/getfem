/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
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
#include <dal_bit_vector.h>
#include <deque>
typedef size_t size_type;

bool quick = false;

void test_speed_part(dal::bit_vector& bv, std::vector<bool>& vb, std::deque<bool>& db) {
  double t0;
  t0=dal::uclock_sec();
  size_type cnt=0, oldcnt;
  for (size_type i=0; i < vb.size(); ++i) {
    cnt += bv[i] ? 1 : 0;
  }
  cerr << "bit_vector   [random_access ]: " << dal::uclock_sec()-t0 << endl; t0 = dal::uclock_sec();

  oldcnt=cnt; cnt = 0;
  for (size_type i=0; i < vb.size(); ++i) {
    cnt += bv.is_in(i) ? 1 : 0;
  }
  cerr << "bit_vector   [is_in         ]: " << dal::uclock_sec()-t0 << endl; t0 = dal::uclock_sec(); assert(oldcnt == cnt);

  oldcnt=cnt; cnt = 0;
  for (size_type i=0; i < vb.size(); ++i) {
    cnt += vb[i] ? 1 : 0;
  }
  cerr << "vector<bool> [random_access ]: " << dal::uclock_sec()-t0 << endl; t0 = dal::uclock_sec(); assert(oldcnt == cnt);

  oldcnt=cnt; cnt = 0;
  for (size_type i=0; i < db.size(); ++i) {
    cnt += db[i] ? 1 : 0;
  }
  cerr << "deque <bool> [random_access ]: " << dal::uclock_sec()-t0 << endl; t0 = dal::uclock_sec(); assert(oldcnt == cnt);

  oldcnt=cnt; cnt = 0;
  for (dal::bit_vector::const_iterator it = bv.begin(); it != bv.end(); ++it) {
    cnt += (*it) ? 1 : 0;
  }
  cerr << "bit_vector   [const_iterator]: " << dal::uclock_sec()-t0 << endl; t0 = dal::uclock_sec(); assert(oldcnt == cnt);

  oldcnt=cnt; cnt = 0;
  for (dal::bv_visitor i(bv); !i.finished(); ++i) {
    cnt ++; 
  }
  cerr << "bit_vector   [bv_visitor    ]: " << dal::uclock_sec()-t0 << endl; t0 = dal::uclock_sec();  assert(oldcnt == cnt);

  oldcnt=cnt; cnt = 0;
  for (size_type i = bv.take_first(); i != size_type(-1); i << bv) {
    cnt ++; 
  }
  cerr << "bit_vector   [operator <<   ]: " << dal::uclock_sec()-t0 << endl; t0 = dal::uclock_sec(); assert(oldcnt == cnt);

  oldcnt=cnt; cnt = 0;
  for (std::vector<bool>::const_iterator it = vb.begin(); it != vb.end(); ++it) {
    cnt += (*it) ? 1 : 0;
  }
  cerr << "vector<bool> [const_iterator]: " << dal::uclock_sec()-t0 << endl; t0 = dal::uclock_sec(); assert(oldcnt == cnt);

  oldcnt=cnt; cnt = 0;
  for (std::deque<bool>::const_iterator it = db.begin(); it != db.end(); ++it) {
    cnt += (*it) ? 1 : 0;
  }
  cerr << "deque <bool> [const_iterator]: " << dal::uclock_sec()-t0 << endl; t0 = dal::uclock_sec(); assert(oldcnt == cnt);
}

void test_speed() {
  size_type N = quick ? 1000000 : 10000000;
  dal::bit_vector   bv; bv.add(0,N);
  std::vector<bool> vb(N,true);
  std::deque<bool> db(N,true);

  dal::bit_vector   bv2;
  std::vector<bool> vb2;
  std::deque<bool> db2;
  double t0 = dal::uclock_sec();
  for (size_type i=0; i < N; ++i) bv2[i] = 1;
  cerr << "bit_vector   [push_back]: " << dal::uclock_sec()-t0 << endl; t0 = dal::uclock_sec();
  for (size_type i=0; i < N; ++i) vb2.push_back(true);
  cerr << "vector<bool> [push_back]: " << dal::uclock_sec()-t0 << endl; t0 = dal::uclock_sec();
  for (size_type i=0; i < N; ++i) db2.push_back(true);
  cerr << "deque<bool>  [push_back]: " << dal::uclock_sec()-t0 << endl; t0 = dal::uclock_sec();


  cerr << "---- full vector ----" << endl;
  test_speed_part(bv,vb,db);

  bv.sup(0,N);
  std::fill(vb.begin(),vb.end(),false);
  std::fill(db.begin(),db.end(),false);
  
  cerr << "---- empty vector ----" << endl;
  test_speed_part(bv,vb,db);
  
  for (size_type i=0; i < N; ++i) {
    bool b = (rand() % 2);
    bv[i] = vb[i] = db[i] = b;
  }

  cerr << "---- random vector (50% true) ----" << endl;
  test_speed_part(bv,vb,db);

  for (size_type i=0; i < N; ++i) {
    bool b = (rand() % 16) ? 0 : 1;
    bv[i] = vb[i] = db[i] = b;
  }

  cerr << "---- random vector (6% true) ----" << endl;
  test_speed_part(bv,vb,db);

  for (size_type i=0; i < N; ++i) {
    bool b = (rand() % 100) ? 0 : 1;
    bv[i] = vb[i] = db[i] = b;
  }

  cerr << "---- random vector (1% true) ----" << endl;
  test_speed_part(bv,vb,db);

  for (size_type i=0; i < N; ++i) {
    bool b = (rand() % 1000) ? 0 : 1;
    bv[i] = vb[i] = db[i] = b;
  }

  cerr << "---- random vector (0.1% true) ----" << endl;
  test_speed_part(bv,vb,db);


}



int main(int argc, char *argv[]) {
  try {
  dal::bit_vector nn;

  if (argc == 2 && strcmp(argv[1],"-quick")==0) quick = true;

  cout << "cardinal : " << nn.card() << endl; assert(nn.card() == 0);
  cout << "first true : " << nn.first_true() << endl;
  //assert(nn.first_true() == 1);
  assert(nn.first_true() == size_type(-1));
  cout << "last  true : " << nn.last_true()  << endl;
  //assert(nn.last_true() == 0);
  assert(nn.last_true() == size_type(-1));
  cout << "first false : " << nn.first_false() << endl;
  assert(nn.first_false() == 0);
  cout << "last  false : " << nn.last_false()  << endl;
  assert(nn.last_false() == 0);
  cout << "first true : " << nn.first_true() << endl;
  //assert(nn.first_true() == 1);
  assert(nn.first_true() == size_type(-1));
  cout << "last  true : " << nn.last_true()  << endl;
  //assert(nn.last_true() == 0);
  assert(nn.last_true() == size_type(-1));
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

  nn.sup(8,4); nn.add(3);
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

  dal::bit_vector u, v, w, z;
  u.add(0, 216);
  v.add(6, 634);
  u &= v;
  u &= w;
  z &= v; return 0;

  test_speed();

  }
  DAL_STANDARD_CATCH_ERROR;
  return 0;
}
