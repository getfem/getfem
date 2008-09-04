// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2008 Yves Renard.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

#include "gmm/gmm_std.h"
#ifndef GETFEM_VERIFY__
  #define GETFEM_VERIFY__
#endif

#include "getfem/dal_tree_sorted.h"


int main(void)
{
  try {
    dal::dynamic_tree_sorted<int> tsa;
    
    cout << tsa.add(0) << endl; cout << tsa << endl;
    cout << tsa.add(1) << endl; cout << tsa << endl;
    cout << tsa.add(6) << endl; cout << tsa << endl;
    cout << tsa.add(5) << endl; cout << tsa << endl;
    cout << tsa.add(4) << endl; cout << tsa << endl;
    cout << tsa.add(2) << endl; cout << tsa << endl;
    cout << tsa.add(3) << endl; cout << tsa << endl;
    tsa.sup(1); cout << tsa << endl;
    tsa.sup(4); cout << tsa << endl;
    
    {
      dal::dynamic_tree_sorted<int>::iterator it = tsa.begin(), end=tsa.end();
      while (it != end) cout << *it++ << " : ";
      cout << endl << endl;
    }
    
    {
      dal::dynamic_tree_sorted<int>::sorted_iterator it = tsa.sorted_begin(),
	end = tsa.sorted_end();
      while (it != end) cout << *it++ << " : ";
    }
    
    cout << "test random ... \n";
    
    tsa.clear();
    
    for (int i = 0; i < 50; i++)
      { 
	int j = rand() % 1000; 
	cout << " add no " << i << " : " << j << endl;
	
	tsa.add(j);
      }
    
    {
      dal::dynamic_tree_sorted<int>::iterator it= tsa.begin(), end = tsa.end();
      while (it != end) cout << *it++ << " : ";
      cout << endl << endl;
    }
    cout << tsa << endl;
    
    for (unsigned long i = 0; i < 50 /*50000*/; i++)
    {
      if (!(i % 10000)) cout <<" no " << i << " nb_elt " << tsa.card() << endl;
      
      if (((rand() & 1) == 1 || tsa.card() < 30) && (tsa.card() < 100))
	{
	  int j = rand() % 10;
	  tsa.add(j);
	  // if (!(i % 10)) cout << tsa << endl;
	  tsa.verify_balance();
	}
      else
	{
	  
	  // cout << " sup of " << k << endl;
	  int k;
	  for (k = int(rand() % (tsa.index()).last_true());
	       !((tsa.index())[k]); 
	       k = int(rand() % (tsa.index()).last_true()) ) {};
	  
	  tsa.sup(k);
	  //    cout << tsa << endl;
	  tsa.verify_balance();
	}
    }
    
    {
      dal::dynamic_tree_sorted<int>::iterator it= tsa.begin(), end = tsa.end();
      while (it != end) cout << *it++ << " : ";
      cout << endl << endl;
    }

    {
      dal::dynamic_tree_sorted<std::string> v;
      for (unsigned i=0; i < 1000; ++i) {
	v.add("toto"); v.add("hop"); v.add_norepeat("grr");
	std::string s; 
	for (int j=0; j < 5; ++j) s.push_back(char('a' + (rand() % 26)));
	v.add(s);
      }
      {
	dal::dynamic_tree_sorted<std::string>::iterator it= v.begin(), end = v.end();
	while (it != end) cout << *it++ << " : ";
	cout << endl << endl;
      }
    }
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0;

}
