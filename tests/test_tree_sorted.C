
#include <dal_std.h>
#ifndef __GETFEM_VERIFY
  #define __GETFEM_VERIFY
#endif

#include <dal_tree_sorted.h>


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
    
    for (unsigned long i = 0; i < 50000; i++)
    {
      if (!(i % 10000)) cout <<" no " << i << " nb_elt " << tsa.card() << endl;
      
      if ((rand() & 1 == 1 || tsa.card() < 30) && (tsa.card() < 100))
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
	  for (k = rand() % (tsa.index()).last_true(); !((tsa.index())[k]); 
	       k = rand() % (tsa.index()).last_true() ) ;
	  
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
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;

}
