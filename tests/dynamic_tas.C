#include <dal_tas.h>

using namespace std;

int main(void)
{
  dal::dynamic_tas<int> t;
  t.add(6);
  cout << "index = " << t.index() << endl;
  cout << "first free place : " << t.index().first_false() << endl;
  assert(t.index().first_false() == 1);

  t.add(4);
  t.add(2);

  cout << "index = " << t.index() << endl;

  t.sup(1);
  cout << "index = " << t.index() << endl;
  cout << "first free place : " << t.index().first_false() << endl;
  assert(t.index().first_false() == 1);


  t.add(3);

  dal::dynamic_tas<int>::iterator it = t.begin(), end = t.end();

  cout << " card = " << t.card() << endl;

  cout << "index = " << t.index() << endl;

  cout << "index.first() = " << t.index().first() << endl;

  for ( ; it != end; ++it) { cout << " : " << *it; }
  cout << endl;
  
  return 0;

}
