#include <bgeot_abstract_linalg.h>


int main(void)
{
  try {
    bgeot::vsvector<double> v(10), w(10);
    v.fill(1.0);
    w.fill(0.0);
    
    bgeot::copy(v, w);
    
    cout << "w = " << w << endl;

    bgeot::fsvector<double, 10> x;

    bgeot::copy(v, x);

    cout << "x = " << x << endl;

    bgeot::add(v, w, x);

    cout << "x = " << x << endl;

    // cout << "sp = " << bgeot::vect_sp(v, w) << endl;

    bgeot::svector<double> z(10);
    z[4] = 1;z[5] = 1;

    bgeot::add(v, z, x);

    cout << "x = " << x << endl;
    
    
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;
}
