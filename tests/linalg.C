#include <bgeot_abstract_linalg.h>


int main(void)
{
  try {
    bgeot::vsvector<double> v(10), w(10);
    v.fill(1.0);
    w.fill(0.0);
    
    bgeot::copy(v, w);
    
    cout << "w = " << w << endl;

    // cout << "sp = " << bgeot::vect_sp(v, w) << endl;
    
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;
}
