/* -*- c++ -*- (enables emacs c++ mode)                                    */
#ifndef BGEOT_IMBRICATED_BOX
#define BGEOT_IMBRICATED_BOX

#include <bgeot_vector.h>

namespace bgeot {
  inline scalar_type sfloor(scalar_type x)
  { return (x >= 0) ? floor(x) : -floor(-x); }

  /// Sort order for a rapid localization in bgeot::geotrans_inv
  struct imbricated_box_less
    : public std::binary_function<base_node, base_node, int>
  { 
    mutable int exp_max, exp_min;
    mutable scalar_type c_max;
    unsigned base;

    /// comparaison function
    int operator()(const base_node &x, const base_node &y) const;
    
    imbricated_box_less(unsigned ba = 10, int emi = -15, int ema = -2) {
      base = ba; exp_max = ema; exp_min = emi;
      c_max = pow(double(base), double(-exp_max));
    }
  };
}
#endif
