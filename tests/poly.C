#include <bgeot_poly.h>


int main(void)
{
  try {
  bgeot::polynomial<double> P(2,2), Q(2,2);
  P[2] = 1.0;
  Q[3] = 2.0;
  cout << "Le nombre de monomes de P est " << P.size() << endl;
  cout << "P = " << P << endl;
  cout << "Q = " << Q << endl;
  cout << "P + Q = " << P + Q << endl;
  cout << "P * 2.0 * Q = " << P * 2.0 * Q << endl;
  bgeot::polynomial<double> R = P * Q;
  cout << "Le nombre de monomes de R est " << R.size() << endl;
  cout << "Le degre de R est " << R.degree() << endl;
  P.direct_product(Q);
  cout << "Produit direct de P et Q : " << P << endl;

  P = bgeot::polynomial<double>(3,1,1);
  P *= 3.0; P *= bgeot::polynomial<double>(3,1,0);
  P += bgeot::polynomial<double>(3,1,1);
  P *= bgeot::polynomial<double>(3,1,2);
  P += bgeot::polynomial<double>(3,1,0);
  cout << "P = " << P << endl;
  

  double tab[3]; tab[0] = 1.0; tab[1] = 2.0; tab[2] = -1.0;

  cout << "P(1.0, 2.0) = " << P.eval(&(tab[0])) << endl;
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;
}
