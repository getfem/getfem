#include <bgeot_approx_integration.h>


int main(void)
{
  try {
  bgeot::papprox_integration pai;

  for (int i = 1; i < 15; ++i)
  {
    pai = bgeot::Gauss_approx_integration(i);

    cout << "methode a " << i << " points " << endl;

    for (int k = 0; k < pai->nb_points(); ++k)
    {
      cout << "Coeff " << k << " : " << pai->integration_coefficients()[k];
      cout << " point en [0,1] : " << pai->integration_points()[k][0];
      cout << " point en [-1,1] : "
	   << 2.0 * pai->integration_points()[k][0] - 1.0 << endl;

    }

    cout << endl << endl;

  }

  pai = bgeot::convex_product_approx_integration(
      bgeot::Gauss_approx_integration(2), bgeot::Gauss_approx_integration(2));

  cout << "methode produit\n";

  for (int k = 0; k < pai->nb_points(); ++k)
  {
    cout << "Coeff " << k << " : " << pai->integration_coefficients()[k];
    cout << " point : " << pai->integration_points()[k] << endl;

  }

  cout << endl << endl;

  for (int n = 1; n < 6; n++)
  {
    cout << "methode Newton-Cotes en dimension " << n << "\n";
    for (int i = 0; i < 2; ++i)
    {
      cout << "methode d'ordre  " << i << "\n";
      pai = bgeot::Newton_Cotes_approx_integration(n, i);
  
      for (int k = 0; k < pai->nb_points(); ++k)
      {
	cout << "Coeff " << k << " : " << pai->integration_coefficients()[k];
	cout << " point : " << pai->integration_points()[k] << endl; 
      }
      
      cout << endl << endl;
    }
  }
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;

}
