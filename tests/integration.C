#include <getfem_integration.h>

using getfem::size_type;

int main(void)
{
  try {
    char meth[500];
    getfem::papprox_integration pai;
    cout.precision(16);
    
    for (size_type i = 1; i < 15; ++i) {
      sprintf(meth, "IM_GAUSS1D(%d)", 2*(i - 1));
      pai = getfem::int_method_descriptor(meth)->method.pai;
      
      cout << "methode a " << i << " points " << endl;
      
      for (size_type k = 0; k < pai->nb_points(); ++k) {
	cout << "Coeff " << k << " : " << pai->integration_coefficients()[k];
	cout << " point en [0,1] : " << pai->integration_points()[k][0];
	cout << " point en [-1,1] : "
	     << 2.0 * pai->integration_points()[k][0] - 1.0 << endl;
	
      }
      
      cout << endl << endl;
      
    }

    sprintf(meth, "IM_PRODUCT(IM_GAUSS1D(2),IM_GAUSS1D(2))");
    pai = getfem::int_method_descriptor(meth)->method.pai;

    cout << "methode produit\n";
    
    for (size_type k = 0; k < pai->nb_points(); ++k) {
      cout << "Coeff " << k << " : " << pai->integration_coefficients()[k];
      cout << " point : " << pai->integration_points()[k] << endl;
      
    }
    
    cout << endl << endl;
    
    for (size_type n = 1; n < 6; n++) {
      cout << "methode Newton-Cotes en dimension " << n << "\n";
      for (size_type i = 0; i < 3; ++i) {
	cout << "methode d'ordre  " << i << "\n";
	sprintf(meth, "IM_NC(%d,%d)", n, i);
	pai = getfem::int_method_descriptor(meth)->method.pai;
	
	for (size_type k = 0; k < pai->nb_points(); ++k) {
	  cout << "Coeff " << k << " : " << pai->integration_coefficients()[k];
	  cout << " point : " << pai->integration_points()[k] << endl; 
	}
	
	cout << endl << endl;
      }
    }

    sprintf(meth, "IM_NC(2, 2)");
    pai = getfem::int_method_descriptor(meth)->method.pai;

    cout << "methode : " << meth << endl;

    cout << "Nb points on convex " << pai->nb_points_on_convex() << endl;
    for (size_type k = 0; k < pai->structure()->nb_faces(); ++k)
      cout << "Nb points on face " << k << " : "
	   <<  pai->nb_points_on_face(k) << endl;
    for (size_type k = 0; k < pai->nb_points(); ++k) {
      cout << "Coeff " << k << " : " << pai->integration_coefficients()[k];
      cout << "\t point : " << pai->integration_points()[k] << endl;
      
    }

    sprintf(meth, "IM_STRUCTURED_COMPOSITE(IM_NC(2, 2), 1)");
    pai = getfem::int_method_descriptor(meth)->method.pai;

    cout << "methode : " << meth << endl;

    cout << "Nb points on convex " << pai->nb_points_on_convex() << endl;
    for (size_type k = 0; k < pai->structure()->nb_faces(); ++k)
      cout << "Nb points on face " << k << " : "
	   <<  pai->nb_points_on_face(k) << endl;
    for (size_type k = 0; k < pai->nb_points(); ++k) {
      cout << "Coeff " << k << " : " << pai->integration_coefficients()[k];
      cout << "\t point : " << pai->integration_points()[k] << endl;
      
    }

  }
  DAL_STANDARD_CATCH_ERROR;
  
  return 0;

}
