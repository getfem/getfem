#include <getfem_integration.h>

using getfem::size_type;


void print_method(getfem::pintegration_method ppi) {
  cout << "methode : " << getfem::name_of_int_method(ppi) << endl;
  getfem::papprox_integration pai = ppi->method.pai;
  cout << "Nb points on convex " << pai->nb_points_on_convex() << endl;
  for (size_type k = 0; k < pai->structure()->nb_faces(); ++k)
    cout << "Nb points on face " << k << " : "
	 <<  pai->nb_points_on_face(k) << endl;
  for (size_type k = 0; k < pai->nb_points(); ++k) {
    cout << "Coeff " << k << " : " << pai->integration_coefficients()[k];
    cout << "\t point : " << pai->integration_points()[k] << endl;
  }
  cout << endl << endl;
}


int main(void)
{
  try {
    char meth[500];
    cout.precision(16);
    
    for (size_type i = 1; i < 15; ++i) {
      sprintf(meth, "IM_GAUSS1D(%ld)", 2*(i - 1));
      print_method(getfem::int_method_descriptor(meth));
    }

    sprintf(meth, "IM_PRODUCT(IM_GAUSS1D(2),IM_GAUSS1D(2))");
    print_method(getfem::int_method_descriptor(meth));
    
    for (size_type n = 1; n < 6; n++) {
      for (size_type i = 0; i < 3; ++i) {
	sprintf(meth, "IM_NC(%ld,%ld)", n, i);
	print_method(getfem::int_method_descriptor(meth));
      }
    }

    sprintf(meth, "IM_NC(2, 2)");
    print_method(getfem::int_method_descriptor(meth));

    sprintf(meth, "IM_STRUCTURED_COMPOSITE(IM_NC(2, 2), 1)");
    print_method(getfem::int_method_descriptor(meth));

    sprintf(meth, "IM_TRIANGLE(3)");
    print_method(getfem::int_method_descriptor(meth));

    sprintf(meth, "IM_TETRAHEDRON(3)");
    print_method(getfem::int_method_descriptor(meth));

    sprintf(meth, "IM_QUAD(3)");
    print_method(getfem::int_method_descriptor(meth));

  }
  DAL_STANDARD_CATCH_ERROR;
  
  return 0;

}
