/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */
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

class exception_cb : public dal::exception_callback  {
  public:
  virtual void callback(const std::string& msg)
  { cerr << msg << endl; *(int *)(0) = 0; }
};


int main(void)
{
  exception_cb cb;
  dal::exception_callback::set_exception_callback(&cb);

  try {
    char meth[500];
    cout.precision(16);
    
    for (size_type i = 1; i < 15; ++i) {
      sprintf(meth, "IM_GAUSS1D(%d)", int(2*(i - 1)));
      print_method(getfem::int_method_descriptor(meth));
    }

    sprintf(meth, "IM_PRODUCT(IM_GAUSS1D(2),IM_GAUSS1D(2))");
    print_method(getfem::int_method_descriptor(meth));
    
    for (size_type n = 1; n < 6; n++) {
      for (size_type i = 0; i < 3; ++i) {
	sprintf(meth, "IM_NC(%d,%d)", int(n), int(i));
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
