//===========================================================================
//
// Copyright (C) 2006-2009 Yves Renard, Julien Pommier.
//
// This file is a part of GetFEM++
//
// GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

#include "getfem_mesh.h"

using std::cout;
using std::cerr;
using std::endl;
using std::cin;
using getfem::scalar_type;
using getfem::size_type;

int ntheta, nphi, nlayers, degre;
scalar_type Rtheta, Rphimax, Rphimin; 
scalar_type X0, Y0, Z0;

getfem::base_node nodepos(int i, int j, int k) 
{
  scalar_type x,y;
  scalar_type theta = i * 2*M_PI / (ntheta);
  scalar_type phi = j * 2*M_PI / (nphi);
  scalar_type Rp = Rphimin + (k*(Rphimax-Rphimin)/(nlayers));

  x = Rp * sin(phi);
  y = Rtheta + Rp * cos(phi);
  
  getfem::base_node n(3);
  n[0] = X0+x;
  n[1] = Y0+y * cos(theta);
  n[2] = Z0+y * sin(theta);
  return n;
}

int main() {
  // coord du centre
  X0 = 0; Y0 = 0; Z0 = 20;

  // rayons
  Rtheta = 15; Rphimax = 5; Rphimin = 4.;

  // nb de mailles
  cerr << "nombre de cellules ntheta   : "; cin >> ntheta;
  cerr << "nombre de cellules nphi     : "; cin >> nphi;
  cerr << "nombre de couches de mailles: "; cin >> nlayers;

  degre = 2;


  ntheta *= degre; nphi *= degre;
  nlayers *= degre;

  getfem::getfem_mesh m;
  bgeot::pgeometric_trans pgt = bgeot::parallelepiped_geotrans(3,degre);

  std::vector<getfem::base_node> N((degre+1)*(degre+1)*(degre+1));
  for (size_type i=0; i < ntheta; i+=degre) {
    for (size_type j=0; j < nphi; j+=degre) {
      for (size_type k=0; k < nlayers; k+=degre) {
	size_type cnt = 0;
	for (size_type ii=0; ii < degre+1; ++ii) {
	  for (size_type jj=0; jj < degre+1; ++jj) {
	    for (size_type kk=0; kk < degre+1; ++kk) {
	      N[cnt++] = nodepos(i+ii,j+jj,k+kk);
	    }
	  }
	}
	m.add_convex_by_points(pgt, N.begin());
      }
    }
  }
  m.write_to_file("donut_regulier.mesh");
  return 0;
}
