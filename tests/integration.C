/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2004  Yves Renard, Julien Pommier.                   */
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
#include <numeric>
#include <getfem_integration.h>
#include <bgeot_comma_init.h>
#include <getfem_mesh_fem.h>
#include <iomanip>
#include <map>

using getfem::size_type;
using bgeot::base_tensor;
using bgeot::base_matrix;
using bgeot::base_vector;
using bgeot::base_node;
using bgeot::scalar_type;
using bgeot::opt_long_scalar_type;
using bgeot::dim_type;

bool do_heavy_checks = false;

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

class matrix_collection {
public:
  std::vector<base_vector> lst;
  std::vector<std::string> im_names;
};

std::map<std::pair<bgeot::pgeometric_trans,size_type>, matrix_collection> ME;

static void check_method(const std::string& im_name, getfem::pintegration_method ppi, size_type k, bgeot::pgeometric_trans pgt) {
  getfem::getfem_mesh m; 
  getfem::mesh_fem mf1(m);
  getfem::mesh_fem mf2(m);
  assert(ppi); assert(pgt);
  cout << "checking " << im_name << "...\n" << std::flush; 
  m.add_convex_by_points(pgt, pgt->convex_ref()->points().begin());
  mf1.set_finite_element(m.convex_index(),getfem::classical_fem(pgt,k/2),ppi);
  mf2.set_finite_element(m.convex_index(),getfem::classical_fem(pgt,(k-k/2)),ppi);
  matrix_collection &mc = ME[std::make_pair(pgt,k)];
  getfem::pmat_elem_type pme = getfem::mat_elem_product(getfem::mat_elem_base(mf1.fem_of_element(0)),getfem::mat_elem_base(mf2.fem_of_element(0)));
  getfem::pmat_elem_computation pmec = getfem::mat_elem(pme, ppi, pgt);
  base_tensor t;
  pmec->gen_compute(t, m.points_of_convex(0), 0);
  mc.lst.push_back(t);
  mc.im_names.push_back(im_name);
}

static void check_im_order(const std::string& s, size_type expected_pk=size_type(-1), size_type expected_qk=size_type(-1)) {
  getfem::pintegration_method ppi = getfem::int_method_descriptor(s);
  size_type pk = 10000, qk = 10000;
  if (!ppi->is_exact()) {
    size_type dim = ppi->approx_method()->dim();
    for (bgeot::power_index idx(dim); idx.degree() <= pk; ++idx) {
      opt_long_scalar_type sum = 0, realsum = 1.;
      for (size_type i=0; i < ppi->approx_method()->nb_points_on_convex(); ++i) {
	opt_long_scalar_type prod = ppi->approx_method()->coeff(i);
	for (size_type d=0; d < dim; ++d) prod *= pow(ppi->approx_method()->point(i)[d], idx[d]);
	sum += prod;
      }
      if (ppi->structure()->basic_structure() == bgeot::simplex_structure(dim)) {        
        size_type fa = 1;
        for (size_type z = 0; z < dim; z++)
          for (size_type k = 1; k <= idx[z]; ++k, ++fa)
            realsum *= opt_long_scalar_type(scalar_type(k)) / opt_long_scalar_type(scalar_type(fa));
        for (size_type k = 0; k < dim; k++) { realsum /= opt_long_scalar_type(scalar_type(fa)); fa++; }
    /*	for (size_type d=dim-1, c=0; d+1 != 0; --d) { c += idx[d]+1; realsum *= opt_long_scalar_type(c); }
      realsum = opt_long_scalar_type(1.)/realsum;*/
      } else if (ppi->structure()->basic_structure() == bgeot::parallelepiped_structure(dim)) {
	for (size_type d=0; d < dim; ++d) realsum *= opt_long_scalar_type(idx[d]+1);
	realsum = opt_long_scalar_type(1.)/realsum;
      }
      if (dal::abs((realsum - sum)/realsum) > 1e-9) { 
        /*	cout << "degree=" << idx.degree() << ", idx=";
          for (size_type d=0; d < dim; ++d) cout << idx[d] << " "; cout << ", realsum=" << realsum << ", sum = " << sum << "\n";*/
	pk = std::min<size_type>(pk,idx.degree());
	qk = std::min<size_type>(qk, *std::max_element(idx.begin(),idx.end()));
      }
    }
  }
  cout << std::setw(70) << getfem::name_of_int_method(ppi) << ", PK DEGREE=" << std::setw(2) << pk-1 
       << ", QK DEGREE=" << std::setw(2) << qk-1 << "\n";
}

const std::vector<size_type>& TRIANGLE_D() { 
  static std::vector<size_type> i_d;
  if (i_d.size() == 0) bgeot::sc(i_d) += 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13;
  return i_d;
}
const std::vector<size_type>& TETRA_D() { 
  static std::vector<size_type> i_d;
  if (i_d.size() == 0) bgeot::sc(i_d) += 1, 2, 3, 5, 6, 8;
  return i_d;
}
const std::vector<size_type>& SIMPLEX4_D() { 
  static std::vector<size_type> i_d;
  if (i_d.size() == 0) bgeot::sc(i_d) += 3;
  return i_d;
}
const std::vector<size_type>& QUAD_D() { 
  static std::vector<size_type> i_d;
  if (i_d.size() == 0) bgeot::sc(i_d) += 2, 3, 5, 7, 9, 17;
  return i_d;
}
const std::vector<size_type>& HEXA_D() { 
  static std::vector<size_type> i_d;
  if (i_d.size() == 0) bgeot::sc(i_d) += 5,9,11;
  return i_d;
}
const std::vector<size_type>& CUBE4D_D() { 
  static std::vector<size_type> i_d;
  if (i_d.size() == 0) bgeot::sc(i_d) += 5,9;
  return i_d;
}

static void check_orders() {
  char s[512];
  for (std::vector<size_type>::const_iterator it = TRIANGLE_D().begin(); it != TRIANGLE_D().end(); ++it) {
    sprintf(s,"IM_TRIANGLE(%d)",*it); check_im_order(s);
  }
  for (std::vector<size_type>::const_iterator it = TETRA_D().begin(); it != TETRA_D().end(); ++it) {
    sprintf(s,"IM_TETRAHEDRON(%d)",*it); check_im_order(s);
  }
  for (std::vector<size_type>::const_iterator it = QUAD_D().begin(); it != QUAD_D().end(); ++it) {
    sprintf(s,"IM_QUAD(%d)",*it); check_im_order(s);
  }
  for (std::vector<size_type>::const_iterator it = TETRA_D().begin(); it != TETRA_D().end(); ++it) {
    sprintf(s,"IM_TETRAHEDRON(%d)",*it); check_im_order(s);
  }
  for (std::vector<size_type>::const_iterator it = SIMPLEX4_D().begin(); it != SIMPLEX4_D().end(); ++it) {
    sprintf(s,"IM_SIMPLEX4D(%d)",*it); check_im_order(s);
  }
  for (std::vector<size_type>::const_iterator it = HEXA_D().begin(); it != HEXA_D().end(); ++it) {
    sprintf(s,"IM_HEXAHEDRON(%d)",*it); check_im_order(s);
  }
  for (std::vector<size_type>::const_iterator it = CUBE4D_D().begin(); it != CUBE4D_D().end(); ++it) {
    sprintf(s,"IM_CUBE4D(%d)",*it); check_im_order(s);
  }
}

static void check_methods() {
  char s[512];
  getfem::pintegration_method ppi;
  for (size_type k=0; k < 15; ++k) {
    sprintf(s,"IM_GAUSS1D(%d)",k); ppi = getfem::int_method_descriptor(s);
    check_method(s,ppi,k,bgeot::simplex_geotrans(1,1));
    sprintf(s,"IM_NC(1,%d)",k); ppi = getfem::int_method_descriptor(s);
    check_method(s,ppi,k,bgeot::simplex_geotrans(1,1));
    sprintf(s,"IM_EXACT_SIMPLEX(1)"); ppi = getfem::int_method_descriptor(s);
    check_method(s,ppi,k,bgeot::simplex_geotrans(1,1));
  }
  for (size_type d=2; d < 5; ++d) {
    for (size_type k=0; k < 7-d; ++k) {
      sprintf(s,"IM_EXACT_SIMPLEX(%d)",d); ppi = getfem::int_method_descriptor(s);
      check_method(s,ppi,k,bgeot::simplex_geotrans(d,1));
    }
  }
  for (std::vector<size_type>::const_iterator it = TRIANGLE_D().begin(); it != TRIANGLE_D().end(); ++it) {
    sprintf(s,"IM_TRIANGLE(%d)",*it); ppi = getfem::int_method_descriptor(s);
    for (size_type k=1; k <= *it; ++k) { 
      check_method(s,ppi,k,bgeot::simplex_geotrans(2,1));
    }
  }
  for (size_type d=2; d < 5; ++d) {
    for (size_type i=1; i < 8; ++i) {
      for (size_type k=0; k < std::min(i,5-d); ++k) {
	sprintf(s,"IM_NC(%d,%d)",d,i); ppi = getfem::int_method_descriptor(s);
	check_method(s,ppi,k,bgeot::simplex_geotrans(d,1));
      }
    }
  }
  for (std::vector<size_type>::const_iterator it = TETRA_D().begin(); it != TETRA_D().end(); ++it) {
    sprintf(s,"IM_TETRAHEDRON(%d)",*it); ppi = getfem::int_method_descriptor(s);
    for (size_type k=1; k <= *it; ++k) { 
      check_method(s,ppi,k,bgeot::simplex_geotrans(3,1));
    }
  }
  for (std::vector<size_type>::const_iterator it = SIMPLEX4_D().begin(); it != SIMPLEX4_D().end(); ++it) {
    sprintf(s,"IM_SIMPLEX4D(%d)",*it); ppi = getfem::int_method_descriptor(s);
    for (size_type k=1; k <= *it; ++k) { 
      check_method(s,ppi,k,bgeot::simplex_geotrans(4,1));
    }
  }
  
  for (size_type d=1; d < 5; ++d) {
    size_type kmax = 0;
    switch (d) {
    case 1: kmax = 10; break;
    case 2: kmax = 10; break;
    case 3: kmax = 6; break;
    default: kmax = 3; break;
    }
    for (size_type k=0; k < kmax; ++k) {
      sprintf(s,"IM_EXACT_PARALLELEPIPED(%d)",d); ppi = getfem::int_method_descriptor(s);
      check_method(s,ppi,k,bgeot::parallelepiped_linear_geotrans(d));
      sprintf(s,"IM_GAUSS_PARALLELEPIPED(%d,%d)",d,k); ppi = getfem::int_method_descriptor(s);
      check_method(s,ppi,k,bgeot::parallelepiped_linear_geotrans(d));
      sprintf(s,"IM_NC_PARALLELEPIPED(%d,%d)",d,k); ppi = getfem::int_method_descriptor(s);
      check_method(s,ppi,k,bgeot::parallelepiped_linear_geotrans(d));
      if (d>1) {
	sprintf(s,"IM_PRODUCT(IM_GAUSS_PARALLELEPIPED(%d,%d),IM_NC(1,%d))",d-1,k,k); ppi = getfem::int_method_descriptor(s);
	check_method(s,ppi,k,bgeot::parallelepiped_linear_geotrans(d));
      }
    }
  }
  for (std::vector<size_type>::const_iterator it = QUAD_D().begin(); it != QUAD_D().end(); ++it) {
    sprintf(s,"IM_QUAD(%d)",*it); ppi = getfem::int_method_descriptor(s);
    for (size_type k=1; k <= size_type(sqrt(scalar_type(*it))); k++) { 
      check_method(s,ppi,k,bgeot::parallelepiped_linear_geotrans(2));
    }
  }
  for (std::vector<size_type>::const_iterator it = HEXA_D().begin(); it != HEXA_D().end(); ++it) {
    sprintf(s,"IM_HEXAHEDRON(%d)",*it); ppi = getfem::int_method_descriptor(s);
    check_method(s,ppi,size_type(::pow(scalar_type(*it),1./3.)), bgeot::parallelepiped_linear_geotrans(3));
  }
  for (std::vector<size_type>::const_iterator it = CUBE4D_D().begin(); it != CUBE4D_D().end(); ++it) {
    sprintf(s,"IM_CUBE4D(%d)", *it); ppi = getfem::int_method_descriptor(s);
    check_method(s,ppi,1,bgeot::parallelepiped_linear_geotrans(4));
  }

  for (size_type d=2; d < 5; ++d) {
    for (size_type k=0; k < 7-d; ++k) {
      sprintf(s,"IM_EXACT_PRISM(%d)",d);
      ppi = getfem::int_method_descriptor(s);
      
      check_method(s, getfem::int_method_descriptor(s), k,
		   bgeot::prism_linear_geotrans(d));
      sprintf(s,"IM_NC_PRISM(%d,%d)",d,k);
      ppi = getfem::int_method_descriptor(s);
      check_method(s, getfem::int_method_descriptor(s), k,
		   bgeot::prism_geotrans(d, std::max(k, size_type(1))));
    
      if (d == 3) {	
	sprintf(s,"IM_PRODUCT(IM_TRIANGLE(6),IM_GAUSS1D(6))");
	ppi = getfem::int_method_descriptor(s);
	check_method(s, getfem::int_method_descriptor(s), k,
		     bgeot::prism_geotrans(d,k));
      }
    }
  }

  {
    sprintf(s,"IM_STRUCTURED_COMPOSITE(IM_GAUSS1D(3),4)");
    check_method(s, getfem::int_method_descriptor(s), 3, bgeot::simplex_geotrans(1,1));
    sprintf(s,"IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(3),4)");
    check_method(s, getfem::int_method_descriptor(s), 3,bgeot::simplex_geotrans(2,1));
    sprintf(s,"IM_STRUCTURED_COMPOSITE(IM_TETRAHEDRON(5),3)");
    check_method(s, getfem::int_method_descriptor(s), 5,bgeot::simplex_geotrans(3,1));
    /* // not implemented ...
      sprintf(s,"IM_STRUCTURED_COMPOSITE(IM_NC(4,2),3)");
      check_method(s, getfem::int_method_descriptor(s), 2,bgeot::simplex_geotrans(4,1));
    */
    sprintf(s,"IM_STRUCTURED_COMPOSITE(IM_QUAD(5),10)"); // QUAD(5) can't integrate Q5 polynomials, but it is sufficiently refined...
    check_method(s, getfem::int_method_descriptor(s), 5, bgeot::parallelepiped_linear_geotrans(2));
    sprintf(s,"IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(3,2),2)");
    check_method(s, getfem::int_method_descriptor(s), 2, bgeot::parallelepiped_linear_geotrans(3));
    sprintf(s,"IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(4,2),2)");
    check_method(s, getfem::int_method_descriptor(s), 2, bgeot::parallelepiped_linear_geotrans(4));
    cerr << "FIXME:  structured_mesh not implemented for prisms\n";
    /*sprintf(s,"IM_STRUCTURED_COMPOSITE(IM_NC_PRISM(3,3),2)");
      check_method(s, getfem::int_method_descriptor(s), 2, bgeot::prism_geotrans(3,1));*/
  }
}

static int inspect_results() {
  static int failcnt = 0;
  for (std::map< std::pair<bgeot::pgeometric_trans,size_type>, matrix_collection>::const_iterator it = ME.begin();
       it != ME.end(); ++it) {
    size_type K = (*it).first.second;
    bgeot::pgeometric_trans pgt = (*it).first.first;
    const matrix_collection &mc = (*it).second;
    cout << "inspecting " << bgeot::name_of_geometric_trans(pgt) << "/K=" << K << "\n";
    //<< " : " << mc.im_names.size() << " integration results\n";
    scalar_type sumref = std::accumulate(mc.lst[0].begin(), mc.lst[0].end(),0.);
    cout << " reference " << std::setw(70) << mc.im_names[0] << " : sum= " << std::setw(6) << sumref << "\n";    
    for (size_type i = 1; i < mc.im_names.size(); ++i) {
      scalar_type sum = std::accumulate(mc.lst[0].begin(), mc.lst[0].end(),0.);
      scalar_type dist = bgeot::vect_dist2(mc.lst[0],mc.lst[i]);
      bool ok  = (dal::abs(sum-sumref) < 1e-10 && dal::abs(dist) < 1e-6);
      if (ok)  cout << "  [OK]     ";
      else     cout << "  [ERROR!] ";
      cout << std::setw(70) << mc.im_names[i] << " : sum= " << std::setw(6) << sum << ", dist=" << std::setw(9) << dist << "\n";
      if (!ok) {
	// cerr << mc.lst[0] << "\n" << mc.lst[i] << "\n";
	//	cerr << " !!integration: error with " << mc.im_names[i] << " or " << mc.im_names[0] << "\n";
	++failcnt;
      }
    }
  }
  return failcnt;
}


static void print_some_methods() {
  char meth[500];
  cout.precision(60);
    
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

  sprintf(meth, "IM_STRUCTURED_COMPOSITE(IM_QUAD(2),3)");
  print_method(getfem::int_method_descriptor(meth));

  sprintf(meth, "IM_TRIANGLE(3)");
  print_method(getfem::int_method_descriptor(meth));

  sprintf(meth, "IM_TETRAHEDRON(3)");
  print_method(getfem::int_method_descriptor(meth));

  sprintf(meth, "IM_HEXAHEDRON(9)");
  print_method(getfem::int_method_descriptor(meth));
}

int main(int argc, char **argv)
{
  exception_cb cb;
  dal::exception_callback::set_exception_callback(&cb);

  if (argc == 2 && strcmp(argv[1], "-all")) do_heavy_checks = true;
  try {
    print_some_methods();
    check_methods();
    int failcnt = inspect_results();
    cout << "\nOrders of some approximate integration methods:\n";
    check_orders();
    if (failcnt) { cerr << "an error occured with " << failcnt << " integration methods\n"; return 1; }
  }
  DAL_STANDARD_CATCH_ERROR;
  
  return 0;

}
