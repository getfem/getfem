#include <bgeot_generic_solvers.h>


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
    z[4] = 1; z[5] = 1;

    bgeot::add(v, z, x);

    cout << "x = " << x << endl;

    bgeot::copy(z, x);
    
    cout << "x = " << x << endl;

    bgeot::vsmatrix<double> m(10, 10);
    m.clear(); m(3, 2) = 1.0;

    bgeot::transposed(m)(6,8) = 2.0;

    cout << "transposed(m)(2,3) = " <<  (bgeot::transposed(m))(2,3) << endl;
    cout << "transposed(m)(3,2) = " <<  (bgeot::transposed(m))(3,2) << endl;
    
    bgeot::smatrix<double> n(10, 10); 
    
    bgeot::copy(m, bgeot::transposed(n));
    bgeot::copy(bgeot::transposed(m), n);
    cout << "m = \n" << m << endl;
    cout << "n = \n" << n << endl;

    cout << "x = " << x << endl;
    bgeot::add(bgeot::scaled(z, 10.0), x);
    cout << "x = " << x << endl;
    bgeot::add(bgeot::scaled(v, -10.0), x);
    cout << "x = " << x << endl;

    bgeot::mult(bgeot::transposed(m), x, x);

    cout << "x = " << x << endl;

    bgeot::clear(x);
    std::vector<size_t> index(3); index[0] = 2; index[1] = 8; index[2] = 3;
    std::vector<double> zz(3); zz[0] = 1; zz[1] = 2; zz[2] = 3; 
    
    bgeot::copy(zz, bgeot::sub_vector(x, index.begin(), index.end()));
    cout << "x = " << x << endl;

    std::vector<size_t> index2(2); index2[0] = 1; index2[2] = 0;
    std::vector<double> zzz(2); zzz[0] = 5; zzz[1] = 6; 
    
    bgeot::copy(zzz, bgeot::sub_vector(bgeot::sub_vector(x, index.begin(),
							 index.end()),
				       index2.begin(), index2.end()));
    cout << "x = " << x << endl;


    std::vector<size_t> index3(3); index3[0] = 2; index3[1] = 4; index3[2] = 3;
    bgeot::reverse_index ri(index3.begin(), index3.end(), 10, m);
    bgeot::vsmatrix<double> mm(3,3);

    cout << "mm(0,0) = " << bgeot::sub_matrix(m, index3.begin(), index3.end(),
		 index3.begin(), index3.end(),ri,ri)(0,0) << endl;

    bgeot::copy(bgeot::sub_matrix(m, index3.begin(), index3.end(),
				  index3.begin(), index3.end(),ri,ri), mm);
    cout << "mm = \n" << mm << endl;
   

  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;
}
