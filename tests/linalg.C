#include <gmm.h>
#include <bgeot_matrix.h>


int main(void)
{
  try {
    std::vector<double> v(10), w(10);
    gmm::clear(v);
    gmm::clear(w);
    std::fill(v.begin(), v.end(), 1.0);
    gmm::copy(v, w);
    
    cout << "w = "; gmm::write(w, cout); cout << endl;

    bgeot::fsvector<double, 10> x;

    gmm::copy(v, x);

    cout << "x = "; gmm::write(x, cout); cout << endl;

    gmm::add(v, w, x);

    cout << "x = "; gmm::write(x, cout); cout << endl;

    // cout << "sp = " << gmm::vect_sp(v, w) << endl;

    gmm::wsvector<double> z(10);
    z[4] = 1; z[5] = 1;

    gmm::add(v, z, x);

    cout << "x = "; gmm::write(x, cout); cout << endl;

    gmm::copy(z, x);
    
    cout << "x = "; gmm::write(x, cout); cout << endl;

    bgeot::vsmatrix<double> m(10, 10);
    gmm::clear(m); m(3, 2) = 1.0;

    gmm::transposed(m)(6,8) = 2.0;

    cout << "transposed(m)(2,3) = " <<  (gmm::transposed(m))(2,3) << endl;
    cout << "transposed(m)(3,2) = " <<  (gmm::transposed(m))(3,2) << endl;
    
    gmm::row_matrix<gmm::wsvector<double> > n(10, 10); 
    
    gmm::copy(m, gmm::transposed(n));
    gmm::copy(gmm::transposed(m), n);
    cout << "m = "; gmm::write(m, cout); cout << endl;
    cout << "n = "; gmm::write(n, cout); cout << endl;

    cout << "x = "; gmm::write(x, cout); cout << endl;
    gmm::add(gmm::scaled(z, 10.0), x);
    cout << "x = "; gmm::write(x, cout); cout << endl;
    gmm::add(gmm::scaled(v, -10.0), x);
    cout << "x = "; gmm::write(x, cout); cout << endl;

    gmm::mult(gmm::transposed(m), x, x);

    cout << "x = "; gmm::write(x, cout); cout << endl;

    gmm::clear(x);
    std::vector<size_t> index(3); index[0] = 2; index[1] = 8; index[2] = 3;
    std::vector<double> zz(3); zz[0] = 1; zz[1] = 2; zz[2] = 3; 
    
    gmm::copy(zz, gmm::sub_vector(x, index.begin(), index.end()));
    cout << "x = "; gmm::write(x, cout); cout << endl;

    std::vector<size_t> index2(2); index2[0] = 1; index2[2] = 0;
    std::vector<double> zzz(2); zzz[0] = 5; zzz[1] = 6; 
    
    gmm::copy(zzz, gmm::sub_vector(gmm::sub_vector(x, index.begin(),
							 index.end()),
				       index2.begin(), index2.end()));
    cout << "x = "; gmm::write(x, cout); cout << endl;


    std::vector<size_t> index3(3); index3[0] = 2; index3[1] = 4; index3[2] = 3;
    gmm::reverse_index ri(index3.begin(), index3.end(), 10, m);
    bgeot::vsmatrix<double> mm(3,3);

    cout << "mm(0,0) = " << gmm::sub_matrix(m, index3.begin(), index3.end(),
		 index3.begin(), index3.end(),ri,ri)(0,0) << endl;

    gmm::copy(gmm::sub_matrix(m, index3.begin(), index3.end(),
				  index3.begin(), index3.end(),ri,ri), mm);
    cout << "mm = "; gmm::write(mm, cout); cout << endl;
   

  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;
}
