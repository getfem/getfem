#include <gmm.h>


// sclalar product working also for matrices (to be done in GMM++ ...
template<class VAR> 
typename gmm::linalg_taits<VAR>::value_type
local_sp(const VAR &X, const VAR &Y)
{ return local_sp(X, Y, typename gmm::linalg_traits<VAR>::linalg_type); }

template<class VAR> 
typename gmm::linalg_taits<VAR>::value_type
local_sp(const VAR &X, const VAR &Y, gmm::abstract_vector)
{ return gmm::vect_sp(X, Y); }

template<class VAR> 
typename gmm::linalg_taits<VAR>::value_type
local_sp(const VAR &X, const VAR &Y, gmm::abstract_matrix) {
  typename gmm::linalg_taits<VAR>::value_type res(0);
  for (gmm::size_type i = 0; i < gmm::mat_nrows(X); ++i) 
    for (gmm::size_type j = 0; j < gmm::mat_ncols(X); ++j)
      res += X(i, j) * Y(i, j);
  return res;
}



// Make a test of the gradient around X.
template <class FUNC, class GRAD, class VAR> 
test_grad(FUNC f, GRAD grad, const VAR &X) {
  
  typedef typename gmm::linalg_taits<VAR>::value_type T;
  typedef typename number_traits<T>::magnitude_type R;
  VAR Y(X), Z(X), G(X);
  
  grad(X, G);
  T valx = f(X);

  R eps(1);
  gmm::fill_random(Z);
  T derdir = local_sp(G, Z);
  for (int i = 0; i < 10; ++i, eps /= R(2)) {
    gmm::add(gmm::scaled(Z, eps), X, Y);
    
    T valy = f(Y);
    T estimate_derdir = (valy - valx) / eps;

    cout << " " << gmm::abs(derdir - estimate_derdir);

  }
  cout << endl;
  
  
  
  

}


template <typename MAT, typename MAT2> void
squared_Frobenius_condition_number_gradient(const MAT& M, MAT2& G) { 
  typedef typename gmm::linalg_traits<MAT>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  
  size_type n = mat_ncols(M);
  gmm::dense_matrix<T> B(n,n), C(n,n);
  gmm::mult(M,gmm::transposed(M),B);
  R trB = gmm::mat_trace(B);
  gmm::lu_inverse(B);
  R trBinv = gmm::mat_trace(B);
  gmm::mult(B,B,C);
  gmm::mult(gmm::scaled(M, T(-2)*trB), C, G);
  gmm::add(gmm::scaled(M, T(2)*trBinv), G);
}


typedef gmm::dense_matrix<double> DM;

struct func {
  double operator()(const DM &M) { return Frobenius_condition_number_sqr(M); } 
};

struc grad {
  void operator()(const DM &M, DM &G)
    { squared_Frobenius_condition_number_gradient(M, G); }
};


int main(void) {


  DM M(2, 2);
  gmm::fill_random(M);
  cout << "M = " << M << endl;

  test_grad(func(), grad(), M);
  



  return 0;
}
