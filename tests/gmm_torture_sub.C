// SQUARED_MATRIX_PARAM
// VECTOR_PARAM;
// VECTOR_PARAM;
// ENDPARAM;

#include <gmm_kernel.h>
#include <gmm_dense_lu.h>

using gmm::size_type;

template <typename MAT1, typename VECT1, typename VECT2>
void test_procedure(const MAT1 &_m1, const VECT1 &_v1, const VECT2 &_v2) {
  VECT1 &v1 = const_cast<VECT1 &>(_v1);
  VECT2 &v2 = const_cast<VECT2 &>(_v2);
  MAT1  &m1 = const_cast<MAT1  &>(_m1);
  typedef typename gmm::linalg_traits<MAT1>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  double prec = gmm::default_tol(R());
  R error, det;

  size_type m = gmm::vect_size(v1), n = m/2;
  std::vector<T> v3(n);

  det = gmm::abs(gmm::lu_det(gmm::sub_matrix(m1, gmm::sub_interval(0,n))));
  det = std::min(det, R(1));
  if (det > R(prec * 10000.0)) {
    gmm::lu_solve(gmm::sub_matrix(m1, gmm::sub_interval(0,n)), v3,
		  gmm::sub_vector(v2, gmm::sub_interval(0,n)));
    gmm::mult(gmm::sub_matrix(m1, gmm::sub_interval(0,n)), v3,
	      gmm::sub_vector(v1, gmm::sub_interval(0,n)));
    gmm::add(gmm::scaled(gmm::sub_vector(v1, gmm::sub_interval(0,n)), T(-1)),
	     gmm::sub_vector(v2, gmm::sub_interval(0,n)), v3);
    if ((error = gmm::vect_norm2(v3)) >= R(prec * 20000.0 / det))
      DAL_THROW(gmm::failure_error, "Error too large: "<< error);
  }

  det = gmm::abs(gmm::lu_det(gmm::sub_matrix(m1, gmm::sub_slice(0,n,1))));
  det = std::min(det, R(1));
  if (det > R(prec * 10000.0)) {
    gmm::lu_solve(gmm::sub_matrix(m1, gmm::sub_slice(0,n,1)), v3,
		  gmm::sub_vector(v2, gmm::sub_slice(0,n,1)));
    gmm::mult(gmm::sub_matrix(m1, gmm::sub_slice(0,n,1)), v3,
	      gmm::sub_vector(v1, gmm::sub_slice(0,n,1)));
    gmm::add(gmm::scaled(gmm::sub_vector(v1, gmm::sub_slice(0,n,1)), T(-1)),
	     gmm::sub_vector(v2, gmm::sub_slice(0,n,1)), v3);
    if ((error = gmm::vect_norm2(v3)) >= R(prec * 20000.0 / det))
      DAL_THROW(gmm::failure_error, "Error too large: "<< error);
  }
  
  gmm::copy(gmm::identity_matrix(), gmm::sub_matrix(gmm::transposed(m1),
		        gmm::sub_interval(0,n), gmm::sub_interval(0,m)));
  gmm::clear(gmm::sub_vector(v2, gmm::sub_interval(n, m-n)));
  cout << "sub matrix of m1 : "
       << gmm::sub_matrix(gmm::transposed(m1), gmm::sub_interval(0,n)) << endl;

  gmm::mult(gmm::sub_matrix(m1, gmm::sub_interval(0,m), 
			    gmm::sub_interval(0,n)),
	    gmm::sub_vector(v2, gmm::sub_interval(0,n)),
	    gmm::scaled(v2, T(-1)), v1);
  if ((error = gmm::vect_norm2(v1)) >= R(prec * 2000.0))
    DAL_THROW(gmm::failure_error, "Error too large: " << error);
  
  
}
