// SQUARED_MATRIX_PARAM;
// DENSE_VECTOR_PARAM;
// VECTOR_PARAM;
// ENDPARAM;


using gmm::size_type;

template <typename MAT1, typename VECT1, typename VECT2>
void test_procedure(const MAT1 &_m1, const VECT1 &_v1, const VECT2 &_v2) {
  VECT1 &v1 = const_cast<VECT1 &>(_v1);
  VECT2 &v2 = const_cast<VECT2 &>(_v2);
  MAT1  &m1 = const_cast<MAT1  &>(_m1);
  typedef typename gmm::linalg_traits<MAT1>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  double prec = gmm::default_tol(R());

  size_type m = gmm::mat_nrows(m1);
  std::vector<T> v3(m);

  R det = gmm::abs(gmm::lu_det(m1)), error;
  det = std::min(det, R(1));

  if (det > R(prec * 10000.0)) {

    gmm::lu_solve(m1, v1, v2);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);

    error = gmm::vect_norm2(v3);
    if (error >= R(prec * 20000.0 / det))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    gmm::lu_inverse(m1);
    gmm::mult(m1, v2, v1);
    gmm::lu_inverse(m1);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    
    error = gmm::vect_norm2(v3);
    if (error >= R(prec * 20000.0 / det))
      DAL_THROW(gmm::failure_error, "Error too large: "<< error);
  }
}
