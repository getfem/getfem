// VECTOR_PARAM;
// RECTANGULAR_MATRIX_PARAM;
// VECTOR_PARAM;
// VECTOR_PARAM;
// ENDPARAM;


using gmm::size_type;

template <class MAT1, class VECT1, class VECT2, class VECT3>
void test_procedure(const VECT1 &_v1, const MAT1 &_m1, const VECT2 &_v2,
		    const VECT3 &_v3) {
  VECT1 &v1 = const_cast<VECT1 &>(_v1);
  VECT2 &v2 = const_cast<VECT2 &>(_v2);
  VECT3 &v3 = const_cast<VECT3 &>(_v3);
  MAT1  &m1 = const_cast<MAT1  &>(_m1);
  typedef typename gmm::linalg_traits<MAT1>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  double prec = gmm::default_tol(R());

  size_type m = gmm::mat_nrows(m1), n = gmm::mat_ncols(m1), mm = std::min(m,n);
  gmm::sub_interval SUBI(0, m);
  gmm::dense_matrix<T> q(m, mm), r(mm, mm);
  std::vector<T> v4(m);
  
  if (m >= n) {
    gmm::mult(m1, v2, v1);
    gmm::qr_factor(m1, q, r);
    gmm::mult(r, v2, v3);
    gmm::mult(q, v3, gmm::scaled(v1, T(-1)), v4);
  }
  else {
    gmm::mult(gmm::sub_matrix(m1, SUBI), gmm::sub_vector(v2, SUBI), v1);
    gmm::qr_factor(gmm::sub_matrix(m1, SUBI), q, r);
    gmm::mult(r, gmm::sub_vector(v2, SUBI),  gmm::sub_vector(v3, SUBI));
    gmm::mult(q, gmm::sub_vector(v3, SUBI), gmm::scaled(v1, T(-1)), v4);
  }
  if (gmm::vect_norm2(v4) >= R(prec * 1000.0))
    DAL_THROW(dal::failure_error, "Error too large: " << gmm::vect_norm2(v1));
  
}
