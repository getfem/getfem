// SQUARED_MATRIX_PARAM
// VECTOR_PARAM;
// VECTOR_PARAM;
// RECTANGULAR_MATRIX_PARAM;
// VECTOR_PARAM;
// VECTOR_PARAM;
// ENDPARAM;


using gmm::size_type;

template <class MAT1, class MAT2, class VECT1, class VECT2, class VECT3,
	  class VECT4>
void test_procedure(const MAT1 &_m1, const VECT1 &_v1, const VECT2 &_v2, 
		    const MAT2 &_m2, const VECT3 &_v3, const VECT4 &_v4) {
  VECT1 &v1 = const_cast<VECT1 &>(_v1);
  VECT2 &v2 = const_cast<VECT2 &>(_v2);
  VECT3 &v3 = const_cast<VECT3 &>(_v3);
  VECT4 &v4 = const_cast<VECT4 &>(_v4);
  MAT1  &m1 = const_cast<MAT1  &>(_m1);
  MAT2  &m2 = const_cast<MAT2  &>(_m2);
  typedef typename gmm::linalg_traits<MAT1>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  double prec = gmm::default_tol(R());

  size_type m = gmm::mat_nrows(m2), n = gmm::mat_ncols(m2);
  size_type nn = std::min(m,n), mm = std::max(m, n);
  std::vector<T> v6(m);

  gmm::lu_solve(m1, v6, v2);
  gmm::mult(m1, v6, v1);
  gmm::add(gmm::scaled(v1, T(-1)), v2, v6);
  if (gmm::vect_norm2(v6) >= R(prec * 1000.0))
    DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v1));

  gmm::lu_solve(gmm::transposed(m1), v6, v2);
  gmm::mult(gmm::transposed(m1), v6, v1);
  gmm::add(gmm::scaled(v1, T(-1)), v2, v6);
  if (gmm::vect_norm2(v6) >= R(prec * 1000.0))
    DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v1));

  gmm::lu_solve(gmm::conjugated(m1), v6, v2);
  gmm::mult(gmm::conjugated(m1), v6, v1);
  gmm::add(gmm::scaled(v1, T(-1)), v2, v6);
  if (gmm::vect_norm2(v6) >= R(prec * 1000.0))
    DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v1));

  gmm::lu_solve(gmm::transposed(gmm::conjugated(m1)), v6, v2);
  gmm::mult(gmm::transposed(gmm::conjugated(m1)), v6, v1);
  gmm::add(gmm::scaled(v1, T(-1)), v2, v6);
  if (gmm::vect_norm2(v6) >= R(prec * 1000.0))
    DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v1));

  gmm::lu_solve(gmm::transposed(gmm::scaled(m1, T(-6))), v6, v2);
  gmm::mult(gmm::transposed(gmm::scaled(m1, T(-6))), v6, v1);
  gmm::add(gmm::scaled(v1, T(-1)), v2, v6);
  if (gmm::vect_norm2(v6) >= R(prec * 1000.0))
    DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v1));

  gmm::dense_matrix<T> q(mm, nn), r(nn, nn);
  if (m >= n) {
    std::vector<T> v5(m);
    gmm::mult(m2, v3, v2);
    gmm::qr_factor(m2, q, r);
    gmm::mult(r, v3, v4);
    gmm::mult(q, v4, gmm::scaled(v2, T(-1)), v5);
    if (gmm::vect_norm2(v4) >= R(prec * 1000.0))
      DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v1));

  }
  else {
    std::vector<T> v5(n);
    gmm::mult(gmm::conjugated(m2), v2, v3);
    gmm::qr_factor(gmm::conjugated(m2), q, r);
    gmm::mult(r, v1, v2);
    gmm::mult(q, v2, gmm::scaled(v3, T(-1)), v5);
    if (gmm::vect_norm2(v4) >= R(prec * 1000.0))
      DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v1));

  }
  
}
