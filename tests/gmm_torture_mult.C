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
void test_procedure2(const MAT1 &_m1, const VECT1 &_v1, const VECT2 &_v2, 
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

  size_type m = gmm::vect_size(v1), n = gmm::vect_size(v3);
  size_type nn = std::min(m,n), mm = std::max(m, n);
  std::vector<T> v6(m);


  gmm::lu_solve(m1, v6, v2);
  gmm::mult(m1, v6, v1);
  gmm::add(gmm::scaled(v1, T(-1)), v2, v6);
  if (gmm::vect_norm2(v6) >= R(prec * 10000.0))
    DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v6));

  gmm::lu_solve(gmm::transposed(m1), v6, v2);
  gmm::mult(gmm::transposed(m1), v6, v1);
  gmm::add(gmm::scaled(v1, T(-1)), v2, v6);
  if (gmm::vect_norm2(v6) >= R(prec * 10000.0))
    DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v6));

  gmm::lu_solve(gmm::conjugated(m1), v6, v2);
  gmm::mult(gmm::conjugated(m1), v6, v1);
  gmm::add(gmm::scaled(v1, T(-1)), v2, v6);
  if (gmm::vect_norm2(v6) >= R(prec * 10000.0))
    DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v6));

  gmm::lu_solve(gmm::transposed(gmm::conjugated(m1)), v6, v2);
  gmm::mult(gmm::transposed(gmm::conjugated(m1)), v6, v1);
  gmm::add(gmm::scaled(v1, T(-1)), v2, v6);
  if (gmm::vect_norm2(v6) >= R(prec * 10000.0))
    DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v6));

  gmm::lu_solve(gmm::transposed(gmm::scaled(m1, T(-6))), v6, v2);
  gmm::mult(gmm::transposed(gmm::scaled(m1, T(-6))), v6, v1);
  gmm::add(gmm::scaled(v1, T(-1)), v2, v6);
  if (gmm::vect_norm2(v6) >= R(prec * 10000.0))
    DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v6));

  gmm::dense_matrix<T> q(mm, nn), r(nn, nn);
  if (m >= n) {
    std::vector<T> v5(m);
    gmm::mult(m2, v3, v2);
    gmm::qr_factor(m2, q, r);
    gmm::mult(r, v3, v4);
    gmm::mult(q, v4, gmm::scaled(v2, T(-1)), v5);
    if (gmm::vect_norm2(v5) >= R(prec * 10000.0))
      DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v5));

  }
  else {
    std::vector<T> v5(n);
    gmm::mult(gmm::conjugated(m2), v2, v3);
    gmm::qr_factor(gmm::conjugated(m2), q, r);
    gmm::mult(r, v2, v1);
    gmm::mult(q, v1, gmm::scaled(v3, T(-1)), v5);
    if (gmm::vect_norm2(v5) >= R(prec * 10000.0))
      DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v5));

  }
  
}


template <class MAT1, class MAT2, class VECT1, class VECT2, class VECT3,
	  class VECT4>
void test_procedure(const MAT1 &_m1, const VECT1 &v1, const VECT2 &v2, 
		    const MAT2 &_m2, const VECT3 &v3, const VECT4 &v4) {
  MAT1  &m1 = const_cast<MAT1  &>(_m1);
  MAT2  &m2 = const_cast<MAT2  &>(_m2);
  typedef typename gmm::linalg_traits<MAT1>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  double prec = gmm::default_tol(R());

  test_procedure2(m1, v1, v2, m2, v3, v4);

  size_type m = gmm::vect_size(v1), n = gmm::vect_size(v3);
  gmm::csr_matrix<T> mm1(m, m);
  gmm::copy(m1, mm1);
  gmm::csc_matrix<T> mm2(m, n);
  gmm::copy(m2, mm2);
  test_procedure2(mm1, v1, v2, mm2, v3, v4);

  size_type mm = m / 2; nn = n / 2;
  test_procedure2(gmm::sub_matrix(mm1, gmm::sub_interval(0, mm),
				  gmm::sub_interval(0, mm)), v1, v2,
		  gmm::sub_matrix(mm2, gmm::sub_interval(0, mm),
				  gmm::sub_interval(0, nn)), v3, v4);

  gmm::add(gmm::scaled(mm1, T(-1)), m1);
  gmm::add(gmm::scaled(mm2, T(-1)), m2);
  
  R error = gmm::mat_norm2(m1) + gmm::mat_norm2(m2);
  if (error >= R(prec * 10000.0))
    DAL_THROW(dal::failure_error, "Error too large: "<< error);

}


  
