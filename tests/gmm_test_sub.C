// SQUARED_MATRIX_PARAM
// VECTOR_PARAM;
// VECTOR_PARAM;
// ENDPARAM;


using gmm::size_type;

template <class MAT1, class VECT1, class VECT2>
void test_procedure(const MAT1 &_m1, const VECT1 &_v1, const VECT2 &_v2) {
  VECT1 &v1 = const_cast<VECT1 &>(_v1);
  VECT2 &v2 = const_cast<VECT2 &>(_v2);
  MAT1  &m1 = const_cast<MAT1  &>(_m1);
  typedef typename gmm::linalg_traits<MAT1>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  double prec = gmm::default_tol(R());

  size_type m = gmm::vect_size(v1), n = m/2;
  std::vector<T> v3(n);

  gmm::lu_solve(gmm::sub_matrix(m1, gmm::sub_interval(0,n)), v3,
		gmm::sub_vector(v2, gmm::sub_interval(0,n)));
  gmm::mult(gmm::sub_matrix(m1, gmm::sub_interval(0,n)), v3,
	    gmm::sub_vector(v1, gmm::sub_interval(0,n)));
  gmm::add(gmm::scaled(gmm::sub_vector(v1, gmm::sub_interval(0,n)), T(-1)),
	   gmm::sub_vector(v2, gmm::sub_interval(0,n)), v3);
  if (gmm::vect_norm2(v3) >= R(prec * 10000.0))
    DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v3));

  gmm::lu_solve(gmm::sub_matrix(m1, gmm::sub_slice(0,n,1)), v3,
		gmm::sub_vector(v2, gmm::sub_slice(0,n,1)));
  gmm::mult(gmm::sub_matrix(m1, gmm::sub_slice(0,n,1)), v3,
	    gmm::sub_vector(v1, gmm::sub_slice(0,n,1)));
  gmm::add(gmm::scaled(gmm::sub_vector(v1, gmm::sub_slice(0,n,1)), T(-1)),
	   gmm::sub_vector(v2, gmm::sub_slice(0,n,1)), v3);
  if (gmm::vect_norm2(v3) >= R(prec * 10000.0))
    DAL_THROW(dal::failure_error, "Error too large: "<< gmm::vect_norm2(v3));

  
  gmm::copy(gmm::identity_matrix(), gmm::sub_matrix(gmm::transposed(m1),
		        gmm::sub_interval(0,n), gmm::sub_interval(0,m)));
  gmm::clear(gmm::sub_vector(v2, gmm::sub_interval(n, m-n)));
  cout << "sub matrix of m1 : "
       << gmm::sub_matrix(gmm::transposed(m1), gmm::sub_interval(0,n)) << endl;

  gmm::mult(gmm::sub_matrix(m1, gmm::sub_interval(0,m), 
			    gmm::sub_interval(0,n)),
	    gmm::sub_vector(v2, gmm::sub_interval(0,n)),
	    gmm::scaled(v2, T(-1)), v1);
  if (gmm::vect_norm2(v1) >= R(prec * 1000.0))
    DAL_THROW(dal::failure_error, "Error too large: " << gmm::vect_norm2(v1));
  
  
}
