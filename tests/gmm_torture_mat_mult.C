// RECTANGULAR_MATRIX_PARAM
// RECTANGULAR_MATRIX_PARAM;
// RECTANGULAR_MATRIX_PARAM;
// ENDPARAM;


using gmm::size_type;

template <typename MAT1, typename MAT2, typename MAT3>
void test_procedure(const MAT1 &_m1, const MAT2 &_m2, const MAT3 &_m3) {
  MAT1  &m1 = const_cast<MAT1  &>(_m1);
  MAT2  &m2 = const_cast<MAT2  &>(_m2);
  MAT3  &m3 = const_cast<MAT3  &>(_m3);
  typedef typename gmm::linalg_traits<MAT1>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  double prec = gmm::default_tol(R());

  size_type k = gmm::mat_nrows(m1);
  size_type l = std::max(gmm::mat_ncols(m1), gmm::mat_nrows(m2));
  size_type m = std::max(gmm::mat_ncols(m2), gmm::mat_nrows(m3));
  size_type n = gmm::mat_ncols(m3);

  size_type mm = std::min(m, k);
  size_type nn = std::min(n, m);

  gmm::dense_matrix<T> m1bis(mm, l), m2bis(l, nn), m3bis(mm, nn);
  gmm::copy(gmm::sub_matrix(m1, gmm::sub_interval(0,mm),
			    gmm::sub_interval(0,l)), m1bis);
  gmm::copy(gmm::sub_matrix(m2, gmm::sub_interval(0,l),
			    gmm::sub_interval(0,nn)), m2bis);
  gmm::mult(m1bis, m2bis, m3bis);
  gmm::mult(gmm::sub_matrix(m1, gmm::sub_interval(0,mm),
			    gmm::sub_interval(0,l)),
	    gmm::sub_matrix(m2, gmm::sub_interval(0,l),
			    gmm::sub_interval(0,nn)),
	    gmm::sub_matrix(m3, gmm::sub_interval(0,mm),
			    gmm::sub_interval(0,nn)));
  gmm::add(gmm::scaled(m3bis, T(-1)),
	   gmm::sub_matrix(m3, gmm::sub_interval(0,mm),
			   gmm::sub_interval(0,nn)));
  
  R error = gmm::mat_norm2(gmm::sub_matrix(m3, gmm::sub_interval(0,mm),
					   gmm::sub_interval(0,nn)));

  if (error >= R(prec * 10000.0)) 
    DAL_THROW(gmm::failure_error, "Error too large: " << error);
  
  gmm::mult(gmm::scaled(gmm::sub_matrix(m1, gmm::sub_interval(0,mm),
					gmm::sub_interval(0,l)), T(-2)),
	    gmm::sub_matrix(m2, gmm::sub_interval(0,l),
			    gmm::sub_interval(0,nn)),
	    gmm::sub_matrix(gmm::transposed(m3), gmm::sub_interval(0,mm),
			    gmm::sub_interval(0,nn)));
  gmm::add(gmm::scaled(m3bis, T(2)),
	   gmm::transposed(gmm::sub_matrix(m3, gmm::sub_interval(0,nn),
					   gmm::sub_interval(0,mm))));
  
  error = gmm::mat_norm2(gmm::sub_matrix(m3, gmm::sub_interval(0,nn),
					 gmm::sub_interval(0,mm)));

  if (error >= R(prec * 10000.0)) 
    DAL_THROW(gmm::failure_error, "Error too large: " << error);
  
}
