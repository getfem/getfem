// RECTANGULAR_MATRIX_PARAM
// SQUARED_MATRIX_PARAM
// ENDPARAM;


using gmm::size_type;


// template <class MAT, class T> void print_for_matlab(const MAT &m, T) { 
//   cout.precision(16);
//   cout << "[ ";
//   for (size_type i = 0; i < gmm::mat_nrows(m); ++i) {
//     for (size_type j = 0; j < gmm::mat_ncols(m); ++j) cout << " " << m(i,j);
//     if (i != gmm::mat_nrows(m)-1) cout << " ; \n";
//   }
//   cout << " ]" << endl;
// }

// template <class MAT, class T> void print_for_matlab(const MAT &m,
// 						    std::complex<T>) { 
//   cout.precision(16);
//   cout << "[ ";
//   for (size_type i = 0; i < gmm::mat_nrows(m); ++i) {
//     for (size_type j = 0; j < gmm::mat_ncols(m); ++j)
//       cout << " (" << m(i,j).real() << "+" << m(i,j).imag() << "*i)" ;
//     if (i != gmm::mat_nrows(m)-1) cout << " ; \n";
//   }
//   cout << " ]" << endl;
// }

// template <class MAT> inline void print_for_matlab(const MAT &m)
// { print_for_matlab(m, gmm::linalg_traits<MAT>::value_type()); }

template <class T> inline T real_or_complex(double a, double b,T) { return a; }
template <class T> inline
std::complex<T> real_or_complex(double a, double b, std::complex<T>)
{ return std::complex<T>(a, b); }


template <class MAT1, class MAT2>
void test_procedure(const MAT1 &_m1, const MAT2 &_m2) {
  MAT1  &m1 = const_cast<MAT1  &>(_m1);
  MAT2  &m2 = const_cast<MAT2  &>(_m2);
  typedef typename gmm::linalg_traits<MAT1>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  double prec = gmm::default_tol(R());
  R error;

  // gmm::qr_factor(A, Q, R) is tested in test_gmm_mult.C

  //
  // test for gmm::qr_factor(A), apply_house_right and apply_house_left
  //
  size_type m = gmm::mat_nrows(m1), n = gmm::mat_ncols(m1);
  size_type k = size_type(rand() % 50);
    
  gmm::dense_matrix<T> dm1(m, n);
  gmm::copy(m1, dm1);
  if (m >= n) {
    gmm::dense_matrix<T> q(k,m), qaux(k,m), q2(m,k), dm1aux(k,n), m1aux(k,n);
    gmm::fill_random(q); gmm::copy(q, qaux);
    gmm::mult(q, m1, m1aux);

    gmm::qr_factor(dm1);
    gmm::copy(dm1, m1);
    
    gmm::apply_house_right(dm1, q);
    for (size_type i = 0; i < m; ++i)
      for (size_type j = 0; j < i; ++j)
	dm1(i, j) = T(0);
    gmm::mult(q, dm1, dm1aux);
    gmm::add(gmm::scaled(m1aux, T(-1)), dm1aux);
    error = gmm::mat_norm2(dm1aux);
    if (error >= R(prec * 10000.0)) 
      DAL_THROW(dal::failure_error, "Error too large: " << error);

    gmm::copy(gmm::identity_matrix(), q);
    gmm::apply_house_right(m1, q);
    size_type min_km = std::min(k, m);
    gmm::dense_matrix<T> a(min_km, min_km), b(min_km, min_km);
    gmm::copy(gmm::identity_matrix(), b);
    if (k > m) gmm::mult(gmm::conjugated(q), q, a);
    else gmm::mult(q, gmm::conjugated(q), a);
    error = gmm::mat_norm2(a);
    if (error >= R(prec * 10000.0)) 
      DAL_THROW(dal::failure_error, "Error too large: " << error);
      
    gmm::copy(gmm::conjugated(qaux), q2);
    gmm::apply_house_left(m1, q2);
    gmm::mult(gmm::conjugated(q2), dm1, dm1aux);
    gmm::add(gmm::scaled(m1aux, T(-1)), dm1aux);
    error = gmm::mat_norm2(dm1aux);
    if (error >= R(prec * 10000.0)) 
      DAL_THROW(dal::failure_error, "Error too large: " << error);

  }
  else {
    gmm::dense_matrix<T> q(k,n), qaux(k,n), q2(n,k), dm1aux(k,m), m1aux(k,m);
    gmm::fill_random(q); gmm::copy(q, qaux);
    gmm::mult(q, gmm::transposed(m1), m1aux);

    gmm::qr_factor(gmm::transposed(dm1));
    gmm::copy(dm1, m1);
    
    gmm::apply_house_right(gmm::transposed(dm1), q);
    for (size_type i = 0; i < n; ++i)
      for (size_type j = 0; j < i; ++j)
	dm1(j, i) = T(0);
    gmm::mult(q, gmm::transposed(dm1), dm1aux);
    gmm::add(gmm::scaled(m1aux, T(-1)), dm1aux);
    error = gmm::mat_norm2(dm1aux);
    if (error >= R(prec * 10000.0)) 
      DAL_THROW(dal::failure_error, "Error too large: " << error);

    gmm::copy(gmm::identity_matrix(), q);
    gmm::apply_house_right(gmm::transposed(m1), q);
    size_type min_km = std::min(k, n);
    gmm::dense_matrix<T> a(min_km, min_km), b(min_km, min_km);
    gmm::copy(gmm::identity_matrix(), b);
    if (k > n) gmm::mult(gmm::conjugated(q), q, a);
    else gmm::mult(q, gmm::conjugated(q), a);
    error = gmm::mat_norm2(a);
    if (error >= R(prec * 10000.0)) 
      DAL_THROW(dal::failure_error, "Error too large: " << error);
      
    gmm::copy(gmm::conjugated(qaux), q2);
    gmm::apply_house_left(gmm::transposed(m1), q2);
    gmm::mult(gmm::conjugated(q2), gmm::transposed(dm1), dm1aux);
    gmm::add(gmm::scaled(m1aux, T(-1)), dm1aux);
    error = gmm::mat_norm2(dm1aux);
    if (error >= R(prec * 10000.0)) 
      DAL_THROW(dal::failure_error, "Error too large: " << error);

  }
  
  //
  // Test for implicit_qr_algorithm
  //

  m = gmm::mat_nrows(m2);
  gmm::dense_matrix<T> cq(m, m), cr(m, m), ca(m, m);  
  std::vector<T> cv(m), eigc(m);
  gmm::fill_random(cq);
  gmm::copy(cq, cr);
  gmm::lu_inverse(cr);

  gmm::fill_random(cv);
  if (m >  0) cv[ 0] = real_or_complex(     0.0,  0.0, cv[0]);
  if (m >  1) cv[ 1] = real_or_complex(     0.0,  0.0, cv[0]);
  if (m >  2) cv[ 2] = real_or_complex(     0.01,-0.1, cv[0]);
  if (m >  3) cv[ 3] = real_or_complex(     0.01, 0.1, cv[0]);
  if (m >  4) cv[ 4] = real_or_complex(    -2.0,  3.0, cv[0]);
  if (m >  5) cv[ 5] = real_or_complex(    -2.0,  3.0, cv[0]);
  if (m >  6) cv[ 6] = real_or_complex(   -50.0,  3.0, cv[0]);
  if (m >  7) cv[ 7] = real_or_complex(   100.0,  1.0, cv[0]);
  if (m >  8) cv[ 8] = real_or_complex(   300.0,  1.0, cv[0]);
  if (m >  9) cv[ 9] = real_or_complex(   500.0,  1.0, cv[0]);
  if (m > 10) cv[10] = real_or_complex(  1000.0,  1.0, cv[0]);
  if (m > 11) cv[11] = real_or_complex(  4000.0,  1.0, cv[0]);
  if (m > 12) cv[12] = real_or_complex(  5000.0,  1.0, cv[0]);
  if (m > 13) cv[13] = real_or_complex( 10000.0,  1.0, cv[0]);
  if (m > 14) cv[14] = real_or_complex( 80000.0,  1.0, cv[0]);
  if (m > 15) cv[15] = real_or_complex(100000.0,  1.0, cv[0]);
  gmm::clear(m2);
  for (size_type l = 0; l < m; ++l) m2(l, l) = cv[l];
  gmm::mult(cq, m2, ca); 
  gmm::mult(ca, cr, ca);
  
  implicit_qr_algorithm(ca, eigc, cq);

  for (size_type l = 0; l < m; ++l) {
    bool found = false;
     for (size_type k = 0; k < m; ++k)
       if (dal::abs(eigc[l] - cv[k]) < sqrt(sqrt(tol))*(dal::abs(eigc[l])+1.0))
	 { cv[k] = -1.123236; found = true; break; }
     if (found == false) {
       cerr << "Eigenvalue " << l << " not found\n" << std::flush;
       DAL_THROW(dal::failure_error, "Error on QR algorithm.");
     }

     std::vector vy(m);
     gmm::mult(ca, gmm::mat_col(cq, l),
	       gmm::scaled(gmm::mat_col(cq, l), -eigcr[l]), vy);
     error = gmm::vect_norm2(vy);
     if (error >= R(prec * 10000.0)) 
       DAL_THROW(dal::failure_error, "Error too large: " << error);
  }


  //
  // Test for symmetric_qr_algorithm
  //

  m = gmm::mat_nrows(m2);
  gmm::dense_matrix<T> cq(m, m), cr(m, m), ca(m, m);  
  std::vector<R> cvr(m), eigcr(m);
  gmm::fill_random(cr);
  gmm::qr_factor(cr, cq, ca);
  gmm::fill_random(cvr);

  gmm::copy(gmm::identity_matrix(), m2);

  if (m >  0) cvr[ 0] = T(R(     0.0 ));
  if (m >  1) cvr[ 1] = T(R(     0.0 ));
  if (m >  2) cvr[ 2] = T(R(     0.01));
  if (m >  3) cvr[ 3] = T(R(     0.01));
  if (m >  4) cvr[ 4] = T(R(    -2.0 ));
  if (m >  5) cvr[ 5] = T(R(    -2.0 ));
  if (m >  6) cvr[ 6] = T(R(   -50.0 ));
  if (m >  7) cvr[ 7] = T(R(   100.0 ));
  if (m >  8) cvr[ 8] = T(R(   300.0 ));
  if (m >  9) cvr[ 9] = T(R(   500.0 ));
  if (m > 10) cvr[10] = T(R(  1000.0 ));
  if (m > 11) cvr[11] = T(R(  4000.0 ));
  if (m > 12) cvr[12] = T(R(  5000.0 ));
  if (m > 13) cvr[13] = T(R( 10000.0 ));
  if (m > 14) cvr[14] = T(R( 80000.0 ));
  if (m > 15) cvr[15] = T(R(100000.0 ));
  gmm::clear(m2);
  for (size_type l = 0; l < m; ++l) m2(l, l) = cv[l];

  gmm::mult(gmm::conjugated(cq), m2, ca); 
  gmm::mult(ca, cq, ca);
  
  symmetric_qr_algorithm(ca, eigcr, cq);

  for (size_type l = 0; l < m; ++l) {
    bool found = false;
     for (size_type k = 0; k < m; ++k)
       if (dal::abs(eigcr[l]-cvr[k])< sqrt(sqrt(tol))*(dal::abs(eigcr[l])+1.0))
	 { cvr[k] = -1.123236; found = true; break; }
     if (found == false) {
       cerr << "Eigenvalue " << l << " not found\n" << std::flush;
       DAL_THROW(dal::failure_error, "Error on QR algorithm.");
     }

     std::vector vy(m);
     gmm::mult(ca, gmm::mat_col(cq, l),
	       gmm::scaled(gmm::mat_col(cq, l), -eigcr[l]), vy);
     error = gmm::vect_norm2(vy);
     if (error >= R(prec * 10000.0)) 
       DAL_THROW(dal::failure_error, "Error too large: " << error);

  }



}
