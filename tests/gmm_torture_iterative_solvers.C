// SQUARED_MATRIX_PARAM;
// DENSE_VECTOR_PARAM;
// VECTOR_PARAM;
// ENDPARAM;

#include <gmm.h>

using gmm::size_type;

bool print_debug = false;

template <typename MAT1, typename VECT1, typename VECT2>
void test_procedure(const MAT1 &_m1, const VECT1 &_v1, const VECT2 &_v2) {
  VECT1 &v1 = const_cast<VECT1 &>(_v1);
  VECT2 &v2 = const_cast<VECT2 &>(_v2);
  MAT1  &m1 = const_cast<MAT1  &>(_m1);
  typedef typename gmm::linalg_traits<MAT1>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  R prec = gmm::default_tol(R());

  gmm::clean(v1, 0.01);
  for (size_type i = 0; i < vect_size(v1); ++i)
    if (gmm::abs(v1[i]) < R(1) / R(100))
      DAL_THROW(gmm::failure_error, "Error in clean");

  static int nexpe = 0;
  if (print_debug) 
    { cout << "Begin experiment " << ++nexpe << "\n\nwith " << m1 << "\n\n"; }

  size_type m = gmm::mat_nrows(m1);
  std::vector<T> v3(m);

  R det = gmm::abs(gmm::lu_det(m1)), error;
  R cond = gmm::condest(m1);

  if (print_debug)
    cout << "condition number = " << cond << " det = " << det << endl;
  if (det == R(0) && cond < R(1) / prec && cond != R(0))
    DAL_THROW(gmm::failure_error, "Inconsistent condition number: " << cond);

  if (prec * cond < R(1)/R(10000) && det != R(0)) {

    gmm::identity_matrix P1;
    gmm::diagonal_precond<MAT1> P2(m1);
    gmm::mr_approx_inverse_precond<MAT1> P3(m1, 10, prec);
    gmm::ilu_precond<MAT1> P4(m1);
    gmm::ilut_precond<MAT1> P5(m1, 10, prec);

    if (print_debug)
      cout << "\nTest for bicgstab with no preconditionner\n";

    gmm::fill_random(v1);
    gmm::iteration iter((double(prec*cond))*100.0, print_debug ? 1:0, 1000*m);
    gmm::bicgstab(m1, v1, v2, P1, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for bicgstab with diagonal preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::bicgstab(m1, v1, v2, P2, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for bicgstab with mr preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::bicgstab(m1, v1, v2, P3, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for bicgstab with ilu preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::bicgstab(m1, v1, v2, P4, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for bicgstab with ilut preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::bicgstab(m1, v1, v2, P5, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for gmres with no preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::gmres(m1, v1, v2, P1, 50, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for gmres with diagonal preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::gmres(m1, v1, v2, P2, 50, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for gmres with mr preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::gmres(m1, v1, v2, P3, 50, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for gmres with ilu preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::gmres(m1, v1, v2, P4, 50, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for gmres with ilut preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::gmres(m1, v1, v2, P5, 50, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for qmr with no preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::qmr(m1, v1, v2, P1, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for qmr with diagonal preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::qmr(m1, v1, v2, P2, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for qmr with ilu preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::qmr(m1, v1, v2, P4, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for qmr with ilut preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::qmr(m1, v1, v2, P5, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);


    
    gmm::dense_matrix<T> m2(m, m), m3(m, m);
    gmm::mult(gmm::conjugated(m1), m1, m2);
    gmm::copy(m1, m3);
    gmm::add(gmm::conjugated(m1), m3);
    gmm::copy(m2, m1);
    gmm::cholesky_precond<MAT1> P6(m1);
    gmm::choleskyt_precond<MAT1> P7(m1, 10, prec);
    
    if (print_debug)
      cout << "\nTest for cg with no preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    iter.set_resmax(double(prec*cond*cond)*10.0);
    gmm::cg(m1, v1, v2, P1, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);
    if (!is_hermitian(m1))
      DAL_THROW(gmm::failure_error, "The matrix is not hermitian");

    if (print_debug)
      cout << "\nTest for cg with diagonal preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::cg(m1, v1, v2, P2, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for cg with ildlt preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::cg(m1, v1, v2, P6, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    if (print_debug)
      cout << "\nTest for cg with ildltt preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::cg(m1, v1, v2, P7, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);

    
    if (print_debug)
      cout << "\nTest for gmres with ildltt preconditionner\n";
    gmm::copy(m3, m1);
    if (!is_hermitian(m1))
      DAL_THROW(gmm::failure_error, "The matrix is not hermitian");
    cond = gmm::condest(m1);
    if (print_debug)
      cout << "condition number for this experiment= " << cond << endl;
    gmm::choleskyt_precond<MAT1> P8(m1, 10, prec);

    iter.init(); gmm::fill_random(v1);
    iter.set_resmax(double(prec*cond)*10.0);
    gmm::gmres(m1, v1, v2, P8, 50, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    if (!(error <= prec * cond * R(20000)))
      DAL_THROW(gmm::failure_error, "Error too large: " << error);
    
  }
}
