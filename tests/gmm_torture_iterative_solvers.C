// SQUARED_MATRIX_PARAM;
// DENSE_VECTOR_PARAM;
// VECTOR_PARAM;
// ENDPARAM;

#include <gmm.h>

using gmm::size_type;

bool print_debug = false;
int nb_fault_allowed = 20;

template <typename MAT1, typename VECT1, typename VECT2>
void test_procedure(const MAT1 &_m1, const VECT1 &_v1, const VECT2 &_v2) {
  VECT1 &v1 = const_cast<VECT1 &>(_v1);
  VECT2 &v2 = const_cast<VECT2 &>(_v2);
  MAT1  &m1 = const_cast<MAT1  &>(_m1);
  typedef typename gmm::linalg_traits<MAT1>::value_type T;
  typedef typename gmm::number_traits<T>::magnitude_type R;
  R prec = gmm::default_tol(R());
  static int nb_fault = 0;
  static int bicgstab_fault(0), gmres_fault(0), qmr_fault(0), cg_fault(0);
  static int ident_fault(0), diagonal_fault(0), ilu_fault(0), ilut_fault(0);
  static int ildlt_fault(0), ildltt_fault(0), mr_fault(0);

  static int bicgstab_nb(0), gmres_nb(0), qmr_nb(0), cg_nb(0);
  static int ident_nb(0), diagonal_nb(0), ilu_nb(0), ilut_nb(0);
  static int ildlt_nb(0), ildltt_nb(0), mr_nb(0);

  static int bicgstab_nb_iter(0), gmres_nb_iter(0), qmr_nb_iter(0);
  static int cg_nb_iter(0), ident_nb_iter(0), diagonal_nb_iter(0);
  static int ilu_nb_iter(0), ilut_nb_iter(0), ildlt_nb_iter(0);
  static int ildltt_nb_iter(0), mr_nb_iter(0);
  static int nexpe = 0;
  ++nexpe;

  gmm::clean(v1, 0.01);
  for (size_type i = 0; i < gmm::vect_size(v1); ++i)
    if (v1[i] != T(0) && gmm::abs(v1[i]) < R(1) / R(101))
      DAL_THROW(gmm::failure_error, "Error in clean");

  if (print_debug) {
    cout << "Begin experiment " << nexpe << "\n\nwith " << m1 << "\n\n"; 
    dal::set_warning_level(3);
  }

  size_type m = gmm::mat_nrows(m1);
  std::vector<T> v3(m);

  R det = gmm::abs(gmm::lu_det(m1)), error;
  R cond = gmm::condest(m1);

  if (print_debug)
    cout << "condition number = " << cond << " det = " << det << endl;
  if (det == R(0) && cond < R(1) / prec && cond != R(0))
    DAL_THROW(gmm::failure_error, "Inconsistent condition number: " << cond);

  if (sqrt(prec) * cond < R(1)/R(100) && det > sqrt(prec)*R(10)) {

    gmm::identity_matrix P1;
    gmm::diagonal_precond<MAT1> P2(m1);
    gmm::mr_approx_inverse_precond<MAT1> P3(m1, 10, prec);
    gmm::ilu_precond<MAT1> P4(m1);
    gmm::ilut_precond<MAT1> P5(m1, 10, prec);

    R detmr = gmm::abs(gmm::lu_det(P3.approx_inverse()));

    if (print_debug)
      cout << "\nTest for bicgstab with no preconditionner\n";

    gmm::fill_random(v1);
    gmm::iteration iter((double(prec*cond))*100.0, print_debug ? 1:0, 300*m);
    gmm::bicgstab(m1, v1, v2, P1, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++bicgstab_nb; ++ident_nb;
    bicgstab_nb_iter += iter.get_iteration();
    ident_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * R(20000))) {
      ++nb_fault; ++bicgstab_fault; ++ident_fault;
      if (nb_fault > nb_fault_allowed)
      DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    if (print_debug)
      cout << "\nTest for bicgstab with diagonal preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::bicgstab(m1, v1, v2, P2, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++bicgstab_nb; ++diagonal_nb;
    bicgstab_nb_iter += iter.get_iteration();
    diagonal_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * R(20000))) {
      ++nb_fault; ++bicgstab_fault; ++diagonal_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    if (detmr > prec * R(100)) {
      if (print_debug)
	cout << "\nTest for bicgstab with mr preconditionner\n";
      
      iter.init(); gmm::fill_random(v1);
      gmm::bicgstab(m1, v1, v2, P3, iter);
      gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
      error = gmm::vect_norm2(v3);
      ++bicgstab_nb; ++mr_nb;
      bicgstab_nb_iter += iter.get_iteration();
      mr_nb_iter += iter.get_iteration();
      if (!(error <= prec * cond * R(20000))) {
	++nb_fault; ++bicgstab_fault; ++mr_fault;
	if (nb_fault > nb_fault_allowed)
	  DAL_THROW(gmm::failure_error, "Error too large: " << error);
      }
    }

    if (print_debug)
      cout << "\nTest for bicgstab with ilu preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::bicgstab(m1, v1, v2, P4, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++bicgstab_nb;  ++ilu_nb;
    bicgstab_nb_iter += iter.get_iteration();
    ilu_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * R(20000))) {
      ++nb_fault; ++bicgstab_fault;  ++ilu_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    if (print_debug)
      cout << "\nTest for bicgstab with ilut preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::bicgstab(m1, v1, v2, P5, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++bicgstab_nb; ++ilut_nb;
    bicgstab_nb_iter += iter.get_iteration();
    ilut_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * R(20000))) {
      ++nb_fault; ++bicgstab_fault; ++ilut_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    if (print_debug)
      cout << "\nTest for gmres with no preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::gmres(m1, v1, v2, P1, 50, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++gmres_nb; ++ident_nb;
    gmres_nb_iter += iter.get_iteration();
    ident_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * R(20000))) {
      ++nb_fault; ++gmres_fault; ++ident_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    if (print_debug)
      cout << "\nTest for gmres with diagonal preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::gmres(m1, v1, v2, P2, 50, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++gmres_nb; ++diagonal_nb;
    gmres_nb_iter += iter.get_iteration();
    diagonal_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * R(20000))) {
      ++nb_fault; ++gmres_fault; ++diagonal_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    if (detmr > prec * R(100)) {
      if (print_debug)
	cout << "\nTest for gmres with mr preconditionner\n";
      
      iter.init(); gmm::fill_random(v1);
      gmm::gmres(m1, v1, v2, P3, 50, iter);
      gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
      error = gmm::vect_norm2(v3);
      ++gmres_nb; ++mr_nb;
      gmres_nb_iter += iter.get_iteration();
      mr_nb_iter += iter.get_iteration();
      if (!(error <= prec * cond * R(20000))) {
	++nb_fault; ++gmres_fault; ++mr_fault;
	if (nb_fault > nb_fault_allowed)
	  DAL_THROW(gmm::failure_error, "Error too large: " << error);
      }
    }

    if (print_debug)
      cout << "\nTest for gmres with ilu preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::gmres(m1, v1, v2, P4, 50, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++gmres_nb; ++ilu_nb;
    gmres_nb_iter += iter.get_iteration();
    ilu_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * R(20000))) {
      ++nb_fault; ++gmres_fault; ++ilu_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    if (print_debug)
      cout << "\nTest for gmres with ilut preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::gmres(m1, v1, v2, P5, 50, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++gmres_nb; ++ilut_nb;
    gmres_nb_iter += iter.get_iteration();
    ilut_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * R(20000))) {
      ++nb_fault; ++gmres_fault; ++ilut_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    if (print_debug)
      cout << "\nTest for qmr with no preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::qmr(m1, v1, v2, P1, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++qmr_nb; ++ident_nb;
    qmr_nb_iter += iter.get_iteration();
    ident_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * R(20000))) {
      ++nb_fault; ++qmr_fault; ++ident_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    if (print_debug)
      cout << "\nTest for qmr with diagonal preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::qmr(m1, v1, v2, P2, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++qmr_nb; ++diagonal_nb;
    qmr_nb_iter += iter.get_iteration();
    diagonal_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * R(20000))) {
      ++nb_fault; ++qmr_fault; ++diagonal_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    if (print_debug)
      cout << "\nTest for qmr with ilu preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::qmr(m1, v1, v2, P4, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++qmr_nb; ++ilu_nb;
    qmr_nb_iter += iter.get_iteration();
    ilu_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * R(20000))) {
      ++nb_fault; ++qmr_fault; ++ilu_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    if (print_debug)
      cout << "\nTest for qmr with ilut preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::qmr(m1, v1, v2, P5, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++qmr_nb; ++ilut_nb;
    qmr_nb_iter += iter.get_iteration();
    ilut_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * R(20000))) {
      ++nb_fault; ++qmr_fault; ++ilut_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }


    
    gmm::dense_matrix<T> m2(m, m), m3(m, m);
    gmm::mult(gmm::conjugated(m1), m1, m2);
    gmm::copy(m1, m3);
    gmm::add(gmm::conjugated(m1), m3);
    gmm::copy(m2, m1);
    gmm::cholesky_precond<MAT1> P6(m1);
    gmm::choleskyt_precond<MAT1> P7(m1, 10, prec);
    
    if (!is_hermitian(m1))
      DAL_THROW(gmm::failure_error, "The matrix is not hermitian");

    if (print_debug)
      cout << "\nTest for cg with no preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    iter.set_resmax(double(prec*cond*cond)*10.0);
    gmm::cg(m1, v1, v2, P1, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++cg_nb; ++ident_nb;
    cg_nb_iter += iter.get_iteration();
    ident_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * cond * R(20000))) {
      ++nb_fault; ++cg_fault; ++ident_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    if (print_debug)
      cout << "\nTest for cg with diagonal preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::cg(m1, v1, v2, P2, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++cg_nb; ++diagonal_nb;
    cg_nb_iter += iter.get_iteration();
    diagonal_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * cond * R(20000))) {
      ++nb_fault; ++cg_fault; ++diagonal_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    if (print_debug)
      cout << "\nTest for cg with ildlt preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::cg(m1, v1, v2, P6, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++cg_nb; ++ildlt_nb;
    cg_nb_iter += iter.get_iteration();
    ildlt_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * cond * R(20000))) {
      ++nb_fault; ++cg_fault; ++ildlt_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    if (print_debug)
      cout << "\nTest for cg with ildltt preconditionner\n";

    iter.init(); gmm::fill_random(v1);
    gmm::cg(m1, v1, v2, P7, iter);
    gmm::mult(m1, v1, gmm::scaled(v2, T(-1)), v3);
    error = gmm::vect_norm2(v3);
    ++cg_nb; ++ildltt_nb;
    cg_nb_iter += iter.get_iteration();
    ildltt_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * cond * R(20000))) {
      ++nb_fault; ++cg_fault; ++ildltt_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }

    
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
    ++gmres_nb; ++ildltt_nb;
    gmres_nb_iter += iter.get_iteration();
    ildltt_nb_iter += iter.get_iteration();
    if (!(error <= prec * cond * R(20000))) {
      ++nb_fault; ++gmres_fault; ++ildltt_fault;
      if (nb_fault > nb_fault_allowed)
	DAL_THROW(gmm::failure_error, "Error too large: " << error);
    }
    
  }

  if (nb_fault > 0 && (nexpe == 100 || print_debug)) {
    cerr << "Convergence failed for:\n";
    if (bicgstab_fault > 0) cerr << bicgstab_fault << " times for bicgstab\n";
    if (gmres_fault > 0) cerr << gmres_fault << " times for gmres\n";
    if (qmr_fault > 0) cerr << qmr_fault << " times for qmr\n";
    if (cg_fault > 0) cerr << cg_fault << " times for cg\n";
    if (ident_fault > 0)
      cerr << ident_fault << " times with no preconditonner\n";
    if (diagonal_fault > 0)
      cerr << diagonal_fault << " times with diagonal preconditonner\n";
    if (mr_fault > 0) cerr << mr_fault << " times with mr preconditonner\n";
    if (ilu_fault > 0) cerr << ilu_fault << " times with ilu preconditonner\n";
    if (ilut_fault > 0)
      cerr << ilut_fault << " times with ilut preconditonner\n";
    if (ildlt_fault > 0)
      cerr << ildlt_fault << " times with ildlt preconditonner\n";
    if (ildltt_fault > 0)
      cerr << ildltt_fault << " times with ildltt preconditonner\n";
  }

  if (nexpe == 100) {
    cout << "\n mean number of iterations\n";
    if (bicgstab_nb) 
      cout << "bicgstab: " << double(bicgstab_nb_iter) / double(bicgstab_nb)
	   << endl;
    if (gmres_nb) 
      cout << "gmres: " << double(gmres_nb_iter) / double(gmres_nb)
	   << endl;
    if (qmr_nb) 
      cout << "qmr: " << double(qmr_nb_iter) / double(qmr_nb)
	   << endl;
    if (cg_nb) 
      cout << "cg: " << double(cg_nb_iter) / double(cg_nb)
	   << endl;
    if (ident_nb) 
      cout << "no precond: " << double(ident_nb_iter) / double(ident_nb)
	   << endl;
    if (diagonal_nb) 
      cout << "diagonal precond: "
	   << double(diagonal_nb_iter) / double(diagonal_nb) << endl;
    if (mr_nb) 
      cout << "mr precond: " << double(mr_nb_iter) / double(mr_nb)
	   << endl;
    if (ilu_nb) 
      cout << "ilu precond: " << double(ilu_nb_iter) / double(ilu_nb)
	   << endl;
    if (ilut_nb) 
      cout << "ilut precond: " << double(ilut_nb_iter) / double(ilut_nb)
	   << endl;
    if (ildlt_nb) 
      cout << "ildlt precond: " << double(ildlt_nb_iter) / double(ildlt_nb)
	 << endl;
    if (ildltt_nb) 
      cout << "ildltt precond: " << double(ildltt_nb_iter) / double(ildltt_nb)
	   << endl;
  }

}
