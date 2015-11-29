/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2011-2015 Tomas Ligursky, Yves Renard, Konstantinos Poulios

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/** @file getfem_continuation.h
    @author Tomas Ligursky <tomas.ligursky@gmail.com>
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @author Konstantinos Poulios <logari81@googlemail.com>
    @date October 17, 2011.
    @brief Inexact Moore-Penrose continuation method.
*/
#ifndef GETFEM_CONTINUATION_H__
#define GETFEM_CONTINUATION_H__

#include <getfem/getfem_model_solvers.h>

namespace getfem {


  //=========================================================================
  // Abstract Moore-Penrose continuation method
  //=========================================================================

  template <typename VECT, typename MAT>
  class virtual_cont_struct {

  protected:
#ifdef _MSC_VER
    const double tau_bp_init = 1.e6;
    const double diffeps = 1.e-8;
#else
    static constexpr double tau_bp_init = 1.e6;
    static constexpr double diffeps = 1.e-8;
#endif

    int singularities;

  private:

    bool non_smooth;
    double scfac, h_init_, h_max_, h_min_, h_inc_, h_dec_;
    size_type maxit_, thrit_;
    double maxres_, maxdiff_, mincos_, delta_max_, delta_min_, thrvar_;
    size_type nbdir_, nbspan_;
    int noisy_;
    double tau_lp, tau_bp_1, tau_bp_2;

    // stored singularities info
    std::map<double, double> tau_bp_graph;
    VECT alpha_hist, tau_bp_hist;
    std::string sing_label;
    VECT x_sing, x_next;
    double gamma_sing, gamma_next;
    std::vector<VECT> tx_sing, tx_predict;
    std::vector<double> tgamma_sing, tgamma_predict;

    // randomized data
    VECT bb_x_, cc_x_;
    double bb_gamma, cc_gamma, dd;

  public:
    /* Compute a unit tangent at (x, gamma) that is accute to the incoming
       vector. */
    void compute_tangent(const VECT &x, double gamma,
                         VECT &tx, double &tgamma) {
      VECT g(x), y(x);
      F_gamma(x, gamma, g);                     // g = F_gamma(x, gamma)
      solve_grad(x, gamma, y, g);               // y = F_x(x, gamma)^-1 * g
      tgamma = 1. / (tgamma - w_sp(tx, y));
      gmm::copy(gmm::scaled(y, -tgamma), tx);   // tx = -tgamma * y

      scale(tx, tgamma, 1./w_norm(tx, tgamma)); // [tx,tgamma] /=  w_norm(tx,tgamma)

      mult_grad(x, gamma, tx, y);               // y = F_x(x, gamma) * tx
      gmm::add(gmm::scaled(g, tgamma), y);      // y += tgamma * g
      double r = norm(y);
      if (r > 1.e-10)
        GMM_WARNING1("Tangent computed with the residual " << r);
    }

  private:

    /* Calculate a tangent vector at (x, gamma) + h * (tX, tGamma) and test
       whether it is close to (tX, tGamma). Informatively, compare it with
       (tx, tgamma), as well. */
    bool test_tangent(const VECT &x, double gamma,
                      const VECT &tX, double tGamma,
                      const VECT &tx, double tgamma, double h) {
      bool res = false;
      double Gamma1, tGamma1(tgamma);
      VECT X1(x), tX1(tx);

      scaled_add(x, gamma, tX, tGamma, h, X1, Gamma1); // [X1,Gamma1] = [x,gamma] + h * [tX,tGamma]
      compute_tangent(X1, Gamma1, tX1, tGamma1);

      double cang = cosang(tX1, tX, tGamma1, tGamma);
      if (noisy() > 1)
        cout << "cos of the angle with the tested tangent " << cang << endl;
      if (cang >= mincos())
        res = true;
      else {
        cang = cosang(tX1, tx, tGamma1, tGamma);
        if (noisy() > 1)
          cout << "cos of the angle with the initial tangent " << cang << endl;
      }
      return res;
    }

    /* Simple tangent switch. */
    bool switch_tangent(const VECT &x, double gamma,
                        VECT &tx, double &tgamma, double &h) {
      double Gamma, tGamma(tgamma);
      VECT X(x), tX(tx);

      if (noisy() > 0) cout << "Trying a simple tangent switch" << endl;
      if (noisy() > 0) cout << "Starting computing a new tangent" << endl;
      h *= 1.5;
      scaled_add(x, gamma, tx, tgamma, h, X, Gamma);  // [X,Gamma] = [x,gamma] + h * [tx,tgamma]
      compute_tangent(X, Gamma, tX, tGamma);
      // One can test the cosine of the angle between (tX, tGamma) and
      // (tx, tgamma), for sure, and increase h_min if it were greater or
      // equal to mincos(). However, this seems to be superfluous.

      if (noisy() > 0)
        cout << "Starting testing the computed tangent" << endl;
      double h_test = -0.9 * h_min();
      bool accepted(false);
      while (!accepted && (h_test > -h_max())) {
        h_test = -h_test
                 + pow(10., floor(log10(-h_test / h_min()))) * h_min();
        accepted = test_tangent(x, gamma, tX, tGamma, tx, tgamma, h_test);
        if (!accepted) {
          h_test *= -1.;
          accepted = test_tangent(x, gamma, tX, tGamma, tx, tgamma, h_test);
        }
      }

      if (accepted) {
        if (h_test < 0) {
          gmm::scale(tX, -1.);
          tGamma *= -1.;
          h_test *= -1.;
        }
        if (noisy() > 0)
          cout << "Tangent direction switched, "
               << "starting computing a suitable step size" << endl;
        h = h_init();
        bool h_adapted = false;
        while (!h_adapted && (h > h_test)) {
          h_adapted = test_tangent(x, gamma, tX, tGamma, tx, tgamma, h);
          h *= h_dec();
        }
        h = h_adapted ? h / h_dec() : h_test;
        copy(tX, tGamma, tx, tgamma);
      } else
        if (noisy() > 0) cout << "Simple tangent switch has failed!" << endl;

      return accepted;
    }

    /* Test for limit points (also called folds or turning points). */
    bool test_limit_point(double tgamma) {
      double tau_lp_old = get_tau_lp();
      set_tau_lp(tgamma);
      return (tgamma * tau_lp_old < 0);
    }

    void init_test_functions(const VECT &x, double gamma,
                             const VECT &tx, double tgamma) {
      set_tau_lp(tgamma);
      if (this->singularities > 1) {
        if (noisy() > 0) cout << "Starting computing an initial value of the "
                              << "test function for bifurcations" << endl;
        set_tau_bp_2(test_function_bp(x, gamma, tx, tgamma));
      }
    }

    /* Test function for bifurcation points for a given matrix. The first part
       of the solution of the augmented system is passed in
       (v_x, v_gamma). */
    double test_function_bp(const MAT &A, const VECT &g,
                            const VECT &tx, double tgamma,
                            VECT &v_x, double &v_gamma) {
      VECT y(g), z(g);
      size_type nn = gmm::vect_size(g);

      solve(A, y, z, g, bb_x(nn));              // [y,z] = A^-1 * [g,bb_x]
      v_gamma = (bb_gamma - sp(tx, z)) / (tgamma - sp(tx, y));
      scaled_add(z, y, -v_gamma, v_x);          // v_x = y - v_gamma*z
      double tau = 1. / (dd - sp(cc_x(nn), v_x) - cc_gamma * v_gamma);
      scale(v_x, v_gamma, -tau);                // [v_x,v_gamma] *= -tau

      // control of the norm of the residual
      mult(A, v_x, y);
      gmm::add(gmm::scaled(g, v_gamma), y);     // y += v_gamma*g
      gmm::add(gmm::scaled(bb_x(nn), tau), y);  // y += bb_x*tau
      double r = sp(tx, v_x) + tgamma * v_gamma + bb_gamma * tau;
      double q = sp(cc_x(nn), v_x) + cc_gamma * v_gamma + dd * tau - 1.;
      r = sqrt(sp(y, y) + r * r + q * q);
      if (r > 1.e-10)
        GMM_WARNING1("Test function evaluated with the residual " << r);

      return tau;
    }

    double test_function_bp(const MAT &A, const VECT &g,
                            const VECT &tx, double tgamma) {
      VECT v_x(g); double v_gamma;
      return test_function_bp(A, g, tx, tgamma, v_x, v_gamma);
    }

    /* Test function for bifurcation points for the gradient computed at
       (x, gamma). */
    double test_function_bp(const VECT &x, double gamma,
                            const VECT &tx, double tgamma,
                            VECT &v_x, double &v_gamma) {
      MAT A; VECT g(x);
      F_x(x, gamma, A);
      F_gamma(x, gamma, g);
      return test_function_bp(A, g, tx, tgamma, v_x, v_gamma);
    }

    double test_function_bp(const VECT &x, double gamma,
                            const VECT &tx, double tgamma) {
      VECT v_x(x); double v_gamma;
      return test_function_bp(x, gamma, tx, tgamma, v_x, v_gamma);
    }

  public:
    /* Test for smooth bifurcation points. */
    bool test_smooth_bifurcation(const VECT &x, double gamma,
                                 const VECT &tx, double tgamma) {
      double tau_bp_1_old = get_tau_bp_1();
      set_tau_bp_1(get_tau_bp_2());
      set_tau_bp_2(test_function_bp(x, gamma, tx, tgamma));
      return (get_tau_bp_2() * get_tau_bp_1() < 0) &&
             (gmm::abs(get_tau_bp_1()) < gmm::abs(tau_bp_1_old));
    }

    /* Test for non-smooth bifurcation points. */
    bool test_nonsmooth_bifurcation(const VECT &x1, double gamma1,
                                    const VECT &tx1, double tgamma1,
                                    const VECT &x2, double gamma2,
                                    const VECT &tx2, double tgamma2) {
      VECT g1(x1), g2(x1), g(x1), tx(x1);

      // compute gradients at the two given points
      MAT A1, A2;
      F_x(x2, gamma2, A2);
      F_gamma(x2, gamma2, g2);
      F_x(x1, gamma1, A1);
      F_gamma(x1, gamma1, g1);
      double tau1 = test_function_bp(A1, g1, tx1, tgamma1);
      double tau2 = test_function_bp(A2, g2, tx2, tgamma2);
      double tau_var_ref = std::max(gmm::abs(tau2 - tau1), 1.e-8);
      set_tau_bp_2(tau1);
      init_tau_bp_graph();
      MAT A(A2);

      // monitor the sign changes of the test function on the convex
      // combination
      size_type nb_changes = 0;
      double delta = delta_min(), tau0 = tau_bp_init, tgamma;
      for (double alpha=0.; alpha < 1.; ) {
        alpha = std::min(alpha + delta, 1.);
        scaled_add(A1, 1.-alpha, A2, alpha, A); // A = (1-alpha)*A1 + alpha*A2
        scaled_add(g1, 1.-alpha, g2, alpha, g); // g = (1-alpha)*g1 + alpha*g2
        scaled_add(tx1, tgamma1, 1.-alpha, tx2, tgamma2, alpha, tx, tgamma);
        //[tx,tgamma] = (1-alpha)*[tx1,tgamma1] + alpha*[tx2,tgamma2]

        tau2 = test_function_bp(A, g, tx, tgamma);
        if ((tau2 * tau1 < 0) && (gmm::abs(tau1) < gmm::abs(tau0)))
          ++nb_changes;
        insert_tau_bp_graph(alpha, tau2);

        if (gmm::abs(tau2 - tau1) < 0.5 * thrvar() * tau_var_ref)
          delta = std::min(2 * delta, delta_max());
        else if (gmm::abs(tau2 - tau1) > thrvar() * tau_var_ref)
          delta = std::max(0.1 * delta, delta_min());
        tau0 = tau1;
        tau1 = tau2;
      }

      set_tau_bp_1(tau_bp_init);
      set_tau_bp_2(tau2);
      return nb_changes % 2;
    }

  private:
    /* Newton-type corrections for the couple ((X, Gamma), (tX, tGamma)).
       The current direction of (tX, tGamma) is informatively compared with
       (tx, tgamma). */
    bool newton_corr(VECT &X, double &Gamma, VECT &tX,
                     double &tGamma, const VECT &tx, double tgamma,
                     size_type &it) {
      bool converged = false;
      double Delta_Gamma, res(0), diff;
      VECT f(X), g(X), Delta_X(X), y(X);

      if (noisy() == 1) cout << "Starting correction" << endl;
      F(X, Gamma, f);                                          // f = F(X, Gamma) = -rhs(X, Gamma)
//CHANGE 1: line search
//double res0 = norm(f);

      for (it=0; it < maxit() && res < 1.e8; ++it) {
        F_gamma(X, Gamma, f, g);                               // g = F_gamma(X, Gamma)
        solve_grad(X, Gamma, Delta_X, y, f, g);                // y = F_x(X, Gamma)^-1 * g
                                                               // Delta_X = F_x(X, Gamma)^-1 * f
        Delta_Gamma = sp(tX, Delta_X) / (sp(tX, y) - tGamma);  // Delta_Gamma = tX.Delta_X / (tX.y - tGamma)
        if (isnan(Delta_Gamma)) {
          if (noisy() > 0) cout << "Newton correction failed with NaN" << endl;
          return false;
        }
        gmm::add(gmm::scaled(y, -Delta_Gamma), Delta_X);       // Delta_X -= Delta_Gamma * y
        scaled_add(X, Gamma, Delta_X, Delta_Gamma, -1,
                   X, Gamma);                                  // [X,Gamma] -= [Delta_X,Delta_Gamma]
        F(X, Gamma, f);                                        // f = F(X, gamma) = -rhs(X, gamma)
        res = norm(f);

//CHANGE 1: line search
//for (size_type ii=0; ii < 4 && (isnan(res) || res > res0); ++ii) { // some basic linesearch
//  scale(Delta_X, Delta_Gamma, 0.5);
//  scaled_add(X, Gamma, Delta_X, Delta_Gamma, 1, X, Gamma);     // [X,Gamma] += [Delta_X,Delta_Gamma]
//  F(X, Gamma, f);                                              // f = F(X, gamma) = -rhs(X, gamma)
//  res = norm(f);
//}

        tGamma = 1. / (tGamma - w_sp(tX, y));                  // tGamma = 1 / (tGamma - k*tX.y)
        gmm::copy(gmm::scaled(y, -tGamma), tX);                // tX = -tGamma * y
        scale(tX, tGamma, 1./w_norm(tX, tGamma));              // [tX,tGamma] /= w_norm(tX,tGamma)

        diff = w_norm(Delta_X, Delta_Gamma);
        if (noisy() > 1)
          cout << " Correction " << std::setw(3) << it << ":"
          << " Gamma = " << std::fixed << std::setprecision(6) << Gamma
          << " residual = " << std::scientific << std::setprecision(3) << res
          << " difference = " << std::scientific << std::setprecision(3) << diff
          << " cosang = " << std::fixed << std::setprecision(6)
                          << cosang(tX, tx, tGamma, tgamma) << endl;

        if (res <= maxres() && diff <= maxdiff()) {
          converged = true;
          // recalculate the final tangent, for sure
          compute_tangent(X, Gamma, tX, tGamma);
          break;
        }
      }
      if (noisy() == 1) cout << "Correction finished with Gamma = "
                             << Gamma << endl;
      return converged;
    }

    bool newton_corr(VECT &X, double &Gamma, VECT &tX,
                     double &tGamma, const VECT &tx, double tgamma) {
      size_type it;
      return newton_corr(X, Gamma, tX, tGamma, tx, tgamma, it);
    }

    /* Try to perform one predictor-corrector step starting from the couple
       ((x, gamma), (tx, tgamma)). Return the resulting couple in the case of
       convergence. */
    bool test_predict_dir(VECT &x, double &gamma,
                          VECT &tx, double &tgamma) {
      bool converged = false;
      double h = h_init(), Gamma, tGamma;
      VECT X(x), tX(x);
      while (!converged) { //step control
        // prediction
        scaled_add(x, gamma, tx, tgamma, h, X, Gamma);   // [X,Gamma] = [x,gamma] + h * [tx,tgamma]
        if (noisy() > 0)
          cout << "(TPD) Prediction   : Gamma = " << Gamma
               << " (for h = " << h << ", tgamma = " << tgamma << ")" << endl;
        copy(tx, tgamma, tX, tGamma);
        //correction
        converged = newton_corr(X, Gamma, tX, tGamma, tx, tgamma);

        if (h > h_min())
          h = std::max(0.199 * h_dec() * h, h_min());
        else
          break;
      }
      if (converged) {
        // check the direction of the tangent found
        scaled_add(X, Gamma, x, gamma, -1., tx, tgamma); // [tx,tgamma] = [X,Gamma] - [x,gamma]
        if (sp(tX, tx, tGamma, tgamma) < 0)
          scale(tX, tGamma, -1.);                        // [tX,tGamma] *= -1
        copy(X, Gamma, x, gamma);
        copy(tX, tGamma, tx, tgamma);
      }
      return converged;
    }

    /* A tool for approximating a smooth bifurcation point close to (x, gamma)
       and locating the two branches emanating from there. */
    void treat_smooth_bif_point(const VECT &x, double gamma,
                                const VECT &tx, double tgamma, double h) {
      double tau0(get_tau_bp_1()), tau1(get_tau_bp_2());
      double gamma0(gamma), Gamma,
             tgamma0(tgamma), tGamma(tgamma), v_gamma;
      VECT x0(x), X(x), tx0(tx), tX(tx), v_x(tx);

      if (noisy() > 0)
        cout  << "Starting locating the bifurcation point" << endl;

      // predictor-corrector steps with a secant-type step-length adaptation
      h *= tau1 / (tau0 - tau1);
      for (size_type i=0; i < 10 && (gmm::abs(h) >= h_min()); ++i) {
        scaled_add(x0, gamma0, tx0, tgamma0, h, X, Gamma); // [X,Gamma] = [x0,gamma0] + h * [tx0,tgamma0]
        if (noisy() > 0)
          cout << "(TSBP) Prediction   : Gamma = " << Gamma
               << " (for h = " << h << ", tgamma = " << tgamma << ")" << endl;
        if (newton_corr(X, Gamma, tX, tGamma, tx0, tgamma0)) {
          copy(X, Gamma, x0, gamma0);
          if (cosang(tX, tx0, tGamma, tgamma0) >= mincos())
            copy(tX, tGamma, tx0, tgamma0);
          tau0 = tau1;
          tau1 = test_function_bp(X, Gamma, tx0, tgamma0, v_x, v_gamma);
          h *= tau1 / (tau0 - tau1);
        } else {
          scaled_add(x0, gamma0, tx0, tgamma0, h, x0, gamma0); // [x0,gamma0] += h*[tx0,tgamma0]
          test_function_bp(x0, gamma0, tx0, tgamma0, v_x, v_gamma);
          break;
        }
      }
      if (noisy() > 0)
        cout  << "Bifurcation point located" << endl;
      set_sing_point(x0, gamma0);
      insert_tangent_sing(tx0, tgamma0);

      if (noisy() > 0)
        cout  << "Starting searching for the second branch" << endl;
      double no = w_norm(v_x, v_gamma);
      scale(v_x, v_gamma, 1./no);                               // [v_x,v_gamma] /= no
      if (test_predict_dir(x0, gamma0, v_x, v_gamma)
          && insert_tangent_sing(v_x, v_gamma))
        { if (noisy() > 0) cout << "Second branch found" << endl; }
      else if (noisy() > 0) cout << "Second branch not found!" << endl;
    }

  public:

    /* A tool for approximating a non-smooth point close to (x, gamma) and
       locating one-sided smooth solution branches emanating from there. It is
       supposed that (x, gamma) is a point on a smooth solution branch within
       the distance of h_min() from the end point of this branch and
       (tx, tgamma) is the corresponding tangent that is directed towards the
       end point. The boolean set_next determines whether the first new
       branch found (if any) is to be chosen for further continuation. */
    void treat_nonsmooth_point(const VECT &x, double gamma,
                               const VECT &tx, double tgamma, bool set_next) {
      double gamma0(gamma), Gamma(gamma);
      double tgamma0(tgamma), tGamma(tgamma);
      double h = h_min(), mcos = mincos();
      VECT x0(x), X(x), tx0(tx), tX(tx);

      // approximate the end point more precisely by a bisection-like algorithm
      if (noisy() > 0)
        cout  << "Starting locating a non-smooth point" << endl;

      scaled_add(x0, gamma0, tx0, tgamma0, h, X, Gamma);      // [X,Gamma] = [x0,gamma0] + h*[tx0,tgamma0]
      if (newton_corr(X, Gamma, tX, tGamma, tx0, tgamma0)) {  // --> X, Gamma, tX, tGamma
        double cang = cosang(tX, tx0, tGamma, tgamma0);
        if (cang >= mcos) mcos = (cang + 1.) / 2.;
      }

      copy(tx0, tgamma0, tX, tGamma);
      h /= 2.;
      for (size_type i = 0; i < 15; i++) {
        scaled_add(x0, gamma0, tx0, tgamma0, h, X, Gamma);    // [X,Gamma] = [x0,gamma0] + h*[tx0,tgamma0]
        if (noisy() > 0)
          cout << "(TNSBP) Prediction   : Gamma = " << Gamma
               << " (for h = " << h << ", tgamma = " << tgamma << ")" << endl;
        if (newton_corr(X, Gamma, tX, tGamma, tx0, tgamma0)
            && (cosang(tX, tx0, tGamma, tgamma0) >= mcos)) {
          copy(X, Gamma, x0, gamma0);
          copy(tX, tGamma, tx0, tgamma0);
        } else {
          copy(tx0, tgamma0, tX, tGamma);
        }
        h /= 2.;
      }
      if (noisy() > 0)
        cout  << "A non-smooth point located" << endl;
      set_sing_point(x0, gamma0);

      // take two reference vectors to span a subspace of directions emanating
      // from the end point
      if (noisy() > 0)
        cout << "Starting a thorough search for other branches" << endl;
      double tgamma1 = tgamma0, tgamma2 = tgamma0;
      VECT tx1(tx0), tx2(tx0);
      scale(tx1, tgamma1, -1.);                                     // [tx1,tgamma1] *= -1
      insert_tangent_sing(tx1, tgamma1);
      h = h_min();
      scaled_add(x0, gamma0, tx0, tgamma0, h, X, Gamma);            // [X,Gamma] = [x0,gamma0] + h*[tx0,tgamma0]
      compute_tangent(X, Gamma, tx2, tgamma2);

      // try systematically the directions of linear combinations of the couple
      // of the reference vectors for finding new possible tangent predictions
      // emanating from the end point
      size_type i1 = 0, i2 = 0, nspan = 0;
      double a, a1, a2, no;

      do {
        for (size_type i = 0; i < nbdir(); i++) {
          a = (2 * M_PI * double(i)) / double(nbdir());
          a1 = sin(a);
          a2 = cos(a);
          scaled_add(tx1, tgamma1, a1, tx2, tgamma2, a2, tX, tGamma); // [tX,tGamma] = a1*[tx1,tgamma1] + a2*[tx2,tgamma2]
          no = w_norm(tX, tGamma);
          scaled_add(x0, gamma0, tX, tGamma, h/no, X, Gamma);         // [X,Gamma] = [x0,gamma0] + h/no * [tX,tGamma]
          compute_tangent(X, Gamma, tX, tGamma);

          if (gmm::abs(cosang(tX, tx0, tGamma, tgamma0)) < mincos()
              || (i == 0 && nspan == 0)) {
            copy(tX, tGamma, tx0, tgamma0);
            if (insert_tangent_predict(tX, tGamma)) {
              if (noisy() > 0)
                cout << "New potential tangent vector found, "
                     << "trying one predictor-corrector step" << endl;
              copy(x0, gamma0, X, Gamma);

              if (test_predict_dir(X, Gamma, tX, tGamma)) {
                if (insert_tangent_sing(tX, tGamma)) {
                  if ((i == 0) && (nspan == 0)
                      // => (tX, tGamma) = (tx2, tgamma2)
                      && (gmm::abs(cosang(tX, tx0, tGamma, tgamma0))
                          >= mincos())) { i2 = 1; }
                  if (set_next) set_next_point(X, Gamma);

                }
                copy(x0, gamma0, X, Gamma);
                copy(tx0, tgamma0, tX, tGamma);
              }

              scale(tX, tGamma, -1.);                               // [tX,tGamma] *= -1
              if (test_predict_dir(X, Gamma, tX, tGamma)
                  && insert_tangent_sing(tX, tGamma) && set_next)
                set_next_point(X, Gamma);
            }
          }
        }

        // heuristics for varying the reference vectors
        bool perturb = true;
        if (i1 + 1 < i2) { ++i1; perturb = false; }
        else if(i2 + 1 < nb_tangent_sing())
          { ++i2; i1 = 0; perturb = false; }
        if (!perturb) {
          copy(get_tx_sing(i1), get_tgamma_sing(i1), tx1, tgamma1);
          copy(get_tx_sing(i2), get_tgamma_sing(i2), tx2, tgamma2);
        } else {
          gmm::fill_random(tX);
          tGamma = gmm::random(1.);
          no = w_norm(tX, tGamma);
          scaled_add(tx2, tgamma2, tX, tGamma, 0.1/no, tx2, tgamma2);
          // [tx2,tgamma2] += 0.1/no * [tX,tGamma]
          scaled_add(x0, gamma0, tx2, tgamma2, h, X, Gamma);        // [X,Gamma] = [x0,gamma0] + h*[tx2,tgamma2]
          compute_tangent(X, Gamma, tx2, tgamma2);
        }
      } while (++nspan < nbspan());

      if (noisy() > 0)
        cout << "Located branches " << nb_tangent_sing() << endl;
    }


    void init_Moore_Penrose_continuation(const VECT &x,
                                         double gamma, VECT &tx,
                                         double &tgamma, double &h) {
      gmm::clear(tx);
      tgamma = (tgamma >= 0) ? 1. : -1.;
      if (noisy() > 0)
        cout << "Starting computing an initial tangent" << endl;
      compute_tangent(x, gamma, tx, tgamma);
      h = h_init();
      if (this->singularities > 0)
        init_test_functions(x, gamma, tx, tgamma);
    }


    /* Perform one step of the (non-smooth) Moore-Penrose continuation.
       NOTE: The new point need not to be saved in the model in the end! */
    void Moore_Penrose_continuation(VECT &x, double &gamma,
                                    VECT &tx, double &tgamma,
                                    double &h, double &h0) {
      bool converged, new_point = false, tangent_switched = false;
      size_type it, step_dec = 0;
      double tgamma0 = tgamma, Gamma, tGamma;
      VECT tx0(tx), X(x), tX(x);

      clear_tau_bp_currentstep();
      clear_sing_data();

      do {
        h0 = h;
        // prediction
        scaled_add(x, gamma, tx, tgamma, h, X, Gamma);              // [X,Gamma] = [x,gamma] + h*[tx,tgamma]
        if (noisy() > 0)
          cout << " Prediction    : Gamma = " << Gamma
               << " (for h = " << std::scientific << std::setprecision(3) << h
               << ", tgamma = " << tgamma << ")" << endl;
        copy(tx, tgamma, tX, tGamma);

        // correction
        converged = newton_corr(X, Gamma, tX, tGamma, tx, tgamma, it);
        double cang(converged ? cosang(tX, tx, tGamma, tgamma) : 0);

        if (converged && cang >= mincos()) {
          new_point = true;
          if (this->singularities > 0) {
            if (test_limit_point(tGamma)) {
              set_sing_label("limit point");
              if (noisy() > 0) cout << "Limit point detected!" << endl;
            }
            if (this->singularities > 1) { // Treat bifurcations
              if (noisy() > 0)
                cout << "New point found, starting computing a test function "
                     << "for bifurcations" << endl;
              if (!tangent_switched) {
                if (test_smooth_bifurcation(X, Gamma, tX, tGamma)) {
                  set_sing_label("smooth bifurcation point");
                  if (noisy() > 0)
                    cout << "Smooth bifurcation point detected!" << endl;
                  treat_smooth_bif_point(X, Gamma, tX, tGamma, h);
                }
              } else if (test_nonsmooth_bifurcation(x, gamma, tx0,
                                                    tgamma0, X, Gamma, tX,
                                                    tGamma)) {
                set_sing_label("non-smooth bifurcation point");
                if (noisy() > 0)
                  cout << "Non-smooth bifurcation point detected!" << endl;
                treat_nonsmooth_point(x, gamma, tx0, tgamma0, false);
              }
            }
          }

//CHANGE 2: avoid false step increases
//if (step_dec == 0 && it < thrit() && h_inc()*(1-cang) < (1-mincos()))
          if (step_dec == 0 && it < thrit())
            h = std::min(h_inc() * h, h_max());
        } else if (h > h_min()) {
          h = std::max(h_dec() * h, h_min());
          step_dec++;
        } else if (this->non_smooth && !tangent_switched) {
          if (noisy() > 0)
            cout << "Classical continuation has failed!" << endl;
          if (switch_tangent(x, gamma, tx, tgamma, h)) {
            tangent_switched = true;
            step_dec = (h >= h_init()) ? 0 : 1;
            if (noisy() > 0)
              cout << "Restarting the classical continuation" << endl;
          } else break;
        } else new_point = true;
      } while (!new_point);

      if (new_point) {
        copy(X, Gamma, x, gamma);
        copy(tX, tGamma, tx, tgamma);
      } else if (this->non_smooth && this->singularities > 1) {
        h0 = h_min();
        treat_nonsmooth_point(x, gamma, tx0, tgamma0, true);
        if (gmm::vect_size(get_x_next()) > 0) {
          if (test_limit_point(tGamma)) {
            set_sing_label("limit point");
            if (noisy() > 0) cout << "Limit point detected!" << endl;
          }
          if (noisy() > 0)
            cout << "Starting computing a test function for bifurcations"
                 << endl;
          bool bifurcation_detected = (nb_tangent_sing() > 2);
          if (bifurcation_detected) {
            // update the stored values of the test function only
            set_tau_bp_1(tau_bp_init);
            set_tau_bp_2(test_function_bp(get_x_next(),
                                          get_gamma_next(),
                                          get_tx_sing(1),
                                          get_tgamma_sing(1)));
          } else
            bifurcation_detected
              = test_nonsmooth_bifurcation(x, gamma, tx, tgamma,
                                           get_x_next(),
                                           get_gamma_next(),
                                           get_tx_sing(1),
                                           get_tgamma_sing(1));
          if (bifurcation_detected) {
            set_sing_label("non-smooth bifurcation point");
            if (noisy() > 0)
              cout << "Non-smooth bifurcation point detected!" << endl;
          }

          copy(get_x_next(), get_gamma_next(), x, gamma);
          copy(get_tx_sing(1), get_tgamma_sing(1), tx, tgamma);
          h = h_init();
          new_point = true;
        }
      }

      if (!new_point) {
        cout << "Continuation has failed!" << endl;
        h0 = h = 0;
      }
    }

    void Moore_Penrose_continuation(VECT &x, double &gamma,
                                    VECT &tx, double &tgamma, double &h) {
      double h0;
      Moore_Penrose_continuation(x, gamma, tx, tgamma, h, h0);
    }

  protected:
    // Linear algebra functions
    void copy(const VECT &v1, const double &a1, VECT &v, double &a) const
    { gmm::copy(v1, v); a = a1; }
    void scale(VECT &v, double &a, double c) const { gmm::scale(v, c); a *= c; }
    void scaled_add(const VECT &v1, const VECT &v2, double c2, VECT &v) const
    { gmm::add(v1, gmm::scaled(v2, c2), v); }
    void scaled_add(const VECT &v1, double c1,
                    const VECT &v2, double c2, VECT &v) const
    { gmm::add(gmm::scaled(v1, c1), gmm::scaled(v2, c2), v); }
    void scaled_add(const VECT &v1, const double &a1,
                    const VECT &v2, const double &a2, double c2,
                    VECT &v, double &a) const
    { gmm::add(v1, gmm::scaled(v2, c2), v); a = a1 + c2*a2; }
    void scaled_add(const VECT &v1, const double &a1, double c1,
                    const VECT &v2, const double &a2, double c2,
                    VECT &v, double &a) const {
     gmm::add(gmm::scaled(v1, c1), gmm::scaled(v2, c2), v);
     a = c1*a1 + c2*a2;
    }
    void scaled_add(const MAT &m1, double c1,
                    const MAT &m2, double c2, MAT &m) const
    { gmm::add(gmm::scaled(m1, c1), gmm::scaled(m2, c2), m); }
    void mult(const MAT &A, const VECT &v1, VECT &v) const
    { gmm::mult(A, v1, v); }

    double norm(const VECT &v) const
    { return gmm::vect_norm2(v); }

    double sp(const VECT &v1, const VECT &v2) const
    { return gmm::vect_sp(v1, v2); }
    double sp(const VECT &v1, const VECT &v2, double w1, double w2) const
    { return sp(v1, v2) + w1 * w2; }

    virtual double intrv_sp(const VECT &v1, const VECT &v2) const = 0;

    double w_sp(const VECT &v1, const VECT &v2) const
    { return scfac * intrv_sp(v1, v2); }
    double w_sp(const VECT &v1, const VECT &v2, double w1, double w2) const
    { return w_sp(v1, v2) + w1 * w2; }
    double w_norm(const VECT &v, double w) const
    { return sqrt(w_sp(v, v) + w * w); }

    double cosang(const VECT &v1, const VECT &v2) const {
      double no = sqrt(intrv_sp(v1, v1) * intrv_sp(v2, v2));
      return (no == 0) ? 0: intrv_sp(v1, v2) / no;
    }
    double cosang(const VECT &v1, const VECT &v2, double w1, double w2) const {
//CHANGE 3: new definition of cosang
//double wgamma(0.1);
//double no = w_norm(v1, wgamma*w1) * w_norm(v2, wgamma*w2);
//return (no == 0) ? 0 : w_sp(v1, v2, wgamma*w1, wgamma*w2) / no;
      double no = sqrt((intrv_sp(v1, v1) + w1*w1)*
                       (intrv_sp(v2, v2) + w2*w2));
      return (no == 0) ? 0 : (intrv_sp(v1, v2) + w1*w2) / no;
    }

  public:

    // Misc. for accessing private data
    int noisy(void) const { return noisy_; }
    double h_init(void) const { return h_init_; }
    double h_min(void) const { return h_min_; }
    double h_max(void) const { return h_max_; }
    double h_dec(void) const { return h_dec_; }
    double h_inc(void) const { return h_inc_; }
    size_type maxit(void) const { return maxit_; }
    size_type thrit(void) const { return thrit_; }
    double maxres(void) const { return maxres_; }
    double maxdiff(void) const { return maxdiff_; }
    double mincos(void) const { return mincos_; }
    double delta_max(void) const { return delta_max_; }
    double delta_min(void) const { return delta_min_; }
    double thrvar(void) const { return thrvar_; }
    size_type nbdir(void) const { return nbdir_; }
    size_type nbspan(void) const { return nbspan_; }

    void set_tau_lp(double tau) { tau_lp = tau; }
    double get_tau_lp(void) const { return tau_lp; }
    void set_tau_bp_1(double tau) { tau_bp_1 = tau; }
    double get_tau_bp_1(void) const { return tau_bp_1; }
    void set_tau_bp_2(double tau) { tau_bp_2 = tau; }
    double get_tau_bp_2(void) const { return tau_bp_2; }
    void clear_tau_bp_currentstep(void) {
      tau_bp_graph.clear();
    }
    void init_tau_bp_graph(void) { tau_bp_graph[0.] = tau_bp_2; }
    void insert_tau_bp_graph(double alpha, double tau) {
      tau_bp_graph[alpha] = tau;
      gmm::resize(alpha_hist, 0);
      gmm::resize(tau_bp_hist, 0);
    }
    const VECT &get_alpha_hist(void) {
      size_type i = 0;
      gmm::resize(alpha_hist, tau_bp_graph.size());
      for (std::map<double, double>::iterator it = tau_bp_graph.begin();
           it != tau_bp_graph.end(); it++) {
        alpha_hist[i] = (*it).first; i++;
      }
      return alpha_hist;
    }
    const VECT &get_tau_bp_hist(void) {
      size_type i = 0;
      gmm::resize(tau_bp_hist, tau_bp_graph.size());
      for (std::map<double, double>::iterator it = tau_bp_graph.begin();
           it != tau_bp_graph.end(); it++) {
        tau_bp_hist[i] = (*it).second; i++;
      }
      return tau_bp_hist;
    }

    void clear_sing_data(void) {
      sing_label = "";
      gmm::resize(x_sing, 0);
      gmm::resize(x_next, 0);
      tx_sing.clear();
      tgamma_sing.clear();
      tx_predict.clear();
      tgamma_predict.clear();
    }
    void set_sing_label(std::string label) { sing_label = label; }
    const std::string get_sing_label(void) const { return sing_label; }
    void set_sing_point(const VECT &x, double gamma) {
      gmm::resize(x_sing, gmm::vect_size(x));
      copy(x, gamma, x_sing, gamma_sing);
    }
    const VECT &get_x_sing(void) const { return x_sing; }
    double get_gamma_sing(void) const { return gamma_sing; }
    size_type nb_tangent_sing(void) const { return tx_sing.size(); }
    bool insert_tangent_sing(const VECT &tx, double tgamma){
      bool is_included = false;
      for (size_type i = 0; (i < tx_sing.size()) && (!is_included); ++i) {
        double cang = cosang(tx_sing[i], tx, tgamma_sing[i], tgamma);
        is_included = (cang >= mincos_);
      }
      if (!is_included) {
        tx_sing.push_back(tx);
        tgamma_sing.push_back(tgamma);
      }
      return !is_included;
    }
    const VECT &get_tx_sing(size_type i) const { return tx_sing[i]; }
    double get_tgamma_sing(size_type i) const { return tgamma_sing[i]; }
    const std::vector<VECT> &get_tx_sing(void) const { return tx_sing; }
    const std::vector<double> &get_tgamma_sing(void) const { return tgamma_sing; }

    void set_next_point(const VECT &x, double gamma) {
      if (gmm::vect_size(x_next) == 0) {
        gmm::resize(x_next, gmm::vect_size(x));
        copy(x, gamma, x_next, gamma_next);
      }
    }
    const VECT &get_x_next(void) const { return x_next; }
    double get_gamma_next(void) const { return gamma_next; }

    bool insert_tangent_predict(const VECT &tx, double tgamma) {
      bool is_included = false;
      for (size_type i = 0; (i < tx_predict.size()) && (!is_included); ++i) {
        double cang = gmm::abs(cosang(tx_predict[i], tx, tgamma_predict[i], tgamma));
        is_included = (cang >= mincos_);
      }
      if (!is_included) {
        tx_predict.push_back(tx);
        tgamma_predict.push_back(tgamma);
      }
      return !is_included;
    }

    void init_border(size_type nbdof) {
      srand(unsigned(time(NULL)));
      gmm::resize(bb_x_, nbdof); gmm::fill_random(bb_x_);
      gmm::resize(cc_x_, nbdof); gmm::fill_random(cc_x_);
      bb_gamma = gmm::random(1.)/scalar_type(nbdof);
      cc_gamma = gmm::random(1.)/scalar_type(nbdof);
      dd = gmm::random(1.)/scalar_type(nbdof);
      gmm::scale(bb_x_, scalar_type(1)/scalar_type(nbdof));
      gmm::scale(cc_x_, scalar_type(1)/scalar_type(nbdof));
    }

  protected:

    const VECT &bb_x(size_type nbdof)
    { if (gmm::vect_size(bb_x_) != nbdof) init_border(nbdof); return bb_x_; }
    const VECT &cc_x(size_type nbdof)
    { if (gmm::vect_size(cc_x_) != nbdof) init_border(nbdof); return cc_x_; }

    size_type estimated_memsize(void) {
      size_type szd = sizeof(double);
      return (this->singularities == 0) ? 0
             : (2 * gmm::vect_size(bb_x_) * szd
                + 4 * gmm::vect_size(get_tau_bp_hist()) * szd
                + (1 + nb_tangent_sing()) * gmm::vect_size(get_x_sing()) * szd);
    }

    // virtual methods

    // solve A * g = L
    virtual void solve(const MAT &A, VECT &g, const VECT &L) const = 0;
    // solve A * (g1|g2) = (L1|L2)
    virtual void solve(const MAT &A, VECT &g1, VECT &g2,
                       const VECT &L1, const VECT &L2) const = 0;
    // F(x, gamma) --> f
    virtual void F(const VECT &x, double gamma, VECT &f) const = 0;
    // (F(x, gamma + eps) - f0) / eps --> g
    virtual void F_gamma(const VECT &x, double gamma, const VECT &f0,
                         VECT &g) const = 0;
    // (F(x, gamma + eps) - F(x, gamma)) / eps --> g
    virtual void F_gamma(const VECT &x, double gamma, VECT &g) const = 0;
    // F_x(x, gamma) --> A
    virtual void F_x(const VECT &x, double gamma, MAT &A) const = 0;
    // solve F_x(x, gamma) * g = L
    virtual void solve_grad(const VECT &x, double gamma,
                            VECT &g, const VECT &L) const = 0;
    // solve F_x(x, gamma) * (g1|g2) = (L1|L2)
    virtual void solve_grad(const VECT &x, double gamma, VECT &g1,
                            VECT &g2, const VECT &L1, const VECT &L2) const = 0;
    // F_x(x, gamma) * w --> y
    virtual void mult_grad(const VECT &x, double gamma,
                           const VECT &w, VECT &y) const = 0;

  public:

    virtual_cont_struct
    (int sing = 0, bool nonsm = false, double sfac=0.,
     double hin = 1.e-2, double hmax = 1.e-1, double hmin = 1.e-5,
     double hinc = 1.3, double hdec = 0.5,
     size_type mit = 10, size_type tit = 4,
     double mres = 1.e-6, double mdiff = 1.e-6, double mcos = 0.9,
     double dmax = 0.005, double dmin = 0.00012, double tvar = 0.02,
     size_type ndir = 40, size_type nspan = 1, int noi = 0)
      : singularities(sing), non_smooth(nonsm), scfac(sfac),
        h_init_(hin), h_max_(hmax), h_min_(hmin), h_inc_(hinc), h_dec_(hdec),
        maxit_(mit), thrit_(tit), maxres_(mres), maxdiff_(mdiff), mincos_(mcos),
        delta_max_(dmax), delta_min_(dmin), thrvar_(tvar),
        nbdir_(ndir), nbspan_(nspan), noisy_(noi),
        tau_lp(0.), tau_bp_1(tau_bp_init), tau_bp_2(tau_bp_init),
        gamma_sing(0.), gamma_next(0.)
    {}
    virtual ~virtual_cont_struct() {}

  };






  //=========================================================================
  // Moore-Penrose continuation method for Getfem models
  //=========================================================================


#ifdef GETFEM_MODELS_H__

  class cont_struct_getfem_model
  : public virtual_cont_struct<base_vector, model_real_sparse_matrix> {

  private:
    mutable model *md;
    std::string parameter_name;
    std::string initdata_name, finaldata_name, currentdata_name;
    gmm::sub_interval I; // for continuation based on a subset of model variables
    rmodel_plsolver_type lsolver;
    double maxres_solve;

    void set_variables(const base_vector &x, double gamma) const;
    void update_matrix(const base_vector &x, double gamma) const;

    // implemented virtual methods

    double intrv_sp(const base_vector &v1, const base_vector &v2) const {
      return (I.size() > 0) ? gmm::vect_sp(gmm::sub_vector(v1,I),
                                           gmm::sub_vector(v2,I))
                            : gmm::vect_sp(v1, v2);
    }

    // solve A * g = L
    void solve(const model_real_sparse_matrix &A, base_vector &g, const base_vector &L) const;
    // solve A * (g1|g2) = (L1|L2)
    void solve(const model_real_sparse_matrix &A, base_vector &g1, base_vector &g2,
               const base_vector &L1, const base_vector &L2) const;
    // F(x, gamma) --> f
    void F(const base_vector &x, double gamma, base_vector &f) const;
    // (F(x, gamma + eps) - f0) / eps --> g
    void F_gamma(const base_vector &x, double gamma, const base_vector &f0,
                 base_vector &g) const;
    // (F(x, gamma + eps) - F(x, gamma)) / eps --> g
    void F_gamma(const base_vector &x, double gamma, base_vector &g) const;

    // F_x(x, gamma) --> A
    void F_x(const base_vector &x, double gamma, model_real_sparse_matrix &A) const;
    // solve F_x(x, gamma) * g = L
    void solve_grad(const base_vector &x, double gamma,
                    base_vector &g, const base_vector &L) const;
    // solve F_x(x, gamma) * (g1|g2) = (L1|L2)
    void solve_grad(const base_vector &x, double gamma, base_vector &g1,
                    base_vector &g2,
                    const base_vector &L1, const base_vector &L2) const;
    // F_x(x, gamma) * w --> y
    void mult_grad(const base_vector &x, double gamma,
                   const base_vector &w, base_vector &y) const;

  public:
    size_type estimated_memsize(void);
    const model &linked_model(void) { return *md; }

    void set_parametrised_data_names
    (const std::string &in, const std::string &fn, const std::string &cn) {
      initdata_name = in;
      finaldata_name = fn;
      currentdata_name = cn;
    }

    void set_interval_from_variable_name(const std::string &varname) {
      if (varname == "") I = gmm::sub_interval(0,0);
      else I = md->interval_of_variable(varname);
    }

    cont_struct_getfem_model
    (model &md_, const std::string &pn, double sfac, rmodel_plsolver_type ls,
     double hin = 1.e-2, double hmax = 1.e-1, double hmin = 1.e-5,
     double hinc = 1.3, double hdec = 0.5, size_type mit = 10,
     size_type tit = 4, double mres = 1.e-6, double mdiff = 1.e-6,
     double mcos = 0.9, double mress = 1.e-8, int noi = 0, int sing = 0,
     bool nonsm = false, double dmax = 0.005, double dmin = 0.00012,
     double tvar = 0.02, size_type ndir = 40, size_type nspan = 1)
      : virtual_cont_struct(sing, nonsm, sfac, hin, hmax, hmin, hinc, hdec,
                            mit, tit, mres, mdiff, mcos, dmax, dmin, tvar,
                            ndir, nspan, noi),
        md(&md_), parameter_name(pn),
        initdata_name(""), finaldata_name(""), currentdata_name(""),
        I(0,0), lsolver(ls), maxres_solve(mress)
    {
      GMM_ASSERT1(!md->is_complex(),
                  "Continuation has only a real version, sorry.");
    }

  };

#endif


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONTINUATION_H__ */
