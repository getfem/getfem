// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2011-2011 Tomas Ligursky, Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// As a special exception, you  may use  this file  as it is a part of a free
// software  library  without  restriction.  Specifically,  if   other  files
// instantiate  templates  or  use macros or inline functions from this file,
// or  you compile this  file  and  link  it  with other files  to produce an
// executable, this file  does  not  by itself cause the resulting executable
// to be covered  by the GNU Lesser General Public License.  This   exception
// does not  however  invalidate  any  other  reasons why the executable file
// might be covered by the GNU Lesser General Public License.
//
//===========================================================================

/** @file getfem_continuation.h
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @author Tomas Ligursky <tomas.ligursky@gmail.com>
    @date October 17, 2011.
    @brief Moore-Penrose  (or Gauss-Newton ?) continuation method.
 */
#ifndef GETFEM_CONTINUATION_H__
#define GETFEM_CONTINUATION_H__


namespace getfem {


  //=========================================================================
  // Abstract Moore-Penrose continuation method
  //=========================================================================

  
  template <typename S, typename VECT> 
  void compute_gamma_derivative_(const S &s, const VECT &y0,
				 double gamma0, VECT &v) {
    if (s.has_gamma_derivative()) {
      s.gamma_derivative(y0, gamma0, v);
    } else {
      VECT F0(y0), F1(y0);
      s.F(y0, gamma0, F0);
      s.F(y0, gamma0 + s.epsilon(), F1);
      s.sub(F1, F0, v);
      s.scale(v, 1./s.epsilon());
    }
  }


  template <typename S, typename VECT>
    void init_continuation(const S &s, const VECT &y0,
			   double gamma0, double step,
			   VECT &t_y, double &t_gamma) {
    
    VECT v(y0);

    // Compute the derivative with respect to the continuation parameter
    compute_gamma_derivative_(s, y0, gamma0, v);
    
    // computation of the tangent
    s.solve_grad(y0, gamma0, v, t_y);
    double no = s.norm(t_y); no = sqrt(no * no + 1.0);
    if (step > 0.) no = -no; 
    s.scale(t_y, 1/no);
    t_gamma = 1./no;
  }



  template <typename S, typename VECT>
    void Moore_Penrose_continuation(const S &s, VECT &y0,
				    double &gamma0, double step,
				    VECT &t_y, double &t_gamma) {
    
    VECT Y0(y0), Y(y0), T0_y(y0), T_y(y0), w(y0), z(y0), F(y0), v(y0), X_y(y0);
    double Gamma0, T0_gamma, Gamma;
    step = abs(step);
    

    for (;;) { // step control

      // prediction
      s.scaled_add(y0, 1., t_y, step, Y0);
      Gamma0 = gamma0 + step*t_gamma;
      s.copy(t_y, T0_y);
      T0_gamma = t_gamma;
      
      // correction
      
      for (;;) { // Newton iterations
	
	
	compute_gamma_derivative_(s, Y0, Gamma0, v);
	
	
	s.solve_grad(Y0, Gamma0, v, T_y);
	s.copy(T_y, w);
	
	double T_gamma = 1. / (t_gamma - s.sp(t_y, T_y));
	s.scale(T_y, -T_gamma);
	
	double no = s.norm(T_y); no = sqrt(no * no + T_gamma*T_gamma);
	s.scale(T_y, 1./no); T_gamma /= no;
	
	
	s.F(y0, gamma0, F);
	s.solve_grad(Y0, Gamma0, F, z);
	
	double X_gamma = -s.sp(t_y, z) / (t_gamma - s.sp(t_y, w));
	s.scaled_add(z, 1., w, -X_gamma, X_y);
	
	s.add(Y0, X_y, Y);
	Gamma = Gamma0 + X_gamma;
	
	....;
	
	s.copy(Y, Y0); Gamma0 = Gamma;
	
      }
      

    }

    s.copy(Y, y0); gamma0 = Gamma; s.copy(T_y, t_y); t_gamma = T_gamma;

  }



  //=========================================================================
  // Moore-Penrose continuation method for Getfem models
  //=========================================================================


#ifdef GETFEM_MODELS_H__
 
  struct S_getfem_model {

    double epsilon_;
    model &md;  // for real models only
    std::string parameter_name;
    rmodel_plsolver_type lsolver;

    typedef std::vector<double> VECT;
    

    // Linear algebra functions
    void sub(const VECT &v1, const VECT &v2, VECT &v)
    { gmm::add(v1, gmm::scaled(v2, -1.), v); }
    void add(const VECT &v1, const VECT &v2, VECT &v)
    { gmm::add(v1, v2, v); }
    void scale(VECT &v, double alpha)
    { gmm::scale(v, alpha); }
    void copy(const VECT &v1, VECT &v)
    { gmm::copy(v1, v); }
    void scaled_add(VECT &v1, double a, const VECT &v2, double b, VECT &v)
    { gmm::add(gmm::scaled(v1, a), gmm::scaled(v2, b), v); }


    // Evaluation of  ...
    void F(const VECT &y, double gamma, VECT &f) {
      md.set_real_variable(parameter_name)[0] = gamma;
      md.to_variables(y);
      md.assembly(model::BUILD_RHS);
      gmm::copy(gmm::scaled(md.real_rhs(), -1.), f);
    }

    void solve_grad(const VECT &y, double gamma, const VECT &L, VECT &g) {
      md.set_real_variable(parameter_name)[0] = gamma;
      md.to_variables(y);
      md.assembly(model::BUILD_MATRIX);

      gmm::iteration iter(...);

      lsolver(md.real_tangent_matrix(), g, L, iter);
    }

    
    // Misc.

    double epsilon(void) { return epsilon_; }
    bool has_gamma_derivative(void) { return false; }




  };

#endif

















}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONTINUATION_H__ */
