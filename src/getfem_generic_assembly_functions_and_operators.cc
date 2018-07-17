/*===========================================================================

 Copyright (C) 2013-2018 Yves Renard

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

===========================================================================*/

#include "getfem/getfem_generic_assembly_semantic.h"
#include "getfem/getfem_generic_assembly_functions_and_operators.h"
#include "getfem/getfem_generic_assembly_compile_and_exec.h"

/**
   Providing for special Math functions unavailable on Intel or MSVS C++
   compilers
*/

#if defined(_MSC_VER) && _MSC_VER < 1800
#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/special_functions/erf.hpp>
typedef double (*BoostMathFunction)(double);
BoostMathFunction const acosh = boost::math::acosh<double>;
BoostMathFunction const asinh = boost::math::asinh<double>;
BoostMathFunction const atanh = boost::math::atanh<double>;
BoostMathFunction const erf = boost::math::erf<double>;
BoostMathFunction const erfc = boost::math::erfc<double>;
#endif

namespace getfem {

  base_matrix& __mat_aux1()
  {
    DEFINE_STATIC_THREAD_LOCAL(base_matrix, m);
    return m;
  }



  //=========================================================================
  // Structure dealing with predefined scalar functions.
  //=========================================================================

  scalar_type ga_predef_function::operator()(scalar_type t_,
					     scalar_type u_) const {
    switch(ftype_) {
    case 0:
      if (nbargs_ == 2)
	return (*f2_)(t_, u_);
      else
	return (*f1_)(t_);
      break;
    case 1:
      t.thrd_cast()[0] = t_; u.thrd_cast()[0] = u_;
      workspace.thrd_cast().assembled_potential() = scalar_type(0);
      ga_function_exec(*gis);
      return workspace.thrd_cast().assembled_potential();
      break;
    }
    return 0.;
  }

  bool ga_predef_function::is_affine(const std::string &varname) const {
    if (ftype_ == 1) {
      for (size_type i = 0; i < workspace.thrd_cast().nb_trees(); ++i) {
	const ga_workspace::tree_description &
	  td = workspace.thrd_cast().tree_info(i);
	if (!(ga_is_affine(*(td.ptree), workspace, varname, "")))
	  return false;
      }
      return true;
    }
    return false;
  }
  
  
  static scalar_type ga_Heaviside(scalar_type t) { return (t >= 0.) ? 1.: 0.; }
  static scalar_type ga_pos_part(scalar_type t) { return (t >= 0.) ? t : 0.; }
  static scalar_type ga_reg_pos_part(scalar_type t, scalar_type eps)
  { return (t >= eps) ? t-eps/2. : ((t <= 0) ? 0. : t*t/(2.*eps)); }
  static scalar_type ga_der_reg_pos_part(scalar_type t, scalar_type eps)
  { return (t >= eps) ? 1. : ((t <= 0) ? 0. : t/eps); }
  static scalar_type ga_der2_reg_pos_part(scalar_type t, scalar_type eps)
  { return (t >= eps) ? 0. : ((t <= 0) ? 0. : 1./eps); }


  static scalar_type ga_half_sqr_pos_part(scalar_type t)
  { return (t >= 0.) ? 0.5*t*t : 0.; }
  static scalar_type ga_neg_part(scalar_type t) { return (t >= 0.) ? 0. : -t; }
  static scalar_type ga_half_sqr_neg_part(scalar_type t)
  { return (t >= 0.) ? 0. : 0.5*t*t; }
  static scalar_type ga_sinc(scalar_type t) {// cardinal sine function sin(t)/t
    if (gmm::abs(t) < 1E-4) {
      scalar_type t2 = t*t;
      return 1-t2/6.+ t2*t2/120.;
    } else {
      return sin(t)/t;
    }
  }
  static scalar_type ga_sqr(scalar_type t) { return t*t; }
  static scalar_type ga_max(scalar_type t, scalar_type u)
  { return std::max(t,u); }
  static scalar_type ga_min(scalar_type t, scalar_type u)
  { return std::min(t,u); }
  static scalar_type ga_abs(scalar_type t) { return gmm::abs(t); }
  static scalar_type ga_sign(scalar_type t) { return (t >= 0.) ? 1.: -1.; }

  // Derivatives of predefined functions
  static scalar_type ga_der_sinc(scalar_type t) {
    if (gmm::abs(t) < 1E-4) {
      scalar_type t2 = t*t;
      return  -t/3. + t*t2/30. -t*t2*t2/840.;
    } else {
      return (t*cos(t) - sin(t))/(t*t);
    }
  }
  static scalar_type ga_der2_sinc(scalar_type t) {
    if (gmm::abs(t) < 1E-4) {
      scalar_type t2 = t*t;
      return  -1./3. + t2/10. -t2*t2/168.;
    } else {
      return ((2. - t*t)*sin(t) - 2.*t*cos(t))/(t*t*t);
    }
  }
  static scalar_type ga_der_sqrt(scalar_type t) { return 0.5/sqrt(t); }
  // static scalar_type ga_der_sqr(scalar_type t) { return 2*t; }
  static scalar_type ga_der_t_pow(scalar_type t, scalar_type u)
  { return u*pow(t,u-1.); }
  static scalar_type ga_der_u_pow(scalar_type t, scalar_type u)
  { return pow(t,u)*log(gmm::abs(t)); }
  static scalar_type ga_der_log(scalar_type t) { return 1./t; }
  static scalar_type ga_der_log10(scalar_type t) { return 1./(t*log(10.)); }
  static scalar_type ga_der_tanh(scalar_type t)
  { return 1.-gmm::sqr(tanh(t)); }
  static scalar_type ga_der_asinh(scalar_type t)
  { return 1./(sqrt(t*t+1.)); }
  static scalar_type ga_der_acosh(scalar_type t)
  { return 1./(sqrt(t*t-1.)); }
  static scalar_type ga_der_atanh(scalar_type t)
  { return 1./(1.-t*t); }
  static scalar_type ga_der_cos(scalar_type t)
  { return -sin(t); }
  static scalar_type ga_der_tan(scalar_type t)
  { return 1.+gmm::sqr(tan(t)); }
  static scalar_type ga_der_asin(scalar_type t)
  { return 1./(sqrt(1.-t*t)); }
  static scalar_type ga_der_acos(scalar_type t)
  { return -1./(sqrt(1.-t*t)); }
  static scalar_type ga_der_atan(scalar_type t)
  { return 1./(1.+t*t); }
  static scalar_type ga_der_t_atan2(scalar_type t, scalar_type u)
  { return u/(t*t+u*u); }
  static scalar_type ga_der_u_atan2(scalar_type t, scalar_type u)
  { return -t/(t*t+u*u); }
  static scalar_type ga_der_erf(scalar_type t)
  { return exp(-t*t)*2./sqrt(M_PI); }
  static scalar_type ga_der_erfc(scalar_type t)
  { return -exp(-t*t)*2./sqrt(M_PI); }
  static scalar_type ga_der_neg_part(scalar_type t)
  { return (t >= 0) ? 0. : -1.; }
  static scalar_type ga_der_t_max(scalar_type t, scalar_type u)
  { return (t-u >= 0) ? 1. : 0.; }
  static scalar_type ga_der_u_max(scalar_type t, scalar_type u)
  { return (u-t >= 0) ? 1. : 0.; }



  //=========================================================================
  // Structure dealing with predefined operators.
  //=========================================================================

  static void ga_init_scalar(bgeot::multi_index &mi) { mi.resize(0); }
  static void ga_init_square_matrix(bgeot::multi_index &mi, size_type N)
  { mi.resize(2); mi[0] = mi[1] = N; }

  // Norm Operator
  struct norm_operator : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 1 || args[0]->sizes().size() > 2) return false;
      ga_init_scalar(sizes);
      return true;
    }

    void value(const arg_list &args, base_tensor &result) const
    { result[0] = gmm::vect_norm2(args[0]->as_vector()); }

    // Derivative : u/|u|
    void derivative(const arg_list &args, size_type,
                    base_tensor &result) const {
      scalar_type no = gmm::vect_norm2(args[0]->as_vector());
      if (no == scalar_type(0))
        gmm::clear(result.as_vector());
      else
        gmm::copy(gmm::scaled(args[0]->as_vector(), scalar_type(1)/no),
                  result.as_vector());
    }

    // Second derivative : (|u|^2 Id - u x u)/|u|^3
    void second_derivative(const arg_list &args, size_type, size_type,
                           base_tensor &result) const {
      const base_tensor &t = *args[0];
      size_type N = t.size();
      scalar_type no = gmm::vect_norm2(t.as_vector());
      scalar_type no3 = no*no*no;

      if (no < 1E-25) no = 1E-25; // In order to avoid infinite values

      for (size_type i = 0; i < N; ++i)
        for (size_type j = 0; j < N; ++j) {
          result[j*N+i] = - t[i]*t[j] / no3;
          if (i == j) result[j*N+i] += scalar_type(1)/no;
        }
    }
  };

  // Norm_sqr Operator
  struct norm_sqr_operator : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 1 || args[0]->sizes().size() > 2) return false;
      ga_init_scalar(sizes);
      return true;
    }

    void value(const arg_list &args, base_tensor &result) const
    { result[0] = gmm::vect_norm2_sqr(args[0]->as_vector()); }

    // Derivative : 2*u
    void derivative(const arg_list &args, size_type,
                    base_tensor &result) const {
      gmm::copy(gmm::scaled(args[0]->as_vector(), scalar_type(2)),
                result.as_vector());
    }

    // Second derivative : Id
    void second_derivative(const arg_list &args, size_type, size_type,
                           base_tensor &result) const {
      const base_tensor &t = *args[0];
      size_type N = t.size();
      gmm::clear(result.as_vector());
      for (size_type i = 0; i < N; ++i)
        result[i*N+i] = scalar_type(2);
    }
  };

  // Det Operator
  struct det_operator : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 1 || args[0]->sizes().size() != 2
          || args[0]->sizes()[0] != args[0]->sizes()[1]) return false;
      ga_init_scalar(sizes);
      return true;
    }

    void value(const arg_list &args, base_tensor &result) const {
      size_type N = args[0]->sizes()[0];
      // base_matrix M(N, N);
      // gmm::copy(args[0]->as_vector(), M.as_vector());
      // result[0] = gmm::lu_det(M);
      result[0] = bgeot::lu_det(&(*(args[0]->begin())), N);
    }

    // Derivative : det(M)M^{-T}
    void derivative(const arg_list &args, size_type,
                    base_tensor &result) const {
      size_type N = args[0]->sizes()[0];
      if (N) {
        __mat_aux1().base_resize(N, N);
        gmm::copy(args[0]->as_vector(), __mat_aux1().as_vector());
        scalar_type det = bgeot::lu_inverse(__mat_aux1());
        if (det == scalar_type(0))
          gmm::clear(result.as_vector());
        else {
          auto it = result.begin();
          auto ita = __mat_aux1().begin();
          for (size_type j = 0; j < N; ++j, ++ita) {
            auto itaa = ita;
            *it = (*itaa) * det; ++it;
            for (size_type i = 1; i < N; ++i, ++it)
              { itaa += N; *it = (*itaa) * det; }
          }
          GA_DEBUG_ASSERT(it == result.end(), "Internal error");
        }
      }
    }

    // Second derivative : det(M)(M^{-T}@M^{-T} - M^{-T}_{jk}M^{-T}_{li})
    void second_derivative(const arg_list &args, size_type, size_type,
                           base_tensor &result) const {
      size_type N = args[0]->sizes()[0];
      __mat_aux1().base_resize(N, N);
      gmm::copy(args[0]->as_vector(), __mat_aux1().as_vector());
      scalar_type det = bgeot::lu_inverse(__mat_aux1());
      if (det == scalar_type(0))
        gmm::clear(result.as_vector());
      else {
        auto it = result.begin();
        auto ita = __mat_aux1().begin(), ita_l = ita;
        for (size_type l = 0; l < N; ++l, ++ita_l) {
          auto ita_lk = ita_l, ita_jk = ita;
          for (size_type k = 0; k < N; ++k, ita_lk += N, ita_jk += N) {
            auto ita_j = ita;
            for (size_type j = 0; j < N; ++j, ++ita_j, ++ita_jk) {
              auto ita_ji = ita_j, ita_li = ita_l;
              for (size_type i = 0; i < N; ++i, ++it, ita_ji += N, ita_li += N)
                *it = ((*ita_ji) * (*ita_lk) - (*ita_jk) * (*ita_li)) * det;
            }
          }
        }
        GA_DEBUG_ASSERT(it == result.end(), "Internal error");
      }
    }
  };

  // Inverse Operator (for square matrices)
  struct inverse_operator : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 1 || args[0]->sizes().size() != 2
          || args[0]->sizes()[0] != args[0]->sizes()[1]) return false;
      ga_init_square_matrix(sizes, args[0]->sizes()[0]);
      return true;
    }

    // Value : M^{-1}
    void value(const arg_list &args, base_tensor &result) const {
      size_type N = args[0]->sizes()[0];
      __mat_aux1().base_resize(N, N);
      gmm::copy(args[0]->as_vector(), __mat_aux1().as_vector());
      bgeot::lu_inverse(__mat_aux1());
      gmm::copy(__mat_aux1().as_vector(), result.as_vector());
    }

    // Derivative : -M^{-1}{ik}M^{-1}{lj}  (comes from H -> M^{-1}HM^{-1})
    void derivative(const arg_list &args, size_type,
                    base_tensor &result) const { // to be verified
      size_type N = args[0]->sizes()[0];
      __mat_aux1().base_resize(N, N);
      gmm::copy(args[0]->as_vector(), __mat_aux1().as_vector());
      bgeot::lu_inverse(__mat_aux1());
      auto it = result.begin();
      auto ita = __mat_aux1().begin(), ita_l = ita;
      for (size_type l = 0; l < N; ++l, ++ita_l) {
        auto ita_k = ita;
        for (size_type k = 0; k < N; ++k, ita_k += N) {
          auto ita_lj = ita_l;
          for (size_type j = 0; j < N; ++j, ita_lj += N) {
            auto ita_ik = ita_k;
            for (size_type i = 0; i < N; ++i, ++it, ++ita_ik)
              *it = -(*ita_ik) * (*ita_lj);
          }
        }
      }
      GA_DEBUG_ASSERT(it == result.end(), "Internal error");
    }

    // Second derivative :
    // M^{-1}{ik}M^{-1}{lm}M^{-1}{nj} + M^{-1}{im}M^{-1}{mk}M^{-1}{lj}
    // comes from (H,K) -> M^{-1}HM^{-1}KM^{-1} + M^{-1}KM^{-1}HM^{-1}
    void second_derivative(const arg_list &args, size_type, size_type,
                           base_tensor &result) const { // to be verified
      size_type N = args[0]->sizes()[0];
      __mat_aux1().base_resize(N, N);
      gmm::copy(args[0]->as_vector(), __mat_aux1().as_vector());
      bgeot::lu_inverse(__mat_aux1());
      base_tensor::iterator it = result.begin();
      for (size_type n = 0; n < N; ++n)
        for (size_type m = 0; m < N; ++m)
          for (size_type l = 0; l < N; ++l)
            for (size_type k = 0; k < N; ++k)
              for (size_type j = 0; j < N; ++j)
                for (size_type i = 0; i < N; ++i, ++it)
                  *it = __mat_aux1()(i,k)*__mat_aux1()(l,m)*__mat_aux1()(n,j)
                    + __mat_aux1()(i,m)*__mat_aux1()(m,k)*__mat_aux1()(l,j);
      GA_DEBUG_ASSERT(it == result.end(), "Internal error");
    }
  };

  //=========================================================================
  // Initialization of predefined functions and operators.
  //=========================================================================

  ga_predef_function::ga_predef_function()
    : expr_(""), derivative1_(""), derivative2_(""), gis(nullptr) {}

  ga_predef_function::ga_predef_function(pscalar_func_onearg f,
					 size_type dtype__,
					 const std::string &der)
    : ftype_(0), dtype_(dtype__), nbargs_(1), f1_(f), expr_(""),
        derivative1_(der), derivative2_("") {}
  ga_predef_function::ga_predef_function(pscalar_func_twoargs f,
					 size_type dtype__,
					 const std::string &der1,
					 const std::string &der2)
    : ftype_(0), dtype_(dtype__), nbargs_(2), f2_(f),
      expr_(""), derivative1_(der1), derivative2_(der2), gis(nullptr) {}
  ga_predef_function::ga_predef_function(const std::string &expr__)
    : ftype_(1), dtype_(3), nbargs_(1), expr_(expr__),
      derivative1_(""), derivative2_(""), t(1, 0.), u(1, 0.), gis(nullptr) {}

  
  ga_predef_function_tab::ga_predef_function_tab() {
    
    ga_predef_function_tab &PREDEF_FUNCTIONS = *this;

    // Power functions and their derivatives
    PREDEF_FUNCTIONS["sqrt"] = ga_predef_function(sqrt, 1, "DER_PDFUNC_SQRT");
    PREDEF_FUNCTIONS["sqr"] = ga_predef_function(ga_sqr, 2, "2*t");
    PREDEF_FUNCTIONS["pow"] = ga_predef_function(pow, 1, "DER_PDFUNC1_POW",
                                                 "DER_PDFUNC2_POW");
    
    PREDEF_FUNCTIONS["DER_PDFUNC_SQRT"] =
      ga_predef_function(ga_der_sqrt, 2, "-0.25/(t*sqrt(t))");
    PREDEF_FUNCTIONS["DER_PDFUNC1_POW"] =
      ga_predef_function(ga_der_t_pow, 2, "u*(u-1)*pow(t,u-2)",
                         "pow(t,u-1)*(u*log(t)+1)");
    PREDEF_FUNCTIONS["DER_PDFUNC2_POW"] =
      ga_predef_function(ga_der_u_pow, 2, "pow(t,u-1)*(u*log(t)+1)",
                         "pow(t,u)*sqr(log(t))");

    // Hyperbolic functions
    PREDEF_FUNCTIONS["exp"] = ga_predef_function(exp, 1, "exp");
    PREDEF_FUNCTIONS["log"] = ga_predef_function(log, 1, "DER_PDFUNC_LOG");
    PREDEF_FUNCTIONS["log10"] =
      ga_predef_function(log10, 1, "DER_PDFUNC_LOG10");
    PREDEF_FUNCTIONS["sinh"] = ga_predef_function(sinh, 1, "cosh");
    PREDEF_FUNCTIONS["cosh"] = ga_predef_function(cosh, 1, "sinh");
    PREDEF_FUNCTIONS["tanh"] = ga_predef_function(tanh, 1, "DER_PDFUNC_TANH");
    PREDEF_FUNCTIONS["asinh"] =
      ga_predef_function(asinh, 1, "DER_PDFUNC_ASINH");
    PREDEF_FUNCTIONS["acosh"] =
      ga_predef_function(acosh, 1, "DER_PDFUNC_ACOSH");
    PREDEF_FUNCTIONS["atanh"] =
      ga_predef_function(atanh, 1, "DER_PDFUNC_ATANH");


    PREDEF_FUNCTIONS["DER_PDFUNC_LOG"] =
      ga_predef_function(ga_der_log, 2, "-1/sqr(t)");
    PREDEF_FUNCTIONS["DER_PDFUNC_LOG10"] =
      ga_predef_function(ga_der_log10, 2, "-1/(sqr(t)*log(10))");
    PREDEF_FUNCTIONS["DER_PDFUNC_TANH"] =
      ga_predef_function(ga_der_tanh, 2, "2*tanh(t)*(sqr(tanh(t))-1)");
    PREDEF_FUNCTIONS["DER_PDFUNC_ASINH"] =
      ga_predef_function(ga_der_asinh, 2, "-t/(pow(t*t+1,1.5))");
    PREDEF_FUNCTIONS["DER_PDFUNC_ACOSH"] =
      ga_predef_function(ga_der_acosh, 2, "-t/(pow(t*t-1,1.5))");
    PREDEF_FUNCTIONS["DER_PDFUNC_ATANH"] =
      ga_predef_function(ga_der_atanh, 2, "2*t/sqr(1-t*t)");


    // Trigonometric functions
    PREDEF_FUNCTIONS["sin"] = ga_predef_function(sin, 1, "cos");
    PREDEF_FUNCTIONS["cos"] = ga_predef_function(cos, 1, "DER_PDFUNC_COS");
    PREDEF_FUNCTIONS["tan"] = ga_predef_function(tan, 1, "DER_PDFUNC_TAN");
    PREDEF_FUNCTIONS["asin"] = ga_predef_function(asin, 1, "DER_PDFUNC_ASIN");
    PREDEF_FUNCTIONS["acos"] = ga_predef_function(acos, 1, "DER_PDFUNC_ACOS");
    PREDEF_FUNCTIONS["atan"] = ga_predef_function(atan, 1, "DER_PDFUNC_ATAN");
    PREDEF_FUNCTIONS["atan2"] = ga_predef_function(atan2,1,"DER_PDFUNC1_ATAN2",
                                                           "DER_PDFUNC2_ATAN2");
    PREDEF_FUNCTIONS["sinc"] = ga_predef_function(ga_sinc, 1,
                                                  "DER_PDFUNC_SINC");
    PREDEF_FUNCTIONS["DER_PDFUNC_SINC"] =ga_predef_function(ga_der_sinc, 1,
                                                            "DER2_PDFUNC_SINC");
    PREDEF_FUNCTIONS["DER2_PDFUNC_SINC"] = ga_predef_function(ga_der2_sinc);


    PREDEF_FUNCTIONS["DER_PDFUNC_COS"] =
      ga_predef_function(ga_der_cos, 2, "-cos(t)");
    PREDEF_FUNCTIONS["DER_PDFUNC_TAN"] =
      ga_predef_function(ga_der_tan, 2, "2*tan(t)/sqr(cos(t))");
    // PREDEF_FUNCTIONS["DER_PDFUNC_TAN"] =
    //  ga_predef_function(ga_der_tan, 2, "2*tan(t)*(1+sqr(tan(t)))");
    PREDEF_FUNCTIONS["DER_PDFUNC_ASIN"] =
      ga_predef_function(ga_der_asin, 2, "t/(pow(1-t*t,1.5))");
    PREDEF_FUNCTIONS["DER_PDFUNC_ACOS"] =
      ga_predef_function(ga_der_acos, 2, "-t/(pow(1-t*t,1.5))");
    PREDEF_FUNCTIONS["DER_PDFUNC_ATAN"] =
      ga_predef_function(ga_der_atan, 2, "-2*t/sqr(1+t*t)");
    PREDEF_FUNCTIONS["DER_PDFUNC1_ATAN2"] =
      ga_predef_function(ga_der_t_atan2, 2, "-2*t*u/sqr(sqr(u)+sqr(t))",
                         "(sqrt(t)-sqr(u))/sqr(sqr(u)+sqr(t))");
    PREDEF_FUNCTIONS["DER_PDFUNC2_ATAN2"] =
      ga_predef_function(ga_der_u_atan2, 2,
                         "(sqrt(t)-sqr(u))/sqr(sqr(u)+sqr(t))",
                         "2*t*u/sqr(sqr(u)+sqr(t))");


    // Error functions
    PREDEF_FUNCTIONS["erf"]
      = ga_predef_function(erf, 1, "DER_PDFUNC_ERF");
    PREDEF_FUNCTIONS["erfc"]
      = ga_predef_function(erfc, 1, "DER_PDFUNC_ERFC");

    PREDEF_FUNCTIONS["DER_PDFUNC_ERF"] =
      ga_predef_function(ga_der_erf, 2, "exp(-t*t)*2/sqrt(pi)");
    PREDEF_FUNCTIONS["DER_PDFUNC_ERFC"] =
      ga_predef_function(ga_der_erfc, 2, "-exp(-t*t)*2/sqrt(pi)");



    // Miscellaneous functions
    PREDEF_FUNCTIONS["Heaviside"] = ga_predef_function(ga_Heaviside); // ga_predef_function(ga_Heaviside, 2, "(0)");
    PREDEF_FUNCTIONS["sign"] = ga_predef_function(ga_sign);
    PREDEF_FUNCTIONS["abs"] = ga_predef_function(ga_abs, 1, "sign");
    PREDEF_FUNCTIONS["pos_part"]
      = ga_predef_function(ga_pos_part, 1, "Heaviside");
    PREDEF_FUNCTIONS["half_sqr_pos_part"]
      = ga_predef_function(ga_half_sqr_pos_part, 1, "pos_part");
    PREDEF_FUNCTIONS["neg_part"]
      = ga_predef_function(ga_neg_part, 1, "DER_PDFUNC_NEG_PART");
    PREDEF_FUNCTIONS["half_sqr_neg_part"]
      = ga_predef_function(ga_half_sqr_neg_part, 2, "-neg_part(t)");
    PREDEF_FUNCTIONS["reg_pos_part"]
      = ga_predef_function(ga_reg_pos_part, 1, "DER_REG_POS_PART", "");
    PREDEF_FUNCTIONS["DER_REG_POS_PART"]
      = ga_predef_function(ga_der_reg_pos_part, 1, "DER2_REG_POS_PART", "");
    PREDEF_FUNCTIONS["DER_REG_POS_PART"]
      = ga_predef_function(ga_der2_reg_pos_part);
    
    PREDEF_FUNCTIONS["max"]
      = ga_predef_function(ga_max, 1, "DER_PDFUNC1_MAX", "DER_PDFUNC2_MAX");
    PREDEF_FUNCTIONS["min"]
      = ga_predef_function(ga_min, 1, "DER_PDFUNC2_MAX", "DER_PDFUNC1_MAX");

    PREDEF_FUNCTIONS["DER_PDFUNC_NEG_PART"]
      = ga_predef_function(ga_der_neg_part, 2, "-Heaviside(-t)");
    PREDEF_FUNCTIONS["DER_PDFUNC1_MAX"] = ga_predef_function(ga_der_t_max);
    PREDEF_FUNCTIONS["DER_PDFUNC2_MAX"] = ga_predef_function(ga_der_u_max);

  }

  ga_spec_function_tab::ga_spec_function_tab() {
    // Predefined special functions
    ga_spec_function_tab &SPEC_FUNCTIONS = *this;
    
    SPEC_FUNCTIONS.insert("pi");
    SPEC_FUNCTIONS.insert("meshdim");
    SPEC_FUNCTIONS.insert("timestep");
    SPEC_FUNCTIONS.insert("qdim");
    SPEC_FUNCTIONS.insert("qdims");
    SPEC_FUNCTIONS.insert("Id");
  }

  ga_spec_op_tab::ga_spec_op_tab() {
    // Predefined special operators
    ga_spec_op_tab &SPEC_OP = *this;
    SPEC_OP.insert("X");
    SPEC_OP.insert("element_size");
    SPEC_OP.insert("element_K");
    SPEC_OP.insert("element_B");
    SPEC_OP.insert("Normal");
    SPEC_OP.insert("Sym");
    SPEC_OP.insert("Skew");
    SPEC_OP.insert("Def");
    SPEC_OP.insert("Trace");
    SPEC_OP.insert("Deviator");
    SPEC_OP.insert("Interpolate");
    SPEC_OP.insert("Interpolate_filter");
    SPEC_OP.insert("Elementary_transformation");
    SPEC_OP.insert("Xfem_plus");
    SPEC_OP.insert("Xfem_minus");
    SPEC_OP.insert("Print");
    SPEC_OP.insert("Reshape");
    SPEC_OP.insert("Swap_indices");
    SPEC_OP.insert("Index_move_last");
    SPEC_OP.insert("Contract");
    SPEC_OP.insert("Diff");
    SPEC_OP.insert("Grad");
  }


  ga_predef_operator_tab::ga_predef_operator_tab(void) {
    // Predefined operators
    ga_predef_operator_tab &PREDEF_OPERATORS = *this;

    PREDEF_OPERATORS.add_method("Norm", std::make_shared<norm_operator>());
    PREDEF_OPERATORS.add_method("Norm_sqr",
                                std::make_shared<norm_sqr_operator>());
    PREDEF_OPERATORS.add_method("Det", std::make_shared<det_operator>());
    PREDEF_OPERATORS.add_method("Inv", std::make_shared<inverse_operator>());
  }



  bool ga_function_exists(const std::string &name) {
    const ga_predef_function_tab &PREDEF_FUNCTIONS
      = dal::singleton<ga_predef_function_tab>::instance(0);
    return PREDEF_FUNCTIONS.find(name) != PREDEF_FUNCTIONS.end();
  }


  void ga_define_function(const std::string &name, size_type nbargs,
                          const std::string &expr, const std::string &der1,
                          const std::string &der2) {
    auto guard = omp_guard{};

    auto &PREDEF_FUNCTIONS = dal::singleton<ga_predef_function_tab>::instance(0);
    if(PREDEF_FUNCTIONS.find(name) != PREDEF_FUNCTIONS.end()) return;

    GMM_ASSERT1(nbargs >= 1 && nbargs <= 2, "Generic assembly only allows "
                "the definition of scalar function with one or two arguments");
    { // Only for syntax analysis
      base_vector t(1);
      ga_workspace workspace;
      workspace.add_fixed_size_variable("t", gmm::sub_interval(0,1), t);
      if (nbargs == 2)
        workspace.add_fixed_size_variable("u", gmm::sub_interval(0,1), t);
      workspace.add_function_expression(expr);
    }

    PREDEF_FUNCTIONS[name] = ga_predef_function(expr);
    ga_predef_function &F = PREDEF_FUNCTIONS[name];
    F.gis = std::make_unique<instruction_set>();
    for (size_type thread = 0; thread < num_threads(); ++thread)
    {
      F.workspace(thread).add_fixed_size_variable("t", gmm::sub_interval(0,1),
                                                  F.t(thread));
      if (nbargs == 2)
        F.workspace(thread).add_fixed_size_variable("u", gmm::sub_interval(0,1),
                                                    F.u(thread));
      F.workspace(thread).add_function_expression(expr);
      ga_compile_function(F.workspace(thread), (*F.gis)(thread), true);
    }
    F.nbargs_ = nbargs;
    if (nbargs == 1) {
      if (der1.size()) { F.derivative1_ = der1; F.dtype_ = 2; }
    } else {
      if (der1.size() && der2.size()) {
        F.derivative1_ = der1;  F.derivative2_ = der2; F.dtype_ = 2;
      }
    }
  }

  void ga_define_function(const std::string &name, pscalar_func_onearg f,
                          const std::string &der) {
    auto guard = omp_guard{};
    ga_predef_function_tab &PREDEF_FUNCTIONS
      = dal::singleton<ga_predef_function_tab>::instance(0);
    PREDEF_FUNCTIONS[name] = ga_predef_function(f, 1, der);
    ga_predef_function &F = PREDEF_FUNCTIONS[name];
    if (der.size() == 0) F.dtype_ = 0;
    else if (!(ga_function_exists(der))) F.dtype_ = 2;
  }

  void ga_define_function(const std::string &name, pscalar_func_twoargs f,
                          const std::string &der1, const std::string &der2) {
    auto guard = omp_guard{};
    ga_predef_function_tab &PREDEF_FUNCTIONS
      = dal::singleton<ga_predef_function_tab>::instance(0);
    PREDEF_FUNCTIONS[name] = ga_predef_function(f, 1, der1, der2);
    ga_predef_function &F = PREDEF_FUNCTIONS[name];
    if (der1.size() == 0 || der2.size() == 0)
      F.dtype_ = 0;
    else if (!(ga_function_exists(der1)) || !(ga_function_exists(der2)))
      F.dtype_ = 2;
  }

  void ga_undefine_function(const std::string &name) {
    auto guard = omp_guard{};
    ga_predef_function_tab &PREDEF_FUNCTIONS
      = dal::singleton<ga_predef_function_tab>::instance(0);
    ga_predef_function_tab::iterator it = PREDEF_FUNCTIONS.find(name);
    if (it != PREDEF_FUNCTIONS.end()) {
      PREDEF_FUNCTIONS.erase(name);
      std::string name0 = "DER_PDFUNC_" + name;
      ga_undefine_function(name0);
      std::string name1 = "DER_PDFUNC1_" + name;
      ga_undefine_function(name1);
      std::string name2 = "DER_PDFUNC2_" + name;
      ga_undefine_function(name2);
    }
  }

  //=========================================================================
  // User defined functions
  //=========================================================================

  ga_function::ga_function(const ga_workspace &workspace_,
                           const std::string &e)
    : local_workspace(true, workspace_), expr(e), gis(0) {}

  ga_function::ga_function(const model &md, const std::string &e)
    : local_workspace(md), expr(e), gis(0) {}

  ga_function::ga_function(const std::string &e)
    : local_workspace(), expr(e), gis(0) {}

  ga_function::ga_function(const ga_function &gaf)
    : local_workspace(gaf.local_workspace), expr(gaf.expr), gis(0)
  { if (gaf.gis) compile(); }

  void ga_function::compile() const {
    if (gis) delete gis;
    gis = new ga_instruction_set;
    local_workspace.clear_expressions();
    local_workspace.add_function_expression(expr);
    ga_compile_function(local_workspace, *gis, false);
  }

  ga_function &ga_function::operator =(const ga_function &gaf) {
    if (gis) delete gis;
    gis = 0;
    local_workspace = gaf.local_workspace;
    expr = gaf.expr;
    if (gaf.gis) compile();
    return *this;
  }

  ga_function::~ga_function() { if (gis) delete gis; gis = 0; }

  const base_tensor &ga_function::eval() const {
    GMM_ASSERT1(gis, "Uncompiled function");
    gmm::clear(local_workspace.assembled_tensor().as_vector());
    ga_function_exec(*gis);
    return local_workspace.assembled_tensor();
  }

  void ga_function::derivative(const std::string &var) {
    GMM_ASSERT1(gis, "Uncompiled function");
    if (local_workspace.nb_trees()) {
      ga_tree tree = *(local_workspace.tree_info(0).ptree);
      ga_derivative(tree, local_workspace, dummy_mesh(), var, "", 1);
      if (tree.root) {
        ga_semantic_analysis(tree, local_workspace, dummy_mesh(),
			     1, false, true);
        // To be improved to suppress test functions in the expression ...
        // ga_replace_test_by_cte do not work in all operations like
        // vector components x(1)
        // ga_replace_test_by_cte(tree.root, false);
        // ga_semantic_analysis(tree, local_workspace, dummy_mesh(), 1,
	//                      false, true);
      }
      expr = ga_tree_to_string(tree);
    }
    if (gis) delete gis;
    gis = 0;
    compile();
  }

} /* end of namespace */
