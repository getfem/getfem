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

===========================================================================*/

/** @file getfem_continuation.cc
    @author Tomas Ligursky <tomas.ligursky@gmail.com>
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @author Konstantinos Poulios <logari81@googlemail.com>
    @date January 12, 2014.
    @brief inexact Moore-Penrose continuation method.
*/

#include "getfem/getfem_continuation.h"

namespace getfem {

  void cont_struct_getfem_model::set_variables
  (const base_vector &x, double gamma) const {
    md->set_real_variable(parameter_name)[0] = gamma;
    if (!currentdata_name.empty()) {
      gmm::add(gmm::scaled(md->real_variable(initdata_name), 1. - gamma),
               gmm::scaled(md->real_variable(finaldata_name), gamma),
               md->set_real_variable(currentdata_name));
    }
    md->to_variables(x);
  }

  void cont_struct_getfem_model::update_matrix
  (const base_vector &x, double gamma) const {
    set_variables(x, gamma);
    if (noisy() > 2) cout << "starting computing tangent matrix" << endl;
    md->assembly(model::BUILD_MATRIX);
  }

  // solve A * g = L
  void cont_struct_getfem_model::solve
  (const model_real_sparse_matrix &A, base_vector &g,
   const base_vector &L) const {
    if (noisy() > 2) cout << "starting linear solver" << endl;
    gmm::iteration iter(maxres_solve, (noisy() >= 2) ? noisy() - 2 : 0,
                        40000);
    (*lsolver)(A, g, L, iter);
    if (noisy() > 2) cout << "linear solver done" << endl;
  }

  // solve A * (g1|g2) = (L1|L2)
  void cont_struct_getfem_model::solve
  (const model_real_sparse_matrix &A, base_vector &g1, base_vector &g2,
   const base_vector &L1, const base_vector &L2) const {
    if (noisy() > 2) cout << "starting linear solver" << endl;
    gmm::iteration iter(maxres_solve, (noisy() >= 2) ? noisy() - 2 : 0,
                        40000);
    (*lsolver)(A, g1, L1, iter);
    iter.init(); (*lsolver)(A, g2, L2, iter); // (can be optimised)
    if (noisy() > 2) cout << "linear solver done" << endl;
  }

  // F(x, gamma) --> f
  void cont_struct_getfem_model::F
  (const base_vector &x, double gamma, base_vector &f) const {
    set_variables(x, gamma);
    md->assembly(model::BUILD_RHS);
    gmm::copy(gmm::scaled(md->real_rhs(), -1.), f);
  }

  // (F(x, gamma + eps) - f0) / eps --> g
  void cont_struct_getfem_model::F_gamma
  (const base_vector &x, double gamma, const base_vector &f0,
   base_vector &g) const {
    const double eps = diffeps;
    F(x, gamma + eps, g);
    gmm::add(gmm::scaled(f0, -1.), g);
    gmm::scale(g, 1./eps);
  }

  // (F(x, gamma + eps) - F(x, gamma)) / eps --> g
  void cont_struct_getfem_model::F_gamma
  (const base_vector &x, double gamma, base_vector &g) const {
    base_vector f0(x);
    F(x, gamma, f0);
    F_gamma(x, gamma, f0, g);
  }

  // F_x(x, gamma) --> A
  void cont_struct_getfem_model::F_x
  (const base_vector &x, double gamma, model_real_sparse_matrix &A) const {
    update_matrix(x, gamma);
    size_type nbdof = md->nb_dof();
    gmm::resize(A, nbdof, nbdof);
    gmm::copy(md->real_tangent_matrix(), A);
  }

  // solve F_x(x, gamma) * g = L
  void cont_struct_getfem_model::solve_grad
    (const base_vector &x, double gamma, base_vector &g,
     const base_vector &L) const {
    update_matrix(x, gamma);
    solve(md->real_tangent_matrix(), g, L);
//int exp;
//cout<<"det="<<MUMPS_determinant(md->real_tangent_matrix(),exp)<<"*2^"<<exp<<endl;
  }

  // solve F_x(x, gamma) * (g1|g2) = (L1|L2)
  void cont_struct_getfem_model::solve_grad
  (const base_vector &x, double gamma, base_vector &g1, base_vector &g2,
   const base_vector &L1, const base_vector &L2) const {
    update_matrix(x, gamma);
    solve(md->real_tangent_matrix(), g1, g2, L1, L2);
  }

  // F_x(x, gamma) * w --> y
  void cont_struct_getfem_model::mult_grad
  (const base_vector &x, double gamma,
   const base_vector &w, base_vector &y) const {
    update_matrix(x, gamma);
    mult(md->real_tangent_matrix(), w, y);
  }

  size_type cont_struct_getfem_model::estimated_memsize(void) {
    return sizeof(cont_struct_getfem_model)
           + virtual_cont_struct::estimated_memsize();
  }

}  /* end of namespace getfem.                                            */
