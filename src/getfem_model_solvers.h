// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_model_solvers.h : Standard solvers for the brick
//           system.
// Date    : June 15, 2004.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2004-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//


/**
   @file getfem_model_solvers.h
   @brief Standard solver for model bricks
   @see getfem_modeling.h
*/

#ifndef GETFEM_MODEL_SOLVERS_H__
#define GETFEM_MODEL_SOLVERS_H__

#include <getfem_modeling.h>
#include <getfem_Newton.h>
#include <gmm_MUMPS_interface.h>

namespace getfem {

  template <typename T> struct sort_abs_val_
  { bool operator()(T x, T y) { return (gmm::abs(x) < gmm::abs(y)); } };

  /** A default solver for the model brick system.  
      
  Of course it could be not very well suited for a particular
  problem, so it could be copied and adapted to change solvers,
  add a special traitement on the problem, etc ...  This is in
  fact a model for your own solver.

  For small problems, a direct solver is used
  (getfem::SuperLU_solve), for larger problems, a conjugate
  gradient gmm::cg (if the problem is coercive) or a gmm::gmres is
  used (preconditionned with an incomplete factorization).

  When MPI/METIS is enabled, a partition is done via METIS, and the
  gmm::additive_schwarz is used as the parallel solver.

  @ingroup bricks
  */
  template <typename MODEL_STATE> void
  standard_solve(MODEL_STATE &MS, mdbrick_abstract<MODEL_STATE> &problem,
		 gmm::iteration &iter) {

    TYPEDEF_MODEL_STATE_TYPES;

    // TODO : take iter into account for the Newton. compute a consistent 
    //        max residu.

    size_type ndof = problem.nb_dof();

    bool is_linear = problem.is_linear();
    std::auto_ptr<abstract_newton_line_search>
      pline_search(new basic_newton_line_search2());

    R alpha(0);
    VECTOR stateinit;

    dal::bit_vector mixvar;
    gmm::iteration iter_linsolv0 = iter;
    iter_linsolv0.set_maxiter(10000);
    if (!is_linear) { 
      iter_linsolv0.reduce_noisy();
      iter_linsolv0.set_resmax(iter.get_resmax()/100.0);
    }


#if GETFEM_PARA_LEVEL > 1
    double t_init = MPI_Wtime();
#endif
    MS.adapt_sizes(problem); // to be sure it is ok, but should be done before
    problem.compute_residual(MS);
#if GETFEM_PARA_LEVEL > 1
    cout << "comput  residual  time = " << MPI_Wtime() - t_init << endl;
    t_init = MPI_Wtime();
#endif
    problem.compute_tangent_matrix(MS);
#if GETFEM_PARA_LEVEL > 1
    cout << "comput tangent time = " << MPI_Wtime() - t_init << endl;
    // cout << "CM = " << MS.constraints_matrix() << endl;
    // cout << "TM = " << MS.tangent_matrix() << endl;
    t_init = MPI_Wtime();
#endif
    MS.compute_reduced_system();

#if GETFEM_PARA_LEVEL > 1
    cout << "comput reduced system time = " << MPI_Wtime() - t_init << endl;
#endif

#if GETFEM_PARA_LEVEL > 1 && GETFEM_PARA_SOLVER == SCHWARZADD_PARA_SOLVER
    /* use a domain partition ? */
    double t_ref = MPI_Wtime();
    std::set<const mesh *> mesh_set;
    for (size_type i = 0; i < problem.nb_mesh_fems(); ++i) {
      mesh_set.insert(&(problem.get_mesh_fem(i).linked_mesh()));
    }

    cout << "You have " << mesh_set.size() << " meshes\n";

    std::vector< std::vector<int> > eparts(mesh_set.size());
    size_type nset = 0;
    int nparts = 64;//(ndof / 1000)+1; // number of sub-domains.

    // Quand il y a plusieurs maillages, on découpe tous les maillages en
    // autant de parties
    // et on regroupe les ddl de chaque partition par numÃ©ro de sous-partie.
    for (std::set<const mesh *>::iterator it = mesh_set.begin();
	 it != mesh_set.end(); ++it, ++nset) {
      int ne = int((*it)->nb_convex());
      int nn = int((*it)->nb_points()), k = 0, etype = 0, numflag = 0;
      int edgecut;
      
      bgeot::pconvex_structure cvs
	= (*it)->structure_of_convex((*it)->convex_index().first())->basic_structure();
      
      if (cvs == bgeot::simplex_structure(2)) { k = 3; etype = 1; }
      else if (cvs == bgeot::simplex_structure(3)) { k = 4; etype = 2; }
      else if (cvs == bgeot::parallelepiped_structure(2)) { k = 4; etype = 4; }
      else if (cvs == bgeot::parallelepiped_structure(3)) { k = 8; etype = 3; }
      else DAL_THROW(failure_error,
		     "This kind of element is not taken into account");

      
      std::vector<int> elmnts(ne*k), npart(nn);
      eparts[nset].resize(ne);
      int j = 0;
      // should be adapted for high order geotrans
      for (dal::bv_visitor i((*it)->convex_index()); !i.finished(); ++i, ++j)
	for (int l = 0; l < k; ++l)
	  elmnts[j*k+l] = (*it)->ind_points_of_convex(i)[l];

      METIS_PartMeshNodal(&ne, &nn, &(elmnts[0]), &etype, &numflag, &nparts,
			  &edgecut, &(eparts[nset][0]), &(npart[0]));
    }

    std::vector<dal::bit_vector> Bidof(nparts);
    size_type apos = 0;
    for (size_type i = 0; i < problem.nb_mesh_fems(); ++i) {
      const mesh_fem &mf = problem.get_mesh_fem(i);
      nset = 0;
      for (std::set<const mesh *>::iterator it = mesh_set.begin();
	   it != mesh_set.end(); ++it, ++nset)
	if (*it == &(mf.linked_mesh())) break; 
      size_type pos = problem.get_mesh_fem_position(i);
      if (pos != apos)
	DAL_THROW(failure_error, "Multiplicators are not taken into account");
      size_type length = mf.nb_dof();
      apos += length;
      for (dal::bv_visitor j(mf.convex_index()); !j.finished(); ++j) {
	size_type k = eparts[nset][j];
	for (size_type l = 0; l < mf.nb_dof_of_element(i); ++l)
	  Bidof[k].add(mf.ind_dof_of_element(j)[l] + pos);
      }
      if (apos != ndof)
	DAL_THROW(failure_error, "Multiplicators are not taken into account");
    }

    std::vector< gmm::row_matrix< gmm::rsvector<value_type> > > Bi(nparts);    
    for (size_type i = 0; i < size_type(nparts); ++i) {
      gmm::resize(Bi[i], ndof, Bidof[i].card());
      size_type k = 0;
      for (dal::bv_visitor j(Bidof[i]); !j.finished(); ++j, ++k)
	Bi[i](j, k) = value_type(1);
    }

    std::vector< gmm::row_matrix< gmm::rsvector<value_type> > > Bib(nparts);
    gmm::col_matrix< gmm::rsvector<value_type> > Bitemp;
    if (problem.nb_constraints() > 0) {
      for (size_type i = 0; i < size_type(nparts); ++i) {
	gmm::resize(Bib[i], gmm::mat_ncols(MS.nullspace_matrix()),
		    gmm::mat_ncols(Bi[i]));
       	gmm::mult(gmm::transposed(MS.nullspace_matrix()), Bi[i], Bib[i]);
	
	gmm::resize(Bitemp, gmm::mat_nrows(Bib[i]), gmm::mat_ncols(Bib[i]));
	gmm::copy(Bib[i], Bitemp);
	R EPS = gmm::mat_norminf(Bitemp) * gmm::default_tol(R());
	size_type k = 0;
	for (size_type j = 0; j < gmm::mat_ncols(Bitemp); ++j, ++k) {
	  // should do Schmidt orthogonalization for most sofisticated cases 
	  if (k != j) Bitemp.swap_col(j, k);
	  if (gmm::vect_norm2(gmm::mat_col(Bitemp, k)) < EPS) --k;
	}
	gmm::resize(Bitemp, gmm::mat_nrows(Bib[i]), k);
	gmm::resize(Bib[i], gmm::mat_nrows(Bib[i]), k);
	gmm::copy(Bitemp, Bib[i]);
	// cout << "Bib[" << i << "] = " <<  Bib[i] << endl;
	// cout << "Bi[" << i << "] = " <<  Bi[i] << endl;
      }
    } else std::swap(Bi, Bib);    
    cout << "METIS time = " << MPI_Wtime() - t_ref << endl;
#endif

    // cout << "RTM = " << MS.reduced_tangent_matrix() << endl;
    
    
    R act_res = MS.reduced_residual_norm(), act_res_new(0);

    while (is_linear || !iter.finished(act_res)) {
      
      size_type nreddof = gmm::vect_size(MS.reduced_residual());
      gmm::iteration iter_linsolv = iter_linsolv0;
      VECTOR d(ndof), dr(nreddof);
      
      if (!(iter.first())) {
	problem.compute_tangent_matrix(MS);
	MS.compute_reduced_system();
#if GETFEM_PARA_LEVEL > 1
        DAL_THROW(failure_error, "oups ...");
#endif
      }
      
      //       if (iter.get_noisy())
      // cout << "tangent matrix " << MS.tangent_matrix() << endl;
      //      cout << "tangent matrix is "
      // 	   << (gmm::is_symmetric(MS.tangent_matrix(),
      //          1E-6 * gmm::mat_maxnorm(MS.tangent_matrix())) ? "" : "not ")
      // 	   <<  "symmetric. ";

      //      cout << "MM = " << MS.reduced_tangent_matrix() << endl;

//       gmm::dense_matrix<value_type> MM(nreddof,nreddof), Q(nreddof,nreddof);
//       std::vector<value_type> eigval(nreddof);
//       gmm::copy(MS.reduced_tangent_matrix(), MM);
//       // gmm::symmetric_qr_algorithm(MM, eigval, Q);
//       gmm::implicit_qr_algorithm(MM, eigval, Q);
//       std::sort(eigval.begin(), eigval.end(), sort_abs_val_<value_type>());
//       cout << "eival = " << eigval << endl;
//       cout << "vectp : " << gmm::mat_col(Q, nreddof-1) << endl;
//       cout << "vectp : " << gmm::mat_col(Q, nreddof-2) << endl;
//       double emax, emin;
//       cout << "condition number" << condition_number(MM,emax,emin) << endl;
      //       cout << "emin = " << emin << endl;
      //       cout << "emax = " << emax << endl;

      if (iter.get_noisy()) {
	cout << "there are " << gmm::mat_nrows(MS.constraints_matrix())
	     << " constraints\n";
	mixvar = problem.mixed_variables();
	cout << "there are " << mixvar.card() << " mixed variables\n";
      }
      size_type dim = problem.dim();

#if (GETFEM_PARA_LEVEL > 1) && (GETFEM_PARA_SOLVER == SCHWARZADD_PARA_SOLVER)
      double t_ref,t_final;
      t_ref=MPI_Wtime();
      cout<<"begin Seq AS"<<endl; // à remmetre en fonction : utilisation
      // de mpi_distributed_matrix a revoir : mettre une référence.
      additive_schwarz(gmm::mpi_distributed_matrix<..>(MS.reduced_tangent_matrix()), dr,
		       gmm::scaled(MS.reduced_residual(), value_type(-1)),
		       gmm::identity_matrix(), Bib, iter_linsolv, gmm::using_superlu(),
		       gmm::using_cg());
      t_final=MPI_Wtime();
      cout<<"temps Seq AS "<< t_final-t_ref<<endl;
#elif GETFEM_PARA_LEVEL == 1 && GETFEM_PARA_SOLVER == MUMPS_PARA_SOLVER
      MUMPS_solve(MS.reduced_tangent_matrix(), dr,
		  gmm::scaled(MS.reduced_residual(), value_type(-1)));
#elif GETFEM_PARA_LEVEL > 1 && GETFEM_PARA_SOLVER == MUMPS_PARA_SOLVER
      double tt_ref=MPI_Wtime();
      MUMPS_distributed_matrix_solve(MS.reduced_tangent_matrix(), dr,
				     gmm::scaled(MS.reduced_residual(),
						 value_type(-1)));
      cout<<"temps MUMPS "<< MPI_Wtime() - tt_ref<<endl;
#else
      // if (0) {
      if (
#ifdef GMM_USES_MUMPS
	  (ndof < 200000 && dim <= 2) || (ndof < 100000 && dim <= 3)
	  || (ndof < 1000)
#else  
	  (ndof < 200000 && dim <= 2) || (ndof < 15000 && dim <= 3)
	  || (ndof < 1000)
#endif
	  ) {
	
	// cout << "M = " << MS.reduced_tangent_matrix() << endl;
	// cout << "L = " << MS.reduced_residual() << endl;
	
	double t = dal::uclock_sec();
#ifdef GMM_USES_MUMPS
	DAL_TRACE2("Solving with MUMPS\n");
	MUMPS_solve(MS.reduced_tangent_matrix(), dr,
		    gmm::scaled(MS.reduced_residual(), value_type(-1)));
#else
	double rcond;
	SuperLU_solve(MS.reduced_tangent_matrix(), dr,
		      gmm::scaled(MS.reduced_residual(), value_type(-1)),
		      rcond);
	if (iter.get_noisy()) cout << "condition number: " << 1.0/rcond<< endl;
#endif
	cout << "resolution time = " << dal::uclock_sec() - t << endl;
      }
      else {
	if (problem.is_coercive()) {
	  gmm::ildlt_precond<T_MATRIX> P(MS.reduced_tangent_matrix());
	  gmm::cg(MS.reduced_tangent_matrix(), dr, 
		  gmm::scaled(MS.reduced_residual(), value_type(-1)),
		  P, iter_linsolv);
	  if (!iter_linsolv.converged()) DAL_WARNING2("cg did not converge!");
	} else {
	  if (mixvar.card() == 0) {
	    gmm::ilu_precond<T_MATRIX> P(MS.reduced_tangent_matrix());
	  
	    gmm::gmres(MS.reduced_tangent_matrix(), dr, 
		       gmm::scaled(MS.reduced_residual(),  value_type(-1)), P,
		       300, iter_linsolv);
	  }
	  else {
	    gmm::ilut_precond<T_MATRIX> P(MS.reduced_tangent_matrix(),
					  20, 1E-10);
	    // gmm::identity_matrix P;
	    gmm::gmres(MS.reduced_tangent_matrix(), dr, 
		       gmm::scaled(MS.reduced_residual(),  value_type(-1)),
		       P, 300, iter_linsolv);
	  }
	  if (!iter_linsolv.converged())
	    DAL_WARNING2("gmres did not converge!");
	}
      } 
#endif // GETFEM_PARA_LEVEL < 2

      MS.unreduced_solution(dr,d);
    
      if (is_linear) {
	gmm::add(d, MS.state());
	return;
      }
      else { // line search for the non-linear case.
	gmm::resize(stateinit, ndof);
	gmm::copy(MS.state(), stateinit);
      
	pline_search->init_search(act_res);
	do {
	  alpha = pline_search->next_try();
	  gmm::add(stateinit, gmm::scaled(d, alpha), MS.state());
	  problem.compute_residual(MS);
	  MS.compute_reduced_residual();
	  act_res_new = MS.reduced_residual_norm();
	} while (!pline_search->is_converged(act_res_new));

	if (alpha != pline_search->converged_value()) {
	  alpha = pline_search->converged_value();
	  gmm::add(stateinit, gmm::scaled(d, alpha), MS.state());
	  act_res_new = pline_search->converged_residual();
	}
      }
      act_res = act_res_new; ++iter;
    
      if (iter.get_noisy()) cout << "alpha = " << alpha << " ";
    }
    
  }
  

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MODEL_SOLVERS_H__  */
