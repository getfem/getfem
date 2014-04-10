/*===========================================================================
 
 Copyright (C) 2009-2012 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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


#include "getfem/getfem_fourth_order.h"
#include "getfem/getfem_omp.h"


namespace getfem {

  // ----------------------------------------------------------------------
  // Bilaplacian brick
  // ----------------------------------------------------------------------

  struct bilap_brick : public virtual_brick {

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
					const model::varnamelist &vl,
					const model::varnamelist &dl,
					const model::mimlist &mims,
					model::real_matlist &matl,
					model::real_veclist &,
					model::real_veclist &,
					size_type region,
					build_version version) const {
      GMM_ASSERT1(matl.size() == 1,
		  "Bilaplacian brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
		  "Bilaplacian brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && (dl.size() == 1 || dl.size() == 2),
		  "Wrong number of variables for bilaplacian brick");

      bool KL = (dl.size() == 2);

      bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0)
	|| md.is_var_newer_than_brick(dl[0], ib);


      if (recompute_matrix) {

	const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
	const mesh &m = mf_u.linked_mesh();
	size_type Q = mf_u.get_qdim();
	GMM_ASSERT1(Q == 1, "Bilaplacian brick is only for a scalar field");
	const mesh_im &mim = *mims[0];
	mesh_region rg(region);
	m.intersect_with_mpi_region(rg);
	const mesh_fem *mf_data = md.pmesh_fem_of_variable(dl[0]);
	const model_real_plain_vector *data = &(md.real_variable(dl[0]));
	size_type sl = gmm::vect_size(*data);
	if (mf_data) sl = sl * mf_data->get_qdim() / mf_data->nb_dof();
	
	GMM_ASSERT1(sl == 1, "Bad format of bilaplacian coefficient");

	const mesh_fem *mf_data2 = 0;
	const model_real_plain_vector *data2 = 0;
	if (KL) {
	  mf_data2 = md.pmesh_fem_of_variable(dl[1]);
	  data2 = &(md.real_variable(dl[1]));
	  size_type sl2 = gmm::vect_size(*data2);
	  if (mf_data2) sl = sl * mf_data2->get_qdim() / mf_data2->nb_dof();
	  GMM_ASSERT1(sl2 == 1, "Bad format of bilaplacian coefficient");
	}
	
	if (KL) {
	  GMM_TRACE2("Stiffness matrix assembly of a bilaplacian term for a "
		     "Kirchhoff-Love plate");
	}
	else {
	  GMM_TRACE2("Stiffness matrix assembly of a bilaplacian term");
	}

	gmm::clear(matl[0]);
	if (mf_data) {
	  if (KL)
	    asm_stiffness_matrix_for_bilaplacian_KL
	      (matl[0], mim, mf_u, *mf_data,  *data, *data2, rg);
	  else
	    asm_stiffness_matrix_for_bilaplacian
	      (matl[0], mim, mf_u, *mf_data,  *data, rg);
	} else {
	  if (KL) {
	    asm_stiffness_matrix_for_homogeneous_bilaplacian_KL
	      (matl[0], mim, mf_u,  *data, *data2, rg);
	  }
	  else
	    asm_stiffness_matrix_for_homogeneous_bilaplacian
	      (matl[0], mim, mf_u,  *data, rg);
	}
      }
    }

    virtual scalar_type asm_real_pseudo_potential(const model &md, size_type,
						  const model::varnamelist &vl,
						  const model::varnamelist &,
						  const model::mimlist &,
						  model::real_matlist &matl,
						  model::real_veclist &,
						  model::real_veclist &,
						  size_type) const {
      const model_real_plain_vector &u = md.real_variable(vl[0]);
      return gmm::vect_sp(matl[0], u, u) / scalar_type(2);
    }


    bilap_brick(void) {
      set_flags("Bilaplacian operator", true /* is linear*/,
		true /* is symmetric */, true /* is coercive */,
		true /* is real */, false /* is complex */);
    }

  };

  size_type add_bilaplacian_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname,
   size_type region) {
    pbrick pbr = new bilap_brick;
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist dl(1, dataname);
    return md.add_brick(pbr, model::varnamelist(1, varname), dl, tl,
			model::mimlist(1, &mim), region);
  }

  size_type add_bilaplacian_brick_KL
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname1, const std::string &dataname2,
   size_type region) {
    pbrick pbr = new bilap_brick;
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist dl(1, dataname1);
    dl.push_back(dataname2);
    return md.add_brick(pbr, model::varnamelist(1, varname), dl, tl,
			model::mimlist(1, &mim), region);
  }



  // ----------------------------------------------------------------------
  //
  // Normal derivative source term brick
  //
  // ----------------------------------------------------------------------

  struct normal_derivative_source_term_brick : public virtual_brick {

    virtual void asm_real_tangent_terms(const model &md, size_type,
					const model::varnamelist &vl,
					const model::varnamelist &dl,
					const model::mimlist &mims,
					model::real_matlist &,
					model::real_veclist &vecl,
					model::real_veclist &,
					size_type region,
					build_version) const {
      GMM_ASSERT1(vecl.size() == 1,
		  "Normal derivative source term brick has one and only "
		  "one term");
      GMM_ASSERT1(mims.size() == 1,
		  "Normal derivative source term brick need one and only "
		  "one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() == 1,
		  "Wrong number of variables for normal derivative "
		  "source term brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_im &mim = *mims[0];
      const model_real_plain_vector &A = md.real_variable(dl[0]);
      const mesh_fem *mf_data = md.pmesh_fem_of_variable(dl[0]);
      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      size_type s = gmm::vect_size(A);
      if (mf_data) s = s * mf_data->get_qdim() / mf_data->nb_dof();

      GMM_ASSERT1(s == mf_u.get_qdim()
		  || s == size_type(mf_u.get_qdim()*gmm::sqr(mf_u.linked_mesh().dim())),
		  dl[0] << ": bad format of normal derivative source term "
		  "data. Detected dimension is " << s << " should be "
		  << size_type(mf_u.get_qdim()));

      GMM_TRACE2("Normal derivative source term assembly");
      if (mf_data)
	asm_normal_derivative_source_term(vecl[0], mim, mf_u, *mf_data, A, rg);
      else
	asm_homogeneous_normal_derivative_source_term(vecl[0], mim, mf_u, A, rg);

    }

    virtual scalar_type asm_real_pseudo_potential(const model &md, size_type,
						  const model::varnamelist &vl,
						  const model::varnamelist &,
						  const model::mimlist &,
						  model::real_matlist &,
						  model::real_veclist &vecl,
						  model::real_veclist &,
						  size_type) const {
      const model_real_plain_vector &u = md.real_variable(vl[0]);
      return -gmm::vect_sp(vecl[0], u);
    }

    virtual void asm_complex_tangent_terms(const model &md, size_type,
					   const model::varnamelist &vl,
					   const model::varnamelist &dl,
					   const model::mimlist &mims,
					   model::complex_matlist &,
					   model::complex_veclist &vecl,
					   model::complex_veclist &,
					   size_type region,
					   build_version) const {
      GMM_ASSERT1(vecl.size() == 1,
		  "Normal derivative source term brick has one and only "
		  "one term");
      GMM_ASSERT1(mims.size() == 1,
		  "Normal derivative source term brick need one and only "
		  "one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() == 1,
		  "Wrong number of variables for normal derivative "
		  "source term brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_im &mim = *mims[0];
      const model_complex_plain_vector &A = md.complex_variable(dl[0]);
      const mesh_fem *mf_data = md.pmesh_fem_of_variable(dl[0]);
      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      size_type s = gmm::vect_size(A);
      if (mf_data) s = s * mf_data->get_qdim() / mf_data->nb_dof();

      GMM_ASSERT1(mf_u.get_qdim() == s,
		  dl[0] << ": bad format of normal derivative source term "
		  "data. Detected dimension is " << s << " should be "
		  << size_type(mf_u.get_qdim()));

      GMM_TRACE2("Normal derivative source term assembly");
      if (mf_data)
	asm_normal_derivative_source_term(vecl[0], mim, mf_u, *mf_data, A, rg);
      else
	asm_homogeneous_normal_derivative_source_term(vecl[0], mim, mf_u, A, rg);

    }

    virtual scalar_type asm_complex_pseudo_potential(const model &md,size_type,
						 const model::varnamelist &vl,
						 const model::varnamelist &,
						 const model::mimlist &,
						 model::complex_matlist &,
						 model::complex_veclist &vecl,
						 model::complex_veclist &,
						 size_type) const {
      const model_complex_plain_vector &u = md.complex_variable(vl[0]);
      return -gmm::real(gmm::vect_hp(vecl[0], u)); /* ? */
    }

    normal_derivative_source_term_brick(void) {
      set_flags("Normal derivative source term", true /* is linear*/,
		true /* is symmetric */, true /* is coercive */,
		true /* is real */, true /* is complex */,
		false /* compute each time */, false /* has a Neumann term */);
    }


  };

  size_type add_normal_derivative_source_term_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname, size_type region) {
    pbrick pbr = new normal_derivative_source_term_brick;
    model::termlist tl;
    tl.push_back(model::term_description(varname));
    model::varnamelist vdata(1, dataname);
    return md.add_brick(pbr, model::varnamelist(1, varname),
			vdata, tl, model::mimlist(1, &mim), region);
  }



  // ----------------------------------------------------------------------
  //
  // K.L. source term brick
  //
  // ----------------------------------------------------------------------

  struct KL_source_term_brick : public virtual_brick {

    virtual void asm_real_tangent_terms(const model &md, size_type,
					const model::varnamelist &vl,
					const model::varnamelist &dl,
					const model::mimlist &mims,
					model::real_matlist &,
					model::real_veclist &vecl,
					model::real_veclist &,
					size_type region,
					build_version) const {
      GMM_ASSERT1(vecl.size() == 1,
		  "Kirchoff Love source term brick has one and only "
		  "one term");
      GMM_ASSERT1(mims.size() == 1,
		  "Kirchoff Love source term brick need one and only "
		  "one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() == 2,
		  "Wrong number of variables for Kirchoff Love "
		  "source term brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_im &mim = *mims[0];
      const model_real_plain_vector &A = md.real_variable(dl[0]);
      const mesh_fem *mf_dataA = md.pmesh_fem_of_variable(dl[0]);
      const model_real_plain_vector &B = md.real_variable(dl[1]);
      const mesh_fem *mf_dataB = md.pmesh_fem_of_variable(dl[1]);
      size_type N = mf_u.linked_mesh().dim();
      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      size_type s = gmm::vect_size(A);
      if (mf_dataA) s = s * mf_dataA->get_qdim() / mf_dataA->nb_dof();
      GMM_ASSERT1(mf_u.get_qdim() == 1 && s == N*N,
		  dl[0] << ": bad format of Kirchoff Love Neumann term "
		  "data. Detected dimension is " << s << " should be "
		  << N*N);

      s = gmm::vect_size(B);
      if (mf_dataB) s = s * mf_dataB->get_qdim() / mf_dataB->nb_dof();
      GMM_ASSERT1(s == N,
		  dl[0] << ": bad format of Kirchoff Love Neumann term "
		  "data. Detected dimension is " << s << " should be "
		  << N);


      GMM_TRACE2("Kirchoff Love Neumann term assembly");
      if (mf_dataA)
	asm_neumann_KL_term(vecl[0], mim, mf_u, *mf_dataA, A, B, rg);
      else
	asm_neumann_KL_homogeneous_term(vecl[0], mim, mf_u, A, B, rg);

    }

    virtual scalar_type asm_real_pseudo_potential(const model &md, size_type,
						  const model::varnamelist &vl,
						  const model::varnamelist &,
						  const model::mimlist &,
						  model::real_matlist &,
						  model::real_veclist &vecl,
						  model::real_veclist &,
						  size_type) const {
      const model_real_plain_vector &u = md.real_variable(vl[0]);
      return -gmm::vect_sp(vecl[0], u);
    }

    KL_source_term_brick(void) {
      set_flags("Kirchoff Love Neumann term", true /* is linear*/,
		true /* is symmetric */, true /* is coercive */,
		true /* is real */, false /* is complex */,
		false /* compute each time */, false /* has a Neumann term */);
    }


  };

  size_type add_Kirchoff_Love_Neumann_term_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname1, const std::string &dataname2,
   size_type region) {
    pbrick pbr = new KL_source_term_brick;
    model::termlist tl;
    tl.push_back(model::term_description(varname));
    model::varnamelist vdata(1, dataname1);
    vdata.push_back(dataname2);
    return md.add_brick(pbr, model::varnamelist(1, varname),
			vdata, tl, model::mimlist(1, &mim), region);
  }
























  // ----------------------------------------------------------------------
  //
  // Normal derivative Dirichlet condition brick
  //
  // ----------------------------------------------------------------------
  // Two variables : with multipliers
  // One variable : penalization

  struct normal_derivative_Dirichlet_condition_brick : public virtual_brick {

    bool R_must_be_derivated;
    mutable getfem::omp_distribute<model_real_sparse_matrix> rB_th;
    mutable getfem::omp_distribute<model_real_plain_vector> rV_th;
    mutable getfem::omp_distribute<model_complex_sparse_matrix> cB_th;
    mutable getfem::omp_distribute<model_complex_plain_vector> cV_th;

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
					const model::varnamelist &vl,
					const model::varnamelist &dl,
					const model::mimlist &mims,
					model::real_matlist &matl,
					model::real_veclist &vecl,
					model::real_veclist &,
					size_type region,
					build_version version) const {
      GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
		  "Normal derivative Dirichlet condition brick has one and "
		  "only one term");
      GMM_ASSERT1(mims.size() == 1,
		  "Normal derivative Dirichlet condition brick need one and "
		  "only one mesh_im");
      GMM_ASSERT1(vl.size() >= 1 && vl.size() <= 2 && dl.size() <= 2,
		  "Wrong number of variables for normal derivative Dirichlet "
		  "condition brick");

      model_real_sparse_matrix& rB = rB_th;
      model_real_plain_vector&  rV = rV_th;

      bool penalized = (vl.size() == 1);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_fem &mf_mult = md.mesh_fem_of_variable(vl[vl.size()-1]);
      const mesh_im &mim = *mims[0];
      const model_real_plain_vector *A = 0, *COEFF = 0;
      const mesh_fem *mf_data = 0;
      bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0)
	|| (penalized && md.is_var_newer_than_brick(dl[0], ib));

      if (penalized) {
	COEFF = &(md.real_variable(dl[0]));
	GMM_ASSERT1(gmm::vect_size(*COEFF) == 1,
		    "Data for coefficient should be a scalar");
      }

      size_type s = 0, ind = (penalized ? 1 : 0);
      if (dl.size() > ind) {
	A = &(md.real_variable(dl[ind]));
	mf_data = md.pmesh_fem_of_variable(dl[ind]);
	s = gmm::vect_size(*A);
	if (mf_data) s = s * mf_data->get_qdim() / mf_data->nb_dof();
	GMM_ASSERT1(s == mf_u.get_qdim() || s == mf_u.linked_mesh().dim(),
		    dl[ind] << ": bad format of normal derivative Dirichlet "
		    "data. Detected dimension is " << s << " should be "
		    << size_type(mf_u.get_qdim()) << " or "
		    << mf_u.linked_mesh().dim());
      }

      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      if (recompute_matrix) {
	GMM_TRACE2("Mass term assembly for normal derivative Dirichlet "
		   "condition");
	if (penalized) {
	  gmm::resize(rB, mf_mult.nb_dof(), mf_u.nb_dof());
	  gmm::clear(rB);
	  asm_normal_derivative_dirichlet_constraints
	  (rB, vecl[0], mim, mf_u, mf_mult,
	   *mf_data, *A, rg, R_must_be_derivated, ASMDIR_BUILDH);
	} else {
	  gmm::clear(matl[0]);
	  asm_normal_derivative_dirichlet_constraints
	    (matl[0], vecl[0], mim, mf_u, mf_mult,
	     *mf_data, *A, rg, R_must_be_derivated, ASMDIR_BUILDH);
	}

	if (penalized) {
	  gmm::mult(gmm::transposed(rB), rB, matl[0]);
	  gmm::scale(matl[0], gmm::abs((*COEFF)[0]));
	}
      }

      if (dl.size() > ind) {
	GMM_TRACE2("Source term assembly for normal derivative Dirichlet "
		   "condition");
	model_real_plain_vector *R = penalized ? &rV : &(vecl[0]);
	if (penalized) { gmm::resize(rV, mf_mult.nb_dof()); gmm::clear(rV); }
	  
	if (mf_data) {
	  if (!R_must_be_derivated) {
	    if (s == mf_u.linked_mesh().dim())
	      asm_normal_source_term(*R, mim, mf_mult, *mf_data, *A, rg);
	    else
	      asm_source_term(*R, mim, mf_mult, *mf_data, *A, rg);
	  }
	  else {
	    asm_real_or_complex_1_param
	      (*R, mim, mf_mult, *mf_data, *A, rg,
	       "R=data(#2); V(#1)+=comp(Base(#1).Grad(#2).Normal())(:,i,j,j).R(i)");
	  }
	} else {
	  GMM_ASSERT1(!R_must_be_derivated, "Incoherent situation");
	  if (s == mf_u.linked_mesh().dim())
	    asm_homogeneous_normal_source_term(*R, mim, mf_mult, *A, rg);
	  else
	    asm_homogeneous_source_term(*R, mim, mf_mult, *A, rg);
	}
	if (penalized) {
	  gmm::mult(gmm::transposed(rB), rV, vecl[0]);
	  gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
	  rV = model_real_plain_vector();
	}
      }


    }

    virtual void asm_complex_tangent_terms(const model &md, size_type ib,
					   const model::varnamelist &vl,
					   const model::varnamelist &dl,
					   const model::mimlist &mims,
					   model::complex_matlist &matl,
					   model::complex_veclist &vecl,
					   model::complex_veclist &,
					   size_type region,
					   build_version version) const {
      GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
		  "Normal derivative Dirichlet condition brick has one and "
		  "only one term");
      GMM_ASSERT1(mims.size() == 1,
		  "Normal derivative Dirichlet condition brick need one and "
		  "only one mesh_im");
      GMM_ASSERT1(vl.size() >= 1 && vl.size() <= 2 && dl.size() <= 2,
		  "Wrong number of variables for normal derivative Dirichlet "
		  "condition brick");

      model_complex_sparse_matrix& cB = cB_th;
      model_complex_plain_vector&  cV = cV_th;

      bool penalized = (vl.size() == 1);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_fem &mf_mult = md.mesh_fem_of_variable(vl[vl.size()-1]);
      const mesh_im &mim = *mims[0];
      const model_complex_plain_vector *A = 0, *COEFF = 0;
      const mesh_fem *mf_data = 0;
      bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0)
	|| (penalized && md.is_var_newer_than_brick(dl[0], ib));

      if (penalized) {
	COEFF = &(md.complex_variable(dl[0]));
	GMM_ASSERT1(gmm::vect_size(*COEFF) == 1,
		    "Data for coefficient should be a scalar");
      }

      size_type s = 0, ind = (penalized ? 1 : 0);
      if (dl.size() > ind) {
	A = &(md.complex_variable(dl[ind]));
	mf_data = md.pmesh_fem_of_variable(dl[ind]);
	s = gmm::vect_size(*A);
	if (mf_data) s = s * mf_data->get_qdim() / mf_data->nb_dof();
	GMM_ASSERT1(s == mf_u.get_qdim() || s == mf_u.linked_mesh().dim(),
		    dl[ind] << ": bad format of normal derivative Dirichlet "
		    "data. Detected dimension is " << s << " should be "
		    << size_type(mf_u.get_qdim()) << " or "
		    << mf_u.linked_mesh().dim());
      }

      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      if (recompute_matrix) {
	GMM_TRACE2("Mass term assembly for normal derivative Dirichlet "
		   "condition");
	if (penalized) {
	  gmm::resize(cB, mf_mult.nb_dof(), mf_u.nb_dof());
	  gmm::clear(cB);
	  asm_normal_derivative_dirichlet_constraints
	  (cB, vecl[0], mim, mf_u, mf_mult,
	   *mf_data, *A, rg, R_must_be_derivated, ASMDIR_BUILDH);
	} else {
	  gmm::clear(matl[0]);
	  asm_normal_derivative_dirichlet_constraints
	    (matl[0], vecl[0], mim, mf_u, mf_mult,
	     *mf_data, *A, rg, R_must_be_derivated, ASMDIR_BUILDH);
	}

	if (penalized) {
	  gmm::mult(gmm::transposed(cB), cB, matl[0]);
	  gmm::scale(matl[0], gmm::abs((*COEFF)[0]));
	}
      }

      if (dl.size() > ind) {
	GMM_TRACE2("Source term assembly for normal derivative Dirichlet "
		   "condition");
	model_complex_plain_vector *R = penalized ? &cV : &(vecl[0]);
	if (penalized) { gmm::resize(cV, mf_mult.nb_dof()); gmm::clear(cV); }
	  
	if (mf_data) {
	  if (!R_must_be_derivated) {
	    if (s == mf_u.linked_mesh().dim())
	      asm_normal_source_term(*R, mim, mf_mult, *mf_data, *A, rg);
	    else
	      asm_source_term(*R, mim, mf_mult, *mf_data, *A, rg);
	  }
	  else {
	    asm_real_or_complex_1_param
	      (*R, mim, mf_mult, *mf_data, *A, rg,
	       "R=data(#2); V(#1)+=comp(Base(#1).Grad(#2).Normal())(:,i,j,j).R(i)");
	  }
	} else {
	  GMM_ASSERT1(!R_must_be_derivated, "Incoherent situation");
	  if (s == mf_u.linked_mesh().dim())
	    asm_homogeneous_normal_source_term(*R, mim, mf_mult, *A, rg);
	  else
	    asm_homogeneous_source_term(*R, mim, mf_mult, *A, rg);
	}
	if (penalized) {
	  gmm::mult(gmm::transposed(cB), cV, vecl[0]);
	  gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
	  cV = model_complex_plain_vector();
	}
      }
    }

    normal_derivative_Dirichlet_condition_brick
    (bool penalized, bool R_must_be_derivated_) {
      R_must_be_derivated = R_must_be_derivated_;
      set_flags(penalized ?
		"Normal derivative Dirichlet with penalization brick"
		: "Normal derivative Dirichlet with multipliers brick",
		true /* is linear*/,
		true /* is symmetric */, penalized /* is coercive */,
		true /* is real */, true /* is complex */,
		false /* compute each time */, false /* has a Neumann term */);
    }
  };

  size_type add_normal_derivative_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &multname, size_type region,
   const std::string &dataname, bool R_must_be_derivated) {
    pbrick pbr = new normal_derivative_Dirichlet_condition_brick(false, R_must_be_derivated);
    model::termlist tl;
    tl.push_back(model::term_description(multname, varname, true));
    model::varnamelist vl(1, varname);
    vl.push_back(multname);
    model::varnamelist dl;
    if (dataname.size()) dl.push_back(dataname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }

  size_type add_normal_derivative_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   const mesh_fem &mf_mult, size_type region,
   const std::string &dataname, bool R_must_be_derivated) {
    std::string multname = md.new_name("mult_on_" + varname);
    md.add_multiplier(multname, mf_mult, varname);
    return add_normal_derivative_Dirichlet_condition_with_multipliers
      (md, mim, varname, multname, region, dataname, R_must_be_derivated);
  }

  size_type add_normal_derivative_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   dim_type degree, size_type region,
   const std::string &dataname, bool R_must_be_derivated) {
    const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
    const mesh_fem &mf_mult = classical_mesh_fem(mf_u.linked_mesh(),
						 degree, mf_u.get_qdim());
    return add_normal_derivative_Dirichlet_condition_with_multipliers
      (md, mim, varname, mf_mult, region, dataname, R_must_be_derivated);
  }


  size_type add_normal_derivative_Dirichlet_condition_with_penalization
  (model &md, const mesh_im &mim, const std::string &varname,
   scalar_type penalisation_coeff, size_type region, 
   const std::string &dataname, bool R_must_be_derivated) {
    std::string coeffname = md.new_name("penalization_on_" + varname);
    md.add_fixed_size_data(coeffname, 1);
    if (md.is_complex())
      md.set_complex_variable(coeffname)[0] = penalisation_coeff;
    else
      md.set_real_variable(coeffname)[0] = penalisation_coeff;
    pbrick pbr = new normal_derivative_Dirichlet_condition_brick(true, R_must_be_derivated);
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist vl(1, varname);
    model::varnamelist dl(1, coeffname);
    if (dataname.size()) dl.push_back(dataname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }

  // + changement du coeff de penalisation avec la fonction de Dirichelt standard.
  


}  /* end of namespace getfem.                                             */

