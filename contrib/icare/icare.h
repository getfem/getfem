// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2007-2008 Yves Renard, Julien Pommier.
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
// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
//  templates or use macros or inline functions from this file, or you compile
//  this file and link it with other files to produce an executable, this
//  file does not by itself cause the resulting executable to be covered by
//  the GNU General Public License.  This exception does not however
//  invalidate any other reasons why the executable file might be covered by
//  the GNU General Public License.
// 
//===========================================================================

#ifndef NAVIER_STOKES_H_
#define NAVIER_STOKES_H_

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_modeling.h"
#include "gmm/gmm.h"

namespace getfem {

  /*
   * div(uuT) term. 
   */
  
  template<typename MAT, typename VEC>
  void asm_NS_uuT(const MAT &M, const mesh_im &mim,
		  const mesh_fem &mf, const VEC &U0,
		  const mesh_region &rg = mesh_region::all_convexes()) {
    generic_assembly assem;    
    assem.set("u=data(#1); "
	      "t = comp(vBase(#1).vBase(#1).vGrad(#1));"
	      "M(#1, #1)+=u(i).t(i,k,:,j,:,j,k)"
	      // ";M(#1, #1)+=u(i).t(i,k,:,k,:,j,j)" // optional (or *0.5)
	      );
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(U0);
    assem.push_mat(const_cast<MAT&>(M));
    assem.assembly(rg);
  }


  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_NS_uuT : public mdbrick_abstract<MODEL_STATE> {

    TYPEDEF_MODEL_STATE_TYPES;
    
    mdbrick_abstract<MODEL_STATE> &sub_problem;
    VECTOR U0;
    size_type num_fem, i1, nbd;
    T_MATRIX K;

    virtual void proper_update(void) {}

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					size_type) {
      const mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      gmm::sub_interval SUBI(i0+this->mesh_fem_positions[num_fem],
			     mf_u.nb_dof());
      gmm::resize(K, mf_u.nb_dof(), mf_u.nb_dof());
      gmm::clear(K);
      asm_NS_uuT(K, *(this->mesh_ims[0]), mf_u, U0, 
		 mf_u.linked_mesh().get_mpi_region());
      gmm::add(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				size_type) {
      const mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      gmm::sub_interval SUBI(i0+this->mesh_fem_positions[num_fem],
			     mf_u.nb_dof());
      typename gmm::sub_vector_type<VECTOR *, gmm::sub_interval>::vector_type
	SUBV = gmm::sub_vector(MS.residual(), SUBI);
      gmm::mult_add(K, gmm::sub_vector(MS.state(), SUBI), SUBV);
    }

    template <class VEC> void set_U0(const VEC &q) {
      gmm::resize(U0, gmm::vect_size(q));
      gmm::copy(q, U0);
    }

    void init_(void) {
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = false;
      this->proper_is_symmetric_ = false;
      this->force_update();
    }

    // Constructor which homogeneous diagonal Q
    mdbrick_NS_uuT(mdbrick_abstract<MODEL_STATE> &problem,
		   size_type num_fem_=0) 
      : sub_problem(problem), num_fem(num_fem_) { init_(); }

  };



  /* ******************************************************************** */
  /*		A non-reflective condition.                               */
  /* ******************************************************************** */
  // idée: construire la condition non-réflective comme une condition de dirichlet en explicitant tous les termes: Unp1= Un+dt*Un.N*lamda




// initialisation des CL
// résolution du pb -> nouveau lambda<--| 
//mise à jour de CL de sortie        -->|  



# define MDBRICK_NS_NONREF1 3212435

 template<typename VECTOR, typename VECT>
  void asm_nonref_right_hand_side(VECTOR &R,
				  const mesh_im &mim,
				  const mesh_fem &mf_u,
				  const mesh_fem &mf_mult,
				  const VECT &Un,
				  const mesh_region &rg) {
    generic_assembly assem;
    // construction du terme de droite dans [M]*Unp1=F
    
    // mise en place de Un.N + Un.N*(dUn/dn).N

    std::stringstream ss;
    ss << "u=data$1(#1); "
    "sel=data$1(#2);"
     "V(#1)+=comp(vBase(#1).vBase(#2))(j,i,:,i).u(j)";
 
    assem.set(ss.str());
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_mult);
    assem.push_data(Un);
    assem.push_vec(R);
    assem.assembly(rg);
  }

    

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_NS_nonref1 : public mdbrick_constraint<MODEL_STATE> {
    
    TYPEDEF_MODEL_STATE_TYPES;

    size_type boundary;
    bool mfdata_set, B_to_be_computed, Un_is_modified;
    gmm::sub_index SUB_CT;
    const mesh_fem *mf_mult;
    VECTOR Un;
    scalar_type dt;
    
    const mesh_fem &mf_u() { return *(this->mesh_fems[this->num_fem]); }
    const mesh_im  &mim() { return *(this->mesh_ims[0]); }

    void compute_constraints(unsigned version) {
      size_type ndu = mf_u().nb_dof(), ndm = mf_mult->nb_dof();
      gmm::row_matrix<gmm::rsvector<value_type> > M(ndm, ndu);
      VECTOR V(ndm);
      GMM_TRACE2("Assembling Dirichlet constraints, version " << version);
      if ((version | ASMDIR_BUILDH)) {
	asm_mass_matrix(M, mim(), mf_u(), *mf_mult,
			mf_u().linked_mesh().get_mpi_region());
      }
      if ((version | ASMDIR_BUILDR)) {
	asm_nonref_right_hand_side(V, mim(), mf_u(), *mf_mult, Un, mf_u().linked_mesh().get_mpi_sub_region(boundary));
      }
      if (version & ASMDIR_BUILDH)
	gmm::copy(gmm::sub_matrix(M, SUB_CT, gmm::sub_interval(0, ndu)), 
		  this->B);
      gmm::copy(gmm::sub_vector(V, SUB_CT), this->CRHS);
    }

    virtual void recompute_B_sizes(void) {
      size_type nd = mf_u().nb_dof();
      dal::bit_vector dof_on_bound = mf_mult->dof_on_set(boundary);
      size_type nb_const = dof_on_bound.card();
      std::vector<size_type> ind_ct;
      for (dal::bv_visitor i(dof_on_bound); !i.finished(); ++i)
	ind_ct.push_back(i);
      SUB_CT = gmm::sub_index(ind_ct);
      gmm::resize(this->B, nb_const, nd);
      gmm::resize(this->CRHS, nb_const);
      B_to_be_computed = true;
    }

    virtual void recompute_B(void) {
      unsigned version = 0;
      if (Un_is_modified) { version = ASMDIR_BUILDR; }
      if (B_to_be_computed) { version = ASMDIR_BUILDR | ASMDIR_BUILDH; }
      if (version) { 
	compute_constraints(version);
	this->parameters_set_uptodate();
	B_to_be_computed = false; Un_is_modified = false;
      }
    }

  public :

    void set_Un(const VECTOR &U) {
      Un = U;
      Un_is_modified = true;
    }

    mdbrick_NS_nonref1(mdbrick_abstract<MODEL_STATE> &problem,
		       size_type bound, scalar_type dt_,
		       const mesh_fem &mf_mult_ = dummy_mesh_fem(), 
		       size_type num_fem_=0)
      : mdbrick_constraint<MODEL_STATE>(problem, num_fem_), 
	boundary(bound) {
      dt = dt_;
      mf_mult = (&mf_mult_ == &dummy_mesh_fem()) ? &(mf_u()) : &mf_mult_;
      GMM_ASSERT1(mf_mult->get_qdim() == mf_u().get_qdim(),
		  "The lagrange multipliers mesh fem "
		  "for the mdbrick_NS_nonref1 brick should have the same "
		  "Qdim as the main mesh_fem");

      this->add_proper_boundary_info(this->num_fem, boundary, 
				     MDBRICK_DIRICHLET);
      this->add_dependency(*mf_mult);
      mfdata_set = false; B_to_be_computed = true; Un_is_modified = false;
      this->force_update();
    }
  };

}
#endif
