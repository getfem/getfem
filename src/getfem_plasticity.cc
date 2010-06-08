// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2000-2010 Yves Renard
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
//===========================================================================


#include "getfem/getfem_models.h"
#include "getfem/getfem_plasticity.h"
#if 0
namespace getfem {



  //=========================================================================
  //
  //  Plasticity Brick
  //
  //=========================================================================

  struct plasticity_brick : public virtual_brick {



    virtual void asm_real_tangent_terms(const model &md, size_type /* ib */,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {
      GMM_ASSERT1(mims.size() == 1,
		  "Plasticity brick need a single mesh_im");
      GMM_ASSERT1(vl.size() == 2,
		  "Plasticity brick need two variables"); /** vl[0] = u */
      GMM_ASSERT1(dl.size() == 3,
		  "Wrong number of data for plasticity brick, "
                  << dl.size() << " should be 1 (2 order tensor).");
      GMM_ASSERT1(matl.size() == 1,  "Wrong number of terms for "
		  "plasticity brick");

      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const mesh_fem &mf_u = *(md.pmesh_fem_of_variable(vl[0]));

      const model_real_plain_vector &sigma = md.real_variable(vl[1]);
      const mesh_fem &mf_sigma = *(md.pmesh_fem_of_variable(vl[1]));
      

  /*    const mesh_fem &mf_params =*(md.pmesh_fem_of_variable(dl[0]));
      const model_real_plain_vector &lambda = md.real_variable(dl[0]);
      const model_real_plain_vector &mu = md.real_variable(dl[1]);
      const model_real_plain_vector &stress_threshold = md.real_variable(dl[2]);
      const model_real_plain_vector &sigma = md.real_variable(dl[3]);	
      const mesh_im &mim = *mims[0];
*/
	size-type N = mf_sigma.linked_mesh().dim());
      std::vector<std::vector<scalar_type> > sigma_bar ;
	for(dim_type i = 0; i<N; ++i)
	sigma_bar[mf_sigma.convex_index()][] = 0;
 
  //  std::vector<std::vector<scalar_type> > &saved_proj ;
	size_type flag_hyp;
     //   const abstract_constraints_projection *t_proj;
    //  mesh_region rg(region);
    //  mf_u.linked_mesh().intersect_with_mpi_region(rg);

    
	
      

      if (version & model::BUILD_MATRIX) {
	gmm::clear(matl[0]);
	GMM_TRACE2("Plasticity stiffness matrix assembly");

	
        flag_hyp = 1;
	//VM_projection gradproj(flag_hyp);
	//plasticity_projection gradproj(mim, mf_u, mf_params, u,stress_threshold, lambda, mu, t_proj, sigma_bar,  saved_proj, flag_hyp , false);
	//asm_lhs_for_plasticity(matl[0], mim, mf_u, mf_params,lambda, mu, &gradproj );
      }


      if (version & model::BUILD_RHS) {
        flag_hyp = 0;
      //  VM_projection proj(flag_hyp);
	//plasticity_projection proj(mim, mf_u, mf_params, u,stress_threshold, lambda, mu, t_proj, sigma_bar,  saved_proj, flag_hyp, false);
	//asm_rhs_for_plasticity(vecl[0], mim, mf_u, mf_params, &proj );
	gmm::scale(vecl[0], scalar_type(-1));
      }

    }



    plasticity_brick(void){
      set_flags("Plasticity brick", false /* is linear*/,
                true /* is symmetric */, false /* is coercive */,
		true /* is real */, false /* is complex */);
    }

  };
  
  //=========================================================================
  //  Add a plasticity brick
  //=========================================================================

  size_type add_plasticity_brick
  (model &md, const mesh_im &mim, const std::string &varname,const std::string &dataname,
   size_type region) {
    pbrick pbr = new plasticity_brick();

    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist dl(1, dataname);
    model::varnamelist vl(1, varname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1,&mim), region);
  }

}  /* end of namespace getfem.  */
#endif
