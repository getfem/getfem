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

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/** @file   getfem_generic_assembly_tree.h
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @date   November 18, 2013.
    @brief  Compilation and execution operations. */


#ifndef GETFEM_GENERIC_ASSEMBLY_COMP_EXEC_H__
#define GETFEM_GENERIC_ASSEMBLY_COMP_EXEC_H__

#include "getfem/getfem_generic_assembly_semantic.h"

namespace getfem {

  class ga_if_hierarchy : public std::vector<size_type> {
    
  public:
    void increment() { (back())++; }
    void child_of(const ga_if_hierarchy &gih)
    { *this = gih; push_back(0); }
    bool is_compatible(const std::list<ga_if_hierarchy> &gihl) {

      std::list<ga_if_hierarchy>::const_iterator it = gihl.begin();
      for (; it != gihl.end(); ++it) {
        if (it->size() <= size()) {
          bool ok = true;
          for (size_type i = 0; i+1 < it->size(); ++i)
            if ((*it)[i] != (*this)[i]) { ok = false; break; }
          if (it->back() > (*this)[it->size()-1]) { ok = false; break; }
          if (ok) return true;
        }
      }
      return false;
    }

    ga_if_hierarchy() : std::vector<size_type>(1) { (*this)[0] = 0; }
  };


  struct ga_instruction {
    virtual int exec() = 0;
    virtual ~ga_instruction() {};
  };

  typedef std::shared_ptr<ga_instruction> pga_instruction;
  typedef std::vector<pga_instruction> ga_instruction_list;

  
  struct gauss_pt_corresp { // For neighbour interpolation transformation
    bgeot::pgeometric_trans pgt1, pgt2;
    papprox_integration pai;
    std::vector<size_type> nodes;
  };

  bool operator <(const gauss_pt_corresp &gpc1,
                  const gauss_pt_corresp &gpc2);

  struct ga_instruction_set {

    papprox_integration pai;       // Current approximation method
    fem_interpolation_context ctx; // Current fem interpolation context.
    base_small_vector Normal;      // Outward unit normal vector to the
                                   // boundary in case of boundary integration
    scalar_type elt_size;          // Estimate of the diameter of the element
                                   // if needed.
    bool need_elt_size;
    scalar_type coeff;             // Coefficient for the Gauss point
    size_type nbpt, ipt;           // Number and index of Gauss point
    bgeot::geotrans_precomp_pool gp_pool;
    fem_precomp_pool fp_pool;
    std::map<gauss_pt_corresp, bgeot::pstored_point_tab> neighbour_corresp;

    using region_mim_tuple = std::tuple<const mesh_im *, const mesh_region *, psecondary_domain>;
    struct region_mim : public region_mim_tuple {
      const mesh_im* mim() const {return std::get<0>(static_cast<region_mim_tuple>(*this));}
      const mesh_region* region() const {return std::get<1>(static_cast<region_mim_tuple>(*this));}
      psecondary_domain psd() const {return std::get<2>(static_cast<region_mim_tuple>(*this));}

      region_mim(const mesh_im *mim_, const mesh_region *region_, psecondary_domain psd)
        : region_mim_tuple(mim_, region_, psd) {}
    };

    std::map<std::string, const base_vector *> extended_vars;
    std::map<std::string, base_vector> really_extended_vars;
    std::map<std::string, gmm::sub_interval> var_intervals;
    size_type nb_dof, max_dof;

    struct variable_group_info {
      const mesh_fem *mf;
      gmm::sub_interval Ir, In;
      scalar_type alpha;
      const base_vector *U;
      const std::string *varname;
      variable_group_info() : mf(0), U(0), varname(0) {}
    };

    struct interpolate_info {
      size_type pt_type;
      bool has_ctx;
      const mesh *m;
      fem_interpolation_context ctx;
      base_node pt_y;
      base_small_vector Normal;
      base_matrix G;
      std::map<std::string, variable_group_info> groups_info;
      std::map<var_trans_pair, base_tensor> derivatives;
      std::map<const mesh_fem *, pfem_precomp> pfps;
    };


    struct secondary_domain_info {
      // const mesh *m;
      papprox_integration pai;
      fem_interpolation_context ctx;
      base_small_vector Normal;
      
      std::map<std::string, base_vector> local_dofs;
      std::map<const mesh_fem *, pfem_precomp> pfps;
      std::map<const mesh_fem *, base_tensor> base;
      std::map<const mesh_fem *, base_tensor> grad;
      std::map<const mesh_fem *, base_tensor> hess;
    };


    struct elementary_trans_info {
      base_matrix M;
      const mesh_fem *mf;
      size_type icv;
    };

    std::set<std::string> transformations;

    struct region_mim_instructions {

      const mesh *m;
      const mesh_im *im;
      ga_if_hierarchy current_hierarchy;
      std::map<std::string, base_vector> local_dofs;
      std::map<const mesh_fem *, pfem_precomp> pfps;
      std::map<const mesh_fem *, base_tensor> base;
      std::map<const mesh_fem *, std::list<ga_if_hierarchy>> base_hierarchy;
      std::map<const mesh_fem *, base_tensor> grad;
      std::map<const mesh_fem *, std::list<ga_if_hierarchy>> grad_hierarchy;
      std::map<const mesh_fem *, base_tensor> hess;
      std::map<const mesh_fem *, std::list<ga_if_hierarchy>> hess_hierarchy;

      std::map<const mesh_fem *, base_tensor>
        xfem_plus_base,  xfem_plus_grad,  xfem_plus_hess,
        xfem_minus_base, xfem_minus_grad, xfem_minus_hess;
      std::map<const mesh_fem *, std::list<ga_if_hierarchy>>
        xfem_plus_base_hierarchy,  xfem_plus_grad_hierarchy,
        xfem_plus_hess_hierarchy,  xfem_minus_base_hierarchy,
        xfem_minus_grad_hierarchy, xfem_minus_hess_hierarchy;

      std::map<std::string, std::set<std::string>> transformations;
      std::set<std::string> transformations_der;
      std::map<std::string, interpolate_info> interpolate_infos;
      std::map<std::string, elementary_trans_info> elementary_trans_infos;
      secondary_domain_info secondary_domain_infos;

      // Instructions being executed at the first Gauss point after
      // a change of integration method only.
      ga_instruction_list begin_instructions;
      // Instructions executed once per element
      ga_instruction_list elt_instructions;
      // Instructions executed on each integration/interpolation point
      ga_instruction_list instructions;
      std::map<scalar_type, std::list<pga_tree_node> > node_list;

    region_mim_instructions(): m(0), im(0) {}
    };

    std::list<ga_tree> trees; // The trees are stored mainly because they
                              // contain the intermediary tensors.
    std::list<ga_tree> interpolation_trees;

    typedef std::map<region_mim, region_mim_instructions> instructions_set;

    instructions_set  whole_instructions;

    ga_instruction_set() { max_dof = nb_dof = 0; need_elt_size = false; ipt=0; }
  };

  
  void ga_exec(ga_instruction_set &gis, ga_workspace &workspace);
  void ga_function_exec(ga_instruction_set &gis);
  void ga_compile(ga_workspace &workspace, ga_instruction_set &gis,
                         size_type order);
  void ga_compile_function(ga_workspace &workspace,
                                  ga_instruction_set &gis, bool scalar);
  void ga_compile_interpolation(ga_workspace &workspace,
				ga_instruction_set &gis);
  void ga_interpolation_exec(ga_instruction_set &gis,
			     ga_workspace &workspace,
			     ga_interpolation_context &gic);
  void ga_interpolation_single_point_exec
    (ga_instruction_set &gis, ga_workspace &workspace,
     const fem_interpolation_context &ctx_x, const base_small_vector &Normal,
     const mesh &interp_mesh);
  
} /* end of namespace */


#endif /* GETFEM_GENERIC_ASSEMBLY_COMP_EXEC_H__  */
