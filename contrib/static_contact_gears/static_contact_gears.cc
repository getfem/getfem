// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2009-2009 Konstantinos Poulios.
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

#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_import.h"   /* import functions (load a mesh from file)    */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_Coulomb_friction.h"
#include "gmm/gmm.h"
#include <map>

/* some Getfem++ types that we will be using */
using bgeot::dim_type;
using bgeot::size_type;   /* = unsigned long */
using bgeot::scalar_type; /* = double */
using bgeot::base_node;   /* geometrical nodes(derived from base_small_vector)*/
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::model_real_plain_vector  plain_vector;
typedef getfem::model_real_sparse_matrix sparse_matrix;

struct contact_node {
  size_type dof;              // dof id of the node in mf_rhs
  size_type master;           // cn id of its master node
  std::vector<size_type> cvs; // list of ids of neigbouring convexes
  std::vector<size_type> fcs; // list of local ids of neigbouring faces
  scalar_type dist;           // Distance from master point
  bool is_active;             // 
  
  contact_node() {
    is_active = false;
    dist = 10.;               // This initial value represents a threshold
  }
};


/*
  structure for the elastostatic problem
*/
struct elastostatic_problem {

  enum { DIRICHLET_BOUNDARY = 0,
         DIRICHLET_BOUNDARY_1 = 1, DIRICHLET_BOUNDARY_2 = 2,
         CONTACT_BOUNDARY_1 = 3, CONTACT_BOUNDARY_2 = 4,
         BLOCK_1 = 5, BLOCK_2 = 6};
  getfem::mesh mesh;         /* the mesh */
  getfem::mesh_im mim;       /* the integration methods */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the elastostatic solution */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  scalar_type lambda, mu;    /* elastic coefficients.                        */

  scalar_type residual;      /* max residual for the iterative solvers       */
  scalar_type rot_angle;     /* rotation angle of the pinion gear            */
  scalar_type threshold;     /* threshold distance for contact finding       */

  size_type N;               /* dimension of the problem                     */

  bool solve(plain_vector &, plain_vector &, plain_vector &, plain_vector &);
  void init(void);
  elastostatic_problem(void) : mim(mesh), mf_u(mesh), mf_rhs(mesh) {}

};


/* Define the problem parameters, import the mesh, set finite element
   and integration methods and detect the boundaries.
*/
void elastostatic_problem::init(void) {

  std::string FEM_TYPE  = "FEM_QK(3,1)";
  std::string INTEGRATION = "IM_HEXAHEDRON(5)";

  /* First step : import the mesh */
  for (size_type ii = 0; ii <= 1; ii++) {
    getfem::mesh tmpmesh;
    getfem::mesh_region block, contact_boundary, dirichlet_boundary;
    if (ii == 0) {
      block = mesh.region(BLOCK_1);
      contact_boundary = mesh.region(CONTACT_BOUNDARY_1);
      dirichlet_boundary = mesh.region(DIRICHLET_BOUNDARY_1);
      getfem::import_mesh("gmsh:./gear1.msh",tmpmesh);
    } else {
      block = mesh.region(BLOCK_2);
      contact_boundary = mesh.region(CONTACT_BOUNDARY_2);
      dirichlet_boundary = mesh.region(DIRICHLET_BOUNDARY_2);
      getfem::import_mesh("gmsh:./gear2.msh",tmpmesh);
//bgeot::base_small_vector T(N);
//T[2] = 0.1;
//tmpmesh.translation(T);
    }
    for (dal::bv_visitor cv(tmpmesh.convex_index()); !cv.finished(); ++cv) {
      size_type newcv =
        mesh.add_convex_by_points(tmpmesh.trans_of_convex(cv),
                                  tmpmesh.points_of_convex(cv).begin());
      block.add(newcv);
      for (size_type f = 0; f < tmpmesh.structure_of_convex(cv)->nb_faces(); f++) {
        if (tmpmesh.region(113).is_in(cv,f)) {
          contact_boundary.add(newcv,f);
        } else if (tmpmesh.region(133).is_in(cv,f) ||
                   tmpmesh.region(142).is_in(cv,f) ||
                   tmpmesh.region(143).is_in(cv,f) ||
                   tmpmesh.region(173).is_in(cv,f) ||
                   tmpmesh.region(182).is_in(cv,f) ||
                   tmpmesh.region(183).is_in(cv,f)  ) {
          dirichlet_boundary.add(newcv,f);
//          mesh.region(DIRICHLET_BOUNDARY).add(newcv,f);
        }
      }
    }
  }

  N = mesh.dim();

  residual = 1e-6;
  rot_angle = -1.5e-2;
  threshold = 10.;

  lambda = 1.18e+5;
  mu = 0.83e+5;

  mf_u.set_qdim(dim_type(N));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  mf_u.set_finite_element(pf_u);

  getfem::pintegration_method ppi = getfem::int_method_descriptor(INTEGRATION);
  mim.set_integration_method(ppi);

  /* set the finite element on mf_rhs */
  std::string data_fem_name = "FEM_QK(3,1)";
  mf_rhs.set_finite_element(mesh.convex_index(), 
                            getfem::fem_descriptor(data_fem_name));

}

/*  Construction and solution of the Model.
*/
bool elastostatic_problem::solve(plain_vector &U, plain_vector &RHS,
                                 plain_vector &Forces, plain_vector &CForces) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  dim_type qdim = mf_u.get_qdim();

  cout << "mf_rhs.nb_dof() :" << mf_rhs.nb_dof() << endl;
  cout << "mf_u.nb_dof()   :" << mf_u.nb_dof() << endl;

  getfem::model md;
  md.add_fem_variable("u", mf_u);

  // Linearized elasticity brick.
//  getfem::mdbrick_isotropic_linearized_elasticity<> ELAS(mim, mf_u, lambda,mu);
  md.add_initialized_scalar_data("lambda", lambda);
  md.add_initialized_scalar_data("mu", mu);
  getfem::add_isotropic_linearized_elasticity_brick(md, mim, "u", "lambda", "mu");

  // Nonlinear elasticity brick.
//  base_vector p(2); p[0] = lambda; p[1] = mu;
//  getfem::SaintVenant_Kirchhoff_hyperelastic_law pl;
//  getfem::mdbrick_nonlinear_elasticity<>  ELAS(pl, mim, mf_u, p);
  
  // Defining the contact condition.
  size_type no_cn = 0;
  sparse_matrix BN;
  plain_vector gap;

  std::map<int,int> dof_to_cnid;
  std::vector<contact_node> cns;
  dal::bit_vector cn1, cn2;
  size_type cns_size=0;
  for (size_type swap = 0; swap <= 1; swap++) {
    size_type CONTACT_BOUNDARY = swap ? CONTACT_BOUNDARY_2 : CONTACT_BOUNDARY_1;
    for (getfem::mr_visitor face(mesh.region(CONTACT_BOUNDARY));
         !face.finished(); ++face) {
      assert(face.is_face());
      getfem::mesh_fem::ind_dof_face_ct
        face_dofs = mf_rhs.ind_basic_dof_of_face_of_element(face.cv(),face.f());
      for (size_type it=0; it < face_dofs.size(); ++it) {
        size_type cnid;
        if (dof_to_cnid.count(face_dofs[it]) > 0) {
          cnid = dof_to_cnid[face_dofs[it]];
          cns[cnid].cvs.push_back(face.cv());
          cns[cnid].fcs.push_back(face.f());
        } else {
          contact_node new_cn;
          new_cn.dof = face_dofs[it];
          new_cn.cvs.push_back(face.cv());
          new_cn.fcs.push_back(face.f());
          cns.push_back(new_cn);
          dof_to_cnid[face_dofs[it]] = cns_size;
          cnid = cns_size;
          cns_size++;
        }
        if (CONTACT_BOUNDARY == CONTACT_BOUNDARY_1) {
          cn1.add(cnid);
        } else if (CONTACT_BOUNDARY == CONTACT_BOUNDARY_2) {
          cn2.add(cnid);
        }
      }
    }
  }

  for (dal::bv_visitor i1(cn1); !i1.finished(); ++i1) { //find minimum distance node pairs
    base_node node1 = mf_rhs.point_of_basic_dof(cns[i1].dof);
    for (dal::bv_visitor i2(cn2); !i2.finished(); ++i2) {
      base_node node2 = mf_rhs.point_of_basic_dof(cns[i2].dof);
      scalar_type dist = gmm::vect_norm2(node1-node2);
      if (dist < cns[i1].dist) {
        cns[i1].dist = dist;
        cns[i1].master = i2;
        cns[i1].is_active = true;
      }
      if (dist < cns[i2].dist) {
        cns[i2].dist = dist;
        cns[i2].master = i1;
        cns[i2].is_active = true;
      }
    }
  }

  for (std::vector<contact_node>::iterator cn = cns.begin(); cn < cns.end(); cn++) {
    if (cn->is_active) {
      contact_node cn_s = *cn;              //slave contact node
      contact_node cn_m = cns[cn->master];  //master contact node
      base_node slave_node = mf_rhs.point_of_basic_dof(cn_s.dof);
      base_node master_node = mf_rhs.point_of_basic_dof(cn_m.dof);
      base_node un_sel(3), proj_node_sel(3), proj_node_ref_sel(3);
      scalar_type is_in_min = 1e5;
      size_type cv_sel, fc_sel;
      std::vector<size_type>::iterator cv,fc;
      for (cv = cn_m.cvs.begin(), fc = cn_m.fcs.begin();
           cv < cn_m.cvs.end() && fc < cn_m.fcs.end(); cv++, fc++) {
        base_node un = mesh.normal_of_face_of_convex(*cv,*fc);
        un /= gmm::vect_norm2(un);
        base_node proj_node(3), proj_node_ref(3);
        //proj_node = slave_node - [(slave_node-master_node)*n] * n
        gmm::add(master_node, gmm::scaled(slave_node, -1.), proj_node);
        gmm::copy(gmm::scaled(un, gmm::vect_sp(proj_node, un)), proj_node);
        gmm::add(slave_node, proj_node);

        bgeot::pgeometric_trans pgt = mesh.trans_of_convex(*cv);
        bgeot::geotrans_inv_convex gic;
        gic.init(mesh.points_of_convex(*cv), pgt);
        gic.invert(proj_node, proj_node_ref);
        scalar_type is_in = pgt->convex_ref()->is_in(proj_node_ref);
        if (is_in < is_in_min) {
          is_in_min = is_in;
          cv_sel = *cv;
          fc_sel = *fc;
          un_sel = un;
          proj_node_sel = proj_node;
          proj_node_ref_sel = proj_node_ref;
        }
      }
      if (is_in_min < 0.05) {
        no_cn++;
        gmm::resize(BN, no_cn, mf_u.nb_dof());
        BN(no_cn-1, cn_s.dof*N)   += -un_sel[0];
        BN(no_cn-1, cn_s.dof*N+1) += -un_sel[1];
        BN(no_cn-1, cn_s.dof*N+2) += -un_sel[2];

        base_matrix G;
        base_matrix M(qdim, mf_u.nb_basic_dof_of_element(cv_sel));
        bgeot::vectors_to_base_matrix(G, mesh.points_of_convex(cv_sel));
        getfem::pfem pf = mf_u.fem_of_element(cv_sel);
        bgeot::pgeometric_trans pgt = mesh.trans_of_convex(cv_sel);
        getfem::fem_interpolation_context
          ctx(pgt, pf, proj_node_ref_sel, G, cv_sel, fc_sel);
        pf->interpolation (ctx, M, qdim);

        plain_vector tmpvec(mf_u.nb_basic_dof_of_element(cv_sel));
        gmm::mult(gmm::transposed(M), un_sel, tmpvec);
        size_type j = 0;
        getfem::mesh_fem::ind_dof_ct master_dof = mf_u.ind_basic_dof_of_element(cv_sel);
        getfem::mesh_fem::ind_dof_ct::const_iterator it ;
        for (it = master_dof.begin(); it != master_dof.end(); it++) {
          BN(no_cn-1, *it) += tmpvec[j];
          ++j;
        }
        gap.push_back(gmm::vect_sp(slave_node-proj_node_sel, un_sel));
      }
    }
  }
//  getfem::mdbrick_Coulomb_friction<> FRICTION(ELAS, BN, gap);
//  FRICTION.set_r(5e1);
  md.add_fixed_size_variable("contact_multiplier", no_cn);
  md.add_initialized_scalar_data("r", 1.);
  md.add_initialized_fixed_size_data("gap", gap);
  getfem::add_basic_contact_brick(md, "u", "contact_multiplier", "r",
                                  BN, "gap");

cout << "no_cn: " << no_cn << endl;

  // Defining the DIRICHLET condition.
  plain_vector F(nb_dof_rhs * N);
  dal::bit_vector
    cn = mf_rhs.basic_dof_on_region(mesh.region(DIRICHLET_BOUNDARY_1));
  for (dal::bv_visitor i(cn); !i.finished(); ++i) {
    base_node node = mf_rhs.point_of_basic_dof(i);
    F[i*N] = -node[1] * rot_angle;
    F[i*N+1] = node[0] * rot_angle;
  }
//  getfem::mdbrick_Dirichlet<> DIRICHLET(FRICTION, DIRICHLET_BOUNDARY_1);
//  DIRICHLET.set_constraints_type(getfem::PENALIZED_CONSTRAINTS); //getfem::ELIMINATED_CONSTRAINTS); //
//  DIRICHLET.rhs().set(mf_rhs, F);
  md.add_initialized_fem_data("DirichletData1", mf_rhs, F);
//  getfem::add_Dirichlet_condition_with_penalization
//    (md, mim, "u", 1e6, DIRICHLET_BOUNDARY_1, "DirichletData1");
  getfem::add_Dirichlet_condition_with_multipliers
    (md, mim, "u", mf_u, DIRICHLET_BOUNDARY_1, "DirichletData1");

  gmm::clear(F);
//  getfem::mdbrick_Dirichlet<> FINAL_MODEL(DIRICHLET, DIRICHLET_BOUNDARY_2);
//  FINAL_MODEL.set_constraints_type(getfem::PENALIZED_CONSTRAINTS); //getfem::ELIMINATED_CONSTRAINTS); //
//  FINAL_MODEL.rhs().set(mf_rhs, F);
  md.add_initialized_fem_data("DirichletData2", mf_rhs, F);
//  getfem::add_Dirichlet_condition_with_penalization
//    (md, mim, "u", 1e6, DIRICHLET_BOUNDARY_2, "DirichletData2");
  getfem::add_Dirichlet_condition_with_multipliers
    (md, mim, "u", mf_u, DIRICHLET_BOUNDARY_2, "DirichletData2");

  // Defining the surface pressure term for the NEUMANN boundary.
//  base_vector f(N);
//  for (size_type i = 0; i < nb_dof_rhs; ++i) {
//    gmm::copy(f, gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
//  }
//  getfem::mdbrick_source_term<> NEUMANN(FRICTION, mf_rhs, F, NEUMANN_BOUNDARY_NUM);

//  getfem::standard_model_state MS(FINAL_MODEL);
  gmm::iteration iter(residual, 1, 40000);

//  getfem::standard_solve(MS, FINAL_MODEL, iter);
//  plinsolv.reset(new getfem::linear_solver_superlu<sparse_matrix,plain_vector>);
  gmm::default_newton_line_search ls(size_t(-1), 5.0/3.0,
                                     1.0/1000.0, 3.0/5.0, 1.6);
//  getfem::model_problem<MODEL_STATE> mdpb(MS, FINAL_MODEL, ls);
//  MS.adapt_sizes(FINAL_MODEL); // to be sure it is ok, but should be done before
//  getfem::classical_Newton(mdpb, iter, *plinsolv);
  getfem::standard_solve(md, iter, getfem::rselect_linear_solver(md,"superlu"), ls);
                      
  gmm::resize(U, mf_u.nb_dof());
  gmm::resize(RHS, md.nb_dof());
  gmm::resize(CForces, mf_u.nb_dof());

//  gmm::copy(ELAS.get_solution(MS), U);
  gmm::copy(md.real_variable("u"), U);
//  gmm::copy(MS.residual(), RHS);
  gmm::copy(md.real_rhs(), RHS);

  gmm::copy(gmm::sub_vector(RHS, md.interval_of_variable("u")), Forces);
  gmm::scale(Forces, -1.0);

//  gmm::mult(gmm::sub_matrix(MS.tangent_matrix(),
//                            gmm::sub_interval(0, mf_u.nb_dof()),
//                            gmm::sub_interval(mf_u.nb_dof(), no_cn) ),
//            gmm::sub_vector(MS.state(), gmm::sub_interval(mf_u.nb_dof(), no_cn)),
  gmm::mult(gmm::sub_matrix(md.real_tangent_matrix(),
                            md.interval_of_variable("u"),
                            md.interval_of_variable("contact_multiplier") ),
            md.real_variable("contact_multiplier"), 
            CForces);
  gmm::scale(CForces, -1.0);

  return (iter.converged());
}

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  elastostatic_problem p;
  p.init();

  plain_vector U(p.mf_u.nb_dof());
  plain_vector RHS(p.mf_u.nb_dof());
  plain_vector Forces(p.mf_u.nb_dof());
  plain_vector CForces(p.mf_u.nb_dof());
  if (!p.solve(U, RHS, Forces, CForces)) cout << "Solve has failed\n";

  p.mesh.write_to_file("simple_gears.mesh");
  p.mf_u.write_to_file("simple_gears.mf", true);
  p.mf_rhs.write_to_file("simple_gears.mfd", true);
  gmm::vecsave("simple_gears.U", U);
  gmm::vecsave("simple_gears.RHS", RHS);
  getfem::vtk_export exp("simple_gears.vtk", true);
  exp.exporting(p.mf_u);
  exp.write_point_data(p.mf_u, U, "elastostatic_displacement");
  exp.write_point_data(p.mf_u, Forces, "forces");
  exp.write_point_data(p.mf_u, CForces, "contact_forces");

  return 0; 
}

