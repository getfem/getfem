/*===========================================================================

 Copyright (C) 2026-2026 Yves Renard.

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program.  If not, see https://www.gnu.org/licenses/.

===========================================================================*/

/* Test of the element_curvature operator of the generic assembly
   language on high order meshes of a portion of sphere, of cylinder
   and of circle. */

#include "getfem/getfem_generic_assembly.h"
#include "getfem/getfem_mesh_im.h"
#include "getfem/getfem_mesh_fem.h"
#include "getfem/getfem_interpolation.h"

using bgeot::base_node;
using bgeot::scalar_type;
using bgeot::size_type;
using bgeot::short_type;
using std::cout;
using std::endl;

// Builds a consistently oriented mesh of dimension N of an hypersurface
// parametrized by f on [0,1]^(N-1) (NX subdivisions in each direction,
// simplices of degree K).
template <typename F>
static void build_mapped_mesh(getfem::mesh &m, size_type N, size_type K,
                              size_type NX, F f) {
  size_type P = N-1;
  bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(P, short_type(K));
  const bgeot::stored_point_tab &refpts = pgt->geometric_nodes();
  scalar_type h = scalar_type(1) / scalar_type(NX);
  std::vector<std::vector<base_node> > simplices;
  for (size_type i = 0; i < NX; ++i) {
    if (P == 1) {
      base_node p0(1), p1(1); p0[0] = i*h; p1[0] = (i+1)*h;
      simplices.push_back({p0, p1});
    } else {
      for (size_type j = 0; j < NX; ++j) {
        base_node p00(i*h, j*h), p10((i+1)*h, j*h);
        base_node p01(i*h, (j+1)*h), p11((i+1)*h, (j+1)*h);
        simplices.push_back({p00, p10, p01});
        simplices.push_back({p10, p11, p01});
      }
    }
  }
  for (const auto &s : simplices) {
    std::vector<base_node> pts;
    for (const base_node &xi : refpts) {
      base_node p = s[0];
      for (size_type k = 0; k < P; ++k)
        gmm::add(gmm::scaled(s[k+1] - s[0], xi[k]), p);
      pts.push_back(f(p));
    }
    m.add_convex_by_points(pgt, pts.begin());
  }
}

static scalar_type integrate(const getfem::mesh_im &mim,
                             const std::string &expr) {
  getfem::ga_workspace workspace;
  workspace.add_expression(expr, mim);
  workspace.assembly(0);
  return workspace.assembled_potential();
}

static void check(const std::string &name, size_type N, size_type K,
                  const getfem::mesh &m, const std::string &errexpr,
                  scalar_type tol) {
  getfem::mesh_im mim(m);
  mim.set_integration_method(getfem::int_method_descriptor
                             (N == 3 ? "IM_TRIANGLE(10)" : "IM_GAUSS1D(10)"));
  scalar_type area = integrate(mim, "1");
  scalar_type err = sqrt(integrate(mim, errexpr) / area);
  cout << name << " K = " << K << ", relative L2 error = " << err << endl;
  GMM_ASSERT1(err < tol, "Error too large for " << name << ": " << err);
}

int main(int argc, char * argv[]) {

  GETFEM_MPI_INIT(argc, argv); // For parallelized version

  scalar_type R = 2.;

  for (size_type K = 2; K <= 4; ++K) {
    scalar_type tol = (K == 2) ? 2E-1 : 2E-2;

    // Portion of sphere (theta, phi) in [0.5,2]x[0,1.5], outward normal
    // with this parametrization : grad_s n = P/R
    getfem::mesh m1;
    build_mapped_mesh(m1, 3, K, 8, [R](const base_node &p) {
        scalar_type th = 0.5 + 1.5*p[0], ph = 1.5*p[1];
        return base_node(R*sin(th)*cos(ph), R*sin(th)*sin(ph), R*cos(th));
      });
    std::string Rs = std::to_string(R), L = "element_curvature";
    std::string Pr = "(Id(3)-X@X/sqr("+Rs+"))";
    check("Sphere, curvature tensor", 3, K, m1,
          "Norm_sqr("+L+" - "+Pr+"/"+Rs+")*sqr("+Rs+")", tol);
    check("Sphere, mixed components in the local basis", 3, K, m1,
          "Norm_sqr(element_B'*"+L+"*element_K - Id(2)/"+Rs+")*sqr("+Rs+")",
          tol);
    check("Sphere, element_B'*element_K", 3, K, m1,
          "Norm_sqr(element_B'*element_K - Id(2))", 1E-10);
    check("Sphere, mean curvature", 3, K, m1,
          "sqr(Trace("+L+")/2 - 1/"+Rs+")*sqr("+Rs+")", tol);
    check("Sphere, square of curvature tensor", 3, K, m1,
          "Norm_sqr("+L+"*"+L+" - "+Pr+"/sqr("+Rs+"))*sqr(sqr("+Rs+"))", tol);
    check("Sphere, Gauss curvature", 3, K, m1,
          "sqr((sqr(Trace("+L+"))-Trace("+L+"*"+L+"))/2 - 1/sqr("+Rs+"))"
          "*sqr(sqr("+Rs+"))", tol);

    // Portion of cylinder of axis z
    getfem::mesh m2;
    build_mapped_mesh(m2, 3, K, 8, [R](const base_node &p) {
        scalar_type ph = 2.*p[0];
        return base_node(R*cos(ph), R*sin(ph), 3.*p[1]);
      });
    std::string Q = "([X(1);X(2);0]@[X(1);X(2);0])/sqr("+Rs+")";
    check("Cylinder, curvature tensor", 3, K, m2,
          "Norm_sqr("+L+" - (Id(3)-[0;0;1]@[0;0;1]-"+Q+")/"+Rs+")"
          "*sqr("+Rs+")", tol);

    // Arc of circle in 2D
    getfem::mesh m3;
    build_mapped_mesh(m3, 2, K, 8, [R](const base_node &p) {
        scalar_type ph = 3.*p[0];
        return base_node(R*cos(ph), R*sin(ph));
      });
    check("Circle, curvature tensor", 2, K, m3,
          "Norm_sqr("+L+" - (Id(2)-X@X/sqr("+Rs+"))/"+Rs+")*sqr("+Rs+")",
          tol);
  }

  // Linear transformation : curvature should vanish
  getfem::mesh m4;
  build_mapped_mesh(m4, 3, 1, 4, [](const base_node &p) {
      return base_node(p[0], p[1], p[0]+2.*p[1]);
    });
  check("Plane, curvature tensor", 3, 1, m4, "Norm_sqr(element_curvature)",
        1E-12);

  GETFEM_MPI_FINALIZE;

  return 0;
}
