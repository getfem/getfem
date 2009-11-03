// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2008 Yves Renard, Julien Pommier.
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

#include <getfemint_misc.h>
#include <getfem/getfem_derivatives.h>
#include <getfem/getfem_interpolation.h>
#include <getfem/getfem_assembling.h>
#include <getfemint_mesh_slice.h>
#include <getfem/getfem_error_estimate.h>
#include <getfem/getfem_convect.h>
using namespace getfemint;

static void
error_for_non_lagrange_elements(const getfem::mesh_fem &mf,
				bool warning_only = false) {
  size_type cnt=0, total=0, cnt_no_fem=0;
  for (dal::bv_visitor cv(mf.linked_mesh().convex_index());
       !cv.finished(); ++cv) {
    if (!mf.convex_index()[cv]) cnt_no_fem++;
    else if (!mf.fem_of_element(cv)->is_lagrange()) cnt++;
    total++;
  }
  if (cnt) {
    if (!warning_only) {
      THROW_ERROR("Error: " << cnt << " elements on " << total << " are NOT "
		  "lagrange elements -- Unable to compute a derivative");
    } else {
      GFI_WARNING(cnt << " elements on " << total
		  << " are NOT lagrange elements");
    }
  }
  // Test suppressed. If ones want to interpolate on a specific region
  // for instance.
  // if (cnt_no_fem) {
  //     if (!warning_only) {
  //       THROW_ERROR("Error: " << cnt_no_fem << " elements on " << total << " have NO FEM!");
  //     } else {
  //       GFI_WARNING(cnt_no_fem << " elements on " << total << " have NO FEM");
  //     }
  //   }
}



static void
mesh_edges_deformation(const getfem::mesh_fem *mf, darray &U, unsigned N,
                      mexargs_in &in, mexargs_out &out) {
  unsigned mdim = mf->linked_mesh().dim();
  if (mf->get_qdim() != mdim) {
    THROW_BADARG( "Error, the mesh is of dimension " << mdim <<
		  " while its Qdim is " << mf->get_qdim());
  }
  bgeot::edge_list el;
  const getfem::mesh &m = mf->linked_mesh();

  build_edge_list(m, el, in);

  darray w   = out.pop().create_darray(mdim, N, unsigned(el.size()));

  bgeot::edge_list::const_iterator it = el.begin();
  size_type nbe = 0;
  while (it != el.end()) {
    bgeot::edge_list::const_iterator nit = it;
    size_type ecnt = 0;

    /* count consecutives edges of one convex */
    while (nit != el.end() && (*nit).cv == (*it).cv) {
      ecnt++; nit++;
    }
    size_type cv = (*it).cv;
    check_cv_fem(*mf, cv);

    bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);

    std::vector<getfem::base_node> pt(ecnt * N);

    /* for every edge of the convex, push the points of its refined edge
       on the vector 'pt' */
    size_type pcnt = 0;
    for (bgeot::edge_list::const_iterator eit = it; eit != nit; eit++) {
      /* build the list of points on the edge, on the reference element */
      /* get local point numbers in the convex */
      bgeot::size_type iA = m.local_ind_of_convex_point(cv, (*eit).i);
      bgeot::size_type iB = m.local_ind_of_convex_point(cv, (*eit).j);

      getfem::base_node A = pgt->convex_ref()->points()[iA];
      getfem::base_node B = pgt->convex_ref()->points()[iB];
      for (size_type i = 0; i < N; i++) {
	pt[pcnt++] = A +  (B-A)*(double(i)/double(N-1));
      }
    }
    if (pcnt != ecnt * N) THROW_INTERNAL_ERROR;

    /* now, evaluate the field U on every point of pt  once */
    getfem::base_matrix pt_val;
    interpolate_on_convex_ref(mf, cv, pt, U, pt_val);

    if (pt_val.size() != ecnt * N * mdim) THROW_INTERNAL_ERROR;

    /* evaluate the point location on the real mesh, adds it
       the 'deformation' field pt_val interpolated from U,
       and write the result in the destination vector */
    for (ecnt = 0; it != nit; it++, ecnt++, nbe++) {
      for (size_type i = 0; i < N; i++) {
	size_type ii = ecnt*N+i;
	getfem::base_node def_pt = pgt->transform(pt[ii],
						  m.points_of_convex(cv));
	for (size_type k = 0; k < mdim; k++) {
	  def_pt[k] += pt_val(k,ii);
	}
	std::copy(def_pt.begin(), def_pt.end(), &w(0, i, nbe));
      }
    }
  }
}

template <typename T> static void
gf_compute_gradient(getfemint::mexargs_out& out,
		    const getfem::mesh_fem& mf,
		    const getfem::mesh_fem& mf_grad,
		    const garray<T> &U,
		    size_type qm) {
  garray<T> DU;
  unsigned N = mf.linked_mesh().dim();
  array_dimensions dims(N);
  unsigned qqdim = unsigned(dims.push_back(U,0,U.ndim()-1,true));

  if (qm != 1) dims.push_back(unsigned(qm));
  dims.push_back(unsigned(mf_grad.nb_dof()));
  DU = out.pop().create_array(dims, T());
  std::vector<T> tmp(mf_grad.nb_dof() * qm * N);
  for (unsigned qq=0; qq < qqdim; ++qq) {
    // compute_gradient also checks that the meshes are the same
    getfem::compute_gradient(mf, mf_grad, gmm::sub_vector(U, gmm::sub_slice(qq, mf.nb_dof(),qqdim)), tmp);
    for (unsigned j=0, pos=qq*N; j < tmp.size(); j+=N) {
      for (unsigned k=0; k < N; ++k) DU[pos+k] = tmp[j+k];
      pos += qqdim*N;
    }
  }
}



template <typename T> static void
gf_compute_hessian(getfemint::mexargs_out& out,
		   const getfem::mesh_fem& mf,
		   const getfem::mesh_fem& mf_hess,
		   const garray<T> &U,
		   size_type qm) {
  garray<T> D2U;
  unsigned N = mf.linked_mesh().dim();
  array_dimensions dims(N); dims.push_back(N);
  unsigned qqdim = unsigned(dims.push_back(U,0,U.ndim()-1,true));

  if (qm != 1) dims.push_back(unsigned(qm));
  dims.push_back(unsigned(mf_hess.nb_dof()));
  D2U = out.pop().create_array(dims, T());
  std::vector<T> tmp(mf_hess.nb_dof() * qm * N * N);
  for (unsigned qq=0; qq < qqdim; ++qq) {
    // compute_gradient also checks that the meshes are the same
    getfem::compute_hessian(mf, mf_hess, gmm::sub_vector(U, gmm::sub_slice(qq, mf.nb_dof(),qqdim)), tmp);
    //cerr << "tmp = " << tmp << "\n";
    for (unsigned j=0, pos=qq*N*N; j < tmp.size(); j+=N*N) {
      for (unsigned k=0; k < N*N; ++k) D2U[pos+k] = tmp[j+k];
      /* TODO A VERFIEIER !*/
      pos += qqdim*N*N;
    }
  }
}


template <typename T> static void
gf_interpolate(getfemint::mexargs_in& in, getfemint::mexargs_out& out,
	       const getfem::mesh_fem& mf, const garray<T> &U) {
  array_dimensions dims;
  dims.push_back(U,0,U.ndim()-1,true);
  if (in.front().is_mesh_fem()) {
    const getfem::mesh_fem& mf_dest = *in.pop().to_const_mesh_fem();
    error_for_non_lagrange_elements(mf_dest, true);
    size_type qmult = mf.get_qdim() / mf_dest.get_qdim();
    if (qmult == 0)
      THROW_ERROR("Cannot interpolate a mesh_fem with qdim = " <<
		  int(mf.get_qdim()) << " onto a mesh_fem whose qdim is "
		  << int(mf_dest.get_qdim()));
    if (qmult != 1) dims.push_back(unsigned(qmult));
    dims.push_back(unsigned(mf_dest.nb_dof()));
    dims.opt_transform_col_vect_into_row_vect();
    garray<T> V = out.pop().create_array(dims,T());
    getfem::interpolation(mf, mf_dest, U, V);
  }
  else if (in.front().is_mesh_slice()) {
    getfem::stored_mesh_slice *sl =
      &in.pop().to_getfemint_mesh_slice()->mesh_slice();

    for (size_type i=0; i < sl->nb_convex(); ++i)
      if (!mf.linked_mesh().convex_index().is_in(sl->convex_num(i)))
      THROW_BADARG("the slice is not compatible with the mesh_fem "
		   "(cannot find convex " << sl->convex_num(i) << ")");

    if (mf.get_qdim() != 1) dims.push_back(mf.get_qdim());
    dims.push_back(unsigned(sl->nb_points()));
    dims.opt_transform_col_vect_into_row_vect();
    garray<T> V = out.pop().create_array(dims, T());
    sl->interpolate(mf, U, V);
  }
  else THROW_BADARG("expecting a mesh_fem or a mesh_slice for interpolation");
}

bool U_is_a_vector(const rcarray &U, const std::string& cmd) {
  if (U.sizes().size() == U.sizes().dim(-1)) return true;
  else THROW_BADARG("the U argument for the function " << cmd
		    << " must be a one-dimensional array");
  return false;
}

/*MLABCOM
  FUNCTION [x] = gf_compute(meshfem MF, vec U, operation [, args])

  Various computations involving the solution U of the finite element problem.

  @FUNC ::COMPUTE('L2 norm')
  @FUNC ::COMPUTE('H1 semi norm')
  @FUNC ::COMPUTE('H1 norm')
  @FUNC ::COMPUTE('H2 semi norm')
  @FUNC ::COMPUTE('H2 norm')
  @FUNC ::COMPUTE('gradient')
  @FUNC ::COMPUTE('hessian')
  @FUNC ::COMPUTE('interpolate on')
  @FUNC ::COMPUTE('extrapolate on')
  @FUNC ::COMPUTE('error estimate')
  @FUNC ::COMPUTE('convect')

  * [U2[,MF2,[,X[,Y[,Z]]]]] = gf_compute(MF,U,'interpolate on Q1 grid',
                               {'regular h', hxyz | 'regular N',Nxyz |
           			   X[,Y[,Z]]}

  Creates a cartesian Q1 mesh fem and interpolates U on it. The
  returned field U2 is organized in a matrix such that in can be drawn
  via the MATLAB command 'pcolor'. The first dimension is the Qdim of
  MF (i.e.  1 if U is a scalar field)

  example (mf_u is a 2D mesh_fem):
   >> Uq=gf_compute(mf_u, U, 'interpolate on Q1 grid', 'regular h', [.05, .05]);
   >> pcolor(squeeze(Uq(1,:,:)));

  * E = gf_compute(MF, U, 'mesh edges deformation', N [,vec or
                   mat CVLIST])
  [OBSOLETE FUNCTION! will be removed in a future release]
  Evaluates the deformation of the mesh caused by the field U (for a
  2D mesh, U must be a [2 x nb_dof] matrix). N is the refinment level
  (N>=2) of the edges.  CVLIST can be used to restrict the computation
  to the edges of the listed convexes ( if CVLIST is a row vector ),
  or to restrict the computations to certain faces of certain convexes
  when CVLIST is a two-rows matrix, the first row containing convex
  numbers and the second face numbers.

  * UP = gf_compute(MF, U, 'eval on triangulated surface', int Nrefine,
                    [vec CVLIST])
  [OBSOLETE FUNCTION! will be removed in a future release]
  Utility function designed for 2D triangular meshes : returns a list
  of triangles coordinates with interpolated U values. This can be
  used for the accurate visualization of data defined on a
  discontinous high order element. On output, the six first rows of UP
  contains the triangle coordinates, and the others rows contain the
  interpolated values of U (one for each triangle vertex) CVLIST may
  indicate the list of convex number that should be consider, if not
  used then all the mesh convexes will be used. U should be a row
  vector.

  $Id$
MLABCOM*/
/*MLABEXT
  if (nargin>=3 & strcmpi(varargin{3},'interpolate on Q1 grid')),
    [varargout{1:nargout}]=gf_compute_Q1grid_interp(varargin{[1 2 4:nargin]}); return;
  end;
  MLABEXT*/

void gf_compute(getfemint::mexargs_in& in, getfemint::mexargs_out& out) {

  if (in.narg() < 3)  THROW_BADARG( "Wrong number of input arguments");

  const getfem::mesh_fem *mf   = in.pop().to_const_mesh_fem();
  rcarray U              = in.pop().to_rcarray();
  in.last_popped().check_trailing_dimension(int(mf->nb_dof()));
  std::string cmd        = in.pop().to_string();

  if (check_cmd(cmd, "L2 norm", in, out, 1, 2, 0, 1) && U_is_a_vector(U,cmd)) {
    /*@FUNC n = ::COMPUTE('L2 norm',@tmim mim[, @mat CVids])
    Compute the L2 norm of the (real or complex) field `U`.

    If `CVids` is given, the norm will be computed only on the listed
    convexes.@*/
    const getfem::mesh_im *mim = in.pop().to_const_mesh_im();
    dal::bit_vector bv = in.remaining() ?
      in.pop().to_bit_vector(&mf->convex_index()) : mf->convex_index();
    if (!U.is_complex())
      out.pop().from_scalar(getfem::asm_L2_norm(*mim, *mf, U.real(), bv));
    else out.pop().from_scalar(getfem::asm_L2_norm(*mim, *mf, U.cplx(), bv));
  } else if (check_cmd(cmd, "H1 semi norm", in, out, 1, 2, 0, 1) && U_is_a_vector(U,cmd)) {
    /*@FUNC n = ::COMPUTE('H1 semi norm',@tmim mim[, @mat CVids])
    Compute the L2 norm of grad(`U`).

    If `CVids` is given, the norm will be computed only on the listed
    convexes.@*/
    const getfem::mesh_im *mim = in.pop().to_const_mesh_im();
    dal::bit_vector bv = in.remaining() ?
      in.pop().to_bit_vector(&mf->convex_index()) : mf->convex_index();
    if (!U.is_complex())
      out.pop().from_scalar(getfem::asm_H1_semi_norm(*mim, *mf, U.real(), bv));
    else out.pop().from_scalar(getfem::asm_H1_semi_norm(*mim, *mf, U.cplx(), bv));
  } else if (check_cmd(cmd, "H1 norm", in, out, 1, 2, 0, 1) && U_is_a_vector(U,cmd)) {
    /*@FUNC n = ::COMPUTE('H1 norm',@tmim mim[, @mat CVids])
    Compute the H1 norm of `U`.

    If `CVids` is given, the norm will be computed only on the listed
    convexes.@*/
    const getfem::mesh_im *mim = in.pop().to_const_mesh_im();
    dal::bit_vector bv = in.remaining() ?
      in.pop().to_bit_vector(&mf->convex_index()) : mf->convex_index();
    if (!U.is_complex())
      out.pop().from_scalar(getfem::asm_H1_norm(*mim, *mf, U.real(), bv));
    else out.pop().from_scalar(getfem::asm_H1_norm(*mim, *mf, U.cplx(), bv));
  } else if (check_cmd(cmd, "H2 semi norm", in, out, 1, 2, 0, 1) && U_is_a_vector(U,cmd)) {
    /*@FUNC n = ::COMPUTE('H2 semi norm',@tmim mim[, @mat CVids])
    Compute the L2 norm of D^2(`U`).

    If `CVids` is given, the norm will be computed only on the listed
    convexes.@*/
    const getfem::mesh_im *mim = in.pop().to_const_mesh_im();
    dal::bit_vector bv = in.remaining() ?
      in.pop().to_bit_vector(&mf->convex_index()) : mf->convex_index();
    if (!U.is_complex())
      out.pop().from_scalar(getfem::asm_H2_semi_norm(*mim, *mf, U.real(), bv));
    else out.pop().from_scalar(getfem::asm_H2_semi_norm(*mim, *mf, U.cplx(), bv));
  } else if (check_cmd(cmd, "H2 norm", in, out, 1, 2, 0, 1) && U_is_a_vector(U,cmd)) {
    /*@FUNC n = ::COMPUTE('H2 norm',@tmim mim[, @mat CVids])
    Compute the H2 norm of `U`.

    If `CVids` is given, the norm will be computed only on the listed
    convexes.@*/
    const getfem::mesh_im *mim = in.pop().to_const_mesh_im();
    dal::bit_vector bv = in.remaining() ?
      in.pop().to_bit_vector(&mf->convex_index()) : mf->convex_index();
    if (!U.is_complex())
      out.pop().from_scalar(getfem::asm_H2_norm(*mim, *mf, U.real(), bv));
    else out.pop().from_scalar(getfem::asm_H2_norm(*mim, *mf, U.cplx(), bv));
  } else if (check_cmd(cmd, "gradient", in, out, 1, 1, 0, 1)) {
    /*@FUNC DU = ::COMPUTE('gradient',@tmf mf_du)
    Compute the gradient of the field `U` defined on @tmf `mf_du`.

    The gradient is interpolated on the @tmf `mf_du`, and returned in
    `DU`. For example, if `U` is defined on a P2 @tmf, `DU` should be
    evaluated on a P1-discontinuous @tmf. `mf` and `mf_du` should
    share the same mesh.<Par>

    `U` may have any number of dimensions (i.e. this function is not
    restricted to the gradient of scalar fields, but may also be used
    for tensor fields). However the last dimension of `U` has to be
    equal to the number of dof of `mf`. For example, if `U` is a
    [3x3xNmf] array (where Nmf is the number of dof of `mf`), `DU` will
    be a [Nx3x3[xQ]xNmf_du] array, where N is the dimension of the mesh,
    Nmf_du is the number of dof of `mf_du`, and the optional Q dimension
    is inserted if `Qdim_mf != Qdim_mf_du`, where Qdim_mf is the Qdim of
    `mf` and Qdim_mf_du is the Qdim of `mf_du`.@*/
    const getfem::mesh_fem *mf_grad = in.pop().to_const_mesh_fem();
    error_for_non_lagrange_elements(*mf_grad, true);
    size_type qm = (mf_grad->get_qdim() == mf->get_qdim()) ? 1 : mf->get_qdim();
    if (!U.is_complex()) gf_compute_gradient<scalar_type>(out, *mf, *mf_grad, U.real(), qm);
    else                 gf_compute_gradient<complex_type>(out, *mf, *mf_grad, U.cplx(), qm);
  } else if (check_cmd(cmd, "hessian", in, out, 1, 1, 0, 1)) {
    /*@FUNC HU = ::COMPUTE('hessian',@tmf mf_h)
    Compute the hessian of the field `U` defined on @tmf `mf_h`.

    See also ::COMPUTE('gradient', @tmf mf_du).@*/
    const getfem::mesh_fem *mf_hess = in.pop().to_const_mesh_fem();
    error_for_non_lagrange_elements(*mf_hess, true);
    size_type qm = (mf_hess->get_qdim() == mf->get_qdim()) ? 1 : mf->get_qdim();
    if (!U.is_complex()) gf_compute_hessian<scalar_type>(out, *mf, *mf_hess, U.real(), qm);
    else                 gf_compute_hessian<complex_type>(out, *mf, *mf_hess, U.cplx(), qm);
  } else if (check_cmd(cmd, "eval on triangulated surface", in, out, 1, 2, 0, 1)) {
    int Nrefine = in.pop().to_integer(1, 1000);
    std::vector<convex_face> cvf;
    if (in.remaining() && !in.front().is_string()) {
      iarray v = in.pop().to_iarray(-1, -1);
      build_convex_face_lst(mf->linked_mesh(), cvf, &v);
    } else build_convex_face_lst(mf->linked_mesh(), cvf, 0);
    if (U.sizes().getn() != mf->nb_dof()) {
      THROW_BADARG("Wrong number of columns (need transpose?)");
    }
    eval_on_triangulated_surface(&mf->linked_mesh(), Nrefine, cvf, out, mf, U.real());
  } else if (check_cmd(cmd, "interpolate on", in, out, 1, 1, 0, 1)) {
    /*@FUNC Ui = ::COMPUTE('interpolate on',{@tmf mfi | @tsl sli})
    Interpolate a field on another @tmf or a @tsl.

    - Interpolation on another @tmf `mfi`:<par>
       `mfi` has to be Lagrangian. If `mf` and `mfi` share the same<par>
       mesh object, the interpolation will be much faster.<par>
    - Interpolation on a @tsl `sli`:<par>
       this is similar to interpolation on a refined P1-discontinuous<par>
       mesh, but it is much faster. This can also be used with<par>
       SLICE:INIT('points') to obtain field values at a given set of<par>
       points.<Par>

    See also ::ASM('interpolation matrix')
    @*/
    if (!U.is_complex()) gf_interpolate(in, out, *mf, U.real());
    else                 gf_interpolate(in, out, *mf, U.cplx());
  } else if (check_cmd(cmd, "extrapolate on", in, out, 1, 1, 0, 1)) {
    /*@FUNC Ue = ::COMPUTE('extrapolate on',@tmf mfe)
    Extrapolate a field on another @tmf.

    If the mesh of `mfe` is stricly included in the mesh of `mf`, this
    function does stricly the same job as ::COMPUTE('interpolate_on').
    However, if the mesh of `mfe` is not exactly included in `mf`
    (imagine interpolation between a curved refined mesh and a coarse
    mesh), then values which are outside `mf` will be
    extrapolated.<Par>

    See also ::ASM('extrapolation matrix')@*/
    const getfem::mesh_fem *mf_dest = in.pop().to_const_mesh_fem();
    error_for_non_lagrange_elements(*mf_dest, true);
    if (!U.is_complex()) {
      darray V = out.pop().create_darray(1, unsigned(mf_dest->nb_dof()));
      getfem::interpolation(*mf, *mf_dest, U.real(), V, 2);
    } else {
      carray V = out.pop().create_carray(1, unsigned(mf_dest->nb_dof()));
      getfem::interpolation(*mf, *mf_dest, U.cplx(), V, 2);
    }
  } else if (check_cmd(cmd, "mesh edges deformation", in, out, 1, 2, 0, 1)) {
    /* a virer qd gf_plot_mesh aura ete refait a neuf */
    unsigned N = in.pop().to_integer(2,10000);
    mesh_edges_deformation(mf, U.real(), N, in, out);
  } else if (check_cmd(cmd, "error_estimate", in, out, 1, 1, 0, 1)) {
    /*@FUNC E = ::COMPUTE('error estimate', @tmim mim)
    Compute an a posteriori error estimate.

    Currently there is only one which is available: for each convex,
    the jump of the normal derivative is integrated on its faces.@*/
    const getfem::mesh_im &mim = *in.pop().to_const_mesh_im();
    darray err =
      out.pop().create_darray_h(unsigned(mim.linked_mesh().convex_index().last_true()+1));
    if (!U.is_complex())
      getfem::error_estimate(mim, *mf, U.real(), err, mim.convex_index());
    else
      getfem::error_estimate(mim, *mf, U.cplx(), err, mim.convex_index());
  } else if (check_cmd(cmd, "convect", in, out, 4, 5, 0, 0)) {
    /*@FUNC E = ::COMPUTE('convect', @tmf mf_v, @dvec V, @scalar dt, @int nt[, @str option])
    Compute a convection of `U` with regards to a steady state velocity
    field `V` with a Characteristic-Galerkin method. This
    method is restricted to pure Lagrange fems for U. `mf_v` should represent
    a continuous finite element method. `dt` is the integration time and `nt`
    is the number of integration step on the caracteristics. `option` is an
    option for the part of the boundary where there is a re-entrant convection.
    `option = 'extrapolation'` for an extrapolation on the nearest element
    or `option = 'unchanged'` for a constant value on that boundary.
    This method is rather dissipative, but stable.
    @*/

    const getfem::mesh_fem *mf_v = in.pop().to_const_mesh_fem();
    rcarray V              = in.pop().to_rcarray();
    scalar_type dt = in.pop().to_scalar();
    size_type nt = in.pop().to_integer(0,100000);
    std::string option;
    if (in.remaining()) option = in.pop().to_string();
    getfem::convect_boundary_option opt;
    if (option.size() == 0)
      opt = getfem::CONVECT_EXTRAPOLATION;
    else if (cmd_strmatch(option, "extrapolation"))
      opt = getfem::CONVECT_EXTRAPOLATION;
    else if (cmd_strmatch(option, "unchanged"))
      opt = getfem::CONVECT_UNCHANGED;
    else 
      THROW_BADARG("Bad option " << option<< " for convect command. "
		   "should be 'extrapolation' or 'unchanged'");

    if (U.is_complex() || V.is_complex())
      THROW_BADARG("Sorry, complex version of convect to be interfaced");
    getfem::convect(*mf, U.real(), *mf_v, V.real(), dt, nt, opt);
  } else  bad_cmd(cmd);
}
