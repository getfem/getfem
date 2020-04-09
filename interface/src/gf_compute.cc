/*===========================================================================

 Copyright (C) 2006-2020 Yves Renard, Julien Pommier.

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
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/

#include <getfemint_misc.h>
#include <getfem/getfem_derivatives.h>
#include <getfem/getfem_interpolation.h>
#include <getfem/getfem_assembling.h>
#include <getfem/getfem_mesh_slice.h>
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
    getfem::compute_hessian(mf, mf_hess,
			    gmm::sub_vector(U, gmm::sub_slice(qq, mf.nb_dof(),
							      qqdim)), tmp);
    for (unsigned j=0, pos=qq*N*N; j < tmp.size(); j+=N*N) {
      for (unsigned k=0; k < N*N; ++k) D2U[pos+k] = tmp[j+k];
      pos += qqdim*N*N;
    }
  }
}

template <typename T> static void
gf_interpolate(getfemint::mexargs_in& in, getfemint::mexargs_out& out,
	       const getfem::mesh_fem& mf, const garray<T> &U) {
  array_dimensions dims;
  dims.push_back(U,0,U.ndim()-1,true);
  if (is_meshfem_object(in.front())) {
      const getfem::mesh_fem& mf_dest = *to_meshfem_object(in.pop());
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
    std::vector<T> VV(V.size());
    getfem::interpolation(mf, mf_dest, U, VV);
    gmm::copy(VV, V);
  }
  else if (is_slice_object(in.front())) {
    getfem::stored_mesh_slice *sl = to_slice_object(in.pop());

    for (size_type i=0; i < sl->nb_convex(); ++i)
      if (!mf.linked_mesh().convex_index().is_in(sl->convex_num(i)))
      THROW_BADARG("the slice is not compatible with the mesh_fem "
		   "(cannot find convex " << sl->convex_num(i) << ")");

    if (mf.get_qdim() != 1) dims.push_back(mf.get_qdim());
    dims.push_back(unsigned(sl->nb_points()));
    dims.opt_transform_col_vect_into_row_vect();
    garray<T> V = out.pop().create_array(dims, T());
    std::vector<T> VV(V.size());
    sl->interpolate(mf, U, VV);
    gmm::copy(VV, V);
  }
  else {
    size_type N = mf.linked_mesh().dim();
    darray st = in.pop().to_darray();
    std::vector<double> PTS(st.begin(), st.end());
    size_type nbpoints = gmm::vect_size(PTS) / N;
    getfem::base_node p(N);
    getfem::mesh_trans_inv mti(mf.linked_mesh());
    for (size_type i = 0; i < nbpoints; ++i) {
      gmm::copy(gmm::sub_vector(PTS, gmm::sub_interval(i*N, N)), p);
      // cout << "adding point" << p << endl;
      mti.add_point(p);
    }
    
    size_type qmult = mf.get_qdim();
    if (qmult != 1) dims.push_back(unsigned(qmult));
    dims.push_back(unsigned(nbpoints));
    dims.opt_transform_col_vect_into_row_vect();
    garray<T> V = out.pop().create_array(dims,T());
    std::vector<T> VV(V.size());
    getfem::base_matrix Maux;
    // cout << "begin interpolation, qmult = " << qmult << endl;
    getfem::interpolation(mf, mti, U, VV, Maux, 0);
    // cout << "end interpolation" << endl;
    gmm::copy(VV, V);
  }
  // else THROW_BADARG("expecting a mesh_fem or a mesh_slice for interpolation");
}

bool U_is_a_vector(const rcarray &U, const std::string& cmd) {
  if (U.sizes().size() == U.sizes().dim(-1)) return true;
  else THROW_BADARG("the U argument for the function " << cmd
		    << " must be a one-dimensional array");
  return false;
}

/*@GFDOC
  @ARGS{@tmf MF, @vec U}
  Various computations involving the solution U to a finite element problem.
@*/








// Object for the declaration of a new sub-command.

struct sub_gf_compute : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   const getfem::mesh_fem *mf,
		   rcarray U) = 0;
};

typedef std::shared_ptr<sub_gf_compute> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_compute {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       const getfem::mesh_fem *mf, rcarray U)		\
      { dummy_func(in); dummy_func(out); code }				\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }                           




void gf_compute(getfemint::mexargs_in& m_in, getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {


    /*@FUNC n = ('L2 norm', @tmim mim[, @mat CVids])
    Compute the L2 norm of the (real or complex) field `U`.

    If `CVids` is given, the norm will be computed only on the listed
    elements.@*/
    sub_command
      ("L2 norm", 1, 2, 0, 1,
       U_is_a_vector(U, "L2 norm");
       const getfem::mesh_im *mim = to_meshim_object(in.pop());
       dal::bit_vector bv = in.remaining() ?
       in.pop().to_bit_vector(&mf->convex_index()) : mf->convex_index();
       if (!U.is_complex())
	 out.pop().from_scalar(getfem::asm_L2_norm(*mim, *mf, U.real(), bv));
       else out.pop().from_scalar(getfem::asm_L2_norm(*mim, *mf, U.cplx(),bv));
       );

    /*@FUNC n = ('L2 dist', @tmim mim, @tmf mf2, @vec U2[, @mat CVids])
    Compute the L2 distance between `U` and `U2`.

    If `CVids` is given, the norm will be computed only on the listed
    elements.@*/
    sub_command
      ("L2 dist", 3, 4, 0, 1,
       U_is_a_vector(U, "L2 dist");
       const getfem::mesh_im *mim = to_meshim_object(in.pop());
       const getfem::mesh_fem *mf_2 = to_meshfem_object(in.pop());
       if (!U.is_complex()) {
         darray st = in.pop().to_darray();
         std::vector<double> V(st.begin(), st.end());
	 dal::bit_vector bv = in.remaining() ?
	   in.pop().to_bit_vector(&mf->convex_index()) : mf->convex_index();

	 out.pop().from_scalar(getfem::asm_L2_dist(*mim, *mf, U.real(),
						   *mf_2, V, bv));
       } else {	 
	 carray st = in.pop().to_carray();
	 std::vector<std::complex<double> > V(st.begin(), st.end());
 	 dal::bit_vector bv = in.remaining() ?
 	   in.pop().to_bit_vector(&mf->convex_index()) : mf->convex_index();
	 out.pop().from_scalar(getfem::asm_L2_dist(*mim, *mf, U.cplx(),
						   *mf_2, V, bv));
       }
       );


   /*@FUNC n = ('H1 semi norm', @tmim mim[, @mat CVids])
    Compute the L2 norm of grad(`U`).

    If `CVids` is given, the norm will be computed only on the listed
    elements.@*/
    sub_command
      ("H1 semi norm", 1, 2, 0, 1,
       U_is_a_vector(U, "H1 semi norm");
       const getfem::mesh_im *mim = to_meshim_object(in.pop());
       dal::bit_vector bv = in.remaining() ?
       in.pop().to_bit_vector(&mf->convex_index()) : mf->convex_index();
       if (!U.is_complex())
	 out.pop().from_scalar(getfem::asm_H1_semi_norm(*mim, *mf,
							U.real(), bv));
       else out.pop().from_scalar(getfem::asm_H1_semi_norm(*mim,
							   *mf, U.cplx(), bv));
       );


    /*@FUNC n = ('H1 semi dist', @tmim mim, @tmf mf2, @vec U2[, @mat CVids])
    Compute the semi H1 distance between `U` and `U2`.

    If `CVids` is given, the norm will be computed only on the listed
    elements.@*/
    sub_command
      ("H1 semi dist", 3, 4, 0, 1,
       U_is_a_vector(U, "H1 semi dist");
       const getfem::mesh_im *mim = to_meshim_object(in.pop());
       const getfem::mesh_fem *mf_2 = to_meshfem_object(in.pop());
       if (!U.is_complex()) {
         darray st = in.pop().to_darray();
         std::vector<double> V(st.begin(), st.end());
	 dal::bit_vector bv = in.remaining() ?
	   in.pop().to_bit_vector(&mf->convex_index()) : mf->convex_index();

	 out.pop().from_scalar(getfem::asm_H1_semi_dist(*mim, *mf, U.real(),
							*mf_2, V, bv));
       } else {
	 carray st = in.pop().to_carray();
	 std::vector<std::complex<double> > V(st.begin(), st.end());
 	 dal::bit_vector bv = in.remaining() ?
 	   in.pop().to_bit_vector(&mf->convex_index()) : mf->convex_index();
	 
	 out.pop().from_scalar(getfem::asm_H1_semi_dist(*mim, *mf, U.cplx(),
							*mf_2, V, bv));
       }
       );


    /*@FUNC n = ('H1 norm', @tmim mim[, @mat CVids])
    Compute the H1 norm of `U`.

    If `CVids` is given, the norm will be computed only on the listed
    elements.@*/
    sub_command
      ("H1 norm", 1, 2, 0, 1,
       U_is_a_vector(U, "H1 norm");
       const getfem::mesh_im *mim = to_meshim_object(in.pop());
       dal::bit_vector bv = in.remaining() ?
       in.pop().to_bit_vector(&mf->convex_index()) : mf->convex_index();
       if (!U.is_complex())
	 out.pop().from_scalar(getfem::asm_H1_norm(*mim, *mf, U.real(), bv));
       else out.pop().from_scalar(getfem::asm_H1_norm(*mim, *mf,
						      U.cplx(), bv));
       );


    /*@FUNC n = ('H2 semi norm', @tmim mim[, @mat CVids])
    Compute the L2 norm of D^2(`U`).

    If `CVids` is given, the norm will be computed only on the listed
    elements.@*/
    sub_command
      ("H2 semi norm", 1, 2, 0, 1,
       U_is_a_vector(U, "H2 semi norm");
       const getfem::mesh_im *mim = to_meshim_object(in.pop());
       dal::bit_vector bv = in.remaining() ?
       in.pop().to_bit_vector(&mf->convex_index()) : mf->convex_index();
       if (!U.is_complex())
	 out.pop().from_scalar(getfem::asm_H2_semi_norm(*mim, *mf,
							U.real(), bv));
       else out.pop().from_scalar(getfem::asm_H2_semi_norm(*mim, *mf,
							   U.cplx(), bv));
       );


    /*@FUNC n = ('H2 norm', @tmim mim[, @mat CVids])
    Compute the H2 norm of `U`.

    If `CVids` is given, the norm will be computed only on the listed
    elements.@*/
    sub_command
      ("H2 norm", 1, 2, 0, 1,
       U_is_a_vector(U, "H2 norm");
       const getfem::mesh_im *mim = to_meshim_object(in.pop());
       dal::bit_vector bv = in.remaining() ?
       in.pop().to_bit_vector(&mf->convex_index()) : mf->convex_index();
       if (!U.is_complex())
	 out.pop().from_scalar(getfem::asm_H2_norm(*mim, *mf, U.real(), bv));
       else out.pop().from_scalar(getfem::asm_H2_norm(*mim,*mf, U.cplx(), bv));
       );


    /*@FUNC DU = ('gradient', @tmf mf_du)
    Compute the gradient of the field `U` defined on @tmf `mf_du`.

    The gradient is interpolated on the @tmf `mf_du`, and returned in
    `DU`. For example, if `U` is defined on a P2 @tmf, `DU` should be
    evaluated on a P1-discontinuous @tmf. `mf` and `mf_du` should
    share the same mesh.

    `U` may have any number of dimensions (i.e. this function is not
    restricted to the gradient of scalar fields, but may also be used
    for tensor fields). However the last dimension of `U` has to be
    equal to the number of dof of `mf`. For example, if `U` is a
    [3x3xNmf] array (where Nmf is the number of dof of `mf`), `DU` will
    be a [Nx3x3[xQ]xNmf_du] array, where N is the dimension of the mesh,
    Nmf_du is the number of dof of `mf_du`, and the optional Q dimension
    is inserted if `Qdim_mf != Qdim_mf_du`, where Qdim_mf is the Qdim of
    `mf` and Qdim_mf_du is the Qdim of `mf_du`.@*/
    sub_command
      ("gradient", 1, 1, 0, 1,
       const getfem::mesh_fem *mf_grad = to_meshfem_object(in.pop());
       error_for_non_lagrange_elements(*mf_grad, true);
       size_type qm
         = (mf_grad->get_qdim() == mf->get_qdim()) ? 1 : mf->get_qdim();
       if (!U.is_complex())
	 gf_compute_gradient<scalar_type>(out, *mf, *mf_grad, U.real(), qm);
       else
	 gf_compute_gradient<complex_type>(out, *mf, *mf_grad, U.cplx(), qm);
       );

    /*@FUNC HU = ('hessian', @tmf mf_h)
    Compute the hessian of the field `U` defined on @tmf `mf_h`.

    See also ::COMPUTE('gradient', @tmf mf_du).@*/
    sub_command
      ("hessian", 1, 1, 0, 1,
       const getfem::mesh_fem *mf_hess = to_meshfem_object(in.pop());
       error_for_non_lagrange_elements(*mf_hess, true);
       size_type qm = (mf_hess->get_qdim() == mf->get_qdim()) ? 1 : mf->get_qdim();
       if (!U.is_complex())
	 gf_compute_hessian<scalar_type>(out, *mf, *mf_hess, U.real(), qm);
       else
	 gf_compute_hessian<complex_type>(out, *mf, *mf_hess, U.cplx(), qm);
       );


    /*@FUNC UP = ('eval on triangulated surface', @int Nrefine, [@vec CVLIST])
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
    @*/
    sub_command
      ("eval on triangulated surface", 1, 2, 0, 1,
       int Nrefine = in.pop().to_integer(1, 1000);
       std::vector<convex_face> cvf;
       if (in.remaining() && !in.front().is_string()) {
	 iarray v = in.pop().to_iarray(-1, -1);
	 build_convex_face_lst(mf->linked_mesh(), cvf, &v);
       } else build_convex_face_lst(mf->linked_mesh(), cvf, 0);
       if (U.sizes().getn() != mf->nb_dof()) {
	 THROW_BADARG("Wrong number of columns (need transpose ?)");
       }
       eval_on_triangulated_surface(&mf->linked_mesh(), Nrefine, cvf, out,
				    mf, U.real());
       );


    /*@FUNC Ui = ('interpolate on', {@tmf mfi | @tsl sli | @vec pts})
    Interpolate a field on another @tmf or a @tsl or a list of points.

    - Interpolation on another @tmf `mfi`:
       `mfi` has to be Lagrangian. If `mf` and `mfi` share the same
       mesh object, the interpolation will be much faster.
    - Interpolation on a @tsl `sli`:
       this is similar to interpolation on a refined P1-discontinuous
       mesh, but it is much faster. This can also be used with
       SLICE:INIT('points') to obtain field values at a given set of
       points.
    - Interpolation on a set of points `pts`

    See also ::ASM('interpolation matrix')
    @*/
    sub_command
      ("interpolate on", 1, 1, 0, 1,
       if (!U.is_complex()) gf_interpolate(in, out, *mf, U.real());
       else                 gf_interpolate(in, out, *mf, U.cplx());
      );


    /*@FUNC Ue = ('extrapolate on', @tmf mfe)
    Extrapolate a field on another @tmf.

    If the mesh of `mfe` is stricly included in the mesh of `mf`, this
    function does stricly the same job as ::COMPUTE('interpolate_on').
    However, if the mesh of `mfe` is not exactly included in `mf`
    (imagine interpolation between a curved refined mesh and a coarse
    mesh), then values which are outside `mf` will be
    extrapolated.

    See also ::ASM('extrapolation matrix')@*/
    sub_command
      ("extrapolate on", 1, 1, 0, 1,
       const getfem::mesh_fem *mf_dest = to_meshfem_object(in.pop());
       error_for_non_lagrange_elements(*mf_dest, true);
       if (!U.is_complex()) {
	 darray V = out.pop().create_darray(1, unsigned(mf_dest->nb_dof()));
	 getfem::interpolation(*mf, *mf_dest, U.real(), V, 2);
       } else {
	 carray V = out.pop().create_carray(1, unsigned(mf_dest->nb_dof()));
	 getfem::interpolation(*mf, *mf_dest, U.cplx(), V, 2);
       }
       );

    
    /*@FUNC E = ('error estimate', @tmim mim)
    Compute an a posteriori error estimate.

    Currently there is only one which is available: for each convex,
    the jump of the normal derivative is integrated on its faces.@*/
    sub_command
      ("error_estimate", 1, 1, 0, 1,
       const getfem::mesh_im &mim = *(to_meshim_object(in.pop()));
       darray err =
       out.pop().create_darray_h
       (unsigned(mim.linked_mesh().convex_index().last_true()+1));
       if (!U.is_complex())
	 getfem::error_estimate(mim, *mf, U.real(), err, mim.convex_index());
       else {
	 getfem::base_vector err_imag(gmm::vect_size(err));
	 getfem::error_estimate(mim, *mf, gmm::imag_part(U.cplx()), err_imag,
				mim.convex_index());
	 getfem::error_estimate(mim, *mf, gmm::real_part(U.cplx()), err,
				mim.convex_index());
	 gmm::add(err_imag, err);
       }
       );

      
#ifdef EXPERIMENTAL_PURPOSE_ONLY
            
    /*@FUNC E = ('error estimate nitsche', @tmim mim, @int GAMMAC, @int GAMMAN, @scalar lambda_, @scalar mu_, @scalar gamma0, @scalar f_coeff, @scalar vertical_force)
    Compute an a posteriori error estimate in the case of Nitsche method.

    Currently there is only one which is available: for each convex,
    the jump of the normal derivative is integrated on its faces.@*/
    sub_command
      ("error_estimate_nitsche", 8, 8, 0, 1,
       const getfem::mesh_im &mim = *to_meshim_object(in.pop());
       int GAMMAC = in.pop().to_integer();
       int GAMMAN = in.pop().to_integer();
       scalar_type lambda = in.pop().to_scalar();
       scalar_type mu = in.pop().to_scalar();
       scalar_type gamma0 = in.pop().to_scalar();
       scalar_type f_coeff = in.pop().to_scalar();
       scalar_type vertical_force = in.pop().to_scalar();
       unsigned si = unsigned(mim.linked_mesh().convex_index().last_true()+1);
       darray err =
       out.pop().create_darray_h(si);
       getfem::base_vector ERR(si);
       getfem::base_vector UU(U.real().size());
       gmm::copy(U.real(), UU);
       getfem::error_estimate_nitsche(mim, *mf, UU, GAMMAC, GAMMAN, lambda, mu, gamma0, f_coeff,vertical_force, ERR);
       gmm::copy(ERR, err);
       );   
      
#endif
      
      
      
      
    /*@FUNC ('convect', @tmf mf_v, @dvec V, @scalar dt, @int nt[, @str option[, @dvec per_min, @dvec per_max]])
    Compute a convection of `U` with regards to a steady state velocity
    field `V` with a Characteristic-Galerkin method. The result is returned
    in-place in `U`.
    This method is restricted to pure Lagrange fems for U. `mf_v` should
    represent a continuous finite element method. `dt` is the integration time
    and `nt` is the number of integration step on the caracteristics. `option`
    is an option for the part of the boundary where there is a re-entrant
    convection.
    `option = 'extrapolation'` for an extrapolation on the nearest element,
    `option = 'unchanged'` for a constant value on that boundary or
    `option = 'periodicity'` for a peridiodic boundary. For this latter option
    the two vectors per_min, per_max has to be given and represent the limits
    of the periodic domain (on components where per_max[k] < per_min[k]
    no operation is done).
    This method is rather dissipative, but stable.
    @*/
    sub_command
      ("convect", 4, 7, 0, 0,
       const getfem::mesh_fem *mf_v = to_meshfem_object(in.pop());
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
       else if (cmd_strmatch(option, "periodicity"))
	 opt = getfem::CONVECT_PERIODICITY;
       else if (cmd_strmatch(option, "unchanged"))
	 opt = getfem::CONVECT_UNCHANGED;
       else 
	 THROW_BADARG("Bad option " << option<< " for convect command. "
		      "should be 'extrapolation', 'unchanged' or "
                      "'periodicity'");

       getfem::base_node per_min;
       getfem::base_node per_max;
       if (in.remaining()) {
         rcarray pmin = in.pop().to_rcarray();
         rcarray pmax = in.pop().to_rcarray();
         size_type N = mf_v->linked_mesh().dim();
         per_min.resize(N);
         per_max.resize(N); 
         gmm::copy(pmin.real(), per_min);
         gmm::copy(pmax.real(), per_max);
       }

       if (U.is_complex() || V.is_complex())
	 THROW_BADARG("Sorry, complex version of convect to be interfaced");
       getfem::convect(*mf, U.real(), *mf_v, V.real(),
                       dt, nt, opt, per_min, per_max);

       );


  }
  
  
  if (m_in.narg() < 3)  THROW_BADARG( "Wrong number of input arguments");

  const getfem::mesh_fem *mf   = to_meshfem_object(m_in.pop());
  rcarray U              = m_in.pop().to_rcarray();
  m_in.last_popped().check_trailing_dimension(int(mf->nb_dof()));
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out, mf, U);
  }
  else bad_cmd(init_cmd);

}

/*@MATLABFUNC [U2[,MF2,[,X[,Y[,Z]]]]] = ('interpolate on Q1 grid', {'regular h', hxyz | 'regular N', Nxyz | X[,Y[,Z]]})
  
  Creates a cartesian Q1 mesh fem and interpolates U on it. The
  returned field U2 is organized in a matrix such that in can be drawn
  via the MATLAB command 'pcolor'. The first dimension is the Qdim of
  MF (i.e.  1 if U is a scalar field)

  example (mf_u is a 2D mesh_fem):
  >> Uq=gf_compute(mf_u, U, 'interpolate on Q1 grid', 'regular h', [.05, .05]);
  >> pcolor(squeeze(Uq(1,:,:)));
 @*/
/*@MATLABEXT
  if (nargin>=3 & strcmpi(varargin{3}, 'interpolate on Q1 grid')),
    [varargout{1:nargout}]=gf_compute_Q1grid_interp(varargin{[1 2 4:nargin]});
    return;
  end;
  @*/
