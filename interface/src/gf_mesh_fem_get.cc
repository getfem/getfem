// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2009 Julien Pommier.
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
// $Id$
#include <getfemint_misc.h>
#include <getfemint_mesh_fem.h>
#include <getfemint_integ.h>
#include <getfemint_pfem.h>
#include <getfem/getfem_export.h>
#include <getfem/getfem_mesh_fem_level_set.h>
#include <getfemint_mesh_levelset.h>
#include <getfem/getfem_mesh_im.h>
#include <getfem/getfem_partial_mesh_fem.h>
#include <getfemint_gsparse.h>

using namespace getfemint;

typedef enum {IS_LAGRANGE, IS_EQUIVALENT, IS_POLYNOMIAL} test_what;
static void
test_fems(test_what what, const getfem::mesh_fem *mf, mexargs_in& in,
          mexargs_out& out) {
  dal::bit_vector cvlst;
  bool return_bool = false;
  dal::bit_vector islst;
  if (in.remaining())
    cvlst = in.pop().to_bit_vector(&mf->linked_mesh().convex_index());
  else { cvlst = mf->linked_mesh().convex_index(); return_bool = true; }
  for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv) {
    bool it_is = false;
    if (!mf->linked_mesh().convex_index()[cv])
      THROW_ERROR( "convex " << cv+1 << " does not exist");
    check_cv_fem(*mf, cv);
    switch (what) { /* hope the compiler will optimize that */
    case IS_LAGRANGE:   it_is = mf->fem_of_element(cv)->is_lagrange(); break;
    case IS_EQUIVALENT: it_is = mf->fem_of_element(cv)->is_equivalent(); break;
    case IS_POLYNOMIAL: it_is = mf->fem_of_element(cv)->is_polynomial(); break;
    }
    if (it_is) islst.add(cv);
  }
  if (return_bool)
    out.pop().from_integer
      ((!(mf->is_reduced()) &&
          islst.card() == mf->linked_mesh().convex_index().card()) ? 1 : 0);
  else
    out.pop().from_bit_vector(islst);
}

static void get_basic_fem_of_convexes(const getfem::mesh_fem& mf,
                                      mexargs_in& in, mexargs_out& out) {
  dal::bit_vector cvlst;
  if (in.remaining())
    cvlst = in.pop().to_bit_vector(&mf.linked_mesh().convex_index());
  else { cvlst = mf.linked_mesh().convex_index(); }
  std::vector<id_type> ids; ids.reserve(cvlst.card());
  for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv) {
    if (mf.convex_index().is_in(cv)) {
      getfemint_pfem *gfi_pf =
        getfemint_pfem::get_from(mf.fem_of_element(cv));
      ids.push_back(gfi_pf->get_id());
    } else
      ids.push_back(id_type(-1));
  }
  out.return_packed_obj_ids(ids, FEM_CLASS_ID);
}

static dal::bit_vector
get_cv_dof_list(getfem::mesh_fem *mf, mexargs_in& in) {
  std::vector<convex_face> cvf;
  dal::bit_vector dof;
  if (in.remaining()) {
    iarray v = in.pop().to_iarray(-2, -1);
    build_convex_face_lst(mf->linked_mesh(), cvf, &v);
  } else build_convex_face_lst(mf->linked_mesh(), cvf, 0);
  for (size_type j = 0; j < cvf.size(); ++j) {
    size_type cv = cvf[j].cv;
    size_type  f = cvf[j].f;
    if (!mf->convex_index().is_in(cv))
      THROW_ERROR( "convex " << cv+1 << " has no FEM!");
    if (f != dim_type(-1)) {
      getfem::mesh_fem::ind_dof_face_ct
        c = mf->ind_basic_dof_of_face_of_element(cv,short_type(f));
      //std::for_each(c.begin(), c.end(), std::bind1st(std::mem_fun(&dal::bit_vector::add),&dof)); // not SGI STL compliant !?
      for (unsigned i=0; i < c.size(); ++i)
        dof.add(c[i]);
    } else {
      getfem::mesh_fem::ind_dof_ct cvdof = mf->ind_basic_dof_of_element(cv);
      for (unsigned i=0; i < cvdof.size(); ++i)
        dof.add(cvdof[i]);
    }
  }
  return dof;
}

static void
non_conformal_dof(getfem::mesh_fem &mf, mexargs_in &in, mexargs_out &out) {
  dal::bit_vector cvlst;
  const getfem::mesh &m = mf.linked_mesh();

  std::vector<bgeot::short_type> dcnt(mf.nb_basic_dof());

  if (in.remaining()) cvlst = in.pop().to_bit_vector(&m.convex_index());
  else cvlst = m.convex_index();

  for (dal::bv_visitor ic(cvlst); !ic.finished(); ++ic) {
    check_cv_fem(mf, ic);
    for (short_type f = 0; f < m.structure_of_convex(ic)->nb_faces(); f++) {
      bgeot::short_type q;
      if (!m.is_convex_having_neighbour(ic, f)) {
        q = 2;
      } else {
        q = 1;
      }
      for (short_type i = 0; i < mf.ind_basic_dof_of_face_of_element(ic,f).size();
           ++i) {
        size_type ind = mf.ind_basic_dof_of_face_of_element(ic,f)[i];
        dcnt[ind]= short_type(dcnt[ind] + q);
      }
    }
  }
  iarray w = out.pop().create_iarray_h(
    unsigned(std::count_if(dcnt.begin(),
                           dcnt.end(),
                           std::bind2nd(std::equal_to<bgeot::short_type>(),1))));
  size_type i,j=0;
  /*
  std::copy_if(dcnt.begin(), dcnt.end(),
               std::bind2nd(std::less_equal<bgeot::short_type>(),1)));
   */
  for (i=0; i < dcnt.size(); ++i) {
    if (dcnt[i] == 1) w[j++] = int(i+config::base_index());
  }
}

static std::string get_vtk_dataset_name(getfemint::mexargs_in &in, int count) {
  std::string s;
  if (in.remaining() && in.front().is_string()) {
    s = in.pop().to_string();
  } else {
    std::stringstream name; name << "dataset" << count;
    s = name.str();
  }
  for (size_type i=0; i < s.length(); ++i)
    if (!isalnum(s[i])) s[i] = '_';
  return s;
}

static std::string get_dx_dataset_name(getfemint::mexargs_in &in) {
  std::string s;
  if (in.remaining() && in.front().is_string()) {
    s = in.pop().to_string();
  }
  for (size_type i=0; i < s.length(); ++i)
    if (!isalnum(s[i])) s[i] = '_';
  return s;
}

template <typename T> static void
interpolate_convex_data(const getfem::mesh_fem *pmf,
                        const garray<T> &u, getfemint::mexargs_out& out) {
  assert(u.dim(u.ndim()-1) == pmf->linked_mesh().convex_index().last_true()+1);
  array_dimensions ad;
  for (unsigned i=0; i < u.ndim()-1; ++i) ad.push_back(u.dim(i));
  ad.push_back(unsigned(pmf->nb_basic_dof()));
  garray<T> w = out.pop().create_array(ad, T());
  size_type q = u.size() / u.dim(u.ndim()-1);
  assert(w.size() == q * pmf->nb_dof());
  std::vector<unsigned> dofcnt(pmf->nb_basic_dof());

  /* stupid smoothing when the mesh_fem is not discontinuous */
  for (dal::bv_visitor cv(pmf->convex_index()); !cv.finished(); ++cv)
    for (unsigned k=0; k < pmf->nb_basic_dof_of_element(cv); ++k)
      dofcnt[pmf->ind_basic_dof_of_element(cv)[k]]++;


  for (dal::bv_visitor cv(pmf->convex_index()); !cv.finished(); ++cv) {
    for (unsigned k=0; k < pmf->nb_basic_dof_of_element(cv); ++k) {
      size_type d = pmf->ind_basic_dof_of_element(cv)[k];
      for (unsigned j=0; j < q; ++j)
        w[d*q+j] += u[cv*q+j] / T(dofcnt[d]);
    }
  }
}

/*MLABCOM
  FUNCTION [x] = gf_mesh_fem_get(meshfem MF, operation [, args])

  General function for inquiry about mesh_fem objects.

  @RDATTR MESHFEM:GET('nbdof')
  @RDATTR MESHFEM:GET('nb basic dof')
  @GET    MESHFEM:GET('basic dof from cv')
  @GET    MESHFEM:GET('basic dof from cvid')
  @GET    MESHFEM:GET('non conformal basic dof')
  @RDATTR MESHFEM:GET('qdim')
  @GET    MESHFEM:GET('fem')
  Example:
     cvid=gf_mesh_get(mf,'cvid');
     [f,c2f]=gf_mesh_fem_get(mf, 'fem');
     for i=1:size(f), sf{i}=gf_fem_get('char',f(i)); end;
     for i=1:size(c2f),
       disp(sprintf('the fem of convex %d is %s',...
            cvid(i),sf{i}));
     end;
  @GET    MESHFEM:GET('is_lagrangian')
  @GET    MESHFEM:GET('is_equivalent')
  @GET    MESHFEM:GET('is_polynomial')
  @RDATTR MESHFEM:GET('is_reduced')
  @GET    MESHFEM:GET('reduction matrix')
  @GET    MESHFEM:GET('extension matrix')
  @GET    MESHFEM:GET('basic dof on region')
  @GET    MESHFEM:GET('dof on region')
  @GET    MESHFEM:GET('basic dof nodes')
  @GET    MESHFEM:GET('basic dof from im')
  @GET    MESHFEM:GET('dof partition')
  @GET    MESHFEM:GET('interpolate_convex_data')
  @GET    MESHFEM:GET('save')
  @GET    MESHFEM:GET('char')
  @GET    MESHFEM:GET('export to vtk')
  @GET    MESHFEM:GET('export to dx')
  @GET    MESHFEM:GET('export to pos')
  @GET    MESHFEM:GET('linked mesh')
  @GET    MESHFEM:GET('memsize')

  * U=MESHFEM:GET('eval', expr [,DOFLST])

  Call gf_mesh_fem_get_eval. This function interpolates an expression on a
  lagrangian mesh_fem (for all dof except if DOFLST is specified). The expression can be a
  numeric constant, or a cell array containing numeric constants, string
  expressions or function handles. For example:
    U1=gf_mesh_fem_get(mf,'eval',1)
    U2=gf_mesh_fem_get(mf,'eval',[1;0]) % output has two rows
    U3=gf_mesh_fem_get(mf,'eval',[1 0]) % output has one row, only valid if qdim(mf)==2
    U4=gf_mesh_fem_get(mf,'eval',{'x';'y.*z';4;@myfunctionofxyz})

  $Id$
MLABCOM*/
/*MLABEXT
  if (nargin>=2 & strcmpi(varargin{2},'eval')),
    [varargout{1:nargout}]=gf_mesh_fem_get_eval(varargin{[1 3:nargin]}); return;
  end;
  MLABEXT*/


void gf_mesh_fem_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out) {
  if (in.narg() < 2) THROW_BADARG( "Wrong number of input arguments");
  getfemint_mesh_fem *mi_mf = in.pop().to_getfemint_mesh_fem();
  getfem::mesh_fem *mf   = &mi_mf->mesh_fem();
  std::string cmd        = in.pop().to_string();
  if (check_cmd(cmd, "nbdof", in, out, 0, 0, 0, 1)) {
    /*@RDATTR n = MESHFEM:GET('nbdof')
    Return the number of degrees of freedom (dof) of the @tmf.@*/
    out.pop().from_integer(int(mf->nb_dof()));
  } else if (check_cmd(cmd, "nb basic dof", in, out, 0, 0, 0, 1)) {
    /*@RDATTR n = MESHFEM:GET('nb basic dof')
    Return the number of basic degrees of freedom (dof) of the @tmf.@*/
    out.pop().from_integer(int(mf->nb_basic_dof()));
  } else if (check_cmd(cmd, "dof from cv", in, out, 1, 1, 0, 1)) {
    /*@GET DOF = MESHFEM:GET('dof from cv',@mat CVids)
    Deprecated function. Use MESHFEM:GET('basic dof from cv') instead. @*/
    infomsg() << "WARNING : gf_mesh_fem_get('dof from cv', ...) is a "
              << "deprecated command.\n"
              << "          Use gf_mesh_fem_get('basic dof from cv', "
              << "...) instead." << endl;
    dal::bit_vector dof = get_cv_dof_list(mf, in);
    out.pop().from_bit_vector(dof);
  } else if (check_cmd(cmd, "basic dof from cv", in, out, 1, 1, 0, 1)) {
    /*@GET DOF = MESHFEM:GET('basic dof from cv',@mat CVids)
    Return the dof of the convexes listed in `CVids`.

    WARNING: the Degree of Freedom might be returned in ANY order, do
    not use this function in your assembly routines. Use 'basic dof from cvid'
    instead, if you want to be able to map a convex number with its
    associated degrees of freedom.

    One can also get the list of basic dof on a set on convex faces, by
    indicating on the second row of `CVids` the faces numbers (with
    respect to the convex number on the first row).@*/
    dal::bit_vector dof = get_cv_dof_list(mf, in);
    out.pop().from_bit_vector(dof);
  } else if (check_cmd(cmd, "dof from cvid", in, out, 0, 1, 0, 2) ||
             check_cmd(cmd, "basic dof from cvid", in, out, 0, 1, 0, 2)) {
    if (check_cmd(cmd, "dof from cvid", in, out, 0, 1, 0, 2)) {
      /*@GET @CELL{DOFs, IDx} = MESHFEM:GET('dof from cvid'[, @mat CVids])
        Deprecated function. Use MESHFEM:GET('basic dof from cvid') instead.
        @*/
      infomsg() << "WARNING : gf_mesh_fem_get('dof from cvid', ...) is a "
                << "deprecated command.\n          Use gf_mesh_fem_get('basic "
                << "dof from cvid', ...) instead." << endl;
    }
    /*@GET @CELL{DOFs, IDx} = MESHFEM:GET('basic dof from cvid'[, @mat CVids])
    Return the degrees of freedom attached to each convex of the mesh.

    If `CVids` is omitted, all the convexes will be considered (equivalent
    to `CVids = 1 ... MESH:GET('max cvid')`).

    `IDx` is a @MATLAB{row }vector, `length(IDx) = length(CVids)+1`.
    `DOFs` is a @MATLAB{row }vector containing the concatenated list
    of dof of each convex in `CVids`. Each entry of `IDx` is the position
    of the corresponding convex point list in `DOFs`. Hence, for example,
    the list of points of the second convex is @MATLAB{DOFs(IDx(2):IDx(3)-1)}
    @PYTHON{DOFs[IDx(2):IDx(3)]}.

    If `CVids` contains convex #id which do not exist in the mesh, their
    point list will be empty.@*/
    dal::bit_vector cvlst;
    if (in.remaining()) cvlst = in.pop().to_bit_vector();
    else cvlst.add(0, mf->linked_mesh().convex_index().last_true() + 1);

    std::vector<size_type> pids;
    std::vector<size_type> idx;
    size_type pcnt = 0;
    for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv) {
      idx.push_back(size_type(pcnt + config::base_index()));
      if (mf->convex_index().is_in(cv))
        for (size_type i = 0; i< mf->nb_basic_dof_of_element(cv); ++i, ++pcnt)
          pids.push_back(size_type(mf->ind_basic_dof_of_element(cv)[i] + config::base_index()));
    }
    idx.push_back(size_type(pcnt + config::base_index()));

    iarray opids = out.pop().create_iarray_h(pids.size());
    if (pids.size()) std::copy(pids.begin(), pids.end(), &opids[0]);
    if (out.remaining() && idx.size()) {
      iarray oidx = out.pop().create_iarray_h(idx.size());
      std::copy(idx.begin(), idx.end(), &oidx[0]);
    }
  } else if (check_cmd(cmd, "non conformal dof", in, out, 0, 1, 0, 1)) {
    /*@GET MESHFEM:GET('non conformal dof'[, @mat CVids])
      Deprecated function. Use MESHFEM:GET('non conformal basic dof') instead.
      @*/
      infomsg() << "WARNING : gf_mesh_fem_get('non conformal dof', ...) is a "
                << "deprecated command.\n          Use gf_mesh_fem_get('non "
                << "conformal basic dof', ...) instead." << endl;
    non_conformal_dof(*mf, in, out);
  } else if (check_cmd(cmd, "non conformal basic dof", in, out, 0, 1, 0, 1)) {
    /*@GET MESHFEM:GET('non conformal basic dof'[, @mat CVids])
    Return partially linked degrees of freedom.

    Return the basic dof located on the border of a convex and which belong
    to only one convex, except the ones which are located on the border
    of the mesh.  For example, if the convex 'a' and 'b' share a common
    face, 'a' has a P1 FEM, and 'b' has a P2 FEM, then the basic dof on the
    middle of the face will be returned by this function (this can be
    useful when searching the interfaces between classical FEM and
    hierarchical FEM).@*/
    non_conformal_dof(*mf, in, out);
  } else if (check_cmd(cmd, "qdim", in, out, 0, 0, 0, 1)) {
    /*@RDATTR MESHFEM:GET('qdim')
    Return the dimension Q of the field interpolated by the @tmf.

    By default, Q=1 (scalar field). This has an impact on the dof numbering.@*/
    out.pop().from_integer(mf->get_qdim());
  } else if (check_cmd(cmd, "fem", in, out, 0, 1, 0, 2)) {
    /*@GET @CELL{FEMs, CV2F} = MESHFEM:GET('fem'[, @mat CVids])
    Return a list of FEM used by the @tmf.

    `FEMs` is an array of all @tfem objects found in the convexes
    given in `CVids`. If `CV2F` was supplied as an output argument,
    it contains, for each convex listed in `CVids`, the index of its
    correspounding FEM in `FEMs`.

    Convexes which are not part of the mesh, or convexes which do not
    have any FEM have their correspounding entry in `CV2F` set to -1.@*/
    get_basic_fem_of_convexes(*mf, in, out);
  } else if (check_cmd(cmd, "convex_index", in, out, 0, 0, 0, 1)) {
    /*@GET CVs = MESHFEM:GET('convex_index')
    Return the list of convexes who have a FEM.@*/
    out.pop().from_bit_vector(mf->convex_index());
  } else if (check_cmd(cmd, "is_lagrangian", in, out, 0, 1, 0, 1)) {
    /*@GET bB = MESHFEM:GET('is_lagrangian'[, @mat CVids])
    Test if the @tmf is Lagrangian.

    Lagrangian means that each base function Phi[i] is such that
    Phi[i](P[j]) = delta(i,j), where P[j] is the dof location of
    the jth base function, and delta(i,j) = 1 if i==j, else 0.<Par>

    If `CVids` is omitted, it returns 1 if all convexes in the mesh
    are Lagrangian. If `CVids` is used, it returns the convex indices
    (with respect to `CVids`) which are Lagrangian.@*/
    test_fems(IS_LAGRANGE , mf, in, out);
  } else if (check_cmd(cmd, "is_equivalent", in, out, 0, 1, 0, 1)) {
    /*@GET bB = MESHFEM:GET('is_equivalent'[, @mat CVids])
    Test if the @tmf is equivalent.

    See MESHFEM:GET('is_lagrangian')@*/
    test_fems(IS_EQUIVALENT, mf, in, out);
  } else if (check_cmd(cmd, "is_polynomial", in, out, 0, 1, 0, 1)) {
    /*@GET bB = MESHFEM:GET('is_polynomial'[, @mat CVids])
    Test if all base functions are polynomials.

    See MESHFEM:GET('is_lagrangian')@*/
    test_fems(IS_POLYNOMIAL, mf, in, out);
  } else if (check_cmd(cmd, "is_reduced", in, out, 0, 0, 0, 1)) {
    /*@RDATTR bB = MESHFEM:GET('is_reduced')
    Return 1 if the optional reduction matrix is applied to the dofs.@*/
    out.pop().from_integer(mf->is_reduced() ? 1 : 0);
  } else if (check_cmd(cmd, "reduction matrix", in, out, 0, 0, 0, 1)) {
    /*@GET bB = MESHFEM:GET('reduction matrix')
    Return the optional reduction matrix.@*/
    getfemint::gf_real_sparse_by_col
      M(gmm::mat_nrows(mf->reduction_matrix()),
        gmm::mat_ncols(mf->reduction_matrix()));
    gmm::copy(mf->reduction_matrix(), M);
    out.pop().from_sparse(M);
  } else if (check_cmd(cmd, "extension matrix", in, out, 0, 0, 0, 1)) {
    /*@GET bB = MESHFEM:GET('extension matrix')
    Return the optional extension matrix.@*/
    getfemint::gf_real_sparse_by_col
      M(gmm::mat_nrows(mf->extension_matrix()),
        gmm::mat_ncols(mf->extension_matrix()));
    gmm::copy(mf->extension_matrix(), M);
    out.pop().from_sparse(M);
  } else if (check_cmd(cmd, "basic dof on region", in, out, 1, 1, 0, 1)) {
    /*@GET DOFs = MESHFEM:GET('basic dof on region',@mat Rs)
    Return the list of basic dof (before the optional reduction) lying on one
    of the mesh regions listed in `Rs`.

    More precisely, this function returns the basic dof whose support is
    non-null on one of regions whose #ids are listed in `Rs` (note
    that for boundary regions, some dof nodes may not lie exactly
    on the boundary, for example the dof of Pk(n,0) lies on the center
    of the convex, but the base function in not null on the convex
    border).@*/
    iarray bnums = in.pop().to_iarray(-1);
    dal::bit_vector bv;
    for (size_type i=0; i < bnums.size(); ++i)
      bv |= mf->basic_dof_on_region(bnums[i]);
    out.pop().from_bit_vector(bv);
  } else if (check_cmd(cmd, "dof on region", in, out, 1, 1, 0, 1)) {
    /*@GET DOFs = MESHFEM:GET('dof on region',@mat Rs)
    Return the list of dof (after the optional reduction) lying on one
    of the mesh regions listed in `Rs`.

    More precisely, this function returns the basic dof whose support is
    non-null on one of regions whose #ids are listed in `Rs` (note
    that for boundary regions, some dof nodes may not lie exactly
    on the boundary, for example the dof of Pk(n,0) lies on the center
    of the convex, but the base function in not null on the convex
    border).

    For a reduced mesh_fem
    a dof is lying on a region if its potential corresponding shape
    function is nonzero on this region. The extension matrix is used
    to make the correspondance between basic and reduced dofs.@*/
    iarray bnums = in.pop().to_iarray(-1);
    dal::bit_vector bv;
    for (size_type i=0; i < bnums.size(); ++i)
      bv |= mf->dof_on_region(bnums[i]);
    out.pop().from_bit_vector(bv);
  } else if (check_cmd(cmd, "dof nodes", in, out, 0, 1, 0, 2)) {
    /*@GET DOFpts = MESHFEM:GET('dof nodes'[, @mat DOFids])
    Deprecated function. Use MESHFEM:GET('basic dof nodes') instead. @*/
    infomsg() << "WARNING : gf_mesh_fem_get('dof nodes', ...) is a deprecated "
              << "command.\n          Use gf_mesh_fem_get('basic dof nodes', "
              << "...) instead." << endl;

    dal::bit_vector dof_lst; dof_lst.add(0, mf->nb_basic_dof());
    if (in.remaining())
      dof_lst = in.pop().to_bit_vector(&dof_lst);
    darray w = out.pop().create_darray(mf->linked_mesh().dim(),
                                       unsigned(dof_lst.card()));
    size_type j = 0;
    for (dal::bv_visitor dof(dof_lst); !dof.finished(); ++dof, ++j) {
      if (mf->point_of_basic_dof(dof).size() != w.getm() || j >= w.getn())
        THROW_INTERNAL_ERROR;
      for (size_type i=0; i < w.getm(); i++)
        w(i,j)= mf->point_of_basic_dof(dof)[i];
    }
  } else if (check_cmd(cmd, "basic dof nodes", in, out, 0, 1, 0, 2)) {
    /*@GET DOFpts = MESHFEM:GET('basic dof nodes'[, @mat DOFids])
    Get location of basic degrees of freedom.

    Return the list of interpolation points for the specified
    dof #IDs in `DOFids` (if `DOFids` is omitted, all basic dof are
    considered).@*/
    dal::bit_vector dof_lst; dof_lst.add(0, mf->nb_basic_dof());
    if (in.remaining())
      dof_lst = in.pop().to_bit_vector(&dof_lst);
    darray w = out.pop().create_darray(mf->linked_mesh().dim(),
                                       unsigned(dof_lst.card()));
    size_type j = 0;
    for (dal::bv_visitor dof(dof_lst); !dof.finished(); ++dof, ++j) {
      if (mf->point_of_basic_dof(dof).size() != w.getm() || j >= w.getn())
        THROW_INTERNAL_ERROR;
      for (size_type i=0; i < w.getm(); i++)
        w(i,j)= mf->point_of_basic_dof(dof)[i];
        // std::copy(mf->point_of_dof(dof).begin(),mf->point_of_dof(dof).end(),
        //           &w(0,j));
    }
  } else if (check_cmd(cmd, "dof partition", in, out, 0, 0, 0, 1)) {
    /*@GET DOFP = MESHFEM:GET('dof partition')
    Get the 'dof_partition' array.

    Return the array which associates an integer (the partition number)
    to each convex of the @tmf. By default, it is an all-zero array.
    The degrees of freedom of each convex of the @tmf are connected
    only to the dof of neighbouring convexes which have the same
    partition number, hence it is possible to create partially
    discontinuous @tmf very easily.@*/
    iarray v = out.pop().create_iarray_h
      (unsigned(mf->linked_mesh().convex_index().last_true()+1));
    for (unsigned cv=0; cv < v.size(); ++cv) v[cv] = mf->get_dof_partition(cv);
  } else if (check_cmd(cmd, "save", in, out, 1, 2, 0, 0)) {
    /*@GET MESHFEM:GET('save',@str filename[, @str opt])
    Save a @tmf in a text file (and optionaly its linked mesh object
    if `opt` is the string 'with_mesh').@*/
    std::string s = in.pop().to_string();
    bool with_mesh = false;
    if (in.remaining()) {
      if (cmd_strmatch(in.pop().to_string(), "with mesh")) {
        with_mesh = true;
      } else THROW_BADARG("expecting string 'with mesh'");
    }
    std::ofstream o(s.c_str());
    if (!o) THROW_ERROR("impossible to write in file '" << s << "'");
    o << "% GETFEM MESH+FEM FILE " << endl;
    o << "% GETFEM VERSION " << GETFEM_VERSION << endl;
    if (with_mesh) mf->linked_mesh().write_to_file(o);
    mf->write_to_file(o);
    o.close();
  } else if (check_cmd(cmd, "char", in, out, 0, 0, 0, 1)) {
    /*@GET MESHFEM:GET('char'[, @str opt])
    Output a string description of the @tmf.

    By default, it does not include the description of the linked mesh
    object, except if `opt` is 'with_mesh'.@*/
    std::stringstream s;
    if (in.remaining() && cmd_strmatch(in.pop().to_string(),"with mesh"))
      mf->linked_mesh().write_to_file(s);
    mf->write_to_file(s);
    out.pop().from_string(s.str().c_str());
  } else if (check_cmd(cmd, "linked mesh", in, out, 0, 0, 0, 1)) {
    /*@GET m = MESHFEM:GET('linked mesh')
    Return a reference to the @tmesh object linked to `mf`.@*/
    out.pop().from_object_id(mi_mf->linked_mesh_id(), MESH_CLASS_ID);
  } else if (check_cmd(cmd, "export to vtk", in, out, 0, -1, 0, 0)) {
    /*@GET MESHFEM:GET('export to vtk',@str filename, ... ['ascii'], U, 'name'...)
    Export a @tmf and some fields to a vtk file.

    The FEM and geometric transformations will be mapped to order 1
    or 2 isoparametric Pk (or Qk) FEMs (as VTK does not handle higher
    order elements). If you need to represent high-order FEMs or
    high-order geometric transformations, you should consider
    SLICE:GET('export to vtk').@*/
    std::string fname = in.pop().to_string();
    bool ascii = false;
    while (in.remaining() && in.front().is_string()) {
      std::string cmd2 = in.pop().to_string();
      if (cmd_strmatch(cmd2, "ascii"))
        ascii = true;
      else THROW_BADARG("expecting 'ascii', got " << cmd2);
    }
    getfem::vtk_export exp(fname, ascii);
    exp.exporting(*mf);
    exp.write_mesh();
    int count = 1;
    while (in.remaining()) {
      const getfem::mesh_fem *mf2 = mf;
      if (in.remaining() >= 2 && in.front().is_mesh_fem())
        mf2 = in.pop().to_const_mesh_fem();
      darray U = in.pop().to_darray();
      in.last_popped().check_trailing_dimension(int(mf2->nb_dof()));
      exp.write_point_data(*mf2, U, get_vtk_dataset_name(in, count));
      count+=1;
    }
  } else if (check_cmd(cmd, "export to dx", in, out, 0, -1, 0, 0)) {
    /*@GET MESHFEM:GET('export to dx',@str filename, ...['as', @str mesh_name][,'edges']['serie',@str serie_name][,'ascii'][,'append'], U, 'name'...)
    Export a @tmf and some fields to an OpenDX file.

    This function will fail if the @tmf mixes different convex types
    (i.e. quads and triangles), or if OpenDX does not handle a specific
    element type (i.e. prism connections are not known by OpenDX).

    The FEM will be mapped to order 1 Pk (or Qk) FEMs. If you need to
    represent high-order FEMs or high-order geometric transformations,
    you should consider SLICE:GET('export to dx').@*/
    std::string fname = in.pop().to_string();
    bool ascii = false, append = false, edges = false;
    std::string mesh_name;
    std::string serie_name;
    while (in.remaining() && in.front().is_string()) {
      std::string cmd2 = in.pop().to_string();
      if (cmd_strmatch(cmd2, "ascii"))
        ascii = true;
      else if (cmd_strmatch(cmd2, "edges"))
        edges = true;
      else if (cmd_strmatch(cmd2, "as") && in.remaining())
        mesh_name = in.pop().to_string();
      else if (cmd_strmatch(cmd2, "append"))
        append = true;
      else if (cmd_strmatch(cmd2, "serie") && in.remaining())
        serie_name = in.pop().to_string();
      else THROW_BADARG("expecting 'ascii', got " << cmd2);
    }
    getfem::dx_export exp(fname, ascii, append);
    exp.exporting(*mf, mesh_name.c_str());
    exp.write_mesh();
    if (edges) exp.exporting_mesh_edges();
    while (in.remaining()) {
      const getfem::mesh_fem *mf2 = mf;
      if (in.remaining() >= 2 && in.front().is_mesh_fem())
        mf2 = in.pop().to_const_mesh_fem();
      darray U = in.pop().to_darray();
      in.last_popped().check_trailing_dimension(int(mf2->nb_dof()));
      exp.write_point_data(*mf2, U, get_dx_dataset_name(in));
      if (serie_name.size()) exp.serie_add_object(serie_name);
    }
  } else if (check_cmd(cmd, "export to pos", in, out, 1, -1, 0, 0)) {
    /*@GET MESHFEM:GET('export to pos',@str filename[, @str name][[,@tmf mf1], @mat U1, @str nameU1[[,@tmf mf2], @mat U2, @str nameU2,...])
    Export a @tmf and some fields to a pos file.

    The FEM and geometric transformations will be mapped to order 1
    isoparametric Pk (or Qk) FEMs (as GMSH does not handle higher
    order elements).@*/
    std::string fname = in.pop().to_string();
    std::string name = "";
    if (in.remaining() && in.front().is_string())
      name = in.pop().to_string();

    getfem::pos_export exp(fname);
    exp.write(*mf,name);
    while (in.remaining()) {
      const getfem::mesh_fem *mf2 = mf;
      if (in.remaining() >= 2 && in.front().is_mesh_fem()) {
        mf2 = in.pop().to_const_mesh_fem();
      }
      darray U = in.pop().to_darray();
      in.last_popped().check_trailing_dimension(int(mf2->nb_dof()));

      if (in.remaining() >= 1 && in.front().is_string()) {
        name = in.pop().to_string();
      } else THROW_BADARG("expecting string darray_name")

      exp.write(*mf2, U, name);
    }
  } else if (check_cmd(cmd, "dof_from_im", in, out, 1, 2, 0, 1)) {
    /*@GET MESHFEM:GET('dof_from_im',@tmim mim[, @int p])
    Return a selection of dof who contribute significantly to the
    mass-matrix that would be computed with `mf` and the integration
    method `mim`.

    `p` represents the dimension on what the integration method
    operates (default `p = mesh dimension`).

    IMPORTANT: you still have to set a valid integration method on
    the convexes which are not crosses by the levelset!@*/
    const getfem::mesh_im &mim = *in.pop().to_const_mesh_im();
    int P = -1;
    if (&mim.linked_mesh() != &mf->linked_mesh())
      THROW_BADARG("the mesh_im uses a different mesh");
    if (in.remaining())
      P = in.pop().to_integer(1, mim.linked_mesh().dim());
    out.pop().from_bit_vector(getfem::select_dofs_from_im(*mf, mim, P));
  } else if (check_cmd(cmd, "interpolate_convex_data", in, out, 1, 1, 0, 1)) {
    /*@GET U = MESHFEM:GET('interpolate_convex_data',@mat Ucv)

    Interpolate data given on each convex of the mesh to the @tmf dof.
    The @tmf has to be lagrangian, and should be discontinuous (typically
    a FEM_PK(N,0) or FEM_QK(N,0) should be used).

    The last dimension of the input vector Ucv should have
    MESH:GET('max cvid') elements.

    Example of use: MESHFEM:GET('interpolate_convex_data', MESH:GET('quality'))@*/
    in.front().check_trailing_dimension
      (int(mf->linked_mesh().convex_index().last_true()+1));
    if (in.front().is_complex())
      interpolate_convex_data(mf, in.pop().to_darray(), out);
    else interpolate_convex_data(mf, in.pop().to_carray(), out);
  } else if (check_cmd(cmd, "memsize", in, out, 0, 0, 0, 1)) {
    /*@GET z = MESHFEM:GET('memsize')
    Return the amount of memory (in bytes) used by the mesh_fem object.

    The result does not take into account the linked mesh object.@*/
    out.pop().from_integer(int(mf->memsize()));
  } else if (check_cmd(cmd, "has_linked_mesh_levelset", in, out, 0, 0, 0, 1)) {
    /*@GET MESHFEM:GET('has_linked_mesh_levelset')
    ???.@*/
    getfem::mesh_fem_level_set *mfls =
      dynamic_cast<getfem::mesh_fem_level_set*>(mf);
    out.pop().from_integer(mfls ? 1 : 0);
  } else if (check_cmd(cmd, "linked_mesh_levelset", in, out, 0, 0, 0, 1)) {
    /*@GET MESHFEM:GET('linked_mesh_levelset')
    ???@*/
    getfem::mesh_fem_level_set *mfls =
      dynamic_cast<getfem::mesh_fem_level_set*>(mf);
    if (mfls) {
      getfem::mesh_level_set *mls =
        const_cast<getfem::mesh_level_set*>(&mfls->linked_mesh_level_set());
      getfemint_mesh_levelset *gfi_mls =
        getfemint_mesh_levelset::get_from(mls);
      assert(gfi_mls);
      out.pop().from_object_id(gfi_mls->get_id(),
                               MESH_LEVELSET_CLASS_ID);
    } else THROW_BADARG("not a mesh_fem using a mesh_levelset");
  } else bad_cmd(cmd);
}
