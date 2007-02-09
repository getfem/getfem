// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2006-2006 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

#include <getfemint_misc.h>
#include <getfemint_mesh_fem.h>
#include <getfemint_integ.h>
#include <getfemint_pfem.h>
#include <getfem/getfem_export.h>
#include <getfem/getfem_mesh_fem_level_set.h>
#include <getfemint_mesh_levelset.h>
#include <getfem/getfem_mesh_im.h>
#include <getfem/getfem_partial_mesh_fem.h>

/*
  $Id$

  ChangeLog:
  $Log: gf_mesh_fem_get.cc,v $
  Revision 1.7  2006/03/28 10:06:35  pommier
  *** empty log message ***

  Revision 1.6  2006/02/14 17:57:17  pommier
  *** empty log message ***

  Revision 1.5  2006/02/10 17:04:03  pommier
  *** empty log message ***

  Revision 1.4  2006/02/06 14:35:39  pommier
  *** empty log message ***

  Revision 1.3  2006/01/18 11:21:52  pommier
  *** empty log message ***

  Revision 1.2  2005/03/08 16:50:12  pommier
  added meshim, many doc updates

  Revision 1.1  2005/01/21 15:31:37  renard
  *** empty log message ***

  Revision 1.9  2005/01/05 16:33:32  pommier
  *** empty log message ***

  Revision 1.8  2004/08/18 10:16:12  pommier
  many updates

  Revision 1.7  2004/06/21 09:03:46  pommier
  more work on gf_spmat and python interface

  Revision 1.6  2004/06/11 10:10:28  pommier
  commit recent work: handling of integer and complex arrays, interface to gmm preconditioners and sparse matrices

  Revision 1.5  2004/03/03 16:04:05  pommier
  renamed matlabint_* to getfemint_*

  Revision 1.4  2003/10/05 16:03:40  pommier
  change all bit_vector loops to bv_visitor

  Revision 1.3  2003/07/31 13:23:48  renard
  ajout du numero de convexe dans les paramtres de gen_compute

  Revision 1.2  2003/07/25 09:04:30  pommier
  *** empty log message ***

  Revision 1.1  2003/05/22 13:18:03  pommier
  regroupement de tous les fichiers dans ./src , et mise en place des RPC

  Revision 1.28  2003/05/05 17:08:51  pommier
  changement de dal::bit_vector::add(i,j) en dal::bit_vector::add(i,nb)

  Revision 1.27  2003/03/14 15:14:46  pommier
  ajout acx_getfem

  Revision 1.26  2003/03/11 15:00:30  pommier
  creation du changelog

  Revision 1.25  2003/02/27 16:34:17  pommier
  improve friendlyness with gcc2.95

  Revision 1.24  2003/02/18 09:36:28  pommier
  update doc

 */


using namespace getfemint;

typedef enum {IS_LAGRANGE, IS_EQUIVALENT, IS_POLYNOMIAL} test_what;
static void
test_fems(test_what what, const getfem::mesh_fem *mf, mexargs_in& in, mexargs_out& out)
{
  dal::bit_vector cvlst;
  bool return_bool = false;
  dal::bit_vector islst;
  if (in.remaining()) cvlst = in.pop().to_bit_vector(&mf->linked_mesh().convex_index());
  else { cvlst = mf->linked_mesh().convex_index(); return_bool = true; }
  for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv) {
    bool it_is = false;
    if (!mf->linked_mesh().convex_index()[cv]) THROW_ERROR( "convex " << cv+1 << " does not exist");
    check_cv_fem(*mf, cv);
    switch (what) { /* hope the compiler will optimize that */
    case IS_LAGRANGE:   it_is = mf->fem_of_element(cv)->is_lagrange(); break;
    case IS_EQUIVALENT: it_is = mf->fem_of_element(cv)->is_equivalent(); break;
    case IS_POLYNOMIAL: it_is = mf->fem_of_element(cv)->is_polynomial(); break;
    }
    if (it_is) {
      islst.add(cv);
    }
  }
  if (return_bool) {
    out.pop().from_integer(islst.card() == mf->linked_mesh().convex_index().card() ? 1 : 0);
  } else {
    out.pop().from_bit_vector(islst);
  }
}

static void
get_fem_of_convexes(const getfem::mesh_fem& mf, mexargs_in& in, mexargs_out& out)
{  
  dal::bit_vector cvlst;
  if (in.remaining()) cvlst = in.pop().to_bit_vector(&mf.linked_mesh().convex_index());
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
      getfem::mesh_fem::ind_dof_face_ct c = mf->ind_dof_of_face_of_element(cv,f); 
    //std::for_each(c.begin(), c.end(), std::bind1st(std::mem_fun(&dal::bit_vector::add),&dof)); // not SGI STL compliant !?
      for (unsigned i=0; i < c.size(); ++i)
	dof.add(c[i]); 
    } else {
      getfem::mesh_fem::ind_dof_ct cvdof = mf->ind_dof_of_element(cv); 
      for (unsigned i=0; i < cvdof.size(); ++i) 
	dof.add(cvdof[i]);
    }
  } 
  return dof;
}

static void
non_conformal_dof(getfem::mesh_fem &mf, mexargs_in &in, mexargs_out &out)
{
  dal::bit_vector cvlst;
  const getfem::mesh &m = mf.linked_mesh();

  std::vector<bgeot::short_type> dcnt(mf.nb_dof());

  if (in.remaining()) cvlst = in.pop().to_bit_vector(&m.convex_index());
  else cvlst = m.convex_index();
  
  for (dal::bv_visitor ic(cvlst); !ic.finished(); ++ic) {
    check_cv_fem(mf, ic);
    for (size_type f = 0; f < m.structure_of_convex(ic)->nb_faces(); f++) {
      bgeot::short_type q;
      if (!m.is_convex_having_neighbour(ic, f)) {
	q = 2;
      } else {
	q = 1;
      }      
      for (size_type i = 0; i < mf.ind_dof_of_face_of_element(ic,f).size(); ++i) {
	dcnt[mf.ind_dof_of_face_of_element(ic,f)[i]]+=q;
      }
    }
  }
  iarray w = out.pop().create_iarray_h(std::count_if(dcnt.begin(), dcnt.end(), 
		     std::bind2nd(std::equal_to<bgeot::short_type>(),1)));
  size_type i,j=0;
  /*
  std::copy_if(dcnt.begin(), dcnt.end(), 
	       std::bind2nd(std::less_equal<bgeot::short_type>(),1)));
  */
  for (i=0; i < dcnt.size(); ++i) {
    if (dcnt[i] == 1) w[j++] = i+config::base_index();
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
  ad.push_back(pmf->nb_dof());
  garray<T> w = out.pop().create_array(ad, T());
  size_type q = u.size() / u.dim(u.ndim()-1);
  assert(w.size() == q * pmf->nb_dof());
  std::vector<unsigned> dofcnt(pmf->nb_dof());

  /* stupid smoothing when the mesh_fem is not discontinuous */
  for (dal::bv_visitor cv(pmf->convex_index()); !cv.finished(); ++cv)
    for (unsigned k=0; k < pmf->nb_dof_of_element(cv); ++k)
      dofcnt[pmf->ind_dof_of_element(cv)[k]]++;


  for (dal::bv_visitor cv(pmf->convex_index()); !cv.finished(); ++cv) {
    for (unsigned k=0; k < pmf->nb_dof_of_element(cv); ++k) {
      size_type d = pmf->ind_dof_of_element(cv)[k];
      for (unsigned j=0; j < q; ++j)
	w[d*q+j] += u[cv*q+j] / T(dofcnt[d]);
    }
  }
}

/*MLABCOM
  FUNCTION [x] = gf_mesh_fem_get(meshfem MF, operation [, args])

  General function for inquiry about mesh_fem objects.

  @RDATTR MESHFEM:GET('nbdof')
  @GET    MESHFEM:GET('dof from cv')
  @GET    MESHFEM:GET('dof from cvid')
  @GET    MESHFEM:GET('non conformal dof')
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
  @GET    MESHFEM:GET('dof on region')
  @GET    MESHFEM:GET('dof nodes')
  @GET    MESHFEM:GET('dof from im')
  @GET    MESHFEM:GET('dof partition')
  @GET    MESHFEM:GET('interpolate_convex_data')
  @GET    MESHFEM:GET('save')
  @GET    MESHFEM:GET('char')
  @GET    MESHFEM:GET('export to vtk')
  @GET    MESHFEM:GET('export to dx')
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


void gf_mesh_fem_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfemint_mesh_fem *mi_mf = in.pop().to_getfemint_mesh_fem();
  getfem::mesh_fem *mf   = &mi_mf->mesh_fem();
  std::string cmd        = in.pop().to_string();
  if (check_cmd(cmd, "nbdof", in, out, 0, 0, 0, 1)) {
    /*@RDATTR I = MESHFEM:GET('nbdof')
      Return the number of Degrees of Freedom (DoF) of the @tmf MF.
      @*/
    out.pop().from_integer(mf->nb_dof());
  } else if (check_cmd(cmd, "dof from cv", in, out, 1, 1, 0, 1)) {
    /*@GET I = MESHFEM:GET('dof from cv', CVLST)
      Return the DoF of the convexes listed in CVLST. 

      WARNING: the Degree of Freedom might be returned in ANY order,
      do not use this function in your assembly routines. Use 'dof
      from cvid' instead, if you want to be able to map a convex
      number with its associated degrees of freedom.

      One can also get the list of dof on a set on convex faces, by
      indicating on the second row of CVLST the faces numbers (with
      respect to the convex number on the first row).
      @*/
    dal::bit_vector dof = get_cv_dof_list(mf, in);

    out.pop().from_bit_vector(dof);
#if 0
  } else if (check_cmd(cmd, "ordered dof from cv", in, out, 1, 1, 0, 1)) {
    /*@GET I = MESHFEM:GET('ordered dof from cv', @ivec CV)
      Return the dof of the convex CV, in the order used in elementary matrices.
      
      Since 2004/10/24, this function is obsolete and shall be removed (replaced by @MESHFEM:GET('dof from cvid')).
      @*/
    size_type cv = in.pop().to_convex_number(mf->linked_mesh());
    check_cv_fem(*mf, cv);
    getfem::mesh_fem::ind_dof_ct cvdof = mf->ind_dof_of_element(cv);
    iarray w = out.pop().create_iarray_h(std::distance(cvdof.begin(), cvdof.end()));
    for (size_type i=0; i < cvdof.size(); i++) {
      w[i] = cvdof[i] + config::base_index();
    }
#endif
  } else if (check_cmd(cmd, "dof from cvid", in, out, 0, 1, 0, 2)) {
    /*@GET  [DOFS,IDX]=MESHFEM:GET('dof from cvid' [,CVLST])
      Return the degrees of freedom attached to each convex of the mesh.\\\\

      If CVLST is omitted, all the convexes will be considered (equivalent to CVLST = 1 ... @MESH:GET('max cvid')).

      IDX is a @MATLAB{row }vector, length(IDX) = length(CVLIST)+1. DOFS is a @MATLAB{row }vector containing the concatenated list of dof of each convex in
      CVLST. Each entry of IDX is the position of the corresponding convex
      point list in DOFS. Hence, for example, the list of points of the
      second convex is @MATLAB{DOFS(IDX(2):IDX(3)-1)}@PYTHON{DOFS[IDX(2):IDX(3)]}.\\

      If CVLST contains convex #id which do not exist in the mesh, their
      point list will be empty.
      @*/
    dal::bit_vector cvlst;
    if (in.remaining()) cvlst = in.pop().to_bit_vector();
    else cvlst.add(0, mf->linked_mesh().convex_index().last_true() + 1);
    
    size_type pcnt = 0;
    /* phase one: count the total number of dof */
    for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv) {
      if (mf->convex_index().is_in(cv)) {
	pcnt += mf->nb_dof_of_element(cv);
      }
    }
    /* phase two: allocation */
    iarray pid = out.pop().create_iarray_h(pcnt);
    bool fill_idx = out.remaining();
    iarray idx; if (fill_idx) idx = out.pop().create_iarray_h(cvlst.card() + 1);

    pcnt = 0;
    size_type cvcnt = 0;
    /* phase three: build the list */
    for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv) {
      if (fill_idx) idx[cvcnt] = pcnt + config::base_index();
      if (mf->convex_index().is_in(cv)) {
	for (getfem::mesh_fem::ind_dof_ct::const_iterator pit = 
	       mf->ind_dof_of_element(cv).begin();
	     pit != mf->ind_dof_of_element(cv).end(); ++pit) {
	  pid[pcnt++] = (*pit) + config::base_index();
	}
      }
      cvcnt++;
    }
    if (fill_idx) idx[idx.size()-1] = pcnt+config::base_index(); /* for the last convex */
  } else if (check_cmd(cmd, "non conformal dof", in, out, 0, 1, 0, 1)) { 
    /*@GET MESHFEM:GET('non conformal dof' [,CVLST])
      Return partially linked degrees of freedom.

      Return the dof located on the border of a convex and which belong
      to only one convex, except the ones which are located on the border
      of the mesh.  For example, if the convex 'a' and 'b' share a common
      face, 'a' has a P1 FEM, and 'b' has a P2 FEM, then the dof on the
      middle of the face will be returned by this function (this can be
      useful when searching the interfaces between classical FEM and
      hierarchical FEM).
      @*/
    non_conformal_dof(*mf, in, out);
  } else if (check_cmd(cmd, "qdim", in, out, 0, 0, 0, 1)) {
    /*@RDATTR MESHFEM:GET('qdim')
      Return the dimension Q of the field interpolated by the mesh_fem.

      By default, Q=1 (scalar field). This has an impact on the DOF numbering.
      @*/
    out.pop().from_integer(mf->get_qdim());
  } else if (check_cmd(cmd, "fem", in, out, 0, 1, 0, 2)) {
    /*@GET [FEMLST, CV2F] = MESHFEM:GET('fem' [, CVLST])
      Return a list of FEM used by the @tmf. 

      FEMLST is an array of all @tfem objects found in the convexes
      given in CVLST. If CV2F was supplied as an output argument, it
      contains, for each convex listed in CVLST, the index of its
      correspounding FEM in FEMLST. 

      Convexes which are not part of the mesh, or convexes which do
      not have any FEM have their correspounding entry in CV2F set to
      -1.
      @*/
    get_fem_of_convexes(*mf, in, out);
  } else if (check_cmd(cmd, "convex_index", in, out, 0, 0, 0, 1)) {
    /*@GET CVLST = MESHFEM:GET('convex_index')
      Return the list of convexes who have a FEM.
      @*/
    out.pop().from_bit_vector(mf->convex_index());
  } else if (check_cmd(cmd, "is_lagrangian", in, out, 0, 1, 0, 1)) {
    /*@GET I = MESHFEM:GET('is_lagrangian' [,CVLST])
      Test if the @tmf is Lagrangian.

      Lagrangian means that each base function Phi[i] is such that<par>
         Phi[i](P[j]) = delta(i,j),<par>
      where P[j] is the DoF location of the jth base function, and
      delta(i,j) = 1 if i==j, else 0.<par>
      
      If CVLST is omitted, it returns 1 if all convexes
      in the mesh are Lagrangian. If CVLST is used, it returns the convex indices (with respect to CVLST) which are Lagrangian.
      @*/
    test_fems(IS_LAGRANGE , mf, in, out);
  } else if (check_cmd(cmd, "is_equivalent", in, out, 0, 1, 0, 1)) {
    /*@GET I = MESHFEM:GET('is_equivalent' [,CVLST])
      Test if the @tmf is equivalent.

      See MESHFEM:GET('is_lagrangian')@*/
    test_fems(IS_EQUIVALENT, mf, in, out);
  } else if (check_cmd(cmd, "is_polynomial", in, out, 0, 1, 0, 1)) {
    /*@GET I = MESHFEM:GET('is_polynomial' [,CVLST])
      Test if all base functions are polynomials.

      See MESHFEM:GET('is_lagrangian')@*/
    test_fems(IS_POLYNOMIAL, mf, in, out);
  } else if (check_cmd(cmd, "dof on boundary", in, out, 1, 1, 0, 1) ||
	     check_cmd(cmd, "dof on region", in, out, 1, 1, 0, 1)) {
    /*@GET DOFLST = MESHFEM:GET('dof on region', RLIST)
      Return the list of dof lying on one of the mesh regions listed in RLIST.

      More precisely, this function returns the DoF whose support is
      non-null on one of regions whose #ids are listed in RLIST
      (note that for boundary regions, some dof nodes may not lie
      exactly on the boundary, for example the dof of PK(n,0) lies on
      the center of the convex, but the base function in not null on
      the convex border). @*/
    iarray bnums = in.pop().to_iarray(-1);
    dal::bit_vector bv;
    for (size_type i=0; i < bnums.size(); ++i) {
      bv |= mf->dof_on_set(bnums[i]);
    }
    out.pop().from_bit_vector(bv);
  } else if (check_cmd(cmd, "dof nodes", in, out, 0, 1, 0, 2)) {
    /*@GET [DOF_XY] = MESHFEM:GET('dof nodes'[, DOFLST])
      Get location of Degrees of Freedom.

      Return the list of interpolation points for the specified dof #IDs in DOFLST
      (if DOFLST is omitted, all DoF are considered).@*/
    dal::bit_vector dof_lst; dof_lst.add(0, mf->nb_dof());
    if (in.remaining())
      dof_lst = in.pop().to_bit_vector(&dof_lst);
    darray w = out.pop().create_darray(mf->linked_mesh().dim(), dof_lst.card());
    size_type j = 0;
    for (dal::bv_visitor dof(dof_lst); !dof.finished(); ++dof, ++j) {
      if (mf->point_of_dof(dof).size() != w.getm() || j >= w.getn()) THROW_INTERNAL_ERROR;
      for (size_type i=0; i < w.getm(); i++) w(i,j)= mf->point_of_dof(dof)[i];
      //      std::copy(mf->point_of_dof(dof).begin(),mf->point_of_dof(dof).end(), &w(0,j));
    }
  } else if (check_cmd(cmd, "dof partition", in, out, 0, 0, 0, 1)) {
    /*@GET [DOFP] = MESHFEM:GET('dof partition')
      Get the 'dof_partition' array.

      Return the array which associates an integer (the partition
      number) to each convex of the mesh_fem. By default, it is an
      all-zero array. The degrees of freedom of each convex of the
      mesh_fem are connected only to the dof of neighbouring convexes
      which have the same partition number, hence it is possible to
      create partially discontinuous mesh_fem very easily.
      @*/
    iarray v = out.pop().create_iarray_h(mf->linked_mesh().convex_index().last_true()+1);
    for (unsigned cv=0; cv < v.size(); ++cv) v[cv] = mf->get_dof_partition(cv);
  } else if (check_cmd(cmd, "save", in, out, 1, 2, 0, 0)) {
    /*@GET MESHFEM:GET('save', filename [, opt])
      Save a @tmf in a text file (and optionaly its linked mesh object if opt is the string 'with_mesh').
      @*/
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
    /*@GET MESHFEM:GET('char' [, opt])
      Output a string description of the mesh_fem. 

      By default, it does not include the description of the linked
      mesh object, except if opt is 'with_mesh' @*/
    std::stringstream s;
    if (in.remaining() && cmd_strmatch(in.pop().to_string(),"with mesh"))
      mf->linked_mesh().write_to_file(s);
    mf->write_to_file(s);
    out.pop().from_string(s.str().c_str());
  } else if (check_cmd(cmd, "linked mesh", in, out, 0, 0, 0, 1)) {
    /*@GET M=MESHFEM:GET('linked mesh')
      Return a reference to the mesh object linked to MF.
      @*/
    out.pop().from_object_id(mi_mf->linked_mesh_id(), MESH_CLASS_ID);
  } else if (check_cmd(cmd, "export to vtk", in, out, 0, -1, 0, 0)) {
    /*@GET MESHFEM:GET('export to vtk', @str FILENAME, ... ['ascii'], U, 'name'...)
      Export a mesh_fem and some fields to a vtk file.

      The FEM and geometric transformations will be mapped to order 1
      or 2 isoparametric PK (or QK) FEMs (as VTK does not handle
      higher order elements). If you need to represent high-order FEMs
      or high-order geometric transformations, you should consider
      @SLICE:GET('export to vtk'). @*/
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
      darray U = in.pop().to_darray(); in.last_popped().check_trailing_dimension(mf2->nb_dof());
      exp.write_point_data(*mf2, U, get_vtk_dataset_name(in, count));
      count+=1;
    }
  } else if (check_cmd(cmd, "export to dx", in, out, 0, -1, 0, 0)) {
    /*@GET MESHFEM:GET('export to dx', @str FILENAME, ... ['as', @str mesh_name][,'edges']['serie',@str serie_name][,'ascii'][,'append'], U, 'name'...)
      Export a mesh_fem and some fields to an OpenDX file.

      This function will fail if the mesh_fem mixes different
      convex types (i.e. quads and triangles), or if OpenDX does not
      handle a specific element type (i.e. prism connections are not
      known by OpenDX).

      The FEM will be mapped to order 1 PK (or QK) FEMs. If you need
      to represent high-order FEMs or high-order geometric
      transformations, you should consider @SLICE:GET('export to dx').
      @*/
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
      darray U = in.pop().to_darray(); in.last_popped().check_trailing_dimension(mf2->nb_dof());
      exp.write_point_data(*mf2, U, get_dx_dataset_name(in));
      if (serie_name.size()) exp.serie_add_object(serie_name);
    }
  } else if (check_cmd(cmd, "dof_from_im", in, out, 1, 2, 0, 1)) {
    /*@GET MESHFEM:GET('dof_from_im', @tmim mim [, @tint P])
      Return a selection of dof who contribute significantly to the
      mass-matrix that would be computed with mf and the integration
      method mim.

      P represents the dimension on what the integration method
      operates (default P = mesh dimension).

      IMPORTANT: you still have to set a valid integration method on
      the convexes which are not crosses by the levelset!
      @*/
    const getfem::mesh_im &mim = *in.pop().to_const_mesh_im();
    int P = -1;
    if (&mim.linked_mesh() != &mf->linked_mesh())
      THROW_BADARG("the mesh_im uses a different mesh");
    if (in.remaining()) 
      P = in.pop().to_integer(1, mim.linked_mesh().dim());
    out.pop().from_bit_vector(getfem::select_dofs_from_im(*mf, mim, P));
  } else if (check_cmd(cmd, "interpolate_convex_data", in, out, 1, 1, 0, 1)) {
    /*@GET U=MESHFEM:GET('interpolate_convex_data', Ucv)

      Interpolate data given on each convex of the mesh to the @tmf
      dof.  The meshfem has to be lagrangian, and should be
      discontinuous (typically a FEM_PK(N,0) or FEM_QK(N,0) should be
      used).

      The last dimension of the input vector Ucv should have
      MESH:GET('max cvid') elements.

      Example of use: MESHFEM:GET('interpolate_convex_data', MESH:GET('quality'))
      @*/
    in.front().check_trailing_dimension(mf->linked_mesh().convex_index().last_true()+1);
    if (in.front().is_complex()) 
      interpolate_convex_data(mf, in.pop().to_darray(), out);
    else interpolate_convex_data(mf, in.pop().to_carray(), out);    
  } else if (check_cmd(cmd, "memsize", in, out, 0, 0, 0, 1)) {
    /*@GET MESHFEM:GET('memsize')
      Return the amount of memory (in bytes) used by the mesh_fem object.

      The result does not take into account the linked mesh object.
      @*/
    out.pop().from_integer(mf->memsize());
  } else if (check_cmd(cmd, "has_linked_mesh_levelset", in, out, 0, 0, 0, 1)) {
    getfem::mesh_fem_level_set *mfls = 
      dynamic_cast<getfem::mesh_fem_level_set*>(mf);
    out.pop().from_integer(mfls ? 1 : 0);
  } else if (check_cmd(cmd, "linked_mesh_levelset", in, out, 0, 0, 0, 1)) {
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
