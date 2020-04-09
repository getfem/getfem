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

#ifdef MAINTAINER_MODE
# include <sys/types.h>
# include <sys/wait.h>
# include <unistd.h>
# include <stdio.h>
# include <errno.h>
#endif
#include <getfemint_misc.h>
#include <getfem/getfem_mesh_fem.h>
#include <getfem/getfem_nonlinear_elasticity.h>
#include <getfem/getfem_plasticity.h>
#include <algorithm>

namespace getfemint {
  gfi_array* checked_gfi_array_create(int ndim, const int *dims, gfi_type_id type, gfi_complex_flag is_complex) {
    if (dims == NULL && ndim != 0) GMM_ASSERT1(false, "");
    gfi_array *t = gfi_array_create(ndim,const_cast<int*>(dims),type,is_complex);
    GMM_ASSERT1(t != NULL, "allocation of " << ndim
                << "-array of " << gfi_type_id_name(type,is_complex)
                << " failed\n");
    return t;
  }

  gfi_array* checked_gfi_array_create_0(gfi_type_id type, gfi_complex_flag is_complex) {
    return checked_gfi_array_create(0,NULL,type,is_complex);
  }

  gfi_array* checked_gfi_array_create_1(int M, gfi_type_id type, gfi_complex_flag is_complex) {
    gfi_array *t = gfi_array_create_1(M,type,is_complex);
    GMM_ASSERT1(t != NULL, "allocation of vector of " << M << " "
                << gfi_type_id_name(type,is_complex) << " failed\n");
    return t;
  }

  gfi_array* checked_gfi_array_create_2(int M, int N, gfi_type_id type, gfi_complex_flag is_complex) {
    gfi_array *t = gfi_array_create_2(M,N,type,is_complex);
    GMM_ASSERT1(t != NULL, "allocation of a " << M << "x" << N << " matrix of "
                << gfi_type_id_name(type,is_complex) << " failed\n");
    return t;
  }

  gfi_array* checked_gfi_array_from_string(const char*s) {
    gfi_array *t = gfi_array_from_string(s);
    GMM_ASSERT1(t != NULL, "allocation of a string of length "
                << strlen(s) << " failed\n");
    return t;
  }

  gfi_array* checked_gfi_create_sparse(int m, int n, int nzmax, gfi_complex_flag is_complex) {
    gfi_array *t = gfi_create_sparse(m,n,nzmax,is_complex);
    GMM_ASSERT1(t != NULL, "allocation of sparse(m=" << m
                << ", n=" << n << ", nzmax=" << nzmax << ") failed\n");
    return t;
  }

  gfi_array *
  convert_to_gfi_sparse(const gf_real_sparse_by_row& smat, double threshold)
  {
    int ni = int(gmm::mat_nrows(smat));
    int nj = int(gmm::mat_ncols(smat));
    unsigned nnz = 0;
    gfi_array *mxA;
    std::vector<int> ccnt(nj);
    std::fill(ccnt.begin(), ccnt.end(), 0);
    std::vector<double> rowmax(ni, 0.0), colmax(nj, 0.0);
    gmm::linalg_traits<gmm::org_type<gmm::linalg_traits<gf_real_sparse_by_row>
       ::const_sub_row_type>::t>::const_iterator it, ite;

    /* first pass : find the maxima / row & column */
    for (int i = 0; i < ni; ++i) {
      it = gmm::vect_const_begin(gmm::mat_const_row(smat,i));
      ite = gmm::vect_const_end(gmm::mat_const_row(smat,i));
      for (; it != ite; ++it) {
        rowmax[i] = std::max(rowmax[i], gmm::abs(*it));
        colmax[it.index()] = std::max(colmax[it.index()], gmm::abs(*it));
      }
    }

    /* second pass : count the significant terms */
    for (int i = 0; i < ni; i++) {
      /* count significatives values */
      it = gmm::vect_const_begin(gmm::mat_const_row(smat,i));
      ite = gmm::vect_const_end(gmm::mat_const_row(smat,i));
      for (; it != ite; ++it) {
        if (*it != 0. && /* (*p).e != 0 => max(rowmax,colmax) != 0 */
            gmm::abs(*it) > threshold * std::max(rowmax[i], colmax[it.index()]) ) {
          ccnt[it.index()]++;
          nnz++;
        }
      }
    }

    mxA = checked_gfi_create_sparse(ni, nj, nnz, GFI_REAL); assert(mxA != NULL);
    double *pr;
    unsigned *ir, *jc;
    pr = gfi_sparse_get_pr(mxA); assert(pr != NULL);
    ir = gfi_sparse_get_ir(mxA); assert(ir != NULL);
    jc = gfi_sparse_get_jc(mxA); assert(jc != NULL); /* dim == nj+1 */

    jc[0] = 0;
    for (int j=0; j < nj; j++) {
      //    mexPrintf("ccnt[%d]=%d\n", j, ccnt[j]);
      jc[j+1] = jc[j] + ccnt[j];
    }
    assert(nnz == jc[nj]);

    //gmm::wsvector<getfem::scalar_type>::const_sorted_iterator its, itse;
    gmm::linalg_traits<gmm::rsvector<getfem::scalar_type> >::const_iterator its, itse;

    /* third pass: filling pr and ir */
    std::fill(ccnt.begin(), ccnt.end(), 0);
    gmm::rsvector<getfem::scalar_type> sorted(gmm::mat_ncols(smat));
    for (int i=0; i < ni; i++) {
      gmm::copy(gmm::mat_const_row(smat, i), sorted);
      its  = gmm::vect_begin(sorted);
      itse = gmm::vect_end(sorted);
      for (; its != itse; ++its) {
        if ((*its) != 0. && /* (*its) != 0 => max(rowmax,colmax) != 0 */
            gmm::abs(*its)/std::max(rowmax[i], colmax[its.index()]) > threshold) {
          ir[jc[its.index()]+ccnt[its.index()]] = i;
          pr[jc[its.index()]+ccnt[its.index()]] = *its;
          ccnt[its.index()]++;
        }
      }
    }
    return mxA;
  }

  gfi_array *
  convert_to_gfi_sparse(const gf_real_sparse_by_col& smat, double threshold)
  {

    int ni = int(gmm::mat_nrows(smat));
    int nj = int(gmm::mat_ncols(smat));
    unsigned nnz = 0;
    gfi_array *mxA;

    std::vector<int> ccnt(nj);
    std::fill(ccnt.begin(), ccnt.end(), 0);
    std::vector<double> rowmax(ni, 0.0), colmax(nj, 0.0);
    gmm::linalg_traits<gmm::org_type<gmm::linalg_traits<gf_real_sparse_by_col>
       ::const_sub_col_type>::t>::const_iterator it, ite;

    /* first pass : find the maxima / row & column */
    for (int j = 0; j < nj; ++j) {
      it = gmm::vect_const_begin(gmm::mat_const_col(smat,j));
      ite = gmm::vect_const_end(gmm::mat_const_col(smat,j));
      for (; it != ite; ++it) {
        rowmax[it.index()] = std::max(rowmax[it.index()], gmm::abs(*it));
        colmax[j] = std::max(colmax[j], gmm::abs(*it));
      }
    }

   /* second pass : count the significant terms */
    for (int j = 0; j < nj; j++) {
      /* count significatives values */
      it = gmm::vect_const_begin(gmm::mat_const_col(smat,j));
      ite = gmm::vect_const_end(gmm::mat_const_col(smat,j));
      for (; it != ite; ++it) {
        if (*it != 0. && /* (*p).e != 0 => max(rowmax,colmax) != 0 */
            gmm::abs(*it) > threshold * std::max(colmax[j], rowmax[it.index()]) ) {
          ccnt[j]++;
          nnz++;
        }
      }
    }

    mxA = checked_gfi_create_sparse(ni, nj, nnz, GFI_REAL); assert(mxA != NULL);
    double *pr;
    unsigned *ir, *jc;
    pr = gfi_sparse_get_pr(mxA); assert(pr != NULL);
    ir = gfi_sparse_get_ir(mxA); assert(ir != NULL);
    jc = gfi_sparse_get_jc(mxA); assert(jc != NULL); /* dim == nj+1 */

    jc[0] = 0;
    for (int j=0; j < nj; j++) {
      //    mexPrintf("ccnt[%d]=%d\n", j, ccnt[j]);
      jc[j+1] = jc[j] + ccnt[j];
    }
    assert(nnz == jc[nj]);

    gmm::linalg_traits<gmm::rsvector<getfem::scalar_type> >::const_iterator its, itse;

    /* third pass: filling pr and ir */
    std::fill(ccnt.begin(), ccnt.end(), 0);

    gmm::rsvector<getfem::scalar_type> sorted(gmm::mat_nrows(smat));
    for (int j=0; j < nj; j++) {
      gmm::copy(gmm::mat_const_col(smat, j), sorted);
      its  = gmm::vect_const_begin(sorted);
      itse = gmm::vect_const_end(sorted);
      for (; its != itse; ++its) {
        if ((*its) != 0. && /* (*its) != 0 => max(rowmax,colmax) != 0 */
            gmm::abs((*its))/std::max(colmax[j], rowmax[its.index()]) > threshold) {
          ir[jc[j]+ccnt[j]] = unsigned(its.index());
          pr[jc[j]+ccnt[j]] = (*its);
          ccnt[j]++;
        }
      }
    }
    return mxA;
  }


  /* merge the edge list of a convex into 'el',
     filtering the edges of the face 'f' if (f != -1) */
  static
  void mesh_edge_list_merge(const getfem::mesh &m, bgeot::edge_list &el, const bgeot::edge_list &elcv, int cv, int f)
  {
    bgeot::edge_list::const_iterator it = elcv.begin();
    for (it = elcv.begin(); it != elcv.end(); it++) {
      if (f != -1) {
        bgeot::mesh_structure::ind_pt_face_ct
          pt = m.ind_points_of_face_of_convex(cv, short_type(f));
        if (std::find(pt.begin(), pt.end(), (*it).i) == pt.end()) continue;
        if (std::find(pt.begin(), pt.end(), (*it).j) == pt.end()) continue;
      }
      el.add(*it);
    }
  }

  /* used by gf_mesh_get and gf_mesh_fem_get */
  void
  build_edge_list(const getfem::mesh &m, bgeot::edge_list &el, mexargs_in &in) {
    iarray v;
    bool all_cv = true;
    bool merge_convex = false;

    if (in.remaining() && !in.front().is_string()) {
      v = in.pop().to_iarray(-1, -1);
      all_cv = false;
    }
    if (in.remaining() && in.front().is_string()) {
      std::string s = in.pop().to_string();
      if (cmd_strmatch(s,"merge convex") || cmd_strmatch(s,"merge"))
        merge_convex = true;
      else bad_cmd(s);
    }
    if (all_cv) { /* the fastest way : edges off all convexes */
      bgeot::mesh_edge_list(m, el, merge_convex);
    } else {              /* slower (and not very smart) way : edges from
                             selected convexes/faces */
      if (v.getm() != 1 && v.getm() != 2)
        THROW_ERROR("wrong number of rows for CVLIST");

      size_type j = 0;
      /* loop over the rows -- if face numbers are given, the performance will
         be better if the columns are grouped by convex numbers */
      while (j < v.getn()) {
        size_type cv = size_type(v(0,unsigned(j))-config::base_index());
        bgeot::edge_list elcv;

        if (!m.convex_index().is_in(cv))
          THROW_ERROR("can't build edges of convex " << cv + config::base_index()
                       << ": there is no such convex in mesh");
        std::vector<size_type> cvpt(m.ind_points_of_convex(cv).begin(), m.ind_points_of_convex(cv).end());
        bgeot::mesh_edge_list_convex(m.structure_of_convex(cv), cvpt, cv, elcv, merge_convex);
        if (v.getm() == 2) { /* face numbers present */
          /* loop while the convex number do not change */
          do {
            mesh_edge_list_merge(m, el, elcv, int(cv), int(v(1,unsigned(j))-config::base_index()));
            j++;
          } while (j < v.getn() && size_type(v(0,unsigned(j))-config::base_index()) == cv);
        } else { /* merge all edges of the convex */
          mesh_edge_list_merge(m,el, elcv, int(cv), -1);
          j++;
        }
      }
    }
  }

  /* build a vector of couples <convex number, face number> */
  void build_convex_face_lst(const getfem::mesh& m, std::vector<convex_face>& l, const iarray *v) {
    l.resize(0);
    if (v) {
      if (v->getm() != 1 && v->getm()!=2) THROW_ERROR("too much rows (2 max)");
      l.resize(v->getn());
      for (unsigned i=0; i < v->getn(); ++i) {
        l[i].cv = size_type(v->operator()(0,i))-config::base_index();
        if (!m.convex_index()[l[i].cv])
          THROW_ERROR("the convex " << l[i].cv+config::base_index() << " is not part of the mesh");
        if (v->getm() == 2) {
          l[i].f  = dim_type(v->operator()(1,i)-config::base_index());
          if (dim_type(l[i].f+1) && l[i].f >=
              m.structure_of_convex(l[i].cv)->nb_faces())
            THROW_ERROR("face " << l[i].f+config::base_index() << " of convex " << l[i].cv+config::base_index() << "("
                        << bgeot::name_of_geometric_trans(m.trans_of_convex(l[i].cv))
                        << ") does not exist");
        }
        else l[i].f = dim_type(-1);
      }
    } else {
      l.reserve(m.convex_index().card());
      for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv)
        l.push_back(convex_face(cv,dim_type(-1)));
    }
  }

  getfem::mesh_region
  to_mesh_region(const iarray &v) {
    getfem::mesh_region rg;
    if (v.getm() != 1 && v.getm()!=2) THROW_ERROR("too much rows for mesh_region description (2 max)");
    for (unsigned i=0; i < v.getn(); ++i) {
      size_type cv = size_type(v(0,i))-config::base_index();
      short_type f  = short_type(-1);
      if (v.getm() == 2)
        f  = short_type(v(1,i)-config::base_index());
      rg.add(cv,f);
    }
    return rg;
  }

  getfem::mesh_region
  to_mesh_region(const getfem::mesh& m, const iarray *v) {
    if (v) {
      getfem::mesh_region rg = to_mesh_region(*v);
      /* check that the region is ok wrt m */
      for (getfem::mr_visitor i(rg); !i.finished(); ++i) {
        size_type cv = i.cv();
        if (!m.convex_index()[cv])
          THROW_ERROR("the convex " << cv+config::base_index() << " is not part of the mesh");
        if (i.is_face()) {
          short_type f  = i.f();
          if (f >= m.structure_of_convex(cv)->nb_faces())
            THROW_ERROR("face " << f+config::base_index() << " of convex " << cv+config::base_index() << "("
                        << bgeot::name_of_geometric_trans(m.trans_of_convex(cv))
                        << ") does not exist");
        }
      }
      return rg;
    } else {
      return getfem::mesh_region(m.convex_index());
    }
  }

  /*
    interpolate the solution (in vector U) on the element cv,
    evaluated at the points 'pt' (given in the reference convex)
  */
  void interpolate_on_convex_ref(const getfem::mesh_fem *mf,
                                 getfem::size_type cv,
                                 const std::vector<getfem::base_node> &pt,
                                 const darray& U,
                                 getfem::base_matrix &pt_val)
  {
    assert(mf->convex_index().is_in(cv));
    assert(!mf->is_reduced());
    getfem::pfem cv_fem = mf->fem_of_element(cv);
    dim_type qdim = mf->get_qdim();
    /* largely inspired by getfem_export.h */

    /* interpolation of the solution.                                  */
    /* faux dans le cas des éléments non tau-equivalents ou vectoriel. */
    if (cv_fem->target_dim() != 1)
      THROW_ERROR("interpolation on vector fem is still to be done! "
                  "(or at least to be tested...)");
    if (U.getn() != mf->nb_dof())
      THROW_ERROR("wrong nb of columns for U");
    assert(cv_fem->is_equivalent());

    pt_val.resize(qdim*U.getm(), pt.size());

    getfem::base_matrix G;

    if (mf->fem_of_element(cv)->need_G())
      bgeot::vectors_to_base_matrix(G, mf->linked_mesh().points_of_convex(cv));
    getfem::base_vector coeff(mf->nb_basic_dof_of_element(cv));
    getfem::base_vector val(qdim);

    getfem::fem_interpolation_context
      ctx(mf->linked_mesh().trans_of_convex(cv), cv_fem,
          getfem::base_node(), G, cv);
    for (size_type row = 0; row < U.getm(); ++row) {
      for (size_type j = 0; j < coeff.size(); j++)
        coeff[j] = U(unsigned(row), unsigned(mf->ind_basic_dof_of_element(cv)[j]));
      for (size_type j = 0; j < pt.size(); ++j) {
        ctx.set_xref(pt[j]);
        cv_fem->interpolation(ctx, coeff, val, qdim);
        for (size_type q = 0; q < qdim; ++q)
          pt_val(row*qdim + q,j) = val[q];
      }
    }
  }

  /* utility function for eval_on_triangulated_surface */
  static void
  eval_sub_nodes(unsigned N, const std::vector<getfem::base_node>& V,
                 std::vector<getfem::base_node>& spt) {
    assert(N>0);
    spt.resize(((N+1)*(N+2))/2);

    /*
      .      pt[0]
      .        /\
      .       /  \
      . pt[2]/____\pt[1]

      refinment:

      .         0     <- layer 0
      .        111    <- layer 1 etc..
      .       22222

    */

    spt[0] = V[0];
    size_type pcnt = 1;

    /* find the three nodes of the each sub-triangle */
    for (size_type layer = 1; layer <= N; layer++) {
      getfem::base_node A,B;
      getfem::scalar_type c;
      c = ((getfem::scalar_type)layer)/N;

      A = V[0] + (V[2] - V[0]) * c;
      B = V[0] + (V[1] - V[0]) * c;

      for (size_type inode = 0; inode <= layer; inode++, pcnt++) {
        spt[pcnt] = A + (B-A) * (scalar_type(inode)/scalar_type(layer));
      }
    }
    if (!(pcnt == spt.size())) THROW_INTERNAL_ERROR;
  }

  /* U may have zero rows, in that case, only the geometric deformation
     of the points will be computed and stored in w
  */
  static void
  add_refined_tri(const getfem::mesh *mesh, size_type cv,
                  const std::vector<getfem::base_node>& vertices, int Nrefine,
                  darray& w, size_type tri_cnt,
                  const getfem::mesh_fem *pmf, const darray& U)
  {
    unsigned mesh_dim = mesh->dim();
    unsigned qdim = (pmf) ? pmf->get_qdim() : 0;
    bgeot::pgeometric_trans pgt = mesh->trans_of_convex(cv);
    std::vector<getfem::base_node> pt(Nrefine * Nrefine);

    eval_sub_nodes(Nrefine, vertices, pt);

    getfem::base_matrix pt_val;
    if (pmf) {
      if (!pmf->convex_index()[cv])
	THROW_ERROR( "convex " << cv+config::base_index() << " has no FEM");
      pt_val.resize(qdim*U.getm(), pt.size());
      interpolate_on_convex_ref(pmf, cv, pt, U, pt_val);
    }

    /* apply the geometric transformation to the points, in order to
       find their real location on the mesh */
    //getfem::base_node P(mesh_dim);
    for (size_type i=0; i < pt.size(); ++i) {
      pt[i] = pgt->transform(pt[i], mesh->points_of_convex(cv));
      /*      P.fill(0.0);
      for (getfem::size_type j = 0; j < pgt->nb_points(); ++j) {
        P.addmul(pgt->poly_vector()[j].eval(pt[i].begin()),
                 mesh->points_of_convex(cv)[j]);
                 }*/
      //pt[i] = P;
    }

    size_type refined_tri_cnt = 0;
    /* find the three nodes of the each sub-triangle */
    for (int layer = 0; layer < Nrefine; layer++) {
      for (int itri = 0; itri < layer*2+1; itri++, refined_tri_cnt++) {
        getfem::size_type n[3];

        if ((itri & 1) == 0) {
          /*
            .           0
            .          /\
            .       2 /__\ 1
          */
          n[0] = (layer*(layer+1))/2 + itri/2;
          n[1] = n[0] + layer+1;
          n[2] = n[1]+1;
        } else {
          /*
            .       1 ____ 2
            .         \  /
            .          \/
            .           0
          */
          n[1] = (layer*(layer+1))/2 + itri/2;
          n[2] = n[1]+1;
          n[0] = n[2]+layer+1;
        }

        if (!(n[0] < pt.size()) || !(n[1] < pt.size()) || !(n[2] < pt.size()))
          { THROW_INTERNAL_ERROR; }
        for (unsigned ipt = 0; ipt < 3; ++ipt) {
          for (unsigned idim = 0; idim < mesh_dim; ++idim) {
            w(unsigned(ipt * mesh_dim + idim),
              unsigned(tri_cnt+refined_tri_cnt)) = pt[n[ipt]][idim];
          }
          if (pmf) {
            for (size_type row = 0; row < U.getm(); ++row) {
              for (size_type q = 0; q < qdim; ++q) {
                w(unsigned(3*mesh_dim + row*qdim*3 + ipt*qdim + q),
                  unsigned(tri_cnt + refined_tri_cnt)) = pt_val(q,n[ipt]);
              }
            }
          }
        }
      }
    }
  }

/*
   interpolate a field U defined on mf onto a P1 discontiuous mesh should work
   for planar segment/triangle/quadrangle meshes, and faces of volumic meshes
   which are a segment/triangle/quadrangle. In order to be able to represent
   accurately high order FEMs, a refinment factor can be specified.

   the aim is to be fast (the matlab graphic routines don't need any additional
   slowness...)

   This function can still be used if pmf is NULL, in order to obtain a nice
   decomposition of a mesh into triangle (interesting for non-linear geometric
   deformations, used by gf_plot_mesh)
*/  void
  eval_on_triangulated_surface(const getfem::mesh* mesh, int Nrefine,
                               const std::vector<convex_face>& cvf,
                               mexargs_out& out,
                               const getfem::mesh_fem *pmf, const darray& U) {
    unsigned mesh_dim = mesh->dim();
    unsigned qdim = (pmf ? pmf->get_qdim() : 0);
    if (mesh_dim < 2 || mesh_dim > 3) {
      THROW_ERROR( "This function do not handle " <<
                   mesh_dim << "D meshes (only 2D or 3D)");
    }

    /* first pass: count nb of triangles */
    size_type nb_tri = 0;
    for (size_type i=0; i < cvf.size(); ++i) {
      bgeot::pconvex_ref cv_ref=mesh->trans_of_convex(cvf[i].cv)->convex_ref();

      bgeot::pconvex_structure cv_struc=cv_ref->structure();

      if (cvf[i].f != dim_type(-1)) {
        cv_struc = cv_struc->faces_structure()[cvf[i].f];
      }
      if (bgeot::basic_structure(cv_struc)->nb_points() == 2)  /* pas les lignes */
        continue;
      if (cv_struc->dim() > 2)
        THROW_ERROR("cannot draw a 3D convex (convex nb " << cvf[i].cv
                    << "), please specify "
                    "faces (see gf_mesh_get(m,'outer faces') for example)");
      size_type t_inc = 0;
      switch (bgeot::basic_structure(cv_struc)->nb_points()) {
      case 3: t_inc = 1; break;
      case 4: t_inc =2; break;
      }
      if (pmf == NULL && mesh->trans_of_convex(cvf[i].cv)->is_linear()) {
        nb_tri += t_inc;
      } else {
        nb_tri += Nrefine*Nrefine*t_inc;
      }
    }

    darray w = out.pop().create_darray(3*mesh_dim + 3*qdim*U.getm(),
                                       unsigned(nb_tri));

    /* second pass */
    size_type tri_cnt = 0;
    std::vector<size_type> ipts;
    for (size_type i=0; i < cvf.size(); ++i) {
      std::vector<getfem::base_node> vertices(3);
      size_type n_refine;
      bgeot::pconvex_ref cv_ref=mesh->trans_of_convex(cvf[i].cv)->convex_ref();
      bgeot::pconvex_structure cv_struc=cv_ref->structure();

      if (pmf == NULL && mesh->trans_of_convex(cvf[i].cv)->is_linear())
        n_refine = 1;
      else n_refine = Nrefine;

      if (cvf[i].f == dim_type(-1)) {
        ipts.resize(cv_ref->nb_points());
        for (size_type j=0; j < ipts.size(); ++j) ipts[j]=j;
      } else {
        const bgeot::convex_ind_ct&
          pts = cv_struc->ind_points_of_face(short_type(cvf[i].f));
        ipts.resize(pts.size());
        std::copy(pts.begin(), pts.end(), ipts.begin());
        cv_struc = cv_struc->faces_structure()[cvf[i].f];
      }
      if (bgeot::basic_structure(cv_struc)->nb_points() == 2) /* pas les lignes */
        continue;
      /* for each point, count the nb of faces it belongs */
      std::vector<short_type> fcnt(cv_struc->nb_points(), 0);
      for (short_type f = 0; f < cv_struc->nb_faces(); ++f) {
        for (short_type fp = 0; fp < cv_struc->nb_points_of_face(f); ++fp) {
          fcnt[cv_struc->ind_points_of_face(f)[fp]]++;
        }
      }
      /* remove unused points */
      size_type j=0;
      for (size_type k=0; k < ipts.size(); ++k) {
        if (fcnt[k]==2) {
          ipts[j] = ipts[k];
          j++;
        }
      }
      ipts.resize(j);

#define DO_TRIANGLE(a,b,c) {                                                  \
        vertices[0] = cv_ref->points()[ipts[a]];                              \
        vertices[1] = cv_ref->points()[ipts[b]];                              \
        vertices[2] = cv_ref->points()[ipts[c]];                              \
        add_refined_tri(mesh, cvf[i].cv, vertices,                            \
                        int(n_refine),                                        \
                        w, tri_cnt, pmf, U); tri_cnt += n_refine*n_refine; }

      if (ipts.size() == 3) {
        DO_TRIANGLE(0,1,2);
      } else if (ipts.size() == 4) {
        DO_TRIANGLE(0,1,2);
        DO_TRIANGLE(1,3,2);
      } else {
        cerr << "convex not handled by eval_on_triangulated_surface: "
             << bgeot::name_of_geometric_trans(mesh->trans_of_convex(cvf[i].cv)) << endl;
        //cerr << "basic_structure(cv_struc)->nb_points() = "
        //     << bgeot::basic_structure(cv_struc)->nb_points() << endl;
        //        GMM_ASSERT1(false, "");
      }
      assert(tri_cnt <= nb_tri);
    }
    GMM_ASSERT1(tri_cnt == nb_tri, "tri_cnt=" << tri_cnt
                << ", nb_tri=" << nb_tri);
  }


  /* because it's such a pain to launch matlab with a debugger.. */
  void attach_gdb() {
#   ifdef MAINTAINER_MODE
    pid_t pid, parent;
    char *cmd;
    char cmd_pid[1024];
    cmd = ::getenv("GETFEM_DEBUG_CMD"); /* i.e. "xterm -e gdb /path/to/bin/matlab %d"
                                           (xterm because matlab does weird things with its tty (no echo of input..))
                                        */
    if (cmd == NULL) {
      cerr << "the environment variable $GETFEM_DEBUG_CMD is not set, "
           << "cancelling the debugger\n"
           << "You should set it to something like \"xterm -e "
           << "gdb /path/to/bin/matlab %d\" (%d will be replaced by the "
           << "matlab pid)\n";
      return;
    }
    parent = ::getpid();
    ::snprintf(cmd_pid, 1024, cmd, parent);
    std::cerr << "internal error detected\nrunning /bin/sh -c \""
              << cmd_pid << "\" ...\n" << std::endl;
    switch ((pid=::fork())) {
    case -1: std::cerr << "can't fork gdb\n"; return;
    case 0: /* child */
      {
        execlp("/bin/sh", "/bin/sh", "-c", cmd_pid, NULL);
        cerr << "execlp failure .. : " << ::strerror(errno) << endl;
        exit(1); /* in case the exec failed */
      } break;
    default:
      {
        int status;
        std::cout << "When you have finished with gdb, "
                  << "just quit it (type q)\n";
        ::waitpid(pid, &status, 0);
        std::cout << "gdb done!\n";
      } break;
    }
#   endif /*MAINTAINER_MODE*/
  }


  const getfem::phyperelastic_law &
  abstract_hyperelastic_law_from_name(const std::string &lawname,
                                      size_type N) {
    static getfem::phyperelastic_law SVK_AHL =
      std::make_shared<getfem::SaintVenant_Kirchhoff_hyperelastic_law>();
    static getfem::phyperelastic_law IMR_AHL =
      std::make_shared<getfem::Mooney_Rivlin_hyperelastic_law>(false,false);
    static getfem::phyperelastic_law CMR_AHL =
      std::make_shared<getfem::Mooney_Rivlin_hyperelastic_law>(true,false);
    static getfem::phyperelastic_law INH_AHL =
      std::make_shared<getfem::Mooney_Rivlin_hyperelastic_law>(false,true);
    static getfem::phyperelastic_law CNH_AHL =
      std::make_shared<getfem::Mooney_Rivlin_hyperelastic_law>(true,true);
    static getfem::phyperelastic_law NHB_AHL =
      std::make_shared<getfem::Neo_Hookean_hyperelastic_law>(true);  // Bonet
    static getfem::phyperelastic_law NHC_AHL =
      std::make_shared<getfem::Neo_Hookean_hyperelastic_law>(false); // Ciarlet
    static getfem::phyperelastic_law CG_AHL =
      std::make_shared<getfem::Ciarlet_Geymonat_hyperelastic_law>();
    static getfem::phyperelastic_law GBK_AHL =
      std::make_shared<getfem::generalized_Blatz_Ko_hyperelastic_law>();
    static getfem::phyperelastic_law PS_SVK_AHL =
      std::make_shared<getfem::plane_strain_hyperelastic_law>(SVK_AHL);
    static getfem::phyperelastic_law PS_IMR_AHL =
      std::make_shared<getfem::plane_strain_hyperelastic_law>(IMR_AHL);
    static getfem::phyperelastic_law PS_CMR_AHL =
      std::make_shared<getfem::plane_strain_hyperelastic_law>(CMR_AHL);
    static getfem::phyperelastic_law PS_INH_AHL =
      std::make_shared<getfem::plane_strain_hyperelastic_law>(INH_AHL);
    static getfem::phyperelastic_law PS_CNH_AHL =
      std::make_shared<getfem::plane_strain_hyperelastic_law>(CNH_AHL);
    static getfem::phyperelastic_law PS_NHB_AHL =
      std::make_shared<getfem::plane_strain_hyperelastic_law>(NHB_AHL);
    static getfem::phyperelastic_law PS_NHC_AHL =
      std::make_shared<getfem::plane_strain_hyperelastic_law>(NHC_AHL);
    static getfem::phyperelastic_law PS_CG_AHL =
      std::make_shared<getfem::plane_strain_hyperelastic_law>(CG_AHL);
    static getfem::phyperelastic_law PS_GBK_AHL =
      std::make_shared<getfem::plane_strain_hyperelastic_law>(GBK_AHL);

    if (cmd_strmatch(lawname, "SaintVenant Kirchhoff") ||
        cmd_strmatch(lawname, "svk"))
      { if (N == 2) return PS_SVK_AHL; else return SVK_AHL; }

    if (cmd_strmatch(lawname, "Mooney Rivlin") ||
        cmd_strmatch(lawname, "mr") ||
        cmd_strmatch(lawname, "incompressible Mooney Rivlin") ||
        cmd_strmatch(lawname, "imr"))
      { if (N == 2) return PS_IMR_AHL; else return IMR_AHL; }

    if (cmd_strmatch(lawname, "compressible Mooney Rivlin") ||
        cmd_strmatch(lawname, "cmr"))
      { if (N == 2) return PS_CMR_AHL; else return CMR_AHL; }

    if (cmd_strmatch(lawname, "neo Hookean") ||
        cmd_strmatch(lawname, "nh") ||
        cmd_strmatch(lawname, "compressible neo Hookean") ||
        cmd_strmatch(lawname, "cnh"))
      { if (N == 2) return PS_CNH_AHL; else return CNH_AHL; }

    if (cmd_strmatch(lawname, "incompressible neo Hookean") ||
        cmd_strmatch(lawname, "inh"))
      { if (N == 2) return PS_INH_AHL; else return INH_AHL; }

    if (cmd_strmatch(lawname, "neo Hookean Bonet") ||
        cmd_strmatch(lawname, "nhb"))
      { if (N == 2) return PS_NHB_AHL; else return NHB_AHL; }

    if (cmd_strmatch(lawname, "neo Hookean Ciarlet") ||
        cmd_strmatch(lawname, "nhc"))
      { if (N == 2) return PS_NHC_AHL; else return NHC_AHL; }

    if (cmd_strmatch(lawname, "Ciarlet Geymonat") ||
        cmd_strmatch(lawname, "cg"))
      { if (N == 2) return PS_CG_AHL; else return CG_AHL; }

    if (cmd_strmatch(lawname, "generalized Blatz Ko") ||
        cmd_strmatch(lawname, "gbk"))
      { if (N == 2) return PS_GBK_AHL; else return GBK_AHL; }


    THROW_BADARG(lawname <<
                 " is not the name of a known hyperelastic law. \\"
                 "Valid names are: SaintVenant Kirchhoff, Mooney Rivlin, "
                 "neo Hookean or Ciarlet Geymonat");
    return SVK_AHL;
  }



  /** This function return the right projection type chosen which could only
      be for the moment the Von Mises projection. */
  const getfem::pconstraints_projection &
  abstract_constraints_projection_from_name(const std::string &projname) {

    static getfem::pconstraints_projection
      VM_proj = std::make_shared<getfem::VM_projection>();

    if (cmd_strmatch(projname, "Von Mises") ||
        cmd_strmatch(projname, "VM")) return VM_proj;

    THROW_BADARG(projname <<
                 " is not the name of a known constraints projection. \\"
                 "Valid names are: Von mises or VM");

    return VM_proj;
  }




  /*
  bool interruptible_iteration::finished(double nr) {
  if (is_cancel_flag_set()) {
      //GFI_WARNING("Interrupted by Ctrl-C");
      THROW_INTERRUPTED();
    }
    return gmm::iteration::finished(nr);
  }
  */

  static void ctrl_c_iteration_callback(const gmm::iteration &) {
    if (is_cancel_flag_set()) {
      THROW_INTERRUPTED();
    }
  }

  interruptible_iteration::interruptible_iteration(double r) :
    gmm::iteration(r) {
    gmm::iteration::set_callback(ctrl_c_iteration_callback);
  }
}

