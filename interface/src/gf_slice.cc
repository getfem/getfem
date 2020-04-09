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
// $Id$

#include <memory>
#include <getfem/getfem_mesh_level_set.h>
#include <getfem/getfem_mesh_slice.h>
#include <getfem/getfem_mesh_fem.h>
#include <getfemint_misc.h>
#include <getfemint_workspace.h>


using namespace getfemint;

namespace getfem {
  class mesh_slice_streamline : public stored_mesh_slice {
    scalar_type EPS;
  public:
    mesh_slice_streamline(mesh_slice_cv_dof_data_base *mfU,
                          std::vector<base_node>& seeds,
                          bool forward, bool backward) : EPS(1e-10) {

      poriginal_mesh = &mfU->pmf->linked_mesh();
      const mesh &ml = mfU->pmf->linked_mesh();
      bgeot::geotrans_inv gti;
      std::vector<base_node> ref_pts(seeds.size());
      std::vector<size_type> indpts(seeds.size());
      gti.add_points(seeds);
      cv2pos.clear(); cv2pos.resize(ml.convex_index().last_true() + 1, size_type(-1));
      for (dal::bv_visitor cv(ml.convex_index()); !cv.finished(); ++cv) {
        size_type nbpt;
        if ((nbpt=gti.points_in_convex(ml.convex(cv), ml.trans_of_convex(cv), ref_pts, indpts))) {
          //infomsg() << nbpt << " seeds found in convex " << cv << ": " << indpts << endl;
          for (size_type i=0; i < nbpt; ++i) {
            if (forward)
              extract_streamline(mfU, cv, seeds[indpts[i]], ref_pts[i], +1);
            if (backward)
              extract_streamline(mfU, cv, seeds[indpts[i]], ref_pts[i], -1);
          }
        }
      }
    }
  private:
    /*
      computes next node with runge-kutta.
       return 0 if next node is on a face of the convex,
       +1 if it falls outside and -1 if it is inside
    */
    int do_runge_kutta(bgeot::geotrans_inv_convex& gti, size_type cv,
                       const base_matrix& G, pfem pf, bgeot::pgeometric_trans pgt,
                       const base_vector& coeff, const base_node& P0,
                       const base_node& refP0, scalar_type h,
                       base_node& P1, base_node& refP1) {
      fem_interpolation_context ctx(pgt,pf,refP0,G,cv);
      base_node k1(P0.size());
      pf->interpolation(ctx,coeff,k1,refP0.size());
      //cerr << "do_runge_kutta: P0=" << P0 << ", refP0=" << refP0 << ", coeff=" << coeff << ", h=" << h << ", k1=" << k1 << endl;
      P1 = P0+k1*(h/2);
      gti.invert(P1, refP1);
      scalar_type in1 = pgt->convex_ref()->is_in(refP1);

      //cerr << ", P1=" << P1 << ", refP1=" << refP1 << ", in1=" << in1 << endl;
      if (gmm::abs(in1) < EPS) return 0;
      else if (in1 > 0) return +1;
      else {
        base_node k2(P0.size()); ctx.set_xref(refP1);
        pf->interpolation(ctx,coeff,k2,k2.size());
        P1 = P0+k2*h;
        gti.invert(P1, refP1);
        in1 = pgt->convex_ref()->is_in(refP1);
        //cerr << "do_runge_kutta2: P1=" << P1 << ", refP1=" << refP1 << ", in1=" << in1 << endl;
        if (gmm::abs(in1) < EPS) return 0;
        else if (in1 > 0) return +1;
      }
      return -1;
    }

    size_type find_convex_of_point(const mesh& ml, size_type cv, const base_node& P, base_node& refP, bgeot::geotrans_inv_convex& gti) {
      /* find on which face is the point (approximately) */
      dim_type f = dim_type(-1);
      scalar_type best_f = 1e10;
      size_type cnt = 0;
      for (short_type i=0; i < ml.structure_of_convex(cv)->nb_faces(); ++i) {
        scalar_type v = gmm::abs(ml.trans_of_convex(cv)->convex_ref()->is_in_face(i,refP));
        cnt++;
        if (v < best_f || cnt == 0) { best_f = v; f = dim_type(i); }
      }

      /* look for the other convex sharing this face */
      bgeot::mesh_structure::ind_set clst;
      ml.neighbors_of_convex(cv, f, clst);
      size_type best = size_type(-1); scalar_type best_v = 1e10;
      cnt = 0;
      for (bgeot::mesh_structure::ind_set::const_iterator it = clst.begin();
           it != clst.end(); ++it) {
        if (*it != cv && ml.structure_of_convex(*it)->dim() == ml.dim()) {
          ++cnt;
          gti.init(ml.convex(*it).points(), ml.trans_of_convex(*it));
          gti.invert(P, refP);
          scalar_type v = ml.trans_of_convex(*it)->convex_ref()->is_in(refP);
          if (v < best_v || cnt == 0) { best_v = v; best = *it; }
        }
      }
      if (cnt == 0) return size_type(-1); else return best;
    }

    void extract_streamline(mesh_slice_cv_dof_data_base *mfU, size_type cv,
                            const base_node& seed, const base_node& seed_ref, double dir) {
      getfem::mesh_slicer::cs_nodes_ct snodes;
      getfem::mesh_slicer::cs_simplexes_ct ssimplexes;
      bool change_convex = true, first=true;
      dim_type sdim = mfU->pmf->linked_mesh().dim();
      base_node P0(seed), refP0(seed_ref), P1, refP1;
      scalar_type h = 0;
      pfem pf = 0;
      bgeot::pgeometric_trans pgt = 0;
      base_vector coeff;
      base_matrix G;
      bgeot::geotrans_inv_convex gti(mfU->pmf->linked_mesh().convex(cv), mfU->pmf->linked_mesh().trans_of_convex(cv));
      bool store_convex_and_stop = false;
      GMM_ASSERT1(sdim == mfU->pmf->get_qdim(), "won't compute streamline "
                  "of a field whose dimension is not equal to mesh dimension");
      size_type cnt = 0;
      do {
        if (change_convex) {
          //cerr << "changement de convexe, ancien=" << cv << endl;
          /* init convex-dependent data */
          if (!first) {

            cv = find_convex_of_point(mfU->pmf->linked_mesh(), cv, P0, refP0, gti);
            //cerr << "nouveau convexe: " << cv << endl;
            if (cv == size_type(-1)) {
              //cerr << "fin au convex " << cv << endl;
              break;
            }
          }
          first = false;
          pf = mfU->pmf->fem_of_element(cv);
          pgt = mfU->pmf->linked_mesh().trans_of_convex(cv);
          mfU->copy(cv, coeff); if (dir<0) gmm::scale(coeff, -1);
          /* get convex size estimation */
          h = 1e10;
          for (size_type i=0; i < pgt->nb_points(); ++i) {
            for (size_type j=0; j < pgt->nb_points(); ++j) {
              if (j!=i) {
                h = std::min(h, gmm::vect_dist2_sqr(mfU->pmf->linked_mesh().points_of_convex(cv)[i],
                                                      mfU->pmf->linked_mesh().points_of_convex(cv)[j]));
              }
            }
          }
          h = sqrt(h);
          scalar_type z=gmm::vect_norminf(coeff);
          if (z > EPS)
            h = h/(10*z);
          else h=0; /* on va s'arreter */

          snodes.resize(1); ssimplexes.resize(0);
          snodes[0].faces.reset(); snodes[0].pt = P0; snodes[0].pt_ref = refP0;
          change_convex = false;
        }
        int rk = do_runge_kutta(gti, cv, G, pf, pgt, coeff, P0, refP0, h, P1, refP1);

        if (gmm::vect_dist2_sqr(P0,P1) < 1e-16) { store_convex_and_stop = true; } // we don't move anymore

        /* basic dichotomy on h if the step is to big */
        if (rk > 0) {
          //cerr << "->>> on va sortir du convexe, debut dichotomie\n";
          change_convex = true;
          scalar_type h0 = 0, h1 = 1;
          while (h1 - h0 > 1e-7*h) {
            rk = do_runge_kutta(gti, cv, G, pf, pgt, coeff, P0, refP0, (h0+h1)/2, P1, refP1);
            if (rk == 0) break;
            else if (rk > 0) h1 = (h0 + h1)/2;
            else h0 = (h0 + h1)/2;
          }
          //cerr << "fin dichotomie, h0=" << h0 << ", h1=" << h1 << ", rk=" << rk << endl;
        }

        //cerr << "AJOUT: cv = " << cv << ": P0=" << P0 << " - P1=" << P1 << endl;

        /* add the segment */
        snodes.push_back(slice_node(P1,refP1));
        snodes.back().faces.reset();
        ssimplexes.push_back(slice_simplex(2));
        ssimplexes.back().inodes[0] = snodes.size()-2;
        ssimplexes.back().inodes[1] = snodes.size()-1;

        if (++cnt > 3000) { infomsg() << "too much iterations for streamline extraction\n"; store_convex_and_stop = true; }

        if (change_convex || store_convex_and_stop) {
          /* store streamline of previous convex */
          dal::bit_vector splx_in; splx_in.add(0, ssimplexes.size());
          set_convex(cv, pgt->convex_ref(), snodes, ssimplexes,
                     dim_type(pgt->convex_ref()->structure()->nb_faces()),
                     splx_in, false);
        }
        P0 = P1; refP0 = refP1;
      } while (!store_convex_and_stop);
    }
  };
}

static getfem::slicer_action*
build_slicers(const getfem::mesh& m, std::vector<std::unique_ptr<getfem::slicer_action>> &slicers,
              const gfi_array *arg) {
  GMM_ASSERT1(gfi_array_get_class(arg) == GFI_CELL, "slices must be "
              "described as imbricated cell arrays");
  mexargs_in in(1, &arg, true);
  std::string cmd = in.pop().to_string();
  if (check_cmd(cmd, "none", in, 0, 0)) {
    slicers.push_back(std::make_unique<getfem::slicer_none>());
  } else if (check_cmd(cmd, "planar", in, 3, 3)) {
    int orient = in.pop().to_integer(-1,2);
    getfem::base_node x0 = in.pop().to_base_node();
    slicers.push_back(std::make_unique<getfem::slicer_half_space>(x0, in.pop().to_base_node(), orient));
  } else if (check_cmd(cmd, "ball", in, 3, 3)) {
    int orient = in.pop().to_integer(-1,2);
    getfem::base_node x0 = in.pop().to_base_node();
    slicers.push_back(std::make_unique<getfem::slicer_sphere>(x0, in.pop().to_scalar(1e-5), orient));
  } else if (check_cmd(cmd, "cylinder", in, 4, 4)) {
    int orient = in.pop().to_integer(-1,2);
    getfem::base_node x0 = in.pop().to_base_node();
    getfem::base_node x1 = in.pop().to_base_node();
    slicers.push_back(std::make_unique<getfem::slicer_cylinder>(x0, x1, in.pop().to_scalar(1e-5), orient));
  } else if (check_cmd(cmd, "isovalues", in, 4, 4)) {
    int orient = in.pop().to_integer(-1,2);
    const getfem::mesh_fem &mf = *to_meshfem_object(in.pop());
    darray U = in.pop().to_darray(int(mf.nb_dof()));
    slicers.push_back(std::make_unique<getfem::slicer_isovalues>(getfem::mesh_slice_cv_dof_data<darray>(mf,U),
                                                   in.pop().to_scalar(), orient));
  } else if (check_cmd(cmd, "boundary", in, 0, 1)) {
    getfem::slicer_action *s1 = 0;
    if (in.remaining()) {
      s1 = build_slicers(m, slicers, in.pop().arg);
    } else {
      slicers.push_back(std::make_unique<getfem::slicer_none>());
      s1 = slicers.back().get();
    }
    getfem::mesh_region cvflst;
    getfem::outer_faces_of_mesh(m, m.convex_index(), cvflst);
    slicers.push_back(std::make_unique<getfem::slicer_boundary>(m,*s1,cvflst));
  } else if (check_cmd(cmd, "explode", in, 1, 1)) {
    scalar_type c = in.pop().to_scalar();
    slicers.push_back(std::make_unique<getfem::slicer_explode>(c));
  } else if (check_cmd(cmd, "union", in, 1, -1)) {
    getfem::slicer_action *s1 = build_slicers(m, slicers, in.pop().arg);
    while (in.remaining()) {
      getfem::slicer_action *s2 = build_slicers(m, slicers, in.pop().arg);
      slicers.push_back(std::make_unique<getfem::slicer_union>(*s1,*s2));
      s1 = slicers.back().get();
    }
  } else if (check_cmd(cmd, "intersection", in, 1, -1)) {
    getfem::slicer_action *s1 = build_slicers(m, slicers, in.pop().arg);
    while (in.remaining()) {
      getfem::slicer_action *s2 = build_slicers(m, slicers, in.pop().arg);
      slicers.push_back(std::make_unique<getfem::slicer_intersect>(*s1,*s2));
      s1 = slicers.back().get();
    }
  } else if (check_cmd(cmd, "diff", in, 2, 2)) {
    getfem::slicer_action *s1 = build_slicers(m, slicers, in.pop().arg);
    getfem::slicer_action *s2 = build_slicers(m, slicers, in.pop().arg);
    slicers.push_back(std::make_unique<getfem::slicer_complementary>(*s2));
    slicers.push_back(std::make_unique<getfem::slicer_intersect>(*s1, *slicers.back()));
  } else if (check_cmd(cmd, "comp", in, 1, 1)) {
    getfem::slicer_action *s = build_slicers(m, slicers, in.pop().arg);
    slicers.push_back(std::make_unique<getfem::slicer_complementary>(*s));
  } else if (check_cmd(cmd, "mesh", in, 1, 1)) {
    const getfem::mesh &m2 = *extract_mesh_object(in.pop());
    slicers.push_back(std::make_unique<getfem::slicer_mesh_with_mesh>(m2));
  } else bad_cmd(cmd);
  return slicers.back().get();
}


/*@GFDOC
  Creation of a mesh slice. Mesh slices are very similar to a
  P1-discontinuous @tmf on which interpolation is very fast. The slice is
  built from a mesh object, and a description of the slicing operation, for
  example::

    sl = SLICE:INIT(@CELL{'planar',+1,@MATLAB{[0;0],[1;0]}@SCILAB{[0;0],[1;0]}@PYTHON{[[0],[0]],[[0],[1]]}}, m, 5)@MATLAB{;}@SCILAB{;}

  cuts the original mesh with the half space {y>0}. Each convex of the
  original @tmesh `m` is simplexified (for example a quadrangle is splitted
  into 2 triangles), and each simplex is refined 5 times.

  Slicing operations can be:

  * cutting with a plane, a sphere or a cylinder
  * intersection or union of slices
  * isovalues surfaces/volumes
  * "points", "streamlines" (see below)

  If the first argument is a @tmf `mf` instead of a @tmesh, and if it is
  followed by a `mf`-field `u`@MATLAB{ (with size(u,1) ==
  MESH_FEM:GET('nbdof'))}, then the deformation `u` will be applied to the
  mesh before the slicing operation.

  The first argument can also be a slice.
@*/

/*@INIT sl = ('.op',sliceop, {@tsl sl|{@tmesh m| @tmf mf, vec U}, int refine}[, @mat CVfids])
  Create a @tsl using `sliceop` operation.

  `sliceop` operation is specified with @MATLAB{Matlab CELL arrays (i.e.
  with braces)}@SCILAB{Scilab CELL arrays (i.e. with braces)} @PYTHON{Tuple
  or List, do not forget the extra parentheses!}. The first element is the
  name of the operation, followed the slicing options:

  * @CELL{'none'} :
    Does not cut the mesh.

  * @CELL{'planar', @int orient, @vec p, @vec n} :
    Planar cut. `p` and `n` define a half-space, `p` being a point belong to
    the boundary of the half-space, and `n` being its normal. If `orient` is
    equal to -1 (resp. 0, +1), then the slicing operation will cut the mesh
    with the "interior" (resp. "boundary", "exterior") of the half-space.
    `orient` may also be set to +2 which means that the mesh will be sliced,
    but both the outer and inner parts will be kept.

  * @CELL{'ball', @int orient, @vec c, @scalar r} :
    Cut with a ball of center `c` and radius `r`.

  * @CELL{'cylinder', @int orient, @vec p1, @vec p2, @scalar r} :
    Cut with a cylinder whose axis is the line `(p1, p2)` and whose radius
    is `r`.

  * @CELL{'isovalues', @int orient, @tmf mf, @vec U, @scalar s} :
    Cut using the isosurface of the field `U` (defined on the @tmf `mf`).
    The result is the set `{x such that :math:`U(x) \leq s`}` or `{x such that
    `U`(x)=`s`}` or `{x such that `U`(x) >= `s`}` depending on the value of
    `orient`.

  * @CELL{'boundary'[, SLICEOP]} :
    Return the boundary of the result of SLICEOP, where SLICEOP is any
    slicing operation. If SLICEOP is not specified, then the whole mesh is
    considered (i.e. it is equivalent to @CELL{'boundary',{'none'}}).

  * @CELL{'explode', @mat Coef} :
    Build an 'exploded' view of the mesh: each convex is shrinked (:math:`0 <
    \text{Coef} \leq 1`). In the case of 3D convexes, only their faces are kept.

  * @CELL{'union', SLICEOP1, SLICEOP2} :
    Returns the union of slicing operations.

  * @CELL{'intersection', SLICEOP1, SLICEOP2} :
    Returns the intersection of slicing operations, for example::

      sl = SLICE:INIT(@CELL{intersection',@CELL{'planar',+1,@MATLAB{[0;0;0],[0;0;1]}@SCILAB{[0;0;0],[0;0;1]}@PYTHON{[[0],[0],[0]],[[0],[0],[1]]}},
                                 @CELL{'isovalues',-1,mf2,u2,0}},mf,u,5)

  * @CELL{'comp', SLICEOP} :
    Returns the complementary of slicing operations.

  * @CELL{'diff', SLICEOP1, SLICEOP2} :
    Returns the difference of slicing operations.

  * @CELL{'mesh', @tmesh m} :
    Build a slice which is the intersection of the sliced mesh with another
    mesh. The slice is such that all of its simplexes are stricly contained
    into a convex of each mesh.
@*/

void gf_slice(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg()  <  2) THROW_BADARG("Wrong number of input arguments");
  if (!out.narg_in_range(1,1)) THROW_BADARG("Wrong number of output arguments");

  const getfem::mesh *mm = 0;
  getfem::mesh_level_set *pmls = 0;
  auto pstored = std::make_shared<getfem::stored_mesh_slice>();

  /* "normal" slices */
  if (in.front().is_cell()) {

    /* build slicers */
    const gfi_array *arg = in.pop().arg;

    /* check the source argument (mesh/mesh_fem or slice) */
    std::unique_ptr<getfem::mesh_slice_cv_dof_data<darray> > mfdef;
    std::unique_ptr<getfem::slicer_action> slicer_def;
    getfem::stored_mesh_slice *source_slice = 0;
    if (is_meshfem_object(in.front()) && in.remaining()  >=  3) {
      const getfem::mesh_fem& mf = *to_meshfem_object(in.pop());
      mm = &mf.linked_mesh();
      darray Udef = in.pop().to_darray(-2, int(mf.nb_dof()));
      if (!(mf.get_qdim() == mm->dim() && Udef.getm() == 1) &&
          !(mf.get_qdim() == 1 && Udef.getm() == mm->dim())) {
        THROW_BADARG("either the mesh_fem must have a Qdim=" << int(mm->dim()) <<
                     ", either the data must have " << int(mm->dim()) << " rows");
      }
      mfdef.reset(new getfem::mesh_slice_cv_dof_data<darray>(mf,Udef));
      slicer_def.reset(new getfem::slicer_apply_deformation(*mfdef.get()));
    } else if (is_slice_object(in.front())) {
      source_slice = to_slice_object(in.pop());
      mm = &(source_slice->linked_mesh());
    } else if (is_mesh_levelset_object(in.front())) {
      pmls = to_mesh_levelset_object(in.pop());
      mm = &(pmls->linked_mesh());
    } else {
      mm = extract_mesh_object(in.pop());
    }

    std::vector<std::unique_ptr<getfem::slicer_action>> slicers;
    getfem::slicer_action * s = build_slicers(*mm, slicers, arg);

    /* create the slicer and registers the actions */
    getfem::mesh_slicer slicer(*mm);
    if (pmls) slicer.using_mesh_level_set(*pmls);
    getfem::slicer_build_stored_mesh_slice slicer_store(*pstored);
    if (slicer_def.get()) slicer.push_back_action(*slicer_def.get());
    slicer.push_back_action(*s);
    slicer.push_back_action(slicer_store);

    /* not building from slice ? */
    if (source_slice == 0) {
      if (in.remaining() == 0) THROW_BADARG("Not enough input arguments");
      size_type nrefine = in.pop().to_integer(1,1000);
      /*std::vector<convex_face> cvf;
      if (in.remaining()) {
        iarray v = in.pop().to_iarray(-2, -1);
        build_convex_face_lst(*mm, cvf, &v);
        } else build_convex_face_lst(*mm, cvf, 0);*/
      getfem::mesh_region rg = to_mesh_region(*mm, in);
      slicer.exec(nrefine, rg);
    } else {
      slicer.exec(*source_slice);
    }
    if (in.remaining()) THROW_BADARG("too much input arguments");

  } else if (in.front().is_string()) {
    /* "special" slices are handle here.. */
    std::string cmd = in.pop().to_string();
    if (check_cmd(cmd, "streamlines", in, 3, 3)) {
    /*@INIT sl = ('streamlines', @tmf mf, @mat U, @dmat S)
      Compute streamlines of the (vector) field `U`, with seed points given
      by the columns of `S`.@*/
      const getfem::mesh_fem *mf = to_meshfem_object(in.front());
      mm = extract_mesh_object(in.pop());
      darray U = in.pop().to_darray(int(mf->nb_dof()));
      darray v = in.pop().to_darray(mm->dim(), -1);
      std::vector<getfem::base_node> seeds(v.getn());
      for (unsigned j=0; j < v.getn(); ++j)
        seeds[j] = v.col_to_bn(j);
      getfem::mesh_slice_cv_dof_data<darray> mfU(*mf,U);
      pstored = std::make_shared<getfem::mesh_slice_streamline>(&mfU, seeds,
								true, true);
    } else if (check_cmd(cmd, "points", in, 2, 2)) {
    /*@INIT sl = ('points', @tmesh m, @dmat Pts)
      Return the "slice" composed of points given by the columns of `Pts`
      (useful for interpolation on a given set of sparse points, see
      ``::COMPUTE('interpolate on',sl)``).@*/
      mm = extract_mesh_object(in.pop());
      pstored = std::make_shared<getfem::stored_mesh_slice>();
      getfem::mesh_slicer slicer(*mm);
      getfem::slicer_build_stored_mesh_slice slicer_store(*pstored);
      slicer.push_back_action(slicer_store);

      darray w = in.pop().to_darray(mm->dim(), -1);
      std::vector<getfem::base_node> N(w.getn());
      for (unsigned i=0; i < w.getn(); ++i) N[i] = w.col_to_bn(i);
      slicer.exec(N);
    } else if (check_cmd(cmd, "load", in, 1, 2)) {
    /*@INIT sl = ('load', @str filename[, @tmesh m])
      Load the slice (and its linked mesh if it is not given as an argument)
      from a text file.@*/
      std::string fname = in.pop().to_string();
      pstored = std::make_shared<getfem::stored_mesh_slice>();
      if (in.remaining()) mm = extract_mesh_object(in.pop());
      else {
	auto m = std::make_shared<getfem::mesh>();
	m->read_from_file(fname);
	store_mesh_object(m);
	mm = m.get();
	workspace().add_hidden_object(store_slice_object(pstored), m);
      }
      pstored->read_from_file(fname, *mm);
    } else bad_cmd(cmd);
  } else THROW_BADARG("a slicer specification (i.e. cell array) or a string "
                      "was expected as the first argument");
  if (mm == 0 || pstored.get() == 0) THROW_INTERNAL_ERROR;

  id_type id = store_slice_object(pstored);
  out.pop().from_object_id(id, SLICE_CLASS_ID);
  workspace().set_dependence(pstored.get(), mm);
}
