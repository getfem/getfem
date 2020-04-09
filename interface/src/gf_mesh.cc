/*===========================================================================

 Copyright (C) 2005-2020 Julien Pommier.

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
#include <getfemint.h>
#include <getfemint_workspace.h>
#include <getfem/getfem_import.h>
#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_mesher.h>

using namespace getfemint;

static void
cartesian_mesh(getfem::mesh *pmesh, getfemint::mexargs_in &in,
	       bool linear=true) {
  getfemint::size_type dim = in.remaining();
  
  if (dim == 0) THROW_BADARG( "not enough input arguments");

  std::vector<darray> ppos(dim);
  std::vector<size_type> npts(dim);
  size_type grid_npoints=1, grid_nconvex=1;
  for (size_type i = 0; i < dim; i++) {
    ppos[i] = in.pop().to_darray();
    npts[i] = ppos[i].size();
    grid_npoints *= npts[i];
    grid_nconvex *= (npts[i]-1);
  }

  /* add the points in 'fortran style' order */
  getfem::base_node pt(dim);
  for (size_type i=0; i < grid_npoints; i++) {
    size_type k = i;
    for (size_type j = 0; j < dim; j++) {
      pt[j] = ppos[j][k % (npts[j])];
      k /= (npts[j]);
    }

    size_type id_pt = pmesh->add_point(pt);
    if (id_pt != i) {
      THROW_ERROR(
		"something has changed in getfem, you need to reconsider "
		"gf_mesh('cartesian')\nfor point " << i <<
		", the index is " << id_pt << endl);
    }
  }


  std::vector<int> ipt(dim);
  std::vector<getfem::base_node> pts(1 << (dim+1));

  bgeot::pgeometric_trans pgt = linear ? bgeot::parallelepiped_linear_geotrans(dim)
                                       : bgeot::parallelepiped_geotrans(dim, 1);

  /* add the convexes */
  for (size_type i=0; i < grid_nconvex; i++) {
    size_type k = i;

    /* find point location */
    for (size_type j = 0; j < dim; j++) {
      ipt[j] = int(k % (npts[j]-1));
      k /= (npts[j]-1);
    }

    /* build the vertices list */
    for (size_type j = 0; j < (unsigned(1)<<dim); j++) {
      pts[j].resize(dim);
      for (size_type d=0; d < dim; d++) {
	if ((j >> d) & 1) {
	  pts[j][d] = ppos[d][ipt[d]+1];
	} else {
	  pts[j][d] = ppos[d][ipt[d]];
	}
      }
    }

    // we don't use the add_parall since the geometric transformation
    // is linear (the mesh is cartesian)
    //pmesh->add_parallelepiped_by_points(dim, pts.begin());
    pmesh->add_convex_by_points(pgt, pts.begin());
  }
}

static void
pyramidal_mesh(getfem::mesh *pmesh, getfemint::mexargs_in &in) {
  getfemint::size_type dim = 3;

  std::vector<darray> ppos(dim);
  std::vector<size_type> npts(dim);
  size_type grid_npoints=1, grid_nconvex=1;
  for (size_type i = 0; i < dim; i++) {
    ppos[i] = in.pop().to_darray();
    npts[i] = ppos[i].size();
    grid_npoints *= npts[i];
    grid_nconvex *= (npts[i]-1);
  }

  /* add the points in 'fortran style' order */
  getfem::base_node pt(dim);
  for (size_type i=0; i < grid_npoints; i++) {
    size_type k = i;
    for (size_type j = 0; j < dim; j++) {
      pt[j] = ppos[j][k % (npts[j])];
      k /= (npts[j]);
    }

    size_type id_pt = pmesh->add_point(pt);
    if (id_pt != i) {
      THROW_ERROR(
		"something has changed in getfem, you need to reconsider "
		"gf_mesh('cartesian')\nfor point " << i <<
		", the index is " << id_pt << endl);
    }
  }

  std::vector<int> ipt(dim);
  std::vector<getfem::base_node> pts(1 << (dim+1));

  bgeot::pgeometric_trans pgt = bgeot::pyramid_QK_geotrans(1);

  /* add the convexes */
  for (size_type i=0; i < grid_nconvex; i++) {
    size_type k = i;

    /* find point location */
    for (size_type j = 0; j < dim; j++) {
      ipt[j] = int(k % (npts[j]-1));
      k /= (npts[j]-1);
    }

    /* build the vertices list */
    for (size_type j = 0; j < (unsigned(1)<<dim); j++) {
      pts[j].resize(dim);
      for (size_type d=0; d < dim; d++) {
	if ((j >> d) & 1) {
	  pts[j][d] = ppos[d][ipt[d]+1];
	} else {
	  pts[j][d] = ppos[d][ipt[d]];
	}
      }
    }
	
    bgeot::base_node barycenter(3);
    std::vector<size_type> iipts(8);
    for (size_type j = 0; j < 8; j++) {
	barycenter += pts[j];
	iipts[j] = pmesh->add_point(pts[j]);
    }
    barycenter /= 8.;
    size_type ib = pmesh->add_point(barycenter);
    pmesh->add_pyramid(iipts[0],iipts[1],iipts[2],iipts[3],ib);
    pmesh->add_pyramid(iipts[7],iipts[6],iipts[5],iipts[4],ib);
    pmesh->add_pyramid(iipts[0],iipts[4],iipts[1],iipts[5],ib);
    pmesh->add_pyramid(iipts[1],iipts[5],iipts[3],iipts[7],ib);
    pmesh->add_pyramid(iipts[3],iipts[7],iipts[2],iipts[6],ib);
    pmesh->add_pyramid(iipts[2],iipts[6],iipts[0],iipts[4],ib);

  }
}

static void
triangles_grid_mesh(getfem::mesh *pmesh, getfemint::mexargs_in &in)
{
  if (in.remaining() != 2) THROW_BADARG( "not enough input arguments");

  darray X = in.pop().to_darray();
  darray Y = in.pop().to_darray();
  if (X.size() < 1 || Y.size() < 1) THROW_BADARG( "bad dimensions");

  size_type ni = Y.size(), nj = X.size();
  for (size_type i=0; i < ni; i++) {
    for (size_type j=0; j < nj; j++) {
      getfem::base_node pt(2);
      pt[0] = X[j]; pt[1] = Y[i];
      //      cerr << "pt = " << pt << endl;
      pmesh->add_point(pt);
    }
  }
  for (size_type i=0; i < ni-1; i++) {
    for (size_type j=0; j < nj-1; j++) {
      //cerr << "i=" << i <<" j=" << j << endl;
      pmesh->add_triangle(i*nj + j, (i+1)*nj + j  , (i+1)*nj + j+1);
      pmesh->add_triangle(i*nj + j,     i*nj + j+1, (i+1)*nj + j+1);
    }
  }
}

/*static void
regular_simplices_mesh(getfem::mesh *pmesh, getfemint::mexargs_in &in) {
  std::vector<size_type> nsubdiv =
    in.pop().to_iarray(-1).to_vector<std::vector<size_type> >();
  unsigned N = nsubdiv.size();
  bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(N, 1);
  bool noised = false;
  getfem::regular_unit_mesh(*pmesh, nsubdiv, pgt, noised);
  if (in.remaining()) {
    darray len  = in.pop().to_darray(N);
    getfem::base_matrix M(N,N);
    for (unsigned i=0; i < N; ++i) M(i,i) = len[i];
    pmesh->transformation(M);
  }
  }*/

static void
regular_simplices_mesh(getfem::mesh *pmesh, getfemint::mexargs_in &in) {
  std::vector<darray> xyz;
  std::vector<size_type> nsubdiv;
  unsigned K = 1;
  bool noised = false;
  while (in.remaining()) {
    if (in.front().is_string()) {
      std::string s = in.pop().to_string();
      if (cmd_strmatch(s, "degree")) {
	if (!in.remaining()) { THROW_BADARG("missing degree"); }
	else K=in.pop().to_integer(1, 10);
	} else if (cmd_strmatch(s, "noised"))
	noised = true;
    } else {
      xyz.push_back(in.pop().to_darray(-1));
      if (xyz.back().size() <= 1) THROW_BADARG("wrong dimensions");
      nsubdiv.push_back(xyz.back().size() - 1);
    }
  }
  unsigned N = unsigned(nsubdiv.size());

  getfem::base_node org(N);
  std::vector<getfem::base_small_vector> vtab(N);
  for (dim_type i = 0; i < N; i++) {
    vtab[i] = getfem::base_small_vector(N);
    (vtab[i])[i] = 1;
  }

  getfem::mesh *pmesh2 = 0, msh;
  if (K == 1) pmesh2 = pmesh;
  else pmesh2 = &msh;

  getfem::parallelepiped_regular_simplex_mesh
    (*pmesh2, dim_type(N), org, vtab.begin(), nsubdiv.begin());

  bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(N, short_type(K));

  if (K > 1) {
    /* build a mesh with a geotrans of degree K */
    for (dal::bv_visitor cv(msh.convex_index()); !cv.finished(); ++cv) {
      std::vector<getfem::base_node> pts(pgt->nb_points());
      for (size_type i=0; i < pgt->nb_points(); ++i) {
	pts[i] = msh.trans_of_convex(cv)->transform
	  (pgt->convex_ref()->points()[i], msh.points_of_convex(cv));
      }
      pmesh->add_convex_by_points(pgt, pts.begin());
    }
    pmesh->optimize_structure(false);
  }

  getfem::base_small_vector diff(N);
  for (unsigned k=0; k < N; ++k) {
    diff[k] = std::abs(xyz[k][1] - xyz[k][0]);
    for (unsigned i=1; i < xyz[k].size(); ++i)
      diff[k] = std::min(diff[k], std::abs(xyz[k][i] - xyz[k][i-1]));
  }

  for (dal::bv_visitor i(pmesh->points().index()); !i.finished(); ++i) {
    getfem::base_node &p = pmesh->points()[i];
    for (unsigned k=0; k < N; ++k) {
      unsigned ii = unsigned(p[k] + 1e-6);
      assert(ii < xyz[k].size());
      scalar_type a = p[k]-ii;
      if (ii != xyz[k].size()-1)
	p[k] = (1-a) * xyz[k][ii] + a*xyz[k][ii+1];
      else p[k] = xyz[k][ii];
      if (noised && ii != 0 && ii != nsubdiv[k])
	p[k] += diff[k] * gmm::random(double()) * 0.2 / K;
    }
  }
  pmesh->points().resort();
}

static void
curved_mesh(getfem::mesh *dest_mesh, getfemint::mexargs_in &in)
{
  const getfem::mesh *src_mesh = extract_mesh_object(in.pop());
  darray F = in.pop().to_darray(src_mesh->points().index().last()+1);

  int dim = src_mesh->dim();
  bgeot::base_node pt(dim+1);
  dest_mesh->clear();
  for (dal::bv_visitor i(src_mesh->points().index()); !i.finished(); ++i) {
    std::copy(src_mesh->points()[i].begin(), src_mesh->points()[i].end(), pt.begin());
    pt[dim] = F[i];
    size_type k = dest_mesh->add_point(pt);
    if (k != i) dest_mesh->swap_points(i,k); /* ca meriterait d'etre teste sur un maillage a trous ca.. */
  }

  for (dal::bv_visitor cv(src_mesh->convex_index()); !cv.finished(); ++cv) {
    dest_mesh->add_convex(src_mesh->trans_of_convex(cv),
			  src_mesh->ind_points_of_convex(cv).begin());
  }
}

static void
prismatic_mesh(getfem::mesh *dest_mesh, getfemint::mexargs_in &in)
{
  const getfem::mesh *src_mesh = extract_mesh_object(in.pop());
  size_type nblay = in.pop().to_integer(1,2500000);
  short_type degree(1);
  if (in.remaining()) degree = short_type(in.pop().to_integer(1,2500000));
  getfem::extrude(*src_mesh, *dest_mesh, nblay, degree);
}

static void
ptND_mesh(getfem::mesh *mesh, bool is2D, getfemint::mexargs_in &in)
{
  darray P = in.pop().to_darray(-1, -1);
  iarray T = in.pop().to_iarray(-1, -1);
  // cout << "T(" << T.getm() << ", " << T.getn() << "), size=" << T.size() << "\n";
  size_type mdim = P.getm();
  size_type N = is2D ? 2 : T.getm() - 1;
  if (is2D && T.getm() != 3 && T.getm() != 4) {
    THROW_BADARG("wrong nb of rows for t, 3 or 4 rows were expected, got " << T.getm());
  } else if (T.getm() < 1 || N > 10) {
    THROW_BADARG("wrong nb of rows for t (dim = 0 or dim > 10)");
  }
  if (mdim == 0 || mdim < N) {
    THROW_BADARG("cannot build simplexes of dimension " << N << " with points of dimension " << mdim);
  }
  id_type zone = 0;
  if (in.remaining()) zone = in.pop().to_integer(1,65000);

  size_type warn_cnt = 0;
  std::vector<id_type> id_tab(P.getn());
  for (unsigned i = 0; i < P.getn(); ++i) {
    id_tab[i] = id_type(mesh->add_point(P.col_to_bn(i)));

    /* une hypothese bien commode pour la "compatibilite pdetool" */
    if (id_tab[i] != i && warn_cnt++ == 0) {
      GMM_WARNING1("The numbering of mesh points will be different, pt#" <<
		   i+config::base_index() << " gets id#" << id_tab[i] + config::base_index());
    }
  }

  std::vector<size_type> ipts(N+1);
  for (size_type i = 0; i < T.getn(); ++i) {
    for (size_type k = 0; k < N+1; ++k) {
      ipts[k] = T(k,i) - config::base_index();
      if (ipts[k] >= P.size()) THROW_BADARG( "Bad triangulation.");
    }
    if (zone == 0 || (T.getm() == N+2 && zone == id_type(T(N+1,i))))
      mesh->add_convex(bgeot::simplex_geotrans(N,1), gmm::index_ref_iterator(id_tab.begin(), ipts.begin()));
  }
}


// Object for the declaration of a new sub-command.

struct sub_gf_mesh : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   getfem::mesh *pmesh) = 0;
};

typedef std::shared_ptr<sub_gf_mesh> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_mesh {					\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       getfem::mesh *pmesh)				\
      { dummy_func(in); dummy_func(out); code }				\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }


/*@GFDOC
  This object is able to store any element in any dimension even if you mix
  elements with different dimensions.

  @MATLAB{Note that for recent (> 6.0) versions of matlab, you should
  replace the calls to 'gf_mesh' with 'gfMesh' (this will instruct Matlab to
  consider the getfem mesh as a regular matlab object that can be
  manipulated with get() and set() methods).}
@*/

void gf_mesh(getfemint::mexargs_in& m_in,
	     getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {

    /*@INIT M = ('empty', @int dim)
      Create a new empty mesh.@*/
    sub_command
      ("empty", 1, 1, 0, 1,
       size_type dim = in.pop().to_integer(1,255);
       getfem::base_node pt(dim);
       /* just to initialize the dimension of the mesh
	  (this is not very nice, i know) */
       pmesh->sup_point(pmesh->add_point(pt));
       );


    /*@INIT M = ('cartesian', @dvec X[, @dvec Y[, @dvec Z,..]])
      Build quickly a regular mesh of quadrangles, cubes, etc.@*/
    sub_command
      ("cartesian", 1, 32, 0, 1,
       cartesian_mesh(pmesh, in);
       );

    /*@INIT M = ('pyramidal', @dvec X[, @dvec Y[, @dvec Z,..]])
      Build quickly a regular mesh of pyramids, etc.@*/
    sub_command
      ("pyramidal", 1, 32, 0, 1,
       pyramidal_mesh(pmesh, in);
       );

    /*@INIT M = ('cartesian Q1', @dvec X, @dvec Y[, @dvec Z,..])
      Build quickly a regular mesh of quadrangles, cubes, etc. with
      Q1 elements.@*/
    sub_command
      ("cartesian Q1", 2, 32, 0, 1,
       cartesian_mesh(pmesh, in, false);
       );


    /*@INIT M = ('triangles grid', @dvec X, @dvec Y)
      Build quickly a regular mesh of triangles.

      This is a very limited and somehow deprecated function (See also
      ``MESH:INIT('ptND')``, ``MESH:INIT('regular simplices')`` and
      ``MESH:INIT('cartesian')``).@*/
    sub_command
      ("triangles grid", 2, 2, 0, 1,
       triangles_grid_mesh(pmesh, in);
       );


    /*@INIT M = ('regular simplices', @dvec X[, @dvec Y[, @dvec Z,...]]['degree', @int k]['noised'])
      Mesh a n-dimensional parallelepiped with simplices (triangles,
      tetrahedrons etc) .

      The optional degree may be used to build meshes with non linear
      geometric transformations.@*/
    sub_command
      ("regular simplices", 1, 32, 0, 1,
       regular_simplices_mesh(pmesh, in);
       );


    /*@INIT M = ('curved', @tmesh m, @dvec F)
      Build a curved (n+1)-dimensions mesh from a n-dimensions mesh `m`.

      The points of the new mesh have one additional coordinate, given by
      the vector `F`. This can be used to obtain meshes for shells. `m` may
      be a @tmf object, in that case its linked mesh will be used.@*/
    sub_command
      ("curved", 2, 2, 0, 1,
       curved_mesh(pmesh, in);
       );


    /*@INIT M = ('prismatic', @tmesh m, @int nl[, @int degree])
      Extrude a prismatic @tmesh `M` from a @tmesh `m`.

      In the additional dimension there are `nl` layers of elements
      distributed from ``0`` to ``1``.
      If the optional parameter `degree` is provided with a value greater
      than the default value of ``1``, a non-linear transformation of
      corresponding degree is considered in the extrusion direction.@*/
    sub_command
      ("prismatic", 2, 3, 0, 1,
       prismatic_mesh(pmesh, in);
       );


    /*@INIT M = ('pt2D', @dmat P, @imat T[, @int n])
      Build a mesh from a 2D triangulation.

      Each column of `P` contains a point coordinate, and each column of `T`
      contains the point indices of a triangle. `n` is optional and is a
      zone number. If `n` is specified then only the zone number `n` is
      converted (in that case, `T` is expected to have 4 rows, the fourth
      containing these zone numbers).

      @MATLAB{Can be used to Convert a "pdetool" triangulation exported in
      variables P and T into a GETFEM mesh.}@*/
    sub_command
      ("pt2D", 2, 3, 0, 1,
       ptND_mesh(pmesh, true, in);
       );


    /*@INIT M = ('ptND', @dmat P, @imat T)
      Build a mesh from a n-dimensional "triangulation".

      Similar function to 'pt2D', for building simplexes meshes from a
      triangulation given in `T`, and a list of points given in `P`. The
      dimension of the mesh will be the number of rows of `P`, and the
      dimension of the simplexes will be the number of rows of `T`.@*/
    sub_command
      ("ptND", 2, 2, 0, 1,
       ptND_mesh(pmesh, 0, in);
       );


    /*@INIT M = ('load', @str filename)
      Load a mesh from a GetFEM ascii mesh file.

      See also ``MESH:GET('save', @str filename)``.@*/
    sub_command
      ("load", 1, 1, 0, 1,
       std::string fname = in.pop().to_string();
       pmesh->read_from_file(fname);
       );


    /*@INIT M = ('from string', @str s)
      Load a mesh from a string description.

      For example, a string returned by ``MESH:GET('char')``.@*/
    sub_command
      ("from string", 1, 1, 0, 1,
       std::stringstream ss(in.pop().to_string());
       pmesh->read_from_file(ss);
       );


    /*@INIT M = ('import', @str format, @str filename)
      Import a mesh.

      `format` may be:

      - 'gmsh' for a mesh created with `Gmsh`
      - 'gid' for a mesh created with `GiD`
      - 'cdb' for a mesh created with `ANSYS`
      - 'am_fmt' for a mesh created with `EMC2`@*/
    sub_command
      ("import", 2, 2, 0, 1,
       std::string fmt = in.pop().to_string();
       std::string fname = in.pop().to_string();
       getfem::import_mesh(fname, fmt, *pmesh);
       );


    /*@INIT M = ('clone', @tmesh m2)
      Create a copy of a mesh.@*/
    sub_command
      ("clone", 1, 1, 0, 1,
       const getfem::mesh *m2 = extract_mesh_object(in.pop());
       pmesh->copy_from(*m2);
       );

    /*@INIT M = ('generate', @tmo mo, @scalar h[, @int K = 1[, @mat vertices]])
      Call the experimental mesher of Getfem on the geometry
      represented by `mo`. please control the conformity of the produced mesh.
      You can help the mesher by adding a priori vertices in the array
      `vertices` which should be of size ``n x m`` where ``n`` n is the
      dimension of the mesh and ``m`` the number of points. `h` is
      approximate diameter of the elements. `K` is the degree of the
      mesh ( > 1 for curved boundaries).  The mesher try to optimize the
      quality of the elements. This operation may be time consuming.
      Note that if the mesh generation fails, because of some random
      procedure used, it can be run again since it will not give necessarily
      the same result due to random procedures used.
      The messages send to the console by the mesh generation can be
      deactivated using `gf_util('trace level', 2)`. More information
      can be obtained by `gf_util('trace level', 4)`. See ``MESHER_OBJECT:INIT``
      to manipulate geometric primitives in order to describe the geometry.
      @*/
    sub_command
      ("generate", 2, 4, 0, 1,

       getfem::pmesher_signed_distance psd = to_mesher_object(in.pop());
       double h = in.pop().to_scalar();
       int K = 1;
       if (in.remaining()) K = in.pop().to_integer(1,6);
       std::vector<getfem::base_node> fixed;
       if (in.remaining()) {
	 darray v = in.pop().to_darray(-1, -1);
	 for (int j=0; j < int(v.getn()); j++) {
	   getfem::base_node pt(v.getm());
	   gmm::copy(v.col_to_bn(j), pt);
	   fixed.push_back(pt);
	 }
       }
       int prefind = 1;
       int max_iter = 400;

       getfem::build_mesh(*pmesh, psd, h, fixed, K, -1, max_iter, prefind);
       );
  }


  if (m_in.narg() < 1)  THROW_BADARG( "Wrong number of input arguments");

  auto mesh = std::make_shared<getfem::mesh>();
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);


  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out, mesh.get());

    m_out.pop().from_object_id(store_mesh_object(mesh), MESH_CLASS_ID);
  }
  else bad_cmd(init_cmd);

}
