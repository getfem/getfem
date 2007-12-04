// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2005-2008 Julien Pommier.
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

#include <getfemint.h>
#include <getfemint_workspace.h>
#include <getfemint_mesh.h>
#include <getfem/getfem_import.h>
#include <getfem/getfem_regular_meshes.h>

using namespace getfemint;

static void
cartesian_mesh(getfem::mesh *pmesh, getfemint::mexargs_in &in)
{
  getfemint::size_type dim = in.remaining();
  
  if (dim == 0) THROW_BADARG( "not enough input arguments");

  std::vector<darray> ppos(dim);
  std::vector<size_type> npts(dim);
  dal::uint32_type grid_npoints=1, grid_nconvex=1;
  for (size_type i = 0; i < dim; i++) {
    ppos[i] = in.pop().to_darray();
    npts[i] = ppos[i].size();
    grid_npoints *= npts[i];
    grid_nconvex *= (npts[i]-1);
  }
  
  /* add the points in 'fortran style' order */
  getfem::base_node pt(dim);
  for (dal::uint32_type i=0; i < grid_npoints; i++) {
    dal::uint32_type k = i;    
    for (getfemint::size_type j = 0; j < dim; j++) {
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
    /*    if (i == grid_npoints-1) {
      cerr << "nb de pts ajoutes: " << grid_npoints << " id dernier =" << id_pt << endl;
      }*/
  }
  

  std::vector<int> ipt(dim);
  std::vector<getfem::base_node> pts(1 << dim+1);
  
  bgeot::pgeometric_trans pgt = bgeot::parallelepiped_linear_geotrans(dim);

  /* add the convexes */
  for (dal::uint32_type i=0; i < grid_nconvex; i++) {
    dal::uint32_type k = i;

    /* find point location */
    for (getfemint::size_type j = 0; j < dim; j++) {
      ipt[j] = k % (npts[j]-1);
      k /= (npts[j]-1);
    }

    /* build the vertices list */
    for (getfemint::size_type j = 0; j < (unsigned(1)<<dim); j++) {
      pts[j].resize(dim);
      for (dal::uint32_type d=0; d < dim; d++) {
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
  unsigned N = nsubdiv.size();

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
    (*pmesh2, N, org, vtab.begin(), nsubdiv.begin());

  bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(N, K);

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
    pmesh->optimize_structure();
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
}

static void
curved_mesh(getfem::mesh *dest_mesh, getfemint::mexargs_in &in)
{
  const getfem::mesh *src_mesh = in.pop().to_const_mesh();
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
  const getfem::mesh *src_mesh = in.pop().to_const_mesh();
  unsigned nblay = in.pop().to_integer(1,2500000);
  getfem::extrude(*src_mesh, *dest_mesh, nblay);
}

static void
ptND_mesh(getfem::mesh *mesh, bool is2D, getfemint::mexargs_in &in)
{
  darray P = in.pop().to_darray(-1, -1);
  iarray T = in.pop().to_iarray(-1, -1);
  cout << "T(" << T.getm() << ", " << T.getn() << "), size=" << T.size() << "\n";
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
  for (size_type i = 0; i < P.getn(); ++i) {
    id_tab[i] = mesh->add_point(P.col_to_bn(i));

    /* une hypothèse bien commode pour la "compatibilité pdetool" */
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

/*MLABCOM

  FUNCTION M = gf_mesh([operation [, args]])
  General constructor for mesh object. Returns a getfem handle to the
  newly created mesh object. Note that for recent (> 6.0) versions of
  matlab, you should replace the calls to 'gf_mesh' with 'gfMesh'
  (this will instruct Matlab to consider the getfem mesh as a regular
  matlab object that can be manipulated with get() and set() methods).

  
  @INIT MESH:INIT ('empty')
  @INIT MESH:INIT ('cartesian')
  @INIT MESH:INIT ('regular simplices')
  @INIT MESH:INIT ('triangles grid')
  @INIT MESH:INIT ('curved')
  @INIT MESH:INIT ('prismatic')
  @INIT MESH:INIT ('pt2D')

  @INIT MESH:INIT ('ptND')
  @INIT MESH:INIT ('load')
  @INIT MESH:INIT ('from string')
  @INIT MESH:INIT ('import')
  @INIT MESH:INIT ('clone')
 
  $Id$
MLABCOM*/
void gf_mesh(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  getfemint_mesh *mi_mesh = 
    getfemint_mesh::get_from(new getfem::mesh);
  out.pop().from_object_id(mi_mesh->get_id(), MESH_CLASS_ID);
  getfem::mesh *pmesh = &mi_mesh->mesh();


  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  std::string cmd    = in.pop().to_string();
  if (check_cmd(cmd, "empty", in, out, 1, 1, 0, 1)) {
    /*@INIT  M=MESH:INIT('empty', @int dim)
      Create a new empty mesh.
      @*/
    size_type dim = in.pop().to_integer(1,255);
    getfem::base_node pt(dim); 
    /* just to initialize the dimension of the mesh
       (this is not very nice, i know) */
    pmesh->sup_point(pmesh->add_point(pt)); 
  } else if (check_cmd(cmd, "cartesian", in, out, 1, 32, 0, 1)) {
    /*@INIT  M=MESH:INIT('cartesian', @dvec X[, @dvec Y[, @dvec Z,..]])
      Build quickly a regular mesh of quadrangles, cubes, etc.
      @*/
    cartesian_mesh(pmesh, in);
  } else if (check_cmd(cmd, "triangles grid", in, out, 2, 2, 0, 1)) {
    /*@INIT  M=MESH:INIT('triangles grid', @dvec X, @dvec Y)
      Build quickly a regular mesh of triangles. 

      This is a very limited and somehow deprecated function (See also
      MESH:INIT('ptND'), MESH:INIT('regular simplices') and
      MESH:INIT('cartesian')). @*/
    triangles_grid_mesh(pmesh, in);
  } else if (check_cmd(cmd, "regular simplices", in, out, 1, 32, 0, 1)) {
    /*@INIT M=MESH:INIT('regular simplices', @dvec X[, @dvec Y[, @dvec Z,.., ]]['degree', @int K]['noised'])
      Mesh a n-dimensionnal parallelepipeded with simplices
      (triangles, tetrahedrons etc) .

      The optional degree may be used to build meshes with non linear
      geometric transformations. 
      @*/
    regular_simplices_mesh(pmesh, in);
  } else if (check_cmd(cmd, "curved", in, out, 2, 2, 0, 1)) {
    /*@INIT  M=MESH:INIT('curved', M0, @dvec F)
      Build a curved (n+1)-dimensions mesh from a n-dimensions mesh M0. 

      The points of the new mesh have one additional coordinate, given by
      the vector F. This can be used to obtain meshes for shells. M0 may
      be a @tmf object, in that case its linked mesh will be used. 
      @*/
    curved_mesh(pmesh, in);
  } else if (check_cmd(cmd, "prismatic", in, out, 2, 2, 0, 1)) {
    /*@INIT  M=MESH:INIT('prismatic', M0, @int NLAY)
      Extrude a prismatic mesh M from a mesh M0. 
      
      In the additional dimension there are NLAY layers of elements built
      from 0 to 1. @*/
    prismatic_mesh(pmesh, in);
  } else if (check_cmd(cmd, "pt2D", in, out, 2, 3, 0, 1)) {
    /*@INIT  M=MESH:INIT('pt2D', @dmat P, @ivec T[, @int N])
      Build a mesh from a 2D triangulation.

      Each column of P contains a point coordinate, and each column of T contains the point indices of a triangle. N is optional and
      is a zone number. If N is specified then only the zone number
      'N' is converted (in that case, T is expected to have 4 rows,
      the fourth containing these zone numbers).@MATLAB{<Par>Can be used to Convert a "pdetool" triangulation exported in variables P and T into a GETFEM mesh.}
      @*/
    ptND_mesh(pmesh, true, in);
  } else if (check_cmd(cmd, "ptND", in, out, 2, 2, 0, 1)) {
    /*@INIT  M=MESH:INIT('ptND', @dmat P, @imat T)
      Build a mesh from a N-dimensional "triangulation".

      Similar function to 'pt2D', for building simplexes meshes from a
      triangulation given in T, and a list of points given in P. The
      dimension of the mesh will be the number of rows of P, and the
      dimension of the simplexes will be the number of rows of T.
      @*/
    ptND_mesh(pmesh, 0, in);
  } else if (check_cmd(cmd, "load", in, out, 1, 1, 0, 1)) {
    /*@INIT  M=MESH:INIT('load', @str FILENAME)
      Load a mesh from a GETFEM++ ascii mesh file. See also
      MESH:GET('save',FILENAME). @*/
    std::string fname = in.pop().to_string();
    pmesh->read_from_file(fname);
  } else if (check_cmd(cmd, "from string", in, out, 1, 1, 0, 1)) {
    /*@INIT  M=MESH:INIT('from string', @str S)
      Load a mesh from a string description. For example, a string returned
      by MESH:GET('char').
      @*/
    std::stringstream ss(in.pop().to_string());
    pmesh->read_from_file(ss);
  } else if (check_cmd(cmd, "import", in, out, 2, 2, 0, 1)) {
    /*@INIT  M=MESH:INIT('import', @str FORMAT, @str FILENAME)
      Import a mesh, FORMAT may be:<Par>
      - 'gmsh'   for a mesh created with gmsh ( http://www.geuz.org/gmsh )<par>
      - 'gid'    for a mesh created with GiD  ( http://gid.cimne.upc.es )<par>
      - 'am_fmt' for a mesh created with emc2 ( http://pauillac.inria.fr/cdrom/www/emc2/fra.htm )
      @*/
    std::string fmt = in.pop().to_string();
    std::string fname = in.pop().to_string();
    getfem::import_mesh(fname, fmt, *pmesh);
  } else if (check_cmd(cmd, "clone", in, out, 1, 1, 0, 1)) {
    /*@INIT  M=MESH:INIT('clone', @tmesh M2)
      Create a copy of a mesh.
      @*/
    const getfem::mesh *m2 = in.pop().to_const_mesh();
    pmesh->copy_from(*m2);
  } else bad_cmd(cmd);
}
