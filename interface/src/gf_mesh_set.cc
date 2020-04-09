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
#include <getfem/getfem_mesh.h>

using namespace getfemint;

static void check_empty_mesh(const getfem::mesh *pmesh) {
  if (pmesh->dim() == bgeot::dim_type(-1) || pmesh->dim() == 0) {
    THROW_ERROR( "mesh object has an invalid dimension");
  }
}

static void set_region(getfem::mesh &mesh, getfemint::mexargs_in& in,
                       bool do_clear=true) {
  unsigned boundary_num  = in.pop().to_integer(1);
  iarray v               = in.pop().to_iarray();

  getfem::mesh_region &rg = mesh.region(boundary_num);
  if (do_clear) rg.clear();
  
  if (v.getm() < 1 || v.getm() > 2 || v.getp() != 1 || v.getq() != 1)
    THROW_BADARG( "Invalid format for the convex or face list");

  /* loop over the edges of mxEdge */
  for (size_type j=0; j < v.getn(); j++) {
    size_type cv = size_type(v(0,j))-config::base_index();

    short_type f  = v.getm() == 2 ? short_type(v(1,j)-config::base_index())
                                  : short_type(-1);
    if (!mesh.convex_index().is_in(cv)) {
      THROW_BADARG( "Invalid convex number '" << cv+config::base_index() << "' at column " << j+config::base_index());
    }
    if (f != short_type(-1) && f >= mesh.structure_of_convex(cv)->nb_faces()) {
      THROW_BADARG( "Invalid face number '" << f+config::base_index() << "' at column " << j+config::base_index());
    }
    if (f == short_type(-1)) rg.add(cv);
    else                     rg.add(cv, f);
  }
}

static void intersect_regions(getfem::mesh &mesh, getfemint::mexargs_in& in) {
  unsigned ir1 = in.pop().to_integer(1,100000000);
  unsigned ir2 = in.pop().to_integer(1,100000000);
  getfem::mesh_region &r1 = mesh.region(ir1);
  getfem::mesh_region &r2 = mesh.region(ir2);
  r1 = getfem::mesh_region::intersection(r1, r2);
}

static void merge_regions(getfem::mesh &mesh, getfemint::mexargs_in& in) {
  unsigned ir1 = in.pop().to_integer(1,100000000);
  unsigned ir2 = in.pop().to_integer(1,100000000);
  getfem::mesh_region &r1 = mesh.region(ir1);
  getfem::mesh_region &r2 = mesh.region(ir2);
  r1 = getfem::mesh_region::merge(r1, r2);
}

static void subtract_regions(getfem::mesh &mesh, getfemint::mexargs_in& in) {
  unsigned ir1 = in.pop().to_integer(1,100000000);
  unsigned ir2 = in.pop().to_integer(1,100000000);
  getfem::mesh_region &r1 = mesh.region(ir1);
  getfem::mesh_region &r2 = mesh.region(ir2);
  r1 = getfem::mesh_region::subtract(r1, r2);
}





// Object for the declaration of a new sub-command.

struct sub_gf_mesh_set : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
                   getfemint::mexargs_out& out,
                   getfem::mesh *pmesh) = 0;
};

typedef std::shared_ptr<sub_gf_mesh_set> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_mesh_set {				\
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
  General function for modification of a mesh object.
@*/

void gf_mesh_set(getfemint::mexargs_in& m_in,
                 getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {


    /*@SET PIDs = ('pts', @mat PTS)
    Replace the coordinates of the mesh points with those given in `PTS`.@*/
    sub_command
      ("pts", 1, 1, 0, 1,
       darray P = in.pop().to_darray
       (pmesh->dim(),
        int(pmesh->points().index().last_true()+1));
       for (dal::bv_visitor i(pmesh->points().index()); !i.finished(); ++i) {
         for (unsigned k=0; k < pmesh->dim(); ++k)
           pmesh->points()[i][k] = P(k,i);
       }
       );


    /*@SET PIDs = ('add point', @mat PTS)
    Insert new points in the mesh and return their #ids.

    `PTS` should be an ``nxm`` matrix , where ``n`` is the mesh
    dimension, and ``m`` is the number of points that will be
    added to the mesh. On output, `PIDs` contains the point #ids
    of these new points.

    Remark: if some points are already part of the mesh (with a small
    tolerance of approximately ``1e-8``), they won't be inserted again,
    and `PIDs` will contain the previously assigned #ids of these
    points.@*/
    sub_command
      ("add point", 1, 1, 0, 1,
       check_empty_mesh(pmesh);
       darray v = in.pop().to_darray(pmesh->dim(), -1);
       iarray w = out.pop().create_iarray_h(v.getn());
       for (int j=0; j < int(v.getn()); j++) {
         w[j] = unsigned(pmesh->add_point(v.col_to_bn(j))
                         + config::base_index());
       }
       );


    /*@SET ('del point', @ivec PIDs)
    Removes one or more points from the mesh.

    `PIDs` should contain the point #ids, such as the one returned by
    the 'add point' command.@*/
    sub_command
      ("del point", 1, 1, 0, 0,
       check_empty_mesh(pmesh);
       iarray v = in.pop().to_iarray();
       
       for (size_type j=0; j < v.size(); j++) {
         id_type id = v[j]-config::base_index();
         if (pmesh->is_point_valid(id)) {
           THROW_ERROR("Can't remove point " << id+config::base_index()
                       << ": a convex is still attached to it.");
         }
         pmesh->sup_point(id);
       }
       );


    /*@SET CVIDs = ('add convex', @tgt GT, @mat PTS)
    Add a new convex into the mesh.

    The convex structure (triangle, prism,...) is given by `GT`
    (obtained with GEOTRANS:INIT('...')), and its points are given by
    the columns of `PTS`. On return, `CVIDs` contains the convex #ids.
    `PTS` might be a 3-dimensional array in order to insert more than
    one convex (or a two dimensional array correctly shaped according
    to Fortran ordering).@*/
    sub_command
      ("add convex", 2, 2, 0, 1,
       check_empty_mesh(pmesh);
       bgeot::pgeometric_trans pgt = to_geotrans_object(in.pop());
       darray v = in.pop().to_darray(pmesh->dim(), int(pgt->nb_points()), -1);
       iarray w = out.pop().create_iarray_h(v.getp());
       
       std::vector<getfemint::id_type> qp(pgt->nb_points());
       /* loop over convexes */
       for (unsigned k=0; k < v.getp(); k++) {
         /* loop over convex points */
         for (unsigned j=0; j < v.getn(); j++) {
           qp[j] = unsigned(pmesh->add_point(v.col_to_bn(j,k)));
         }
         id_type cv_id = id_type(pmesh->add_convex(pgt, qp.begin()));
         w[k] = cv_id+config::base_index();
       }
       );


    /*@SET ('del convex', @mat CVIDs)
    Remove one or more convexes from the mesh.

    `CVIDs` should contain the convexes #ids, such as the ones
    returned by the 'add convex' command.@*/
    sub_command
      ("del convex", 1, 1, 0, 0,
       check_empty_mesh(pmesh);
       iarray v = in.pop().to_iarray();
       
       for (size_type j=0; j < v.size(); j++) {
         id_type id = v[j]-config::base_index();
         if (pmesh->convex_index().is_in(id)) {
           pmesh->sup_convex(id);
         } else {
           THROW_ERROR("Can't delete convex " << id+config::base_index()
                       << ", it is not part of the mesh");
         }
       }
       );


    /*@SET ('del convex of dim', @ivec DIMs)
    Remove all convexes of dimension listed in `DIMs`.

    For example; ``MESH:SET('del convex of dim', [1,2])`` remove
    all line segments, triangles and quadrangles.@*/
    sub_command
      ("del convex of dim", 1, 1, 0, 0,
       dal::bit_vector bv = in.pop().to_bit_vector(NULL, 0);
       for (dal::bv_visitor_c cv(pmesh->convex_index());
            !cv.finished(); ++cv) {
         if (bv.is_in(pmesh->structure_of_convex(cv)->dim()))
           pmesh->sup_convex(cv);
       }
       );


    /*@SET ('translate', @vec V)
      Translates each point of the mesh from `V`.@*/
    sub_command
      ("translate", 1, 1, 0, 0,
       check_empty_mesh(pmesh);
       darray v = in.pop().to_darray(pmesh->dim(),1);
       pmesh->translation(v.col_to_bn(0));
       );


    /*@SET ('transform', @mat T)
    Applies the matrix `T` to each point of the mesh.

    Note that `T` is not required to be a ``NxN`` matrix (with
    ``N = MESH:GET('dim')``). Hence it is possible to transform
    a 2D mesh into a 3D one (and reciprocally).@*/
    sub_command
      ("transform", 1, 1, 0, 0,
       check_empty_mesh(pmesh);
       darray v = in.pop().to_darray(-1,-1); //pmesh->dim());
       pmesh->transformation(v.row_col_to_bm());
       );


    /*@SET ('boundary', @int rnum, @dmat CVFIDs)
    DEPRECATED FUNCTION. Use 'region' instead.@*/
    sub_command
      ("boundary", 2, 2, 0, 0,
       set_region(*pmesh, in);
       );


    /*@SET ('region', @int rnum, @dmat CVFIDs)
    Assigns the region number `rnum` to the set of convexes or/and convex
    faces provided in the matrix `CVFIDs`.

    The first row of `CVFIDs` contains convex #ids, and the second row
    contains a face number in the convex (or @MATLAB{0}@SCILAB{0}@PYTHON{``-1``}
    for the whole convex (regions are usually used to store a list of
    convex faces, but you may also use them to store a list of convexes).

    If a vector is provided (or a one row matrix) the region will represent
    the corresponding set of convex.@*/
    sub_command
      ("region", 2, 2, 0, 0,
       set_region(*pmesh, in);
       );

    /*@SET ('extend region', @int rnum, @dmat CVFIDs)
    Extends the region identified by the region number `rnum` to include
    the set of convexes or/and convex faces provided in the matrix
    `CVFIDs`, see also ``MESH:SET('set region)``.@*/
    sub_command
      ("extend region", 2, 2, 0, 0,
       set_region(*pmesh, in, false);
       );

    /*@SET ('region intersect', @int r1, @int r2)
    Replace the region number `r1` with its intersection with region number `r2`.@*/
    sub_command
      ("region intersect", 2, 2, 0, 0,
       intersect_regions(*pmesh,in);
       );


    /*@SET ('region merge', @int r1, @int r2)
    Merge region number `r2` into region number `r1`.@*/
    sub_command
      ("region merge", 2, 2, 0, 0,
       merge_regions(*pmesh,in);
       );


    /*@SET ('region subtract', @int r1, @int r2)
    Replace the region number `r1` with its difference with region
    number `r2`.@*/
    sub_command
      ("region subtract", 2, 2, 0, 0,
       subtract_regions(*pmesh,in);
       );


    /*@SET ('delete boundary', @int rnum, @dmat CVFIDs)
    DEPRECATED FUNCTION. Use 'delete region' instead.@*/
    sub_command
      ("delete boundary", 1, 1, 0, 0,
       dal::bit_vector lst = in.pop().to_bit_vector(&pmesh->regions_index(),0);
       for (dal::bv_visitor b(lst); !b.finished(); ++b)
         pmesh->sup_region(b);
       );


    /*@SET ('delete region', @ivec RIDs)
      Remove the regions whose #ids are listed in `RIDs`@*/
    sub_command
      ("delete region", 1, 1, 0, 0,
       dal::bit_vector lst = in.pop().to_bit_vector(&pmesh->regions_index(),0);
       // pmesh->sup_region(1); ??
       for (dal::bv_visitor b(lst); !b.finished(); ++b)
         pmesh->sup_region(b);
       );


    /*@SET ('merge', @tmesh m2[, @scalar  tol])
      Merge with the @tmesh `m2`.
      
      Overlapping points, within a tolerance radius `tol`, will not be
      duplicated. If `m2` is a @tmf object, its linked mesh will be used.@*/
    sub_command
      ("merge", 1, 2, 0, 0,
       const getfem::mesh *pmesh2 = extract_mesh_object(in.pop());
       scalar_type tol(0);
       if (in.remaining()) tol = in.pop().to_scalar();
       for (dal::bv_visitor cv(pmesh2->convex_index()); !cv.finished(); ++cv)
         pmesh->add_convex_by_points(pmesh2->trans_of_convex(cv),
                                     pmesh2->points_of_convex(cv).begin(),
                                     tol);
       );


    /*@SET ('optimize structure'[, @int with_renumbering])
    Reset point and convex numbering.

    After optimisation, the points (resp. convexes) will
    be consecutively numbered from @MATLAB{1 to MESH:GET('max pid')
    (resp. MESH:GET('max cvid'))}@SCILAB{1 to MESH:GET('max pid')
    (resp. MESH:GET('max cvid'))}@PYTHON{``0`` to
    ``MESH:GET('max pid')-1`` (resp. ``MESH:GET('max cvid')-1``)}.@*/
    sub_command
      ("optimize structure", 0, 1, 0, 0,
       bool with_renumbering = true;
       if (in.remaining()) with_renumbering = (in.pop().to_integer(0,1) != 0);
       pmesh->optimize_structure(with_renumbering);
       );


    /*@SET ('refine'[, @ivec CVIDs])
    Use a Bank strategy for mesh refinement.

    If `CVIDs` is not given, the whole mesh is refined. Note
    that the regions, and the finite element methods and
    integration methods of the @tmf and @tmim objects linked
    to this mesh will be automagically refined.@*/
    sub_command
      ("refine", 0, 1, 0, 0,
       dal::bit_vector bv = pmesh->convex_index();
       if (in.remaining()) bv = in.pop().to_bit_vector(&pmesh->convex_index());
       pmesh->Bank_refine(bv);
       );

  }


  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");
  getfem::mesh *pmesh = to_mesh_object(m_in.pop());
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
              it->second->arg_in_max, it->second->arg_out_min,
              it->second->arg_out_max);
    it->second->run(m_in, m_out, pmesh);
  }
  else bad_cmd(init_cmd);

}
