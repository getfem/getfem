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
#include <map>
#include <getfemint_misc.h>
#include <getfem/getfem_export.h>
#include <getfem/getfem_mesh.h>

using namespace getfemint;

static void check_empty_mesh(const getfem::mesh *pmesh)
{
  if (pmesh->dim() == bgeot::dim_type(-1) || pmesh->dim() == 0) {
    THROW_ERROR( "mesh object has an invalid dimension");
  }
}

static void
get_structure_or_geotrans_of_convexes(const getfem::mesh& m,
                                      mexargs_in& in, mexargs_out& out,
                                      id_type class_id)
{
  dal::bit_vector cvlst;
  if (in.remaining()) cvlst = in.pop().to_bit_vector(&m.convex_index());
  else cvlst = m.convex_index();
  std::vector<id_type> ids; ids.reserve(cvlst.card());
  for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv) {
    if (class_id == CVSTRUCT_CLASS_ID)
      ids.push_back(store_cvstruct_object(m.structure_of_convex(cv)));
    else ids.push_back(store_geotrans_object(m.trans_of_convex(cv)));
  }
  out.return_packed_obj_ids(ids, class_id);
}

/* apply geometric transformation to edge lists */
static void
transform_edge_list(const getfem::mesh &m, unsigned N, const bgeot::edge_list &el, darray &w)
{
  bgeot::edge_list::const_iterator it;
  size_type cv = size_type(-1);
  getfem::mesh::ref_convex cv_ref;
  bgeot::pgeometric_trans pgt;
  size_type ecnt = 0;
  for (it = el.begin(); it != el.end(); it++, ecnt++) {
    if (cv != (*it).cv) {
      cv = (*it).cv;
      cv_ref = m.convex(cv);
      pgt = m.trans_of_convex(cv);
    }
    /* build the list of points on the edge, on the reference element */
    /* get local point numbers in the convex */
    bgeot::size_type iA = m.local_ind_of_convex_point(cv, (*it).i);
    bgeot::size_type iB = m.local_ind_of_convex_point(cv, (*it).j);

    if (iA >= cv_ref.nb_points() ||
        iB >= cv_ref.nb_points() || iA == iB) {
      THROW_INTERNAL_ERROR;
    }
    getfem::base_node A = pgt->convex_ref()->points()[iA];
    getfem::base_node B = pgt->convex_ref()->points()[iB];
    for (size_type i = 0; i < N; i++) {
      getfem::base_node pt;

      pt = A +  (B-A)*(double(i)/double(N-1));
      pt = pgt->transform(pt, m.points_of_convex(cv));
      std::copy(pt.begin(), pt.end(), &w(0,i, ecnt));
    }
  }
}

typedef dal::dynamic_tree_sorted<convex_face> mesh_faces_list;

static void
faces_from_pid(const getfem::mesh &m, mexargs_in &in, mexargs_out &out)
{
  mesh_faces_list lst;
  dal::bit_vector convex_tested, pts = in.pop().to_bit_vector(&m.points().index());

  for (dal::bv_visitor ip(pts); !ip.finished(); ++ip) {
    /* iterator over all convexes attached to point ip */
    bgeot::mesh_structure::ind_cv_ct::const_iterator
      cvit = m.convex_to_point(ip).begin(),
      cvit_end = m.convex_to_point(ip).end();

    for ( ; cvit != cvit_end; cvit++) {
      size_type ic = *cvit;
      if (!convex_tested.is_in(ic)) {
        convex_tested.add(ic);
        /* loop over the convex faces */
        for (short_type f = 0; f < m.structure_of_convex(ic)->nb_faces(); f++) {
          bgeot::mesh_structure::ind_pt_face_ct pt = m.ind_points_of_face_of_convex(ic, f);
          bool ok = true;
          for (unsigned ii=0; ii < pt.size(); ++ii) {
            if (!pts.is_in(pt[ii])) {
              ok = false; break;
            }
          }
          if (ok) lst.add(convex_face(ic,f));
        }
      }
    }
  }
  iarray w = out.pop().create_iarray(2,unsigned(lst.size()));
  for (size_type j=0; j < lst.size(); j++) {
    w(0,j) = int(lst[j].cv+config::base_index());
    w(1,j) = int(short_type(lst[j].f+config::base_index()));
  }
}


struct mesh_faces_by_pts_list_elt  {
  std::vector<size_type> ptid; // point numbers of faces
  int cnt; // number of convexes sharing that face
  size_type cv;
  short_type f;
  inline bool operator < (const mesh_faces_by_pts_list_elt &e) const
  {
    if (ptid.size() < e.ptid.size()) return true;
    if (ptid.size() > e.ptid.size()) return false;
    return ptid < e.ptid;
  }

  mesh_faces_by_pts_list_elt(size_type cv_, short_type f_,
                             std::vector<size_type>& p) : cv(cv_), f(f_) {
    cnt = 0;
    if (p.size() == 0) THROW_INTERNAL_ERROR;
    std::sort(p.begin(), p.end());
//     cerr << "ajout cv=" << cv << " f=" << f << " pt=[";
//     std::copy(p.begin(), p.end(), std::ostream_iterator<size_type,char>(cerr, ","));
//     cerr << "]\n";
    ptid = p;
  }
  mesh_faces_by_pts_list_elt() {}
};
typedef dal::dynamic_tree_sorted<mesh_faces_by_pts_list_elt> mesh_faces_by_pts_list;


static void
outer_faces(const getfem::mesh &m, mexargs_in &in, mexargs_out &out,
            const std::string &condition="")
{
  mesh_faces_by_pts_list lst;
  dal::bit_vector cvlst, checked_pids, rejected_pids;

  bool with_normal(condition == "direction");
  bool with_box(condition == "box");
  bool with_ball(condition == "ball");

  darray normal_vector;
  bgeot::base_node un, pmin, pmax, center;
  scalar_type threshold(0), radius(0);
  if (with_normal) {
    normal_vector = in.pop().to_darray();
    scalar_type angle = in.pop().to_scalar();
    threshold = gmm::vect_norm2(normal_vector)*cos(angle);
  } else if (with_box) {
    darray p1 = in.pop().to_darray(m.dim());
    darray p2 = in.pop().to_darray(m.dim());
    pmin.resize(m.dim());
    pmax.resize(m.dim());
    for (size_type k=0; k < m.dim(); ++k) {
      pmin[k] = std::min(p1[k],p2[k]);
      pmax[k] = std::max(p1[k],p2[k]);
    }
  } else if (with_ball) {
    darray p1 = in.pop().to_darray(m.dim());
    center.resize(m.dim());
    for (size_type k=0; k < m.dim(); ++k) {
      center[k] = p1[k];
    }
    radius = in.pop().to_scalar();
  }

  dim_type elm_dim = m.dim();
  if (in.remaining() && in.front().is_integer())
   elm_dim = dim_type(in.pop().to_integer());

  if (in.remaining()) cvlst = in.pop().to_bit_vector(&m.convex_index());
  else cvlst = m.convex_index();

  for (dal::bv_visitor ic(cvlst); !ic.finished(); ++ic) {
    if (m.structure_of_convex(ic)->dim() == elm_dim) {
      for (short_type f = 0; f < m.structure_of_convex(ic)->nb_faces(); f++) {
        bgeot::mesh_structure::ind_pt_face_ct pt
          = m.ind_points_of_face_of_convex(ic, f);
        std::vector<size_type> p(pt.begin(), pt.end());
        size_type idx = lst.add_norepeat(mesh_faces_by_pts_list_elt(ic,f,p));
        lst[idx].cnt++;
      }
    }
//     else { // DEPRECATED
//       bgeot::mesh_structure::ind_cv_ct pt = m.ind_points_of_convex(ic);
//       std::vector<size_type> p(pt.begin(), pt.end());
//       size_type idx =
//         lst.add_norepeat(mesh_faces_by_pts_list_elt(ic, short_type(-1), p));
//       lst[idx].cnt++;
//     }
  }
  size_type fcnt = 0;
  for (size_type i = 0; i < lst.size(); i++) {
    if (lst[i].cnt == 1) {
      if (with_normal) {
        un = m.mean_normal_of_face_of_convex(lst[i].cv, lst[i].f);
        if (gmm::vect_sp(normal_vector, un) < threshold - 1E-8) {
          lst[i].cnt = -1;
        }
      } else if (with_box) {
        for (std::vector<size_type>::const_iterator pid=lst[i].ptid.begin();
             pid != lst[i].ptid.end(); ++pid) {
          if (!checked_pids.is_in(*pid)) {
            checked_pids.add(*pid);
            for (size_type k=0; k < m.dim(); ++k)
              if (m.points()[*pid][k] < pmin[k] ||
                  m.points()[*pid][k] > pmax[k]) {
                rejected_pids.add(*pid);
                break;
              }
          }
          if (rejected_pids.is_in(*pid)) {
            lst[i].cnt = -1;
            break;
          }
        }
      } else if (with_ball) {
        for (std::vector<size_type>::const_iterator pid=lst[i].ptid.begin();
             pid != lst[i].ptid.end(); ++pid) {
          if (!checked_pids.is_in(*pid)) {
            checked_pids.add(*pid);
            scalar_type checked_radius = scalar_type(0.0);
            for (size_type k=0; k < m.dim(); ++k) {
              checked_radius += pow(m.points()[*pid][k] - center[k], 2);
            }
            checked_radius = std::sqrt(checked_radius);
            for (size_type k=0; k < m.dim(); ++k)
              if (checked_radius > radius) {
                rejected_pids.add(*pid);
                break;
              }
          }
          if (rejected_pids.is_in(*pid)) {
            lst[i].cnt = -1;
            break;
          }
        }
      }
      if (lst[i].cnt == 1) fcnt++;
    }

  }

  iarray w = out.pop().create_iarray(2,unsigned(fcnt));
  fcnt = 0;
  for (size_type j=0; j < lst.size(); j++) {
    if (lst[j].cnt == 1) {
      w(0,fcnt) = int(lst[j].cv+config::base_index());
      w(1,fcnt) = int(short_type(lst[j].f+config::base_index()));
      fcnt++;
    }
  }

}

static void 
all_faces(const getfem::mesh &m, mexargs_in &in, mexargs_out &out) {
  dal::bit_vector cvlst;

  if (in.remaining()) cvlst = in.pop().to_bit_vector(&m.convex_index());
  else cvlst = m.convex_index();
  getfem::mesh_region mr;
  for (dal::bv_visitor ic(cvlst); !ic.finished(); ++ic) mr.add(ic);
  getfem::mesh_region mrr =  all_faces_of_mesh(m, mr);

  unsigned fcnt = 0;
  for (getfem::mr_visitor i(mrr); !i.finished(); ++i) ++fcnt;
  iarray w = out.pop().create_iarray(2, fcnt);
  fcnt = 0;
  for (getfem::mr_visitor i(mrr); !i.finished(); ++i) {
    w(0,fcnt) = int(i.cv()+config::base_index());
    w(1,fcnt) = int(short_type(i.f()+config::base_index()));
    fcnt++;
  }
}



static void 
inner_faces(const getfem::mesh &m, mexargs_in &in, mexargs_out &out) {
  dal::bit_vector cvlst;

  if (in.remaining()) cvlst = in.pop().to_bit_vector(&m.convex_index());
  else cvlst = m.convex_index();
  getfem::mesh_region mr;
  for (dal::bv_visitor ic(cvlst); !ic.finished(); ++ic) mr.add(ic);
  getfem::mesh_region mrr =  inner_faces_of_mesh(m, mr);

  unsigned fcnt = 0;
  for (getfem::mr_visitor i(mrr); !i.finished(); ++i) ++fcnt;
  iarray w = out.pop().create_iarray(2, fcnt);
  fcnt = 0;
  for (getfem::mr_visitor i(mrr); !i.finished(); ++i) {
    w(0,fcnt) = int(i.cv()+config::base_index());
    w(1,fcnt) = int(short_type(i.f()+config::base_index()));
    fcnt++;
  }
}


static bgeot::base_node
normal_of_face(const getfem::mesh& mesh, size_type cv, short_type f, size_type node) {
  if (!mesh.convex_index().is_in(cv)) THROW_BADARG("convex " << cv+1 << " not found in mesh");
  if (f >= mesh.structure_of_convex(cv)->nb_faces())
    THROW_BADARG("convex " << cv+1 << " has only " <<
                 mesh.structure_of_convex(cv)->nb_faces() <<
                 ": can't find face " << f+1);
  if (node >= mesh.structure_of_convex(cv)->nb_points_of_face(f))
    THROW_BADARG("invalid node number: " << node);
  bgeot::base_node N = mesh.normal_of_face_of_convex(cv, f, node);
  N /= gmm::vect_norm2(N);
  gmm::clean(N, 1e-14);
  return N;
}



// Object for the declaration of a new sub-command.

struct sub_gf_mesh_set : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
                   getfemint::mexargs_out& out,
                   const getfem::mesh *pmesh) = 0;
};

typedef std::shared_ptr<sub_gf_mesh_set> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_mesh_set {				\
      virtual void run(getfemint::mexargs_in& in,			\
                       getfemint::mexargs_out& out,			\
                       const getfem::mesh *pmesh)			\
      { dummy_func(in); dummy_func(out); code }				\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }


/*@GFDOC
  General mesh inquiry function. All these functions accept also a
  mesh_fem argument instead of a mesh M (in that case, the mesh_fem
  linked mesh will be used). @MATLAB{Note that when your mesh is
  recognized as a Matlab object , you can simply use "get(M, 'dim')"
  instead of "gf_mesh_get(M, 'dim')".}
@*/

void gf_mesh_get(getfemint::mexargs_in& m_in,
                 getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {

    /*@RDATTR d = ('dim')
    Get the dimension of the mesh (2 for a 2D mesh, etc).@*/
    sub_command
      ("dim", 0, 0, 0, 1,
       out.pop().from_integer(pmesh->dim());
       );


    /*@RDATTR np = ('nbpts')
    Get the number of points of the mesh.@*/
    sub_command
      ("nbpts", 0, 0, 0, 1,
       out.pop().from_integer(int(pmesh->nb_points()));
       );


    /*@RDATTR nc = ('nbcvs')
      Get the number of convexes of the mesh.@*/
    sub_command
      ("nbcvs", 0, 0, 0, 1,
       out.pop().from_integer(int(pmesh->nb_convex()));
       );


    /*@GET P = ('pts'[, @ivec PIDs])
    Return the list of point coordinates of the mesh.

    Each column of the returned matrix contains the coordinates of one
    point. If the optional argument `PIDs` was given, only the points
    whose #id is listed in this vector are returned. Otherwise, the
    returned matrix will have MESH:GET('max_pid') columns, which might
    be greater than MESH:GET('nbpts') (if some points of the mesh have
    been destroyed and no call to MESH:SET('optimize structure') have
    been issued). The columns corresponding to deleted points will be
    filled with NaN. You can use MESH:GET('pid') to filter such invalid
    points.@*/
    sub_command
      ("pts", 0, 1, 0, 1,
       double nan = ::nan("");
       if (!in.remaining()) {
         dal::bit_vector bv = pmesh->points().index();
         darray w = out.pop().create_darray(pmesh->dim(), unsigned(bv.last_true()+1));
         for (size_type j = 0; j < bv.last_true()+1; j++) {
           for (size_type i = 0; i < pmesh->dim(); i++) {
             w(i,j) = (bv.is_in(j)) ? (pmesh->points()[j])[i] : nan;
           }
         }
       } else {
         dal::bit_vector pids = in.pop().to_bit_vector();
         darray w = out.pop().create_darray(pmesh->dim(), unsigned(pids.card()));
         size_type cnt=0;
         for (dal::bv_visitor j(pids); !j.finished(); ++j) {
           if (!pmesh->points().index().is_in(j))
             THROW_ERROR("point " << j+config::base_index() << " is not part of the mesh");
           for (size_type i = 0; i < pmesh->dim(); i++) {
             w(i,cnt) = (pmesh->points()[j])[i];
           }
           cnt++;
         }
       }
       );


    /*@GET Pid = ('pid')
    Return the list of points #id of the mesh.

    Note that their numbering is not supposed to be contiguous from
    @MATLAB{1 to MESH:GET('nbpts')}@PYTHON{0 to MESH:GET('nbpts')-1}@SCILAB{1 to MESH:GET('nbpts')},
    especially if some points have been removed from the mesh. You
    can use MESH:SET('optimize_structure') to enforce a contiguous
    numbering. @MATLAB{Pid is a row vector.}@*/
    sub_command
      ("pid", 0, 0, 0, 1,
       out.pop().from_bit_vector(pmesh->points().index());
       );


    /*@GET PIDs = ('pid in faces', @imat CVFIDs)
    Return point #id listed in `CVFIDs`.

    `CVFIDs` is a two-rows matrix, the first row lists convex #ids,
    and the second lists face numbers. On return, `PIDs` is a
    @MATLAB{row }vector containing points #id.@*/
    sub_command
      ("pid in faces", 1, 1, 0, 1,
       check_empty_mesh(pmesh);
       
       iarray v = in.pop().to_iarray(2,-1);
       
       dal::bit_vector pids;
       for (size_type j=0; j < v.getn(); j++) {
         size_type cv = v(0,j) - config::base_index();
         short_type f = short_type(v(1,j) - config::base_index());
         
         if (pmesh->convex_index().is_in(cv)) {
           if (short_type(-1)==f){
             for (unsigned i=0; i < pmesh->nb_points_of_convex(cv); ++i)
               pids.add(pmesh->ind_points_of_convex(cv)[i]);
           } else if (f < pmesh->structure_of_convex(cv)->nb_faces()) {
             for (unsigned i=0; i < pmesh->structure_of_convex(cv)->nb_points_of_face(f); ++i)
               pids.add(pmesh->ind_points_of_face_of_convex(cv,f)[i]);
           }
         }
       }

       out.pop().from_bit_vector(pids);
       );


    /*@GET PIDs = ('pid in cvids', @imat CVIDs)
      Return point #id listed in `CVIDs`.
      
      `PIDs` is a @MATLAB{row }vector containing points #id.@*/
    sub_command
      ("pid in cvids", 1, 1, 0, 1,
       check_empty_mesh(pmesh);
       
       dal::bit_vector cvlst = in.pop().to_bit_vector();
       dal::bit_vector pids;
       
       for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv)
         if (pmesh->convex_index().is_in(cv))
           for (unsigned i=0; i < pmesh->nb_points_of_convex(cv); ++i)
             pids.add(pmesh->ind_points_of_convex(cv)[i]);
       
       out.pop().from_bit_vector(pids);
       );


    /*@GET PIDs = ('pid in regions', @imat RIDs)
    Return point #id listed in `RIDs`.

    `PIDs` is a @MATLAB{row }vector containing points #id.@*/
    sub_command
      ("pid in regions", 1, 1, 0, 1,
       check_empty_mesh(pmesh);
       
       dal::bit_vector rlst = in.pop().to_bit_vector(&pmesh->regions_index(), 0);
       dal::bit_vector pids;
       
       for (dal::bv_visitor r(rlst); !r.finished(); ++r) {
         if (pmesh->regions_index().is_in(r)) {
           for (getfem::mr_visitor i(pmesh->region(r)); !i.finished(); ++i) {
             if (short_type(-1)==i.f()) {
               for (unsigned j=0; j < pmesh->nb_points_of_convex(i.cv()); ++j)
                 pids.add(pmesh->ind_points_of_convex(i.cv())[j]);
             } else {
               for (unsigned j=0; j < pmesh->structure_of_convex(i.cv())->nb_points_of_face(i.f()); ++j)
                 pids.add(pmesh->ind_points_of_face_of_convex(i.cv(),i.f())[j]);
             }
           }
         }
       }
       
       out.pop().from_bit_vector(pids);
       );


    /*@GET PIDs = ('pid from coords', @mat PTS[, @scalar radius=0])
    Return point #id whose coordinates are listed in `PTS`.

    `PTS` is an array containing a list of point coordinates. On
    return, `PIDs` is a @MATLAB{row }vector containing points
    #id for each point found in `eps` range, and -1 for those
    which where not found in the mesh.@*/
    sub_command
      ("pid from coords", 1, 2, 0, 1,
       check_empty_mesh(pmesh);
       darray v = in.pop().to_darray(pmesh->dim(), -1);
       scalar_type radius = 0;
       if (in.remaining()) radius = in.pop().to_scalar(0);
       iarray w = out.pop().create_iarray_h(v.getn());
       for (unsigned j=0; j < v.getn(); j++) {
         getfem::base_node P = v.col_to_bn(j);
         getfem::size_type id = size_type(-1);
         if (!std::isnan(P[0])) id = pmesh->search_point(P,radius);
         if (id == getfem::size_type(-1)) w[j] = -1;
         else w[j] = int(id + config::base_index());
       }
       );


    /*@GET @CELL{Pid, IDx} = ('pid from cvid'[, @imat CVIDs])
    Return the points attached to each convex of the mesh.

    If `CVIDs` is omitted, all the convexes will be considered
    (equivalent to `CVIDs = MESH:GET('max cvid')`). `IDx` is a
    @MATLAB{row }vector, length(IDx) = length(CVIDs)+1. `Pid` is a
    @MATLAB{row }vector containing the concatenated list of #id of
    points of each convex in `CVIDs`. Each entry of `IDx` is the
    position of the corresponding convex point list in `Pid`. Hence,
    for example, the list of #id of points of the second convex is
    @MATLAB{Pid(IDx(2):IDx(3)-1)}@SCILAB{Pid(IDx(2):IDx(3)-1)}@PYTHON{Pid[IDx(2):IDx(3)]}.

    If `CVIDs` contains convex #id which do not exist in the mesh,
    their point list will be empty.@*/
    sub_command
      ("pid from cvid", 0, 1, 0, 2,
       dal::bit_vector cvlst;
       if (in.remaining()) cvlst = in.pop().to_bit_vector();
       else cvlst.add(0, pmesh->convex_index().last_true()+1);
       
       std::vector<size_type> pids;
       std::vector<size_type> idx;
       size_type pcnt = 0;
       for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv) {
         idx.push_back(size_type(pcnt + config::base_index()));
         if (pmesh->convex_index().is_in(cv))
           for (size_type i=0; i < pmesh->nb_points_of_convex(cv); ++i,++pcnt)
             pids.push_back(size_type(pmesh->ind_points_of_convex(cv)[i] + config::base_index()));
       }
       idx.push_back(size_type(pcnt + config::base_index()));
       
       iarray opids = out.pop().create_iarray_h(unsigned(pids.size()));
       if (pids.size()) std::copy(pids.begin(), pids.end(), &opids[0]);
       if (out.remaining() && idx.size()) {
         iarray oidx = out.pop().create_iarray_h(unsigned(idx.size()));
         std::copy(idx.begin(), idx.end(), &oidx[0]);
       }
       );


    /*@GET @CELL{Pts, IDx} = ('pts from cvid'[, @imat CVIDs])
    Search point listed in `CVID`.

    Return `Pts` and `IDx`.
    If `CVIDs` is omitted, all the convexes will be considered
    (equivalent to `CVIDs = MESH:GET('max cvid')`). `IDx` is a
    @MATLAB{row }vector, length(IDx) = length(CVIDs)+1. `Pts` is a
    @MATLAB{row }vector containing the concatenated list of points
    of each convex in `CVIDs`. Each entry of `IDx` is the position
    of the corresponding convex point list in `Pts`. Hence, for
    example, the list of points of the second convex is
    @MATLAB{Pts(:,IDx(2):IDx(3)-1)}@SCILAB{Pts(:,IDx(2):IDx(3)-1)}@PYTHON{Pts[:,IDx[2]:IDx[3]]}.

    If `CVIDs` contains convex #id which do not exist in the mesh,
    their point list will be empty.@*/
    sub_command
      ("pts from cvid", 0, 1, 0, 2,
       check_empty_mesh(pmesh);
       
       dal::bit_vector cvlst;
       if (in.remaining()) cvlst = in.pop().to_bit_vector();
       else cvlst.add(0, pmesh->convex_index().last_true()+1);
       
       getfem::base_vector pts;
       std::vector<size_type> idx;
       size_type pcnt = 0;
       for (dal::bv_visitor cv(cvlst); !cv.finished(); ++cv) {
         idx.push_back(pcnt + config::base_index());
         if (pmesh->convex_index().is_in(cv))
           for (size_type i=0; i < pmesh->nb_points_of_convex(cv); ++i,++pcnt)
             for (size_type k=0; k< pmesh->dim(); ++k)
               pts.push_back(pmesh->points_of_convex(cv)[i][k]);
       }
       idx.push_back(pcnt + config::base_index());
       
       darray opts = out.pop().create_darray(pmesh->dim(),
                                    unsigned(pts.size()/pmesh->dim()));
       if (pts.size()) std::copy(pts.begin(), pts.end(), &opts[0]);
       if (out.remaining() && idx.size()) {
         iarray oidx = out.pop().create_iarray_h(unsigned(idx.size()));
         std::copy(idx.begin(), idx.end(), &oidx[0]);
       }
       );


    /*@GET CVid = ('cvid')
    Return the list of all convex #id.

    Note that their numbering is not supposed to be contiguous from
    @MATLAB{1 to MESH:GET('nbcvs')}@SCILAB{1 to MESH:GET('nbcvs')}@PYTHON{0 to MESH:GET('nbcvs')-1},
    especially if some points have been removed from the mesh. You
    can use MESH:SET('optimize_structure') to enforce a contiguous
    numbering. @MATLAB{CVid is a row vector.}@*/
    sub_command
      ("cvid", 0, 0, 0, 1,
       out.pop().from_bit_vector(pmesh->convex_index());
       );


    /*@GET m = ('max pid')
      Return the maximum #id of all points in the mesh (see 'max cvid').@*/
    sub_command
      ("max pid", 0, 0, 0, 1,
       int i = pmesh->points().index().card() ? int(pmesh->points().index().last_true()) : -1;
       out.pop().from_integer(i + config::base_index());
       );


    /*@GET m = ('max cvid')
      Return the maximum #id of all convexes in the mesh (see 'max pid').@*/
    sub_command
      ("max cvid", 0, 0, 0, 1,
       int i = pmesh->convex_index().card()? int(pmesh->convex_index().last_true()) : -1;
       out.pop().from_integer(i + config::base_index());
       );


    /*@GET [E,C] = ('edges' [, CVLST][, 'merge'])
    [OBSOLETE FUNCTION! will be removed in a future release]

    Return the list of edges of mesh M for the convexes listed in the
    row vector CVLST. E is a 2 x nb_edges matrix containing point
    indices. If CVLST is omitted, then the edges of all convexes are
    returned. If CVLST has two rows then the first row is supposed to
    contain convex numbers, and the second face numbers, of which the
    edges will be returned.  If 'merge' is indicated, all common
    edges of convexes are merged in a single edge.  If the optional
    output argument C is specified, it will contain the convex number
    associated with each edge.@*/
    sub_command
      ("edges", 0, 2, 0, 2,
       bgeot::edge_list el;
       build_edge_list(*pmesh, el, in);
       iarray w   = out.pop().create_iarray(2, unsigned(el.size()));
       /* copy the edge list to the matlab array */
       for (size_type j = 0; j < el.size(); j++) {
         w(0,j) = int(el[j].i + config::base_index());
         w(1,j) = int(el[j].j + config::base_index());
       }
       if (out.remaining()) {
         iarray cv = out.pop().create_iarray_h(unsigned(el.size()));
         for (size_type j = 0; j < el.size(); j++) {
           cv[j] = int(el[j].cv+config::base_index());
         }
       }
       );


    /*@GET [E,C] = ('curved edges', @int N [, CVLST])
    [OBSOLETE FUNCTION! will be removed in a future release]

    Return E and C.
    More sophisticated version of MESH:GET('edges') designed for
    curved elements. This one will return N (N>=2) points of the
    (curved) edges. With N==2, this is equivalent to
    MESH:GET('edges'). Since the points are no more always part of
    the mesh, their coordinates are returned instead of points
    number, in the array E which is a [ mesh_dim x 2 x nb_edges ]
    array.  If the optional output argument C is specified, it will
    contain the convex number associated with each edge.@*/
    sub_command
      ("curved edges", 1, 2, 0, 2,
       bgeot::edge_list el;
       size_type N = in.pop().to_integer(2,10000);
       build_edge_list(*pmesh, el, in);
       darray w   = out.pop().create_darray(pmesh->dim(), unsigned(N),
                                            unsigned(el.size()));
       transform_edge_list(*pmesh, unsigned(N), el, w);
       if (out.remaining()) {
         iarray cv = out.pop().create_iarray_h(unsigned(el.size()));
         for (size_type j = 0; j < el.size(); j++) {
           cv[j] = int(el[j].cv + config::base_index());
         }
       }
       );


    /*@GET PIDs = ('orphaned pid')
      Return point #id which are not linked to a convex.@*/
    sub_command
      ("orphaned pid", 0, 0, 0, 1,
       dal::bit_vector bv = pmesh->points().index();
       for (dal::bv_visitor cv(pmesh->convex_index()); !cv.finished(); ++cv) {
         for (unsigned i=0; i < pmesh->nb_points_of_convex(cv); ++i)
           bv.sup(pmesh->ind_points_of_convex(cv)[i]);
       }
       out.pop().from_bit_vector(bv);
       );


    /*@GET CVIDs = ('cvid from pid', @ivec PIDs[, @bool share=False])
    Return convex #ids related with the point #ids given in `PIDs`.
    
    If `share=False`, search convex whose vertex #ids are in `PIDs`.
    If `share=True`, search convex #ids that share the point #ids
    given in `PIDs`. `CVIDs` is a @MATLAB{row} vector (possibly
    empty).@*/
    sub_command
      ("cvid from pid", 1, 2, 0, 1,
       check_empty_mesh(pmesh);
       dal::bit_vector pts = in.pop().to_bit_vector(&pmesh->points().index());
       dal::bit_vector cvchecked;
       
       bool share = false;
       if (in.remaining() && in.front().is_bool())
         share = in.pop().to_bool();
       
       std::vector<size_type> cvlst;
       
       /* loop over the points */
       for (dal::bv_visitor ip(pts); !ip.finished(); ++ip) {
         /* iterators over the convexes attached to point ip */
         bgeot::mesh_structure::ind_cv_ct::const_iterator
           cvit = pmesh->convex_to_point(ip).begin();
         bgeot::mesh_structure::ind_cv_ct::const_iterator
           cvit_end = pmesh->convex_to_point(ip).end();
         
         for ( ; cvit != cvit_end; cvit++) {
           size_type ic = *cvit;
           
           //cerr << "cv = " << ic+1 << endl;
           
           if (!cvchecked.is_in(ic)) {
             if (share) cvlst.push_back(ic);
             else { /* check that each point of the convex is in the list */
               bool ok = true;
               bgeot::mesh_structure::ind_cv_ct cvpt = pmesh->ind_points_of_convex(ic);
               for (unsigned ii=0; ii < cvpt.size(); ++ii) {
                 if (!pts.is_in(cvpt[ii])) {
                   ok = false; break;
                 }
               }
               if (ok) cvlst.push_back(ic);
             }
             cvchecked.add(ic);
           }
         }
       }
       iarray w = out.pop().create_iarray_h(unsigned(cvlst.size()));
       for (size_type j=0; j < cvlst.size(); j++)
         w[j] = int(cvlst[j]+config::base_index());
       );


    /*@GET CVFIDs = ('faces from pid', @ivec PIDs)
    Return the convex faces whose vertex #ids are in `PIDs`.

    `CVFIDs` is a two-rows matrix, the first row lists convex #ids,
    and the second lists face numbers (local number in the convex).
    For a convex face to be returned, EACH of its points have to be
    listed in `PIDs`.@*/
    sub_command
      ("faces from pid", 1, 1, 0, 1,
       check_empty_mesh(pmesh);
       faces_from_pid(*pmesh, in, out);
       );


    /*@GET CVFIDs = ('outer faces'[, dim][, CVIDs])
    Return the set of faces not shared by two elements.

    The output `CVFIDs` is a two-rows matrix, the first row lists
    convex #ids, and the second one lists face numbers (local number
    in the convex). If `dim` is provided, the function is forced to
    detect faces of elements that have dimension `dim`, e.g. `dim`=2 will
    detect edges of surface elements, even if these belong to a 3D mesh.
    If `CVIDs` is not given, all convexes are considered, and the
    function basically returns the mesh boundary. If `CVIDs`
    is given, it returns the boundary of the convex set whose #ids are
    listed in `CVIDs`.@*/
    sub_command
      ("outer faces", 0, 2, 0, 1,
       check_empty_mesh(pmesh);
       outer_faces(*pmesh, in, out);
       );

    /*@GET CVFIDs = ('inner faces'[, CVIDs])
    Return the set of faces shared at least by two elements in CVIDs.
    Each face is represented only once and is arbitrarily chosen
    between the two neighbor elements. @*/
    sub_command
      ("inner faces", 0, 1, 0, 1,
       check_empty_mesh(pmesh);
       inner_faces(*pmesh, in, out);
       );
    
    /*@GET CVFIDs = ('all faces'[, CVIDs])
    Return the set of faces of the in CVIDs (in all the mesh if CVIDs is
    omitted). Note that the face shared by two neighbor elements will be
    represented twice. @*/
    sub_command
      ("all faces", 0, 1, 0, 1,
       check_empty_mesh(pmesh);
       all_faces(*pmesh, in, out);
       );

    /*@GET CVFIDs = ('outer faces with direction', @vec v, @scalar angle[, dim][, CVIDs])
    Return the set of faces not shared by two convexes and with a mean outward vector lying within an angle `angle` (in radians) from vector `v`.

    The output `CVFIDs` is a two-rows matrix, the first row lists convex
    #ids, and the second one lists face numbers (local number in the
    convex). The argument `dim` works as in outer_faces().
    If `CVIDs` is given, it returns portion of the boundary of
    the convex set defined by the #ids listed in `CVIDs`.@*/
    sub_command
      ("outer faces with direction", 2, 4, 0, 1,
       check_empty_mesh(pmesh);
       outer_faces(*pmesh, in, out, "direction");
       );

    
    /*@GET CVFIDs = ('outer faces in box', @vec pmin, @vec pmax[, dim][, CVIDs])
    Return the set of faces not shared by two convexes and lying within the box defined by the corner points `pmin` and `pmax`.

    The output `CVFIDs` is a two-rows matrix, the first row lists convex
    #ids, and the second one lists face numbers (local number in the
    convex). The argument `dim` works as in outer_faces().
    If `CVIDs` is given, it returns portion of the boundary of
    the convex set defined by the #ids listed in `CVIDs`.@*/
    sub_command
      ("outer faces in box", 2, 4, 0, 1,
       check_empty_mesh(pmesh);
       outer_faces(*pmesh, in, out, "box");
       );

    /*@GET CVFIDs = ('outer faces in ball', @vec center, @scalar radius[, dim][, CVIDs])
    Return the set of faces not shared by two convexes and lying within the ball of corresponding `center` and `radius`.

    The output `CVFIDs` is a two-rows matrix, the first row lists convex
    #ids, and the second one lists face numbers (local number in the
    convex). The argument `dim` works as in outer_faces().
    If `CVIDs` is given, it returns portion of the boundary of
    the convex set defined by the #ids listed in `CVIDs`.@*/
    sub_command
      ("outer faces in ball", 2, 4, 0, 1,
       check_empty_mesh(pmesh);
       outer_faces(*pmesh, in, out, "ball");
       );

    /*@GET CVFIDs = ('adjacent face', @int cvid, @int fid)
    Return convex face of the neighbor element if it exists.
    If the convex have more than one neighbor
    relatively to the face ``f`` (think to bar elements in 3D for instance),
    return the first face found. @*/
    sub_command
      ("adjacent face", 2, 2, 0, 1,
       check_empty_mesh(pmesh);
       size_type cv = in.pop().to_convex_number(*pmesh);
       short_type f = in.pop().to_face_number(pmesh->structure_of_convex(cv)->nb_faces());
       bgeot::convex_face cvf = pmesh->adjacent_face(cv, f);
       getfem::mesh_region flst;
       if (cvf.cv != size_type(-1))
	 flst.add(cvf.cv,cvf.f);
       out.pop().from_mesh_region(flst);
       );

    /*@GET CVFIDs = ('faces from cvid'[, @ivec CVIDs][, 'merge'])
    Return a list of convex faces from a list of convex #id.

    `CVFIDs` is a two-rows matrix, the first row lists convex #ids,
    and the second lists face numbers (local number in the convex).
    If `CVIDs` is not given, all convexes are considered. The optional
    argument 'merge' merges faces shared by the convex of `CVIDs`.@*/
    sub_command
      ("faces from cvid", 0, 2, 0, 1,
       check_empty_mesh(pmesh);
       dal::bit_vector bv;
       if (in.remaining() && !in.front().is_string())
         bv = in.pop().to_bit_vector(&pmesh->convex_index());
       else
         bv = pmesh->convex_index();
       bool merge = false;
       if (in.remaining() && in.front().is_string()) {
         std::string s = in.pop().to_string();
         if (cmd_strmatch(s, "merge")) merge = true;
         else bad_cmd(s);
       }
       getfem::mesh_region flst;
       size_type cnt = 0;
       for (dal::bv_visitor cv(bv); !cv.finished(); ++cv, ++cnt) {
         for (short_type f=0; f < pmesh->structure_of_convex(cv)->nb_faces(); ++f) {
           bool add = true;
           if (merge) {
             bgeot::mesh_structure::ind_set neighbors;
             pmesh->neighbors_of_convex(cv, f, neighbors);
             for (bgeot::mesh_structure::ind_set::const_iterator it = neighbors.begin();
                  it != neighbors.end(); ++it) {
               if (*it < cv) { add = false; break; }
             }
           }
           if (add) flst.add(cv,f);
         }
       }
       out.pop().from_mesh_region(flst);
       );


    /*@GET [@mat T] = ('triangulated surface', @int Nrefine [,CVLIST])
    [DEPRECATED FUNCTION! will be removed in a future release]

    Similar function to MESH:GET('curved edges') : split (if
    necessary, i.e. if the geometric transformation if non-linear)
    each face into sub-triangles and return their coordinates in T
    (see also ::COMPUTE('eval on P1 tri mesh'))@*/
    sub_command
      ("triangulated surface", 1, 2, 0, 1,
       int Nrefine = in.pop().to_integer(1, 1000);
       std::vector<convex_face> cvf;
       if (in.remaining() && !in.front().is_string()) {
         iarray v = in.pop().to_iarray(-1, -1);
         build_convex_face_lst(*pmesh, cvf, &v);
       } else build_convex_face_lst(*pmesh, cvf, 0);
       eval_on_triangulated_surface(pmesh, Nrefine, cvf, out, NULL, darray());
       );


    /*@GET N = ('normal of face', @int cv, @int f[, @int nfpt])
    Return the normal vector of convex `cv`, face `f` at the `nfpt` point of the face.

    If `nfpt` is not specified, then the normal is evaluated at each
    geometrical node of the face.@*/
    sub_command
      ("normal of face", 2, 3, 0, 1,
       size_type cv = in.pop().to_convex_number(*pmesh);
       short_type f  = in.pop().to_face_number(pmesh->structure_of_convex(cv)->nb_faces());
       size_type node = 0;
       if (in.remaining()) node = in.pop().to_integer(config::base_index(),10000)-config::base_index();
       bgeot::base_node N = normal_of_face(*pmesh, cv, f, node);
       out.pop().from_dcvector(N);
       );


    /*@GET N = ('normal of faces', @imat CVFIDs)
    Return matrix of (at face centers) the normal vectors of convexes.

    `CVFIDs` is supposed a two-rows matrix, the first row lists convex
    #ids, and the second lists face numbers (local number in the convex).@*/
    sub_command
      ("normal of faces", 1, 1, 0, 1,
       iarray v = in.pop().to_iarray(2,-1);
       darray w = out.pop().create_darray(pmesh->dim(), v.getn());
       for (size_type j=0; j < v.getn(); j++) {
         size_type cv = v(0,j) - config::base_index();
         short_type f  = short_type(v(1,j) - config::base_index());
         bgeot::base_node N = normal_of_face(*pmesh, cv, f, 0);
         for (size_type i=0; i < pmesh->dim(); ++i) w(i,j)=N[i];
       }
       );


    /*@GET CVIDs = ('convexes in box', @vec pmin, @vec pmax)
    Return the set of convexes lying entirely within the box defined by the corner points `pmin` and `pmax`.

    The output `CVIDs` is a two-rows matrix, the first row lists convex
    #ids, and the second one lists face numbers (local number in the
    convex). If `CVIDs` is given, it returns portion of the boundary of
    the convex set defined by the #ids listed in `CVIDs`.@*/
    sub_command
      ("convexes in box", 2, 2, 0, 1,
       check_empty_mesh(pmesh);
       int dim = pmesh->dim();
       darray p1 = in.pop().to_darray(dim);
       darray p2 = in.pop().to_darray(dim);
       bgeot::base_node pmin(dim);
       bgeot::base_node pmax(dim);
       for (int k=0; k < dim; ++k) {
         pmin[k] = std::min(p1[k],p2[k]);
         pmax[k] = std::max(p1[k],p2[k]);
       }
       getfem::mesh_region mr = select_convexes_in_box(*pmesh, pmin, pmax);
       iarray w = out.pop().create_iarray_h(unsigned(mr.size()));
       size_type j(0);
       for (getfem::mr_visitor i(mr); !i.finished(); ++i, ++j)
         w[j] = int(i.cv()+config::base_index());
       );

    /*@GET Q = ('quality'[, @ivec CVIDs])
    Return an estimation of the quality of each convex (:math:`0 \leq Q \leq 1`).@*/
    sub_command
      ("quality", 0, 1, 0, 1,
       dal::bit_vector bv;
       if (in.remaining()) bv = in.pop().to_bit_vector(&pmesh->convex_index());
       else bv = pmesh->convex_index();
       darray w = out.pop().create_darray_h(unsigned(bv.card()));
       size_type cnt = 0;
       for (dal::bv_visitor cv(bv); !cv.finished(); ++cv, ++cnt)
         w[cnt] = pmesh->convex_quality_estimate(cv);
       );


    /*@GET A = ('convex area'[, @ivec CVIDs])
    Return an estimate of the area of each convex.@*/
    sub_command
      ("convex area", 0, 1, 0, 1,
       dal::bit_vector bv;
       if (in.remaining()) bv = in.pop().to_bit_vector(&pmesh->convex_index());
       else bv = pmesh->convex_index();
       darray w = out.pop().create_darray_h(unsigned(bv.card()));
       size_type cnt = 0;
       for (dal::bv_visitor cv(bv); !cv.finished(); ++cv, ++cnt)
         w[cnt] = pmesh->convex_area_estimate(cv);
       );

    /*@GET A = ('convex radius'[, @ivec CVIDs])
    Return an estimate of the radius of each convex.@*/
    sub_command
      ("convex radius", 0, 1, 0, 1,
       dal::bit_vector bv;
       if (in.remaining()) bv = in.pop().to_bit_vector(&pmesh->convex_index());
       else bv = pmesh->convex_index();
       darray w = out.pop().create_darray_h(unsigned(bv.card()));
       size_type cnt = 0;
       for (dal::bv_visitor cv(bv); !cv.finished(); ++cv, ++cnt)
         w[cnt] = pmesh->convex_radius_estimate(cv);
       );


    /*@GET @CELL{S, CV2S} = ('cvstruct'[, @ivec CVIDs])
    Return an array of the convex structures.

    If `CVIDs` is not given, all convexes are considered. Each convex
    structure is listed once in `S`, and `CV2S` maps the convexes
    indice in `CVIDs` to the indice of its structure in `S`.@*/
    sub_command
      ("cvstruct", 0, 1, 0, 2,
       get_structure_or_geotrans_of_convexes(*pmesh, in, out,
                                             CVSTRUCT_CLASS_ID);
       );


    /*@GET @CELL{GT, CV2GT} = ('geotrans'[, @ivec CVIDs])
    Returns an array of the geometric transformations.

    See also MESH:GET('cvstruct').@*/
    sub_command
      ("geotrans", 0, 1, 0, 2,
       get_structure_or_geotrans_of_convexes(*pmesh, in, out,
                                             GEOTRANS_CLASS_ID);
       );


    /*@GET RIDs = ('boundaries')
    DEPRECATED FUNCTION. Use 'regions' instead. @*/
    sub_command
      ("boundaries", 0, 0, 0, 1,
       iarray w
       = out.pop().create_iarray_h(unsigned(pmesh->regions_index().card()));
       size_type i=0;
       for (dal::bv_visitor k(pmesh->regions_index()); !k.finished();
            ++k, ++i) {
         w[i] = int(k);
       }
       if (i != w.size()) THROW_INTERNAL_ERROR;
       );


    
    /*@GET RIDs = ('regions')
    Return the list of valid regions stored in the mesh.@*/
    sub_command
      ("regions", 0, 0, 0, 1,
       iarray w
       = out.pop().create_iarray_h(unsigned(pmesh->regions_index().card()));
       size_type i=0;
       for (dal::bv_visitor k(pmesh->regions_index());
            !k.finished(); ++k, ++i) {
         w[i] = int(k);
       }
       if (i != w.size()) THROW_INTERNAL_ERROR;
       );


    /*@GET RIDs = ('boundary')
    DEPRECATED FUNCTION. Use 'region' instead. @*/
    sub_command
      ("boundary", 1, 1, 0, 1,
       check_empty_mesh(pmesh);
       
       std::vector<unsigned> cvlst;
       std::vector<short> facelst;
       
       dal::bit_vector rlst = in.pop().to_bit_vector(0,0);
       
       for (dal::bv_visitor rnum(rlst); !rnum.finished(); ++rnum) {
         if (pmesh->regions_index().is_in(rnum)) {
           for (getfem::mr_visitor i(pmesh->region(rnum));
                !i.finished(); ++i) {
             cvlst.push_back(unsigned(i.cv()));
             facelst.push_back(i.f());
           }
         }
       }
       
       iarray w = out.pop().create_iarray(2, unsigned(cvlst.size()));
       for (size_type j=0; j < cvlst.size(); j++) {
         w(0,j) = cvlst[j]+config::base_index();
         w(1,j) = short_type(facelst[j]+config::base_index());
       }
       );

    

    /*@GET CVFIDs = ('region', @ivec RIDs)
    Return the list of convexes/faces on the regions `RIDs`.

    `CVFIDs` is a two-rows matrix, the first row lists convex #ids,
    and the second lists face numbers (local number in the convex).
    (and @MATLAB{0}@SCILAB{0}@PYTHON{-1} when the whole convex is in the
    regions).@*/
    sub_command
      ("region", 1, 1, 0, 1,
       check_empty_mesh(pmesh);
       
       std::vector<unsigned> cvlst;
       std::vector<short> facelst;
       
       dal::bit_vector rlst = in.pop().to_bit_vector(0,0);
       
       for (dal::bv_visitor rnum(rlst); !rnum.finished(); ++rnum) {
         if (pmesh->regions_index().is_in(rnum)) {
           for (getfem::mr_visitor i(pmesh->region(rnum));
                !i.finished(); ++i) {
             cvlst.push_back(unsigned(i.cv()));
             facelst.push_back(i.f());
           }
         }
       }

       iarray w = out.pop().create_iarray(2, unsigned(cvlst.size()));
       for (size_type j=0; j < cvlst.size(); j++) {
         w(0,j) = cvlst[j]+config::base_index();
         w(1,j) = short_type(facelst[j]+config::base_index());
       }
       );


    /*@GET ('save', @str filename)
    Save the mesh object to an ascii file.

    This mesh can be restored with MESH:INIT('load', filename).@*/
    sub_command
      ("save", 1, 1, 0, 0,
       std::string fname = in.pop().to_string();
       pmesh->write_to_file(fname);
       );


    /*@GET s = ('char')
    Output a string description of the mesh.@*/
    sub_command
      ("char", 0, 0, 0, 1,
       std::stringstream s;
       pmesh->write_to_file(s);
       out.pop().from_string(s.str().c_str());
       );


    /*@GET ('export to vtk', @str filename, ... [,'ascii'][,'quality'])
    Exports a mesh to a VTK file .

    If 'quality' is specified, an estimation of the quality of each
    convex will be written to the file.

    See also MESH_FEM:GET('export to vtk'), SLICE:GET('export to vtk').@*/
    sub_command
      ("export to vtk", 1, 3, 0, 1,
       bool write_q = false; bool ascii = false;
       std::string fname = in.pop().to_string();
       while (in.remaining() && in.front().is_string()) {
         std::string cmd2 = in.pop().to_string();
         if (cmd_strmatch(cmd2, "ascii"))
           ascii = true;
         else if (cmd_strmatch(cmd2, "quality"))
           write_q = true;
         else THROW_BADARG("expecting 'ascii' or 'quality', got " << cmd2);
       }
       getfem::vtk_export exp(fname, ascii);
       exp.exporting(*pmesh);
       exp.write_mesh(); if (write_q) exp.write_mesh_quality(*pmesh);
       );


    /*@GET ('export to dx', @str filename, ... [,'ascii'][,'append'][,'as',@str name,[,'serie',@str serie_name]][,'edges'])
    Exports a mesh to an OpenDX file.

    See also MESH_FEM:GET('export to dx'), SLICE:GET('export to dx').@*/
    sub_command
      ("export to dx", 1, 3, 0, 1,
       bool ascii = false; bool append = false; bool edges = false;
       std::string fname = in.pop().to_string();
       std::string mesh_name;
       std::string serie_name;
       while (in.remaining() && in.front().is_string()) {
         std::string cmd2 = in.pop().to_string();
         if (cmd_strmatch(cmd2, "ascii"))
           ascii = true;
         else if (cmd_strmatch(cmd2, "edges"))
           edges = true;
         else if (cmd_strmatch(cmd2, "append"))
           append = true;
         else if (cmd_strmatch(cmd2, "as") && in.remaining())
           mesh_name = in.pop().to_string();
         else if (cmd_strmatch(cmd2, "serie") && in.remaining() && mesh_name.size())
           serie_name = in.pop().to_string();
         else THROW_BADARG("expecting 'ascii' or 'append', 'serie', or 'as' got " << cmd2);
       }
       getfem::dx_export exp(fname, ascii, append);
       exp.exporting(*pmesh, mesh_name);
       exp.write_mesh();
       if (edges) exp.exporting_mesh_edges();
       if (serie_name.size()) exp.serie_add_object(serie_name,mesh_name);
       );


    /*@GET ('export to pos', @str filename[, @str name])
    Exports a mesh to a POS file .

    See also MESH_FEM:GET('export to pos'), SLICE:GET('export to pos').@*/
    sub_command
      ("export to pos", 1, 2, 0, 0,
       std::string fname = in.pop().to_string();
       std::string name = "";
       if (in.remaining()) name = in.pop().to_string();
       
       getfem::pos_export exp(fname);
       exp.write(*pmesh,name);
       );


    /*@GET z = ('memsize')
      Return the amount of memory (in bytes) used by the mesh.@*/
    sub_command
      ("memsize", 0, 0, 0, 1,
       out.pop().from_integer(int(pmesh->memsize()));
       );


    /*@GET ('display')
      displays a short summary for a @tmesh object.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       infomsg() << "gfMesh object in dimension " << int(pmesh->dim())
       << " with " << pmesh->nb_points() << " points and "
       << pmesh->convex_index().card() << " elements\n";
       );

  }

  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");
  const getfem::mesh *pmesh = extract_mesh_object(m_in.pop());
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
