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
#include <getfemint_workspace.h>
#include <getfem/getfem_mesh_slice.h>
#include <getfem/getfem_export.h>

using namespace getfemint;

static void fmt_pt_povray(std::ofstream &f, const getfem::base_node &pt) {
  char s[100];
  if (pt.size() == 0) GMM_THROW(getfemint_error, "empty point");
  sprintf(s, "<%g,%g,%g>",pt[0],pt.size() > 1 ? pt[1] : 0., pt.size() > 2 ? pt[2] : 0.);
  f << s;
}

static void fmt_pt_povray(std::ofstream &f, const getfem::base_node &pt, const getfem::base_node &n_) {
  getfem::base_node n = 1./gmm::vect_norm2(n_) * n_;
  fmt_pt_povray(f,pt); f << ","; fmt_pt_povray(f,n);
}

static void
export_slice_to_povray(std::ofstream &f, const getfem::stored_mesh_slice& sl) {
  //time_t t = time(NULL);
  //f << "# exported for getfem " << ctime(&t) << "\n";
  f << "mesh {\n";
  size_type igncnt = 0;
  const getfem::mesh &m = sl.linked_mesh();
  for (size_type ic=0; ic < sl.nb_convex(); ++ic) {
    for (getfem::mesh_slicer::cs_simplexes_ct::const_iterator it=sl.simplexes(ic).begin();
         it != sl.simplexes(ic).end(); ++it) {
      if (it->dim() == 2) {
        const getfem::slice_node &A = sl.nodes(ic)[it->inodes[0]];
        const getfem::slice_node &B = sl.nodes(ic)[it->inodes[1]];
        const getfem::slice_node &C = sl.nodes(ic)[it->inodes[2]];
        getfem::slice_node::faces_ct common_faces = (A.faces & B.faces & C.faces);
        short_type fnum = 0;
        for (; common_faces.any(); ++fnum) if (common_faces[fnum]) break;
        if (fnum < m.structure_of_convex(sl.convex_num(ic))->nb_faces()) {
          f << "smooth_triangle {";
          fmt_pt_povray(f,A.pt,m.normal_of_face_of_convex(sl.convex_num(ic),fnum,A.pt_ref));
          fmt_pt_povray(f,B.pt,m.normal_of_face_of_convex(sl.convex_num(ic),fnum,B.pt_ref));
          fmt_pt_povray(f,C.pt,m.normal_of_face_of_convex(sl.convex_num(ic),fnum,C.pt_ref));
          f << "}\n";
        } else {
          f << "triangle {";
          fmt_pt_povray(f,A.pt);
          fmt_pt_povray(f,B.pt);
          fmt_pt_povray(f,C.pt);
          f << "}\n";
        }
      } else ++igncnt;
    }
  }
  f << "}\n";
  if (igncnt) cout << igncnt << " simplexes of dim != 2 ignored\n";
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
interpolate_convex_data(const getfem::stored_mesh_slice *sl,
                        const garray<T> &u, getfemint::mexargs_out& out) {
  assert(u.dim(u.ndim()-1) == sl->linked_mesh().convex_index().last_true()+1);
  array_dimensions ad;
  for (unsigned i=0; i < u.ndim()-1; ++i) ad.push_back(u.dim(i));
  ad.push_back(unsigned(sl->nb_points()));
  garray<T> w = out.pop().create_array(ad, T());
  size_type q = u.size() / u.dim(u.ndim()-1);
  size_type pos = 0;
  for (size_type i=0; i < sl->nb_convex(); ++i) {
    for (unsigned j=0; j < q; ++j) {
      T v = u[(sl->convex_num(i))*q + j];
      for (unsigned k=0; k < sl->nodes(i).size(); ++k) {
        w[pos++] = v;
      }
    }
  }
  assert(pos == w.size());
}

/*@GFDOC
  General function for querying information about @tslice objects.
@*/






// Object for the declaration of a new sub-command.

struct sub_gf_slice_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   getfem::stored_mesh_slice *sl) = 0;
};

typedef std::shared_ptr<sub_gf_slice_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_slice_get {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       getfem::stored_mesh_slice *sl)			\
      { dummy_func(in); dummy_func(out); dummy_func(sl); code }		\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }




void gf_slice_get(getfemint::mexargs_in& m_in,
		  getfemint::mexargs_out& m_out) {

  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {

    /*@RDATTR d = ('dim')
      Return the dimension of the slice (2 for a 2D mesh, etc..).@*/
    sub_command
      ("dim", 0, 0, 0, 1,
       out.pop().from_integer(int(sl->dim()));
       );


    /*@GET a = ('area')
      Return the area of the slice.@*/
    sub_command
      ("area", 0, 0, 0, 1,
       getfem::slicer_compute_area s; sl->replay(s);
       out.pop().from_scalar(s.area());
       );


    /*@GET CVids = ('cvs')
    Return the list of convexes of the original mesh contained in the slice.@*/
    sub_command
      ("cvs", 0, 0, 0, 1,
       iarray w = out.pop().create_iarray_h(unsigned(sl->nb_convex()));
       for (size_type i=0; i < sl->nb_convex(); ++i)
	 w[i] = int(sl->convex_num(i) + config::base_index());
       );


    /*@RDATTR n = ('nbpts')
      Return the number of points in the slice.@*/
    sub_command
      ("nbpts", 0, 0, 0, 1,
       out.pop().from_integer(int(sl->nb_points()));
       );


    /*@RDATTR ns = ('nbsplxs'[, @int dim])
    Return the number of simplexes in the slice.

    Since the slice may contain points (simplexes of dim 0), segments
    (simplexes of dimension 1), triangles etc., the result is a vector
    of size SLICE:GET('dim')+1, except if the optional argument `dim`
    is used.@*/
    sub_command
      ("nbsplxs", 0, 1, 0, 1,
       std::vector<size_type> v; sl->nb_simplexes(v);
       if (in.remaining()) {
	 size_type i= in.pop().to_integer(0,100);
	 out.pop().from_integer(int(i < v.size() ? v[i] : 0));
       } else {
	 out.pop().from_ivector(v);
       }
       );


    /*@GET P = ('pts')
      Return the list of point coordinates.@*/
    sub_command
      ("pts", 0, 0, 0, 1,
       darray w = out.pop().create_darray(unsigned(sl->dim()), unsigned(sl->nb_points()));
       for (size_type ic=0, cnt=0; ic < sl->nb_convex(); ++ic) {
	 for (getfem::mesh_slicer::cs_nodes_ct::const_iterator it=sl->nodes(ic).begin();
	      it != sl->nodes(ic).end(); ++it) {
	   for (size_type k=0; k < sl->dim(); ++k)
	     w[cnt++] = it->pt[k];
	 }
       }
       );


    /*@GET @CELL{S, CV2S} = ('splxs',@int dim)
      Return the list of simplexes of dimension `dim`.
      
      On output, S has 'dim+1' rows, each column contains the point
      numbers of a simplex.  The vector `CV2S` can be used to find the
      list of simplexes for any convex stored in the slice. For example
      '@MATLAB{S(:,CV2S(4):CV2S(5)-1)}@SCILAB{S(:,CV2S(4):CV2S(5)-1)}@PYTHON{S[:,CV2S[4]:CV2S[5]]}'
      gives the list of simplexes for the fourth convex.@*/
    sub_command
      ("splxs", 1, 1, 0, 2,
       size_type sdim = in.pop().to_integer(0,int(sl->dim()));
       iarray w = out.pop().create_iarray(unsigned(sdim+1), unsigned(sl->nb_simplexes(sdim)));
       size_type Scnt = size_type(-1);
       iarray cv2splx;
       if (out.remaining()) {
	 cv2splx = out.pop().create_iarray_h(unsigned(sl->nb_convex()+1));
	 Scnt = config::base_index();
       }
       for (size_type ic=0, cnt=0, pcnt=0; ic < sl->nb_convex(); ++ic) {
	 size_type scnt = 0;
	 for (getfem::mesh_slicer::cs_simplexes_ct::const_iterator it=sl->simplexes(ic).begin();
	      it != sl->simplexes(ic).end(); ++it) {
	   if (it->dim() == sdim) {
	     for (size_type k=0; k < sdim+1; ++k)
	       w[cnt++] = int(it->inodes[k] + pcnt + config::base_index());
	     scnt++;
	   }
	 }
	 pcnt += sl->nodes(ic).size();
	 if (Scnt != size_type(-1)) { cv2splx[ic] = int(Scnt); Scnt+=scnt; }
       }
       if (Scnt != size_type(-1)) cv2splx[sl->nb_convex()] = int(Scnt);
       );


    /*@GET @CELL{P, E1, E2} = ('edges')
      Return the edges of the linked mesh contained in the slice.
      
      `P` contains the list of all edge vertices, `E1` contains
      the indices of each mesh edge in `P`, and `E2` contains the
      indices of each "edges" which is on the border of the slice.
      This function is useless except for post-processing purposes.@*/
    sub_command
      ("edges", 0, 0, 3, 3,
       getfem::mesh m;
       dal::bit_vector slice_edges;
       getfem::mesh_slicer slicer(sl->linked_mesh());
       getfem::slicer_build_edges_mesh action(m,slice_edges);
       slicer.push_back_action(action); slicer.exec(*sl);

       /* return a point list, a connectivity array, and optionnaly a list of edges with are part of the slice */
       double nan = ::nan("");
       dal::bit_vector bv = m.points().index();
       darray P = out.pop().create_darray(m.dim(), unsigned(bv.last_true()+1));
       iarray T1 = out.pop().create_iarray(2, unsigned(m.nb_convex() - slice_edges.card()));
       iarray T2 = out.pop().create_iarray(2, unsigned(slice_edges.card()));
       for (size_type j = 0; j < bv.last_true()+1; j++) {
	 for (size_type i = 0; i < m.dim(); i++) {
	   P(i,j) = (bv.is_in(j)) ? (m.points()[j])[i] : nan;
	 }
       }
       iarray::iterator itt1=T1.begin(); iarray::iterator itt2=T2.begin();
       for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv) {
	 if (!slice_edges[cv]) {
	   itt1[0] = int(m.ind_points_of_convex(cv)[0]);
	   itt1[1] = int(m.ind_points_of_convex(cv)[1]);
	   // gmm::copy_n(m.ind_points_of_convex(cv).begin(), 2, itt1);
	   itt1[0] += config::base_index();
	   itt1[1] += config::base_index(); itt1 += 2;
	 } else {
	   itt2[0] = int(m.ind_points_of_convex(cv)[0]);
	   itt2[1] = int(m.ind_points_of_convex(cv)[1]);
	   // gmm::copy_n(m.ind_points_of_convex(cv).begin(), 2, itt2);
	   itt2[0] += config::base_index();
	   itt2[1] += config::base_index(); itt2 += 2;
	 }
       }
       );


    /*@GET Usl = ('interpolate_convex_data', @mat Ucv)
    Interpolate data given on each convex of the mesh to the slice nodes.

    The input array `Ucv` may have any number of dimensions, but its
    last dimension should be equal to MESH:GET('max cvid').

    Example of use: SLICE:GET('interpolate_convex_data', MESH:GET('quality')).@*/
    sub_command
      ("interpolate_convex_data", 1, 1, 0, 1,
       in.front().check_trailing_dimension(int(sl->linked_mesh().convex_index().last_true()+1));
       if (in.front().is_complex())
	 interpolate_convex_data(sl, in.pop().to_darray(), out);
       else interpolate_convex_data(sl, in.pop().to_carray(), out);
       );

    /*@GET m = ('linked mesh')
      Return the mesh on which the slice was taken.@*/
    sub_command
      ("linked mesh", 0, 0, 0, 1,
       id_type id = workspace().object((const void *)(&(sl->linked_mesh())));
       if (id == id_type(-1)) {
	 auto pst = workspace().hidden_object(workspace().object(sl),
					      &sl->linked_mesh());
	 if (!pst.get()) THROW_INTERNAL_ERROR;
	 std::shared_ptr<getfem::mesh> pm = 
	   std::const_pointer_cast<getfem::mesh>
	   (std::dynamic_pointer_cast<const getfem::mesh>(pst));
	 id = store_mesh_object(pm);
       }
       out.pop().from_object_id(id, MESH_CLASS_ID);
       );

    /*@GET m = ('mesh')
      Return the mesh on which the slice was taken
      (identical to 'linked mesh')@*/
    sub_command
      ("mesh", 0, 0, 0, 1,
       id_type id = workspace().object((const void *)(&(sl->linked_mesh())));
       if (id == id_type(-1)) THROW_INTERNAL_ERROR;
       out.pop().from_object_id(id, MESH_CLASS_ID);
       );


    /*@GET z = ('memsize')
    Return the amount of memory (in bytes) used by the slice object.@*/
    sub_command
      ("memsize", 0, 0, 0, 1,
       out.pop().from_integer(int(sl->memsize()));
       );


    /*@GET ('export to vtk', @str filename, ...)
    Export a slice to VTK.

    Following the `filename`, you may use any of the following options:

    - if 'ascii' is not used, the file will contain binary data
      (non portable, but fast).
    - if 'edges' is used, the edges of the original mesh will be
      written instead of the slice content.

    More than one dataset may be written, just list them. Each dataset
    consists of either:

    - a field interpolated on the slice (scalar, vector or tensor),
      followed by an optional name.
    - a mesh_fem and a field, followed by an optional name.

    Examples:

    - SLICE:GET('export to vtk', 'test.vtk', Usl, 'first_dataset', mf,
      U2, 'second_dataset')
    - SLICE:GET('export to vtk', 'test.vtk', 'ascii', mf,U2)
    - SLICE:GET('export to vtk', 'test.vtk', 'edges', 'ascii', Uslice)@*/
    sub_command
      ("export to vtk",1, -1, 0, 0,
       std::string fname = in.pop().to_string();
       bool ascii = false;
       bool edges = false;
       while (in.remaining() && in.front().is_string()) {
	 std::string cmd2 = in.pop().to_string();
	 if (cmd_strmatch(cmd2, "ascii"))
	   ascii = true;
	 else if (cmd_strmatch(cmd2, "edges"))
	   edges = true;
	 else THROW_BADARG("expecting 'ascii' or 'edges', got " << cmd2);
       }
       getfem::vtk_export exp(fname, ascii);
       getfem::stored_mesh_slice sl_edges;
       const getfem::stored_mesh_slice *vtk_slice = sl;
       getfem::mesh m_edges;
       if (edges) {
	 vtk_slice = &sl_edges;
	 dal::bit_vector slice_edges;
	 getfem::mesh_slicer slicer(sl->linked_mesh());
	 getfem::slicer_build_edges_mesh action(m_edges,slice_edges);
	 slicer.push_back_action(action); slicer.exec(*sl);
	 sl_edges.build(m_edges, getfem::slicer_none());
       }
       exp.exporting(*vtk_slice);
       exp.write_mesh();
       int count = 1;
       if (in.remaining()) {
	 do {
	   if (in.remaining() >= 2 && is_meshfem_object(in.front())) {
	     const getfem::mesh_fem &mf = *to_meshfem_object(in.pop());
	     
	     darray U = in.pop().to_darray();
	     in.last_popped().check_trailing_dimension(int(mf.nb_dof()));
	     
	     exp.write_point_data(mf,U,get_vtk_dataset_name(in, count));
	   } else if (in.remaining()) {
	     darray slU = in.pop().to_darray();
	     in.last_popped().check_trailing_dimension(int(vtk_slice->nb_points()));
	     
	     exp.write_sliced_point_data(slU,get_vtk_dataset_name(in, count));
	   } else THROW_BADARG("don't know what to do with this argument")
	       count+=1;
	 } while (in.remaining());
       }
       );


    /*@GET ('export to vtu', @str filename, ...)
    Export a slice to VTK(XML).

    Following the `filename`, you may use any of the following options:

    - if 'ascii' is not used, the file will contain binary data
      (non portable, but fast).
    - if 'edges' is used, the edges of the original mesh will be
      written instead of the slice content.

    More than one dataset may be written, just list them. Each dataset
    consists of either:

    - a field interpolated on the slice (scalar, vector or tensor),
      followed by an optional name.
    - a mesh_fem and a field, followed by an optional name.

    Examples:

    - SLICE:GET('export to vtu', 'test.vtu', Usl, 'first_dataset', mf,
      U2, 'second_dataset')
    - SLICE:GET('export to vtu', 'test.vtu', 'ascii', mf,U2)
    - SLICE:GET('export to vtu', 'test.vtu', 'edges', 'ascii', Uslice)@*/
    sub_command
      ("export to vtu",1, -1, 0, 0,
       std::string fname = in.pop().to_string();
       bool ascii = false;
       bool edges = false;
       while (in.remaining() && in.front().is_string()) {
         std::string cmd2 = in.pop().to_string();
         if (cmd_strmatch(cmd2, "ascii"))
           ascii = true;
         else if (cmd_strmatch(cmd2, "edges"))
           edges = true;
         else THROW_BADARG("expecting 'ascii' or 'edges', got " << cmd2);
       }
       getfem::vtu_export exp(fname, ascii);
       getfem::stored_mesh_slice sl_edges;
       const getfem::stored_mesh_slice *vtu_slice = sl;
       getfem::mesh m_edges;
       if (edges) {
         vtu_slice = &sl_edges;
         dal::bit_vector slice_edges;
         getfem::mesh_slicer slicer(sl->linked_mesh());
         getfem::slicer_build_edges_mesh action(m_edges,slice_edges);
         slicer.push_back_action(action); slicer.exec(*sl);
         sl_edges.build(m_edges, getfem::slicer_none());
       }
       exp.exporting(*vtu_slice);
       exp.write_mesh();
       int count = 1;
       if (in.remaining()) {
         do {
           if (in.remaining() >= 2 && is_meshfem_object(in.front())) {
             const getfem::mesh_fem &mf = *to_meshfem_object(in.pop());

             darray U = in.pop().to_darray();
             in.last_popped().check_trailing_dimension(int(mf.nb_dof()));

             exp.write_point_data(mf,U,get_vtk_dataset_name(in, count));
           } else if (in.remaining()) {
             darray slU = in.pop().to_darray();
             in.last_popped().check_trailing_dimension(int(vtu_slice->nb_points()));

             exp.write_sliced_point_data(slU,get_vtk_dataset_name(in, count));
           } else THROW_BADARG("don't know what to do with this argument")
               count+=1;
         } while (in.remaining());
       }
       );


    /*@GET ('export to pov', @str filename)
      Export a the triangles of the slice to POV-RAY.@*/
    sub_command
      ("export to pov",1, 1, 0, 0,
       std::string fname = in.pop().to_string();
       std::ofstream f(fname.c_str());
       export_slice_to_povray(f, *sl);
       );


    /*@GET ('export to dx', @str filename, ...)
    Export a slice to OpenDX.

    Following the `filename`, you may use any of the following
    options:

    - if 'ascii' is not used, the file will contain binary data
      (non portable, but fast).
    - if 'edges' is used, the edges of the original mesh will be
      written instead of the slice content.
    - if 'append' is used, the opendx file will not be overwritten,
      and the new data will be added at the end of the file.

    More than one dataset may be written, just list them. Each dataset
    consists of either:

    - a field interpolated on the slice (scalar, vector or tensor),
      followed by an optional name.
    - a mesh_fem and a field, followed by an optional name.@*/
    sub_command
      ("export to dx", 1, -1, 0, 0,
       std::string fname = in.pop().to_string();
       bool ascii = false;
       bool edges = false;
       bool append = false;
       std::string mesh_name; std::string serie_name;
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
	 else if (cmd_strmatch(cmd2, "serie") && in.remaining())
	   serie_name = in.pop().to_string();
	 else THROW_BADARG("expecting 'ascii' or 'edges' or 'append' or 'as', got " << cmd2);
       }
       getfem::dx_export exp(fname, ascii, append);
       
       exp.exporting(*sl, mesh_name.c_str());
       exp.write_mesh();
       if (edges) exp.exporting_mesh_edges();
       if (in.remaining()) {
	 do {
	   if (in.remaining() >= 2 && is_meshfem_object(in.front())) {
	     const getfem::mesh_fem &mf = *to_meshfem_object(in.pop());
	     
	     darray U = in.pop().to_darray();
	     in.last_popped().check_trailing_dimension(int(mf.nb_dof()));
	     
	     exp.write_point_data(mf,U,get_dx_dataset_name(in));
	   } else if (in.remaining()) {
	     darray slU = in.pop().to_darray();
	     in.last_popped().check_trailing_dimension(int(sl->nb_points()));
	     
	     exp.write_sliced_point_data(slU,get_dx_dataset_name(in));
	   } else THROW_BADARG("don't know what to do with this argument");
	   if (serie_name.size()) exp.serie_add_object(serie_name);
	 } while (in.remaining());
       }
       );


    /*@GET ('export to pos', @str filename[, @str name][[,@tmf mf1], @mat U1, @str nameU1[[,@tmf mf1], @mat U2, @str nameU2,...])
    Export a slice to Gmsh.

    More than one dataset may be written, just list them.
    Each dataset consists of either:

    - a field interpolated on the slice (scalar, vector or tensor).
    - a mesh_fem and a field.@*/
    sub_command
      ("export to pos", 1, -1, 0, 0,
       std::string fname = in.pop().to_string();
       getfem::pos_export exp(fname);
       
       std::string name = "";
       if (in.remaining() && in.front().is_string())
	 name = in.pop().to_string();
       exp.write(*sl,name);
       while (in.remaining()) {
	 if (in.remaining() >= 3 && is_meshfem_object(in.front())) {
	     const getfem::mesh_fem *mf = to_meshfem_object(in.pop());
	   
	   darray U = in.pop().to_darray();
	   in.last_popped().check_trailing_dimension(int(mf->nb_dof()));
	   
	   if (in.remaining() >= 1 && in.front().is_string())
	     name = in.pop().to_string();
	   else THROW_BADARG("expecting string darray_name")
	     
	     exp.write(*mf, U, name);
	 } else if (in.remaining() >=2) {
	   darray slU = in.pop().to_darray();
	   in.last_popped().check_trailing_dimension(int(sl->nb_points()));
	   
	   if (in.remaining() >= 1 && in.front().is_string())
	     name = in.pop().to_string();
	   else THROW_BADARG("expecting string darray_name")
	     
	     exp.write(*sl, slU, name);
	 }
       }
       );


    /*@GET s = ('char')
      Output a (unique) string representation of the @tslice.

      This can be used to perform comparisons between two
      different @tslice objects.
      This function is to be completed.
      @*/
    sub_command
      ("char", 0, 0, 0, 1,
       GMM_ASSERT1(false, "Sorry, function to be done");
       // std::string s = ...;
       // out.pop().from_string(s.c_str());
       );


    /*@GET ('display')
      displays a short summary for a @tslice object.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       infomsg() << "gfSlice object in dimension " << sl->dim()
       << " and " << sl->nb_points() << " points.\n";
       );

  }
 
  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  getfem::stored_mesh_slice *sl = to_slice_object(m_in.pop());
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out, sl);
  }
  else bad_cmd(init_cmd);

}
