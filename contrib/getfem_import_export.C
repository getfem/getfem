#include <list>
#include <getfem_mesh_fem.h>
#include <getfem_import.h>
#include <getfem_export.h>

using std::cerr;
using std::cout;
using std::cin;
using std::endl;
using getfem::scalar_type;
using getfem::size_type;
using getfem::short_type;
using getfem::dim_type;

namespace getfem {
  size_type simplex_npoints(dim_type dim, int nrefine=1) {
    size_type n=1;
    for (size_type i=1; i <= dim; ++i) { n *= (nrefine+i); n /= i; }
    return n;
  }


  /* sort convexes by their basic_structures/structures/convexes_of_reference */
  class _compare_convex_structures {
    const getfem_mesh& m;
  public:
    _compare_convex_structures(const getfem_mesh& _m) : m(_m) {}
    bool operator()(size_type a, size_type b) {
      bgeot::pconvex_ref pa = m.trans_of_convex(a)->convex_ref();
      bgeot::pconvex_ref pb = m.trans_of_convex(b)->convex_ref();
      if (pa->structure()->basic_structure() < pb->structure()->basic_structure())
	return true;
      else if (pa->structure()->basic_structure() > pb->structure()->basic_structure())
	return false;
      else if (pa->structure() < pb->structure())
	return true;
      else if (pa->structure() > pb->structure())
	return false;
      else return (pa < pb ? true : (a < b ? true : false));
    }
  };

  void ordered_convex_structures(const getfem_mesh& m, const dal::bit_vector& cvlist, 
				 std::vector<size_type>& ordered_cvlist) {
    ordered_cvlist.resize(0); ordered_cvlist.reserve(cvlist.card());
    dal::bit_vector::const_iterator it = cvlist.begin();
    for (; it != cvlist.end(); ++it) {
      if (*it) ordered_cvlist.push_back(it.index());
    }
    std::sort(ordered_cvlist.begin(), ordered_cvlist.end(), _compare_convex_structures(m));
  }

  template <typename IT>
  void transform_edge_list(const getfem::getfem_mesh &m, size_type nrefine, const bgeot::edge_list &el, IT out_it)
  {
    bgeot::edge_list::const_iterator it;
    size_type cv = size_type(-1);
    getfem::getfem_mesh::ref_convex cv_ref;
    bgeot::pgeometric_trans pgt = NULL;
    for (it = el.begin(); it != el.end(); it++) {
      if (cv != (*it).cv) {
	cv = (*it).cv;
	cv_ref = m.convex(cv);
	pgt = m.trans_of_convex(cv);
      }
      /* build the list of points on the edge, on the reference element */
      /* get local point numbers in the convex */
      bgeot::size_type iA = m.local_ind_of_convex_point(cv, (*it).i);
      bgeot::size_type iB = m.local_ind_of_convex_point(cv, (*it).j);    
      getfem::base_node A = pgt->convex_ref()->points()[iA];
      getfem::base_node B = pgt->convex_ref()->points()[iB];
      for (size_type i = 0; i < nrefine+1; i++) {
	getfem::base_node pt;
	
	pt = A +  (B-A)*(i/(double)(nrefine));
	pt = pgt->transform(pt, m.points_of_convex(cv));
	out_it = std::copy(pt.begin(), pt.end(), out_it);
      }
    }
  }



  /**
     This object stores internally a list of pointers to the mesh
     objects which have been written in the opendx file: DO NOT
     modify/destroy the exported mesh objects during the exportation
   */
  class opendx_export {
    struct _exported_mesh {
      size_type  nrefine;
      bool continuous_data;
      size_type  pts_object_num, cvs_object_num;
      size_type  nb_pts, nb_simplexes;
      size_type  sdim; 
      
      const getfem::getfem_mesh *pm;
      std::vector<size_type> cvlist;
      std::vector<const getfem_mesh *> refined_convexes_mesh;
      std::vector<size_type> pts_id; // list of node numbers (used only when the field continuity is enforced)
      _exported_mesh(const getfem::getfem_mesh *_pm, size_type _nrefine, bool _continuous_data) :
	nrefine(_nrefine), continuous_data(_continuous_data), pts_object_num(0), cvs_object_num(0), 
	nb_pts(0), nb_simplexes(0), sdim(0), pm(_pm)
      {}
    };
    std::vector<_exported_mesh*> exported_meshes;

    const getfem::getfem_mesh *current_mesh;
    const getfem::mesh_fem *_mf;

    std::ofstream of;
    std::ostream& o;
    size_type object_cnt;
    bool closed, binary;
  public:
    static const char* indianness() {
      static int i=0x12345678;
      char *p = (char*)&i;
      if (*p == 0x12) return "msb";
      else if (*p == 0x78) return "lsb";
      else return "this is very strange..";
    }
    /**
       open a .dx file and prepares it
    */
    opendx_export(const std::string& fname, bool bin=true) : of(fname.c_str()), o(of), object_cnt(1), closed(false),binary(bin) {
      if (!o) DAL_THROW(file_not_found_error,
			"dx file '" << fname << "' can't be opened");
      write_header();
    }
    opendx_export(std::ostream& os) : o(os), object_cnt(1), closed(false) {
      write_header();
    }    
    ~opendx_export() {
      for (size_type i=0; i < exported_meshes.size(); ++i) delete exported_meshes[i];
      close();
    }
    void close() {
      if (!closed) o << "end" << endl; closed = true;
    }
    /**
       write the edges of the mesh into the opendx file (for mesh visualization)
    */
    void export_mesh_edges(const getfem::getfem_mesh& m, size_type nrefine);
    /**
       export a field into the opendx file 
    */
    template<class VECT> void export_solution(const getfem::mesh_fem& mf, const VECT& U, size_type nrefine, const char *field_name, bool continuous_data);
  private:
    void write_header() {
      o << "# data file for IBM OpenDX, generated by GetFem++ v " << GETFEM_VERSION << endl;
    }
    _exported_mesh* export_mesh_data_no_repeat(const getfem::getfem_mesh *pm, size_type nrefine, bool continuous_data);
    _exported_mesh* export_mesh_data(const getfem::getfem_mesh *pm, size_type nrefine, bool continuous_data);
    void bwritefloat(float f) {
      o.write((char*)&f, sizeof(f));
    }
    void bwriteint(int i) {
      o.write((char*)&i, sizeof(i));
    }
  };

  template <typename CONT> scalar_type volume_of_simplex(CONT pts) {
    if (pts.size() == 0) return 0;
    size_type N = pts.size()-1;
    base_matrix m(N,N);
    typename CONT::const_iterator it0 = pts.begin(), it = pts.begin()+1;
    for (size_type j=0; j < N; ++j, ++it)
      for (size_type i=0; i < N; ++i)
	m(i,j) = (*it)[i] - (*it0)[i];
    return gmm::lu_det(m);
  }

  /**
     write the point list and convex list of a mesh into the opendx data file,
     but only if it has not been already done. 
     @output the filled _exported_mesh structure
  */
  opendx_export::_exported_mesh* opendx_export::export_mesh_data_no_repeat(const getfem::getfem_mesh *pm, 
									   size_type nrefine, bool continuous_data) {
    for (size_type i=0; i < exported_meshes.size(); ++i) 
      if (exported_meshes[i]->pm == pm && exported_meshes[i]->nrefine == nrefine &&
	  exported_meshes[i]->continuous_data == continuous_data) return exported_meshes[i];
    _exported_mesh *em = export_mesh_data(pm,nrefine,continuous_data);
    exported_meshes.push_back(em);
    return em;
  }

  opendx_export::_exported_mesh* opendx_export::export_mesh_data(const getfem::getfem_mesh *pm, 
								 size_type nrefine, bool continuous_data) {
    typedef dal::dynamic_tree_sorted<base_node,dal::lexicographical_less<base_node,dal::approx_less<scalar_type> > > PT_TAB;
    PT_TAB pt_tab;
    _exported_mesh *em = new _exported_mesh(pm,nrefine,continuous_data);
    em->pts_object_num = object_cnt++;
    em->cvs_object_num = object_cnt++;

    const getfem::getfem_mesh& m = *pm;
    ordered_convex_structures(m, m.convex_index(), em->cvlist);
    em->refined_convexes_mesh.reserve(em->cvlist.size());


    /* check the consistency of convex dimensions */
    em->sdim = dim_type(-1); /* common dimension of convexes */
    for (size_type i=0; i < em->cvlist.size(); ++i) {
      size_type cv = em->cvlist[i];
      bgeot::pconvex_structure cvs = m.structure_of_convex(cv);      
      if ((cvs->dim()) != em->sdim) {
	if (em->sdim == dim_type(-1)) em->sdim = cvs->dim();
	else DAL_THROW(failure_error, "opendx won't handle a set of convexes of different dimensions\n");
      }
    }

    /* count total number of points and simplices */
    {
      bgeot::pconvex_ref prev_cvr = 0;
      const getfem_mesh *ref_convex_mesh = 0;
      size_type big_nb_pts = 0; // count all points, even duplicate ones
      for (size_type i=0; i < em->cvlist.size(); ++i) {
	bgeot::pconvex_ref cvr = m.trans_of_convex(em->cvlist[i])->convex_ref();
	if (prev_cvr != cvr) {
	  prev_cvr = cvr;
	  ref_convex_mesh = getfem::refined_simplex_mesh_for_convex(cvr, nrefine);
	}
	em->refined_convexes_mesh[i] = ref_convex_mesh;
	big_nb_pts += em->refined_convexes_mesh[i]->nb_points(); 
	em->nb_simplexes += em->refined_convexes_mesh[i]->nb_convex();
      }
      em->pts_id.reserve(big_nb_pts);
    }

    

    /* 
       build the point list .
       if continuous_data == true, then eliminate duplicate points (which occur on faces
       of adjacent convexes)
    */
    {
      _geotrans_precomp gp;
      bgeot::stored_point_tab pts, cvm_pts;
      for (size_type i=0; i < em->cvlist.size(); ++i) {
	size_type cv = em->cvlist[i];
	const getfem_mesh *cvm = em->refined_convexes_mesh[i];
	if (i == 0 || cvm != em->refined_convexes_mesh[i-1]) {
	  cvm_pts.resize(cvm->nb_points());
	  std::copy(cvm->points().begin(), cvm->points().end(), cvm_pts.begin());
	  geotrans_precomp_not_stored(m.trans_of_convex(cv), 
				      &cvm_pts, gp);
	  pts.resize(cvm->nb_points());
	}
	gp.transform(m.points_of_convex(cv), pts);
	for (size_type j=0; j < pts.size(); ++j) {
	  if (!continuous_data)
	    em->pts_id.push_back(pt_tab.add(pts[j]));
	  else {
	    em->pts_id.push_back(pt_tab.add_norepeat(pts[j]));
	  }
	}
      }
    }

    em->nb_pts = pt_tab.size();
    /* now write the points */
    {
      o << "object " << em->pts_object_num << " class array type float rank 1 shape " 
	<< int(m.dim()) 
	<< " items " << em->nb_pts;
      if (binary) o << " " << indianness() << " binary";
      o << " data follows" << endl;
      for (size_type j=0; j < pt_tab.size(); ++j) {
	for (size_type k=0; k < pt_tab[j].size(); ++k) {
	  if (!binary) 
	    o << pt_tab[j][k] << " ";
	  else bwritefloat(pt_tab[j][k]);
	}
	if (!binary) o << endl;
      }
      o << endl;
    }

    /* write convex list */
    {
      o << "object " << em->cvs_object_num << " class array type int rank 1 shape " 
	<< simplex_npoints(em->sdim) 
	<< " items " << em->nb_simplexes;
      if (binary) o << " " << indianness() << " binary";
      o<< " data follows" << endl;
      
      //bgeot::pconvex_ref prev_cvr = NULL;
      size_type count_pts = 0;
      for (size_type i=0; i < em->cvlist.size(); ++i) {
	//scalar_type v = volume_of_simplex(em->pm->points_of_convex(em->cvlist[i]));
	//cerr << "volume du simplex " << em->cvlist[i] << ": v=" << v << endl;
	//bgeot::pconvex_ref cvr = m.trans_of_convex(cvlist[i])->convex_ref();
	const getfem_mesh *cvm = em->refined_convexes_mesh[i];
	dal::bit_vector bv = cvm->convex_index();
	size_type cv;
	for (cv << bv; cv != size_type(-1); cv << bv) {
	  std::vector<size_type> P(cvm->nb_points_of_convex(cv));
	  std::copy(cvm->ind_points_of_convex(cv).begin(), cvm->ind_points_of_convex(cv).end(), P.begin());
	  //if (v < 0 && P.size() > 1) std::swap(P[0],P[1]);
	  for (size_type j = 0; j < cvm->nb_points_of_convex(cv); ++j) {
	    if (!binary) o << em->pts_id[count_pts + P[j]] << " ";
	    else bwriteint(em->pts_id[count_pts + P[j]]);
	  }
	  if (!binary) o << endl;
	}
	count_pts += cvm->nb_points();
	if (!binary) o << endl;
      }
      assert(count_pts == em->pts_id.size());
      std::string s_elem_type("unknown");
      switch (em->sdim) {
      case 1: s_elem_type = "lines"; break;
      case 2: s_elem_type = "triangles"; break;
      case 3: s_elem_type = "tetrahedra"; break;
      default: DAL_WARNING(1, "OpenDX won't handle element types of dimension " << int(em->sdim));
      }
      o << "attribute \"element type\" string \"" << s_elem_type 
	<< "\"\nattribute \"ref\" string \"positions\"" << endl << endl;
    }
    return em;
  }

  template<class VECT> void
  opendx_export::export_solution(const getfem::mesh_fem& mf, const VECT& U, size_type nrefine, const char *field_name, bool continuous_data) {
    _exported_mesh* em = export_mesh_data_no_repeat(&mf.linked_mesh(),nrefine,continuous_data);
    const getfem::getfem_mesh& m = *em->pm;
    size_type data_object_num = object_cnt++;
    size_type Q = mf.get_qdim();
    size_type R = U.size() / mf.nb_dof();
    cerr << "R=" << R << endl;
    std::vector<scalar_type> UU(em->nb_pts * Q * R,0.);
    std::vector<size_type> Ucnt(em->nb_pts,0);
    /* compute data */
    {
      _fem_precomp fp;
      bgeot::stored_point_tab pts, cv_pts, cvm_pts;
      base_vector coeff;
      base_node val(Q);
      pfem fem = 0, femprec = 0;
      getfem::base_matrix G;
      size_type pcnt = 0;
      for (size_type i=0; i < em->cvlist.size(); ++i) {
	size_type cv = em->cvlist[i];
	const getfem_mesh *cvm = em->refined_convexes_mesh[i];
	getfem::ref_mesh_dof_ind_ct cv_dof = mf.ind_dof_of_element(cv);
	size_type nbd = mf.nb_dof_of_element(cv);
	femprec = fem;
	fem = mf.fem_of_element(cv);

	if (!fem->is_equivalent()) {
	  cv_pts.resize(m.nb_points_of_convex(cv));
	  std::copy(m.points_of_convex(cv).begin(),m.points_of_convex(cv).end(),cv_pts.begin());
	  getfem::transfert_to_G(G, cv_pts);
	}
	else G.clear();

	if (i == 0 || cvm != em->refined_convexes_mesh[i-1] || fem != femprec) {
	  cvm_pts.resize(cvm->nb_points());
	  std::copy(cvm->points().begin(), cvm->points().end(), cvm_pts.begin());
	  fem_precomp_not_stored(fem, &cvm_pts, fp);
	}
	coeff.resize(nbd);
	for (size_type r = 0; r < R; ++r) {
	  for (size_type j = 0; j < coeff.size(); ++j) {
	    coeff[j] = U[R*cv_dof[j]+r];
	  }
	  for (size_type j = pcnt; j < pcnt+cvm->nb_points(); ++j) {
	    fem->interpolation(&fp, j-pcnt, G, m.trans_of_convex(cv), coeff, val, Q);
	    for (size_type q = 0; q < Q; ++q) {
	      //cerr << "out = " << R*(em->pts_id[pcnt]*Q + q)+r << "cv = " << cv << ", pid = " << em->pts_id[pcnt] << ", q=" << q << ", val=" << val[q] << endl;
	      UU[R*(em->pts_id[j]*Q + q)+r] += val[q];
	    }
	    Ucnt[em->pts_id[j]]++;
	  }
	}
	pcnt += cvm->nb_points();
      }

      for (size_type i=0; i < Ucnt.size(); ++i) 
	for (size_type r = 0; r < R; ++r)
	  for (size_type q = 0; q < Q; ++q)
	    UU[R*(i*Q +q)+r] /= scalar_type(R*Ucnt[i]);
    }
    /* write data */

    o << "object " << data_object_num << " class array type float rank ";
    if (Q == 1) o << "0";     /* scalar data */
    else if (R == 1) o << "1 shape " << Q; /* or vector data */
    else o << "2 shape " << R << " " << Q; /* or tensor data */
    o << " items " << em->nb_pts;
    if (binary) o << " " << indianness() << " binary";
    o << " data follows" << endl;
    for (size_type i=0; i < UU.size(); ++i) { 
      if (!binary) {
	o << UU[i]; if (((i+1) % 10) == 0) o << "\n"; else o << " "; 
      } else bwritefloat(UU[i]);
    }

    o << endl << "attribute \"dep\" string \"positions\"" << endl;
    
    /* write footer */
    o << "object \"" << field_name << "\" class field" << endl
      << "component \"positions\" value " << em->pts_object_num << endl
      << "component \"connections\" value " << em->cvs_object_num << endl 
      << "component \"data\" value " << data_object_num << endl << endl;
  }

  void opendx_export::export_mesh_edges(const getfem::getfem_mesh& m, size_type nrefine) {
    bgeot::edge_list el;
    bgeot::mesh_edge_list(m, el, false);
    std::vector<scalar_type> epts(m.dim() * (nrefine+1) * el.size());
    transform_edge_list(m, nrefine, el, epts.begin());
    size_type edges_pts_object_num = object_cnt++;
    size_type edges_conn_object_num = object_cnt++;
    /* write point list */
    {
      o << "object " << edges_pts_object_num << " class array type float rank 1 shape " 
	<< int(m.dim()) 
	<< " items " << (nrefine+1) * el.size();
      if (binary) o << " " << indianness() << " binary";
      o << " data follows" << endl;      
      for (size_type i=0; i < epts.size(); ++i) { 
	if (!binary) { o << epts[i]; if ((i+1) % m.dim()) o << " "; else o << endl; }
	else bwritefloat(epts[i]);
      }
    }

    /* write connections (edge) list */
    {
      o << "object " << edges_conn_object_num << " class array type int rank 1 shape 2 "
	<< " items " << nrefine*el.size(); 
      if (binary) o << " " << indianness() << " binary";
      o << " data follows" << endl;
      for (size_type i=0; i < el.size(); ++i) { 
	for (size_type j=0; j < nrefine; ++j) {
	  if (!binary) { o << i * (nrefine+1) + j << " " << i * (nrefine+1) + j+1 << endl;}
	  else { bwriteint(i * (nrefine+1) + j); bwriteint(i * (nrefine+1) + j+1); }
	}
      }
      o << "attribute \"element type\" string \"lines\"" << endl
	<< "attribute \"ref\" string \"positions\"" << endl << endl;
    }

    o << "object \"mesh\" class field" << endl
      << "component \"positions\"  value " << edges_pts_object_num << endl
      << "component \"connections\"  value " << edges_conn_object_num << endl << endl;
  }
}

class cmd_line_opts {
public:
  typedef std::vector<const char*> argv_type;
  argv_type argv;
  int carg_num;
  std::string carg_name;
public:
  cmd_line_opts(int _argc, char *_argv[]) : carg_num(0) {
    argv.reserve(_argc);
    for (int i=1; i < _argc; ++i) argv.push_back(_argv[i]);
  }
  bool operator()(const char *s) {    
    if (search_arg(s) != -1) { pop(); return true; }
    else return false;
  }
  std::string pop() {
    std::string s;
    if (carg_num < int(argv.size())) s = argv[carg_num];
    else DAL_THROW(dal::failure_error, "not enough arguments for option " << carg_name);
    argv.erase(argv.begin()+carg_num);
    return s;
  }

private:
  int search_arg(const char *s) {
    int i;
    carg_name = s;
    carg_num = -1;
    for (i=0; i < int(argv.size()); ++i)
      if (strcmp(argv[i],s) == 0)
	carg_num = i; 
    return carg_num;
  }
};

void
check_empty(getfem::getfem_mesh *m=0, getfem::mesh_fem *mf=0, std::vector<scalar_type> *U=0) {
  if (m && m->nb_points() == 0) DAL_THROW(dal::failure_error, "a mesh is required !");
  if (mf && mf->nb_dof() == 0) DAL_THROW(dal::failure_error, "a mesh_fem is required !");
  if (U && U->size() == 0) DAL_THROW(dal::failure_error, "a solution is required !");
}

int 
main(int argc, char *argv[]) {  
  cmd_line_opts args(argc,argv);
  if (args("-help")) {
    cerr << "syntax blahblah\n"; return 0;
  }


  try {
    getfem::getfem_mesh m;
    getfem::mesh_fem mf(m);
    std::vector<scalar_type> U;
    size_type Nrefine = 2;
    if (args("-mesh")) {
      std::string s = args.pop();
      cout << "loading getfem mesh from " << s << endl;
      m.read_from_file(s);
    }
    if (args("-mf")) {
      std::string s = args.pop();
      cout << "loading getfem mesh_fem from " << s << endl;
      mf.read_from_file(s);
    }
    if (args("-read_solution")) {
      std::string s = args.pop();
      cout << "loading mesh, mesh_fem and solution from " << s << endl;
      short_type K;
      getfem::load_solution(s,m,mf,U,K);
    }
    if (args("-import")) {
      std::string fmt = args.pop();
      std::string fname = args.pop();
      cout << "importing mesh from " << fmt << " file '" << fname << "'" << endl;
      getfem::import_mesh(fname, fmt, m);
    }
    if (args("-data")) {
      check_empty(&m,&mf);
      std::ifstream f(args.pop().c_str());
      std::deque<scalar_type> u;
      while (f.good()) {
	scalar_type v;
	f >> v; if (f.good()) u.push_back(v);
      }
      U.resize(u.size()); std::copy(u.begin(),u.end(), U.begin());
      if (U.size() % mf.nb_dof()) {
	cerr << "wrong data or wrong mesh_fem : data has " << U.size() << " values which is not a multiple of the number of dof of the mesh_fem (" << mf.nb_dof() << ")" << endl;
	return 1;
      } else {
	cerr << "read " << U.size() << " scalar values (";
	if (U.size()/mf.nb_dof() == 1)
	  if (mf.get_qdim() == 1) cerr << "scalar field"; else cerr << "vector field";
	else
	  if (mf.get_qdim() == 1) cerr << "vector field"; else cerr << "tensor field";
	cerr << ")" << endl;
      }
      for (size_type i = 0; i < mf.nb_dof(); ++i) f >> U[i];      
    }

    if (m.nb_points()) {
      cerr << "Mesh: dimension = " << m.dim() << ", convexes=" << m.nb_convex() 
	   << ", pts = " << m.nb_points() << endl;
    }
    if (mf.nb_dof()) {
      cerr << "Mesh_Fem: nb_dof = " << mf.nb_dof() 
	   << ", Qdim = " << int(mf.get_qdim()) << endl;
      U.resize(mf.nb_dof()); U[0] = 1;
    }
    if (args("-refine")) {
      Nrefine = atoi(args.pop().c_str());
    }
    if (args("-export")) {
      std::string fmt = args.pop();
      std::string fname = args.pop();
      cout << "exporting to " << fmt << " file '" << fname << "'" << endl;
      if (fmt.compare("opendx")==0) {
	check_empty(&m,&mf,&U);
	getfem::opendx_export ex(fname, !args("-ascii"));
	ex.export_mesh_edges(mf.linked_mesh(),Nrefine);
	ex.export_solution(mf,U,Nrefine,"deformation",args("-cont"));
      }
    }
  } DAL_STANDARD_CATCH_ERROR;
}
