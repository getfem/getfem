#include <bgeot_geotrans_inv.h>
#include <getfem_regular_meshes.h>

using bgeot::size_type;
using bgeot::dim_type;
using bgeot::scalar_type;
using bgeot::base_node;
using bgeot::base_vector;

dim_type N, MESH_TYPE;
scalar_type LX, LY, LZ;
size_type NX, NB_POINTS;

ftool::md_param PARAM;

getfem::getfem_mesh mesh;

int main(int argc, char *argv[])
{
  try {

    PARAM.read_command_line(argc, argv);
    N = PARAM.int_value("N", "Domaine dimension");
    NB_POINTS = PARAM.int_value("NB_POINTS", "Nb points");
    LX = PARAM.real_value("LX", "Size in X");
    LY = PARAM.real_value("LY", "Size in Y");
    LZ = PARAM.real_value("LZ", "Size in Y");
    NX = PARAM.int_value("NX", "Nomber of sace steps ");
    MESH_TYPE = PARAM.int_value("MESH_TYPE", "Mesh type ");
 
    cout << "Mesh generation\n";

    base_node org(N); org.fill(0.0);
    std::vector<base_vector> vtab(N);
    std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NX);
    for (dim_type i = 0; i < N; i++) { 
      vtab[i] = base_vector(N); vtab[i].fill(0.0);
      (vtab[i])[i] = ((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NX);
    }
    // if (N > 1) vtab[N-1][0] = incline * LX / scalar_type(NX);
    
    switch (MESH_TYPE) {
    case 0 : getfem::parallelepiped_regular_simplex_mesh
		     (mesh, N, org, vtab.begin(), ref.begin()); break;
    case 1 : getfem::parallelepiped_regular_mesh
		     (mesh, N, org, vtab.begin(), ref.begin()); break;
    case 2 : getfem::parallelepiped_regular_prism_mesh
		     (mesh, N, org, vtab.begin(), ref.begin()); break;
    default : DAL_THROW(dal::internal_error, "Unknown type of mesh");
    }
    
    mesh.optimize_structure();


    scalar_type exectime = ftool::uclock_sec(), total_time = 0.0;

    bgeot::geotrans_inv gti;
    bgeot::base_node pt(N);
    dal::dynamic_array<base_node> ptab;
    dal::dynamic_array<size_type> itab;

    for (size_type i = 0; i < NB_POINTS; ++i) {
      for (dim_type k = 0; k < N; ++k) pt[k] = drand48();
      gti.add_point(pt);
    }

    cout << "Time to sort points : " << ftool::uclock_sec() - exectime << endl;
    total_time += ftool::uclock_sec() - exectime;
    
    dal::bit_vector nn = mesh.convex_index();
    size_type nbtot = 0;
    for (size_type cv = nn.take_first(); cv != size_type(-1); cv << nn) {

      size_type nb = gti.points_in_convex(mesh.convex(cv),
					  mesh.trans_of_convex(cv),
					  ptab, itab);
      nbtot += nb;
    }

    cout << "Time to invert geo trans : " << ftool::uclock_sec() - exectime
	 << endl;
    cout << "Total number : " << nbtot << endl;
    total_time += ftool::uclock_sec() - exectime;
    
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;
}
