/*===========================================================================

 Copyright (C) 2002-2015 Jean-Yves Heddebaut.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

/**
   @file nonlinear_membrane.cc
   @brief Nonlinear membrane problem (large strain).
   Due to Jean-Yves Heddebaut <jyhed@alpino.be>.
*/

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_nonlinear_elasticity.h"
#include "getfem/getfem_mesh_fem.h"
#include "getfem/getfem_superlu.h"
#include "gmm/gmm.h"

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::base_vector;
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/*
  structure for the membrane problem
*/
struct membrane_problem {
  
	enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1};
	getfem::mesh mesh;         /* the mesh */
	getfem::mesh_im mim;       /* the integration methods */
	getfem::mesh_fem mf_u;     /* main mesh_fem, for the membrane solution */
	getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
	getfem::mesh_fem mf_coef;  /* mesh_fem used to represent pde coefficients  */
	getfem::mesh_fem mf_vm;		/* mesh_fem used to represent pde coefficients for tresca/von mises  */
	scalar_type p0, p1, p2, p3;    /* elastic coefficients.                        */
	scalar_type p4,p5;			//pretension on X' and Y' axis
	scalar_type residual;		/* max residual for the iterative solvers         */
	unsigned int print_convexes;	/*if 1, print convex description*/
	int bdy_type;				/*specify type of BL imposed on specific convex faces, 0 for dirichlet, 1 for Neumann, -1 if no BL on specific elements*/
	std::string datafilename;
	bgeot::md_param PARAM;
	
	bool solve  (plain_vector &U,getfem::base_vector &VM);
	void init(void);
	membrane_problem(void) : mim(mesh), mf_u(mesh), mf_rhs(mesh), mf_coef(mesh),mf_vm(mesh) {};
	
};



/* Read parameters from the .param file, build the mesh, set finite element
   and integration methods and selects the boundaries.
 */
 
void membrane_problem::init(void) {
  print_convexes = unsigned(PARAM.int_value("PRINT_CONVEXES", "1 to print convex description"));
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string VONMISES_FEM_TYPE  = PARAM.string_value("VONMISES_FEM_TYPE","Von Mises FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION","Name of integration method");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "VONMISES_FEM_TYPE="  << VONMISES_FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";
  
  /* First step : build the mesh */
  bgeot::pgeometric_trans pgt =bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type NMesh = pgt->dim();
  cout<<"NMesh="<<NMesh<<endl;
  std::vector<size_type> nsubdiv(NMesh);
  nsubdiv[0]=PARAM.int_value("NX", "Nomber of space steps ");
  nsubdiv[1]=PARAM.int_value("NY", "Nomber of space steps ");
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,PARAM.int_value("MESH_NOISED") != 0);
  
  bgeot::base_matrix M(NMesh,NMesh);
  for (size_type i=0; i < NMesh; ++i) {
    static const char *t[] = {"LX","LY"};
    M(i,i) = PARAM.real_value(t[i],t[i]);
  }
  /* scale the unit mesh to [LX,LY,..] and incline it */
  mesh.transformation(M);
  
  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;
  p0 = PARAM.real_value("P0", "Young modulus along X' axis");
  p1 = PARAM.real_value("P1", "Poison coefficient");
  p2 = PARAM.real_value("P2", "Young modulus along Y' axis");
  p3 = PARAM.real_value("P3", "Shear modulus");
  p4 = PARAM.real_value("PRETENSION_X", "Pretension on X' axis");
  p5 = PARAM.real_value("PRETENSION_Y", "Pretension on Y' axis");
  
  mf_u.set_qdim(3);	//3dim def on 2dim mesh
  
  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  getfem::pfem pf_vm = getfem::fem_descriptor(VONMISES_FEM_TYPE);
  getfem::pintegration_method ppi = getfem::int_method_descriptor(INTEGRATION);
  mim.set_integration_method(ppi);
  mf_u.set_finite_element(pf_u);
  mf_vm.set_finite_element(pf_vm);
  
  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }
  
  
  
  /* set the finite element on mf_coef. Here we use a very simple element
   *  since the only function that need to be interpolated on the mesh_fem 
   * is f(x)=1 ... */
  mf_coef.set_finite_element(mesh.convex_index(),getfem::classical_fem(pgt,0));
  
  /* set boundary conditions*/
  
  /* BL on outer faces*/
  
  int  bdyUp = int(PARAM.int_value("bdyUp", "0 if dirichlet, 1 if neumann, -1 if no BL on upper face"));
  int  bdyLow = int(PARAM.int_value("bdyLow", "0 if dirichlet, 1 if neumann, -1 if no BL on lower face"));
  int  bdyLeft = int(PARAM.int_value("bdyLeft", "0 if dirichlet, 1 if neumann, -1 if no BL on left face"));
  int  bdyRight = int(PARAM.int_value("bdyRight", "0 if dirichlet, 1 if neumann, -1 if no BL on right face"));
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::mr_visitor it(border_faces); !it.finished(); ++it) {
    assert(it.is_face());
    base_node un = mesh.normal_of_face_of_convex(it.cv(), it.f());
    un /= gmm::vect_norm2(un);
    if ((gmm::abs(un[NMesh-1] - 1.0) < 1.0E-7) ){//lower faces 
      if (bdyLow==0){
	mesh.region(DIRICHLET_BOUNDARY_NUM).add(it.cv(),it.f());
	cout<<"dirichlet Low,convex n="<<it.cv()<<"face="<<it.f()<<" n="<<un<<endl;
      }else if (bdyLow==1){
	mesh.region(NEUMANN_BOUNDARY_NUM).add(it.cv(),it.f());
	cout<<"neumann Low,convex n="<<it.cv()<<"face="<<it.f()<<" n="<<un<<endl;
      }
    } else if ((gmm::abs(un[NMesh-1] + 1.0) < 1.0E-7)) {//upper faces
      if (bdyUp==0){
	mesh.region(DIRICHLET_BOUNDARY_NUM).add(it.cv(),it.f());
	cout<<"dirichlet Up,convex n="<<it.cv()<<"face="<<it.f()<<" n="<<un<<endl;
      }else if (bdyUp==1){
	mesh.region(NEUMANN_BOUNDARY_NUM).add(it.cv(),it.f());
	cout<<"neumann Up,convex n="<<it.cv()<<"face="<<it.f()<<" n="<<un<<endl;
      }
    } else if ((gmm::abs(un[0] + 1.0) < 1.0E-7) ) {//left side
      if (bdyLeft==0){
	mesh.region(DIRICHLET_BOUNDARY_NUM).add(it.cv(),it.f());
	cout<<"dirichlet Left,convex n="<<it.cv()<<"face="<<it.f()<<" n="<<un<<endl;
      }else if (bdyLeft==1){
	mesh.region(NEUMANN_BOUNDARY_NUM).add(it.cv(),it.f());
	cout<<"neumann Left,convex n="<<it.cv()<<"face="<<it.f()<<" n="<<un<<endl;
      }
    } else if( (gmm::abs(un[0] - 1.0) < 1.0E-7) ) {//right side
      if (bdyRight==0){
	mesh.region(DIRICHLET_BOUNDARY_NUM).add(it.cv(),it.f());
	cout<<"dirichlet Right,convex n="<<it.cv()<<"face="<<it.f()<<" n="<<un<<endl;
      }else if (bdyRight==1){
	mesh.region(NEUMANN_BOUNDARY_NUM).add(it.cv(),it.f());
	cout<<"neumann Right,convex n="<<it.cv()<<"face="<<it.f()<<" n="<<un<<endl;
      }
    }
  }
  
  
  /* BL on specific elements*/
  
  bdy_type= int(PARAM.int_value("bdy_type", "set to 0 to impose dirichlet BL on next convex faces, 1 to set neumann BL, -1 not to set any BL on this specific elements"));
  long int bdy_elements[8][2];
  bdy_elements [0][0] = PARAM.int_value("bdy_element1", "element nbr , set to -1 if not used");
  bdy_elements [0][1]= PARAM.int_value("bdy_face1", "face local nbr of BL element 1 to be  fixed");
  bdy_elements [1][0] = PARAM.int_value("bdy_element2", "element nbr , set to -1 if not used");
  bdy_elements [1][1]= PARAM.int_value("bdy_face2", "face local nbr of BL element 2 to be  fixed");
  bdy_elements [2][0] = PARAM.int_value("bdy_element3", "element nbr , set to -1 if not used");
  bdy_elements [2][1]= PARAM.int_value("bdy_face3", "face local nbr of BL element 3 to be  fixed");
  bdy_elements [3][0] = PARAM.int_value("bdy_element4", "element nbr , set to -1 if not used");
  bdy_elements [3][1]= PARAM.int_value("bdy_face4", "face local nbr of BL element 4 to be  fixed");
  bdy_elements [4][0] = PARAM.int_value("bdy_element5", "element nbr , set to -1 if not used");
  bdy_elements [4][1]= PARAM.int_value("bdy_face5", "face local nbr of BL element 5 to be  fixed");
  bdy_elements [5][0] = PARAM.int_value("bdy_element6", "element nbr , set to -1 if not used");
  bdy_elements [5][1]= PARAM.int_value("bdy_face6", "face local nbr of BL element 6 to be  fixed");
  bdy_elements [6][0] = PARAM.int_value("bdy_element7", "element nbr , set to -1 if not used");
  bdy_elements [6][1]= PARAM.int_value("bdy_face7", "face local nbr of BL element 7 to be  fixed");
  bdy_elements [7][0] = PARAM.int_value("bdy_element8", "element nbr , set to -1 if not used");
  bdy_elements [7][1]= PARAM.int_value("bdy_face8", "face local nbr of BL element 8 to be  fixed");	
  
  if(bdy_type==0 || bdy_type==1)
    for (unsigned int elm = 0; elm <8; ++elm) 
      if(bdy_elements [elm][1]>=0) {
        mesh.region(bdy_type).add(bdy_elements[elm][0],
                                  bgeot::short_type(bdy_elements[elm][1]));
        cout << "BL " << (bdy_type==0 ? "Dirichlet" : "Neumann")
             << " on individual element ,convex n=" << bdy_elements [elm][0]
             << "face=" << bdy_elements [elm][1] << endl;
      }
  
}


/*
  Set initial random displacements on vertical dofs if not fixed. This will avoid tangent matrix singularities along the axis perpendicular to the element.
  This is not needed if there is a pretension.
*/

void setInitialDisp(plain_vector &U,
		    size_type nb_dof,
		    dal::bit_vector dirichlet_dof,
		    float initialDispAmplitude){
  for (size_type dof = 0; dof < nb_dof; ++dof) {
    if(!dirichlet_dof[dof]&& (dof%3==2)){
      U[dof]=(double(rand()%100)*initialDispAmplitude/100)-.5*initialDispAmplitude;
    }
  }
  //	cout<<"init displ="<<U<<endl;
}

/*
  print convex description 
*/


void convexDescription (getfem::mesh_fem &mymesh_fem){
  dal::bit_vector nn = mymesh_fem.linked_mesh().convex_index(); 
  bgeot::size_type i;
  cout<<"Convex description:"<<endl<<endl; 
  for (i << nn; i != bgeot::size_type(-1); i << nn) 
    { 
      cout << "Convex " << i << endl; 
      bgeot::pconvex_structure cvs = mymesh_fem.linked_mesh().structure_of_convex(i); 
      cout << "Number of vertices : " << cvs->nb_points() << endl; 
      cout << "Number of faces : " << cvs->nb_faces() << endl; 
      for (bgeot::short_type f = 0; f < cvs->nb_faces(); ++f) 
	{ 
	  cout << "face " << f << " has " << cvs->nb_points_of_face(f)<<" points,"; 
	  cout << " with global indexes : "<<endl; 
	  for (bgeot::short_type k = 0; k < cvs->nb_points_of_face(f); ++k) 
	    {
	      cout << mymesh_fem.linked_mesh().ind_points_of_convex(i)[cvs->ind_points_of_face(f)[k]] << " "; 
	      cout << " of coordinates "<< mymesh_fem.linked_mesh().points()[mymesh_fem.linked_mesh().ind_points_of_convex(i)[cvs->ind_points_of_face(f)[k]]] << " "<<endl; 
	    }
	} 
    }
  
  GMM_ASSERT1(!mymesh_fem.is_reduced(), "To be adapted");
  
  //print dof coordinate 
  cout<<"Dof of points"<<endl<<endl;
  for  (size_type ii = 0; ii < mymesh_fem.nb_basic_dof(); ++ii){
    cout<<"dof "<< ii<<"coord="<<mymesh_fem.point_of_basic_dof(ii)<<endl;
  }
  
}//end convexDescription




/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/
bool membrane_problem::solve  (plain_vector &U,getfem::base_vector &VM) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type NMesh = mesh.dim();
  float punctualForce1,punctualForce2;
  //  size_type punctualDof1 = PARAM.int_value("PUNCTUAL_DOF1","dof nbr where punctual force is applied");
  // size_type punctualDof2 = PARAM.int_value("PUNCTUAL_DOF2","dof nbr where punctual force is applied");
  float totalPunctualForce1 = float(PARAM.real_value("PUNCTUAL_FORCE1","intensity of punctual force"));
  float totalPunctualForce2 = float(PARAM.real_value("PUNCTUAL_FORCE2","intensity of punctual force"));
  int INITIAL_DISP  = int(PARAM.int_value("INITIAL_DISP", "1 to set init vert displ on free dof"));
  int nb_step = int(PARAM.int_value("NBSTEP"));
  base_vector p(6); p[0] = p0; p[1] = p1; p[2] = p2;p[3] = p3;p[4] = p4;p[5] = p5;
  unsigned int NFem =mf_u.get_qdim();
  
  if(print_convexes==1) convexDescription(mf_u);
  
  
  getfem::abstract_hyperelastic_law *pl = new getfem::membrane_elastic_law();
  
  cout << "law=membrane, parameters=" << p <<endl;

  getfem::model md;
  md.add_fem_variable("u", mf_u);

  md.add_initialized_fixed_size_data("params", p);
  add_nonlinear_elasticity_brick(md,  mim, "u", *pl, "params");  
  
  // Defining the volumic source term.
  unsigned int src_type =  unsigned(PARAM.int_value("src_type","type of source term, 0 for volumic, 1 if src term if applied to neumann bdy"));
  base_vector fsrc(NFem);
  fsrc[0] = PARAM.real_value("FORCEX","Amplitude of the source term");
  fsrc[1] = PARAM.real_value("FORCEY","Amplitude of the source term");
  fsrc[2] = PARAM.real_value("FORCEZ","Amplitude of the source term");
  
  //src term is applied to neumann limit region if src_type=1,otherwize volumic src term on all convexes
  //	jyh:to set punctual force on individual dof, use getfem_modeling_jyh.h iso getfem_modeling.h (set alias to /Users/oracle/fem/getfem_modeling_jyh.h
  //	in /usr/local/include/getfem)
  
  md.add_initialized_fixed_size_data("VolumicData", fsrc);

  if (src_type==1)
    getfem::add_source_term_brick(md, mim, "u", "VolumicData",NEUMANN_BOUNDARY_NUM);
  else
    getfem::add_source_term_brick(md, mim, "u", "VolumicData");
  
  // Dirichlet condition using dirichlet brick
  
  //to allow to impose an expansion of the domain (one border is moved towards positive values of the corresponding axis, the opposite border if moved
  //the same distance towards negative values), set this param to 1
  bool opposite_bdy_reversed = (1== PARAM.int_value("opposite_bdy_reversed","if set to 1,imposed displ applied to dirichlet bdy will be inversed on the opposite bdy"));
  base_vector LH(NFem);	//half dimensions of domain
  LH[0]=PARAM.real_value("LX")/2;
  LH[1]=PARAM.real_value("LY")/2;
  
  base_vector fd(NFem),fdr(NFem);
  fd[0] = PARAM.real_value("dx","Amplitude of the imposed X displacement on dirichlet BL");
  fd[1] = PARAM.real_value("dy","Amplitude of the imposed Y displacement on dirichlet BL");
  fd[2] = PARAM.real_value("dz","Amplitude of the imposed Z displacement on dirichlet BL");
  fdr[2]=fd[2];
  GMM_ASSERT1(!mf_rhs.is_reduced(), "to be adapted");
  plain_vector Frhs(nb_dof_rhs * NFem);
  for (size_type i = 0; i < nb_dof_rhs; ++i) {
    const base_node P = mf_rhs.point_of_basic_dof(i);
    for (size_type j = 0; j < 2; ++j)
      fdr[j]=(opposite_bdy_reversed && P[j]<LH[j] ? int(-1)*fd[j]:fd[j]);
    gmm::copy(fdr, gmm::sub_vector(Frhs, gmm::sub_interval(i*NFem, NFem)));
  }
  
  md.add_initialized_fem_data("DirichletData", mf_rhs, Frhs);
  if (PARAM.int_value("DIRICHLET_VERSION") == 0)
    getfem::add_Dirichlet_condition_with_multipliers
      (md, mim, "u", mf_u, DIRICHLET_BOUNDARY_NUM, "DirichletData");
  else
    getfem::add_Dirichlet_condition_with_penalization
      (md, mim, "u", 1E15, DIRICHLET_BOUNDARY_NUM, "DirichletData");

  
  
  // Generic solver.
  gmm::iteration iter(residual, 1, PARAM.int_value("MAXITER"));
  


  
  //set initial value for MS.state to prevent initial null transverse rigidity
  //not needed if prestressed
  GMM_ASSERT1(!mf_u.is_reduced(), "To be adapted");
  if(INITIAL_DISP==1) {
    setInitialDisp(md.set_real_variable("u"), mf_u.nb_dof(),
		   mf_u.basic_dof_on_region(DIRICHLET_BOUNDARY_NUM),
		   float(PARAM.real_value("initialDispAmplitude",
    "Amplitude of initial displacement set to avoid 3th dim singularities")));
  }
  
  for (int step = 0; step < nb_step; ++step)  {
    //increment external forces
    plain_vector Dfsrc(fsrc);
    gmm::copy(gmm::scaled(fsrc, (step+1.)/(scalar_type)nb_step), Dfsrc);
    gmm::copy(Dfsrc, md.set_real_variable("VolumicData"));
    
    punctualForce1=(totalPunctualForce1)*(float(step+1))/float(nb_step);
    punctualForce2=(totalPunctualForce2)*(float(step+1))/float(nb_step);
    
    /* increment  imposed displacement  */
    plain_vector DFrhs(Frhs);
    gmm::copy(gmm::scaled(Frhs, (step+1.)/(scalar_type)nb_step), DFrhs);
    gmm::copy(DFrhs, md.set_real_variable("DirichletData"));
    
    cout << "step " << step << endl;
    iter = gmm::iteration(residual, int(PARAM.int_value("NOISY")),
                          PARAM.int_value("MAXITER"));
    
    /* let the default non-linear solve (Newton) do its job */
    getfem::standard_solve(md, iter);
    
    
    
    pl->reset_unvalid_flag();
    md.assembly(getfem::model::BUILD_RHS);
    if (pl->get_unvalid_flag()) 
      GMM_WARNING1("The solution is not completely valid, the determinant "
                   "of the transformation is negative on "
                   << pl->get_unvalid_flag() << " gauss points");
    gmm::resize(U, mf_u.nb_dof());
    gmm::copy(md.real_variable("u"), U);		
  }
  
  
  //clean to avoid vtk export error
  gmm::clean(U, 1E-20);
  cout<<"U final="<<U<<endl;
  cout<<"max displ=" << *( std::max_element( U.begin(), U.end() ) )<<endl;
  cout<<"mim displ=" << *( std::min_element( U.begin(), U.end() ) )<<endl;
  
  //compute von mises stress
  // unsigned int PRINT_STRESSES = unsigned(PARAM.int_value("PRINT_STRESSES","1 to print stresses"));
  getfem::compute_Von_Mises_or_Tresca(md, "u", *pl, "params", mf_vm, VM, false);
  
  //clean to avoid vtk export error
  gmm::clean(VM, 1E-20);
  //transform 2D mesh in a 3D mesh to allow correct export of 3D displacements
  bgeot::base_matrix Mt(NFem,NMesh);
  gmm::clear(Mt);
  for (size_type i=0; i < NMesh; ++i) Mt(i,i) = 1;
  mesh.transformation(Mt);
  
  getfem::vtk_export exp(datafilename + ".vtk",PARAM.int_value("VTK_EXPORT")==1);
  //  exp.exporting(mesh);
  exp.exporting(mf_u);
  exp.write_point_data(mf_u, U, "displacement");
  exp.exporting(mf_vm);
  exp.write_point_data(mf_vm,VM, "Von Mises stress");
  cout << "export done, you can view the data file with (for example)\n"
    "mayavi2 -d " << datafilename << ".vtk -f WarpVector -m Surface -m Outline\n";
  
  return (iter.converged());
}

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {
  
  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.
  
  try
    {    
      membrane_problem p;
      p.PARAM.read_command_line(argc, argv);
      p.init();
      p.mesh.write_to_file(p.datafilename + ".mesh");
      p.mf_u.write_to_file(p.datafilename + ".mf", true);
      p.mf_rhs.write_to_file(p.datafilename + ".mfd", true);
      plain_vector U(p.mf_u.nb_dof());
      getfem::base_vector VM(p.mf_vm.nb_dof());
      if (!p.solve(U,VM)) cerr << "Solve has failed\n";
    }
  
  GMM_STANDARD_CATCH_ERROR;	
  return 0; 
}
