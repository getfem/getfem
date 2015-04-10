.. $Id: install.rst 4738 2014-07-27 12:25:54Z renard $

.. include:: ../replaces.txt

.. _tut-thermo_elec_coupling:

Example of Thermo-elastic and electrical coupling (|gf| 5.0)
============================================================

This example aims to present a simple example of a multiphysics problem with a nonlinear coupling of a displacement field, a temperature field and en electric potential field. It also aims to compare the use of the C++ librairy and the different interfaces. The corresponding demo file is present in the test directories of |gf| (`tests/`, `interface/tests/python`, `interface/scr/scilab/demos` and `interface/tests/matlab`).

The problem setting
-------------------


The weak formulation
--------------------


Implementation in C++ and with the interface
--------------------------------------------




Initialisation
**************




========== ================================================
**C++**    .. code-block:: c++                             

             #include "getfem/getfem_model_solvers.h" 
             #include "getfem/getfem_export.h"
             #include "gmm/gmm.h"
             #include "getfem/getfem_mesher.h"
             #include "getfem/getfem_generic_assembly.h"
           
             using bgeot::size_type;
             using bgeot::base_node;
             using bgeot::base_small_vector;
             typedef getfem::model_real_plain_vector plain_vector;

             int main(void) {
---------- ------------------------------------------------
**Python** .. code-block:: python                                     

             import getfem as gf
             import numpy as np
---------- ------------------------------------------------
**Scilab** .. code-block:: matlab                             

             gf_workspace('clear all'); // for multiple runs
---------- ------------------------------------------------
**Matlab** .. code-block:: matlab                                    

             gf_workspace('clear all'); % for multiple runs
========== ================================================


Parameters of the model
***********************

|degre|



=========== ================================================
**C++**     .. code-block:: c++                             

                double epsilon = 1.; // Thickness of the plate (cm)
                double E = 21E6;     // Young Modulus (N/cm^2)
                double nu = 0.3;     // Poisson ratio
                double clambda = E*nu/((1+nu)*(1-2*nu));
                double cmu = E/(2*(1+nu));
                double clambdastar = 2*clambda*cmu/(clambda+2*cmu);
                double F = 100E2;    // Force density at the right boundary (N/cm^2)
                double kappa = 4.;   // Thermal conductivity (W/(cm K))
                double D = 10.;      // Heat transfert coefficient (W/(K cm^2))
                double air_temp = 20;// Temperature of the air in oC.
                double alpha_th = 16.6E-6; // Thermal expansion coefficient (/K)
                double T0 = 20.;     // Reference temperature in oC
                double rho_0 = 1.754E-8; // Resistance temperature coefficient at T0
                double alpha = 0.0039; // Second resistance temperature coefficient

                double h = 2.        // Approximate mesh size
                bgeot::dim_type elements_degree = 2; // Degree of the finite element methods
----------- ------------------------------------------------
**Scripts** .. code-block:: python                                     

                epsilon = 1.; E = 21E6; nu = 0.3;
                clambda = E*nu/((1+nu)*(1-2*nu));
                cmu = E/(2*(1+nu));
                clambdastar = 2*clambda*cmu/(clambda+2*cmu);
                F = 100E2; kappa = 4.; D = 10;
                air_temp = 20; alpha_th = 16.6E-6;
                T0 = 20; rho_0 = 1.754E-8;
                alpha = 0.0039;

                h = 2; elements_degree = 2;     

=========== ================================================




Mesh generation
***************



========== ===========================================================================
**C++**    .. code-block:: c++                             

                getfem::mesh mesh;
                getfem::mesher_rectangle mo1(base_node(0., 0.), base_node(100., 25.));
                getfem::mesher_ball mo2(base_node(25., 12.5), 8.);
                getfem::mesher_ball mo3(base_node(50., 12.5), 8.);
                getfem::mesher_ball mo4(base_node(75., 12.5), 8.);
                getfem::mesher_union mo5(mo2, mo3, mo4);
                getfem::mesher_setminus mo(mo1, mo5);

                std::vector<getfem::base_node> fixed;
                getfem::build_mesh(mesh, mo, h, fixed, 2, -2);
---------- ---------------------------------------------------------------------------
**Python** .. code-block:: python                                     

                mo1 = gf.MesherObject('rectangle', [0., 0.], [100., 25.])
                mo2 = gf.MesherObject('ball', [25., 12.5], 8.)
                mo3 = gf.MesherObject('ball', [50., 12.5], 8.)
                mo4 = gf.MesherObject('ball', [75., 12.5], 8.)
                mo5 = gf.MesherObject('union', mo2, mo3, mo4)
                mo  = gf.MesherObject('set minus', mo1, mo5)

                mesh = gf.Mesh('generate', mo, h, 2)
---------- ---------------------------------------------------------------------------
**Scilab** .. code-block:: matlab                             

                mo1 = gf_mesher_object('rectangle', [0 0], [100 25]);
                mo2 = gf_mesher_object('ball', [25 12.5], 8);
                mo3 = gf_mesher_object('ball', [50 12.5], 8);
                mo4 = gf_mesher_object('ball', [75 12.5], 8);
                mo5 = gf_mesher_object('union', mo2, mo3, mo4);
                mo  = gf_mesher_object('set minus', mo1, mo5);

                mesh = gf_mesh('generate', mo, h, 2);
---------- ---------------------------------------------------------------------------
**Matlab** .. code-block:: matlab

                mo1 = gf_mesher_object('rectangle', [0 0], [100 25]);
                mo2 = gf_mesher_object('ball', [25 12.5], 8);
                mo3 = gf_mesher_object('ball', [50 12.5], 8);
                mo4 = gf_mesher_object('ball', [75 12.5], 8);
                mo5 = gf_mesher_object('union', mo2, mo3, mo4);
                mo  = gf_mesher_object('set minus', mo1, mo5);

                mesh = gf_mesh('generate', mo, h, 2);
========== ===========================================================================



Boundary selection
******************


========== ===========================================================================
**C++**    .. code-block:: c++                             

                getfem::mesh_region border_faces;
                getfem::outer_faces_of_mesh(mesh, border_faces);
                getfem::mesh_region fb1
                  = getfem::select_faces_in_box(mesh, border_faces, base_node(1., 1.),
                                 base_node(99., 24.));
                getfem::mesh_region fb2
                  = getfem::select_faces_of_normal(mesh, border_faces,
                                     base_small_vector( 1., 0.), 0.01);
                getfem::mesh_region fb3
                  = getfem::select_faces_of_normal(mesh, border_faces,
                                     base_small_vector(-1., 0.), 0.01);
                getfem::mesh_region fb4
                  = getfem::select_faces_of_normal(mesh, border_faces,
                                     base_small_vector(0.,  1.), 0.01);
                getfem::mesh_region fb5
                  = getfem::select_faces_of_normal(mesh, border_faces,
                                     base_small_vector(0., -1.), 0.01);

                size_type RIGHT_BOUND=1, LEFT_BOUND=2, TOP_BOUND=3, BOTTOM_BOUND=4;
                mesh.region( RIGHT_BOUND) = getfem::mesh_region::subtract(fb2, fb1);
                mesh.region(  LEFT_BOUND) = getfem::mesh_region::subtract(fb3, fb1);
                mesh.region(   TOP_BOUND) = getfem::mesh_region::subtract(fb4, fb1);
                mesh.region(BOTTOM_BOUND) = getfem::mesh_region::subtract(fb5, fb1);
---------- ---------------------------------------------------------------------------
**Python** .. code-block:: python                                     

                fb1 = mesh.outer_faces_in_box([1., 1.], [99., 24.])
                fb2 = mesh.outer_faces_with_direction([ 1., 0.], 0.01)
                fb3 = mesh.outer_faces_with_direction([-1., 0.], 0.01)
                fb4 = mesh.outer_faces_with_direction([0.,  1.], 0.01)
                fb5 = mesh.outer_faces_with_direction([0., -1.], 0.01)

                RIGHT_BOUND=1; LEFT_BOUND=2; TOP_BOUND=3; BOTTOM_BOUND=4; HOLE_BOUND=5;

                mesh.set_region( RIGHT_BOUND, fb2)
                mesh.set_region(  LEFT_BOUND, fb3)
                mesh.set_region(   TOP_BOUND, fb4)
                mesh.set_region(BOTTOM_BOUND, fb5)
                mesh.set_region(  HOLE_BOUND, fb1)
                mesh.region_subtract( RIGHT_BOUND, HOLE_BOUND)
                mesh.region_subtract(  LEFT_BOUND, HOLE_BOUND)
                mesh.region_subtract(   TOP_BOUND, HOLE_BOUND)
                mesh.region_subtract(BOTTOM_BOUND, HOLE_BOUND)
---------- ---------------------------------------------------------------------------
**Scilab** .. code-block:: matlab                             

                fb1 = gf_mesh_get(mesh, 'outer faces in box', [1 1], [99 24]);
                fb2 = gf_mesh_get(mesh, 'outer faces with direction', [ 1 0], 0.01);
                fb3 = gf_mesh_get(mesh, 'outer faces with direction', [-1 0], 0.01);
                fb4 = gf_mesh_get(mesh, 'outer faces with direction', [0  1], 0.01);
                fb5 = gf_mesh_get(mesh, 'outer faces with direction', [0 -1], 0.01);

                RIGHT_BOUND=1; LEFT_BOUND=2; TOP_BOUND=3; BOTTOM_BOUND=4; HOLE_BOUND=5;
                gf_mesh_set(mesh, 'region',  RIGHT_BOUND, fb2);
                gf_mesh_set(mesh, 'region',   LEFT_BOUND, fb3);
                gf_mesh_set(mesh, 'region',    TOP_BOUND, fb4);
                gf_mesh_set(mesh, 'region', BOTTOM_BOUND, fb5);
                gf_mesh_set(mesh, 'region',   HOLE_BOUND, fb1);
                gf_mesh_set(mesh, 'region subtract',  RIGHT_BOUND, HOLE_BOUND);
                gf_mesh_set(mesh, 'region subtract',   LEFT_BOUND, HOLE_BOUND);
                gf_mesh_set(mesh, 'region subtract',    TOP_BOUND, HOLE_BOUND);
                gf_mesh_set(mesh, 'region subtract', BOTTOM_BOUND, HOLE_BOUND);
---------- ---------------------------------------------------------------------------
**Matlab** .. code-block:: matlab

                fb1 = gf_mesh_get(mesh, 'outer faces in box', [1 1], [99 24]);
                fb2 = gf_mesh_get(mesh, 'outer faces with direction', [ 1 0], 0.01);
                fb3 = gf_mesh_get(mesh, 'outer faces with direction', [-1 0], 0.01);
                fb4 = gf_mesh_get(mesh, 'outer faces with direction', [0  1], 0.01);
                fb5 = gf_mesh_get(mesh, 'outer faces with direction', [0 -1], 0.01);

                RIGHT_BOUND=1; LEFT_BOUND=2; TOP_BOUND=3; BOTTOM_BOUND=4; HOLE_BOUND=5;
                gf_mesh_set(mesh, 'region',  RIGHT_BOUND, fb2);
                gf_mesh_set(mesh, 'region',   LEFT_BOUND, fb3);
                gf_mesh_set(mesh, 'region',    TOP_BOUND, fb4);
                gf_mesh_set(mesh, 'region', BOTTOM_BOUND, fb5);
                gf_mesh_set(mesh, 'region',   HOLE_BOUND, fb1);
                gf_mesh_set(mesh, 'region subtract',  RIGHT_BOUND, HOLE_BOUND);
                gf_mesh_set(mesh, 'region subtract',   LEFT_BOUND, HOLE_BOUND);
                gf_mesh_set(mesh, 'region subtract',    TOP_BOUND, HOLE_BOUND);
                gf_mesh_set(mesh, 'region subtract', BOTTOM_BOUND, HOLE_BOUND);
========== ===========================================================================


Mesh draw
*********

========== ===========================================================================
**C++**    .. code-block:: c++                             

                getfem::vtk_export exp("mesh.vtk", false);
                exp.exporting(mesh);
                exp.write_mesh();
                // You can view the mesh for instance with
                // mayavi2 -d mesh.vtk -f ExtractEdges -m Surface
---------- ---------------------------------------------------------------------------
**Python** .. code-block:: python                                     

                mesh.export_to_vtk('mesh.vtk');
                # You can view the mesh for instance with
                # mayavi2 -d mesh.vtk -f ExtractEdges -m Surface
---------- ---------------------------------------------------------------------------
**Scilab** .. code-block:: matlab

                scf(1);
                gf_plot_mesh(mesh, 'refine', 8, 'curved', 'on', 'regions', ...
                             [RIGHT_BOUND LEFT_BOUND TOP_BOUND BOTTOM_BOUND]);
                title('Mesh');
                sleep(1000);
---------- ---------------------------------------------------------------------------
**Matlab** .. code-block:: matlab

                gf_plot_mesh(mesh, 'refine', 8, 'curved', 'on', 'regions', ...
                             [RIGHT_BOUND LEFT_BOUND TOP_BOUND BOTTOM_BOUND]);
                title('Mesh');
                pause(1);
========== ===========================================================================


