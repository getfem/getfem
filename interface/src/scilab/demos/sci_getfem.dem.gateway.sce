// ====================================================================
// Yann COLLETTE
// Copyright 2009-2020
// This file is released into the public domain
// ====================================================================
demopath = get_absolute_file_path("sci_getfem.dem.gateway.sce");

subdemolist = ["Bilaplacian",                     "demo_bilaplacian.sce"; ..
               "Crack",                           "demo_crack.sce"; ..
	       "Continuation method",             "demo_continuation.sce"; ..
	       "Convection in a rotating cavity", "demo_convection_rotating_cavity.sce"; ..
	       "Fictitious domains",              "demo_fictitious_domains.sce"; ..
	       "Final Laplacian",                 "demo_finallaplacian.sce"; ..
	       "Laplacian 1D",                    "demo_laplacian1D.sce"; ..
	       "Laplacian",                       "demo_laplacian.sce"; ..
	       "Mesh generation",                 "demo_mesh_generation.sce"; ..
	       "Mortar",                          "demo_mortar.sce"; ..
	       "Navier-Stokes",                   "demo_Navier_Stokes.sce"; ..
	       "Nonlinear Elasticity",            "demo_nonlinear_elasticity.sce"; ..
	       "Plasticity",                      "demo_plasticity.sce"; ..
	       "Plate",                           "demo_plate.sce"; ..
	       "Refine",                          "demo_refine.sce"; ..
	       "Slice",                           "demo_slices.sce"; ..
	       "Static Contact",                  "demo_static_contact.sce"; ..
	       "Step by Step",                    "demo_step_by_step.sce"; ..
	       "Stokes 2D Poiseuille Arc",        "demo_stokes_2D_poiseuille_arc.sce"; ..
	       "Stokes 2D Poiseuille",            "demo_stokes_2D_poiseuille.sce"; ..
	       "Stokes 2D tube",                  "demo_stokes_2D_tube.sce"; ..
	       "Stokes 3D tank",                  "demo_stokes_3D_tank.sce"; ..
	       "Structural optimization",         "demo_structural_optimization.sce"; ..
	       "topological optimization",        "demo_topological_optimization.sce"; ..
	       "Tripod",                          "demo_tripod.sce"; ..
	       "Wave 2D",                         "demo_wave2D.sce"; ..
	       "Wave equation",                   "demo_wave_equation.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
