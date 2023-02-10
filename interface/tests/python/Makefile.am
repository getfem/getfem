#  Copyright (C) 1999-2020 Yves Renard
#
#  This file is a part of GetFEM
#
#  GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
#  under  the  terms  of the  GNU  Lesser General Public License as published
#  by  the  Free Software Foundation;  either version 3 of the License,  or
#  (at your option) any later version along with the GCC Runtime Library
#  Exception either version 3.1 or (at your option) any later version.
#  This program  is  distributed  in  the  hope  that it will be useful,  but
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#  License and GCC Runtime Library Exception for more details.
#  You  should  have received a copy of the GNU Lesser General Public License
#  along  with  this program;  if not, write to the Free Software Foundation,
#  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

if QHULL
optpy = check_levelset.py                               \
        demo_thermo_elasticity_electrical_coupling.py
else
optpy =
endif


EXTRA_DIST= 						\
	check_export.py 				\
	check_export_vtu.py 				\
	check_global_functions.py			\
	check_levelset.py				\
	check_asm.py    				\
	check_bspline_mesh_fem.py    			\
	check_secondary_domain.py    			\
	check_mixed_mesh.py    				\
	demo_crack.py 					\
	demo_fictitious_domains.py 			\
	demo_laplacian.py 				\
	demo_parallel_laplacian.py 			\
	demo_laplacian_HHO.py 				\
	demo_laplacian_pyramid.py			\
	demo_laplacian_DG.py 				\
	demo_laplacian_aposteriori.py			\
	demo_mortar.py 					\
	demo_phase_field.py              		\
	demo_plasticity.py              		\
	demo_plate.py 					\
	demo_unit_disk.py 				\
        demo_truss.py                                   \
	demo_static_contact.py 				\
	demo_dynamic_contact_1D.py			\
	demo_Mindlin_Reissner_plate.py 			\
        demo_thermo_elasticity_electrical_coupling.py	\
	demo_step_by_step.py				\
	demo_stokes_3D_tank.py          		\
	demo_stokes_3D_tank_draw.py			\
	demo_finite_strain_plasticity.py                \
	demo_finite_strain_plasticity_3D.py             \
	demo_cracked_thermo_elastic_body.py             \
	demo_nonlinear_elasticity.py		 	\
	demo_elasticity_HHO.py                          \
	demo_elastic_ring_contact.py                    \
	demo_wheel_contact.py                           \
	demo_tripod.py 					\
	demo_tripod_alt.py 				\
	demo_wave.py					\
	demo_wave_equation.py				\
	getfem_tvtk.py

if BUILDPYTHON
TESTS = 						\
	check_export.py					\
	check_export_vtu.py				\
	check_global_functions.py			\
	check_asm.py	  				\
	check_bspline_mesh_fem.py	  		\
	check_secondary_domain.py 			\
	check_mixed_mesh.py  				\
	demo_truss.py                                   \
	demo_wave.py					\
	demo_wave_equation.py				\
	demo_laplacian.py				\
	demo_laplacian_HHO.py				\
	demo_elasticity_HHO.py                          \
	$(optpy)

AM_TESTS_ENVIRONMENT = \
	export PYTHONPATH=$(top_builddir)/interface/src/python; \
	export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(top_builddir)/src/.libs;
LOG_COMPILER = $(PYTHON)
endif

CLEANFILES = *.vtk *.dx *.pyc tank_3D* tripod* plate* *.pos *.dx results/* __pycache__/*
