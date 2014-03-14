EXTRA_DIST= 				\
	check_export.py 		\
	check_global_functions.py	\
	check_levelset.py		\
	demo_crack.py 			\
	demo_fictitious_domains.py 	\
	demo_laplacian.py 		\
	demo_mortar.py 			\
	demo_plasticity.py              \
	demo_plate.py 			\
	demo_static_contact.py 		\
	demo_step_by_step.py		\
	demo_stokes_3D_tank.py          \
	demo_stokes_3D_tank_draw.py	\
	demo_tripod.py 			\
	demo_tripod_alt.py 		\
	demo_wave.py			\
	getfem_tvtk.py

if BUILDPYTHON
TESTS = 				\
	check_export.py			\
	check_global_functions.py	\
	demo_wave.py			\
	demo_laplacian.py		\
	check_levelset.py

## to be replaced with AM_TESTS_ENVIRONMENT in automake versions >= 1.13
TESTS_ENVIRONMENT = \
	PYTHONPATH='$(srcdir):$(srcdir)/../../src/python:.:../../src/python'; \
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(top_builddir)/src/.libs; \
	export PYTHONPATH; \
	export LD_LIBRARY_PATH;
LOG_COMPILER = python
endif

CLEANFILES = *.vtk *.dx *.pyc tank_3D* tripod* plate* *.pos *.dx
