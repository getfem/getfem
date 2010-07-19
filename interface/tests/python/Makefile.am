EXTRA_DIST= 				\
	demo_tripod.py 			\
	demo_tripod_alt.py 		\
	demo_plate.py 			\
	check_export.py 		\
	check_levelset.py		\
	getfem_tvtk.py 			\
	demo_laplacian.py 		\
	demo_plasticity.py              \
	demo_wave.py			\
	check_global_functions.py	\
	demo_crack.py 			\
	demo_fictitious_domains.py 	\
	demo_step_by_step.py		\
	demo_stokes_3D_tank.py          \
	demo_stokes_3D_tank_draw.py	\
	getfem_tvtk.py

if BUILDPYTHON
TESTS = 						\
	$(srcdir)/check_export.py			\
	$(srcdir)/check_global_functions.py		\
	$(srcdir)/demo_wave.py				\
	$(srcdir)/demo_laplacian.py			\
	$(srcdir)/check_levelset.py

TESTS_ENVIRONMENT = \
	PYTHONPATH=$(srcdir):$(srcdir)/../../src/python:.:../../src/python \
	python
endif

CLEANFILES = *.vtk *.dx *.pyc tank_3D* tripod* plate* *.pos *.dx
