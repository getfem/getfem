EXTRA_DIST=demo_tripod.py demo_tripod_alt.py demo_plate.py \
	check_export.py getfem_tvtk.py demo_laplacian.py demo_wave.py

if BUILDPYTHON
TESTS = 						\
	$(srcdir)/check_export.py			\
	$(srcdir)/check_global_functions.py			\
	$(srcdir)/demo_wave.py				\
	$(srcdir)/demo_laplacian.py			\
	$(srcdir)/check_levelset.py

TESTS_ENVIRONMENT = \
	PYTHONPATH=$(srcdir):$(srcdir)/../../src/python:.:../../src/python \
	python
endif

clean: clean-am
	rm -f *.vtk
	rm -f *.dx
	rm -f *.pyc
	rm -f tank_3D*
	rm -f tripod*
	rm -f plate*
