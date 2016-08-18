check_PROGRAMS = 

CLEANFILES = 

if BUILDPYTHON
TESTS = test_plasticity.py

AM_TESTS_ENVIRONMENT = \
	export PYTHONPATH=$(top_builddir)/interface/src/python; \
	export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(top_builddir)/src/.libs;
EXTRA_DIST = \
	test_plasticity.py \
	test_small_strain_plasticity.py

LOG_COMPILER = python
endif

