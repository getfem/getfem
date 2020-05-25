.. _documenting-index:

###############
  Documenting
###############


The GetFEM library has a substantial body of documentation, much of it
contributed by various authors. The markup used for the GetFEM documentation
is `reStructuredText`_, developed by the `docutils`_ project, amended by custom
directives and using a toolset named `Sphinx`_ to postprocess the HTML output.

This document describes the style guide for our documentation, the custom
reStructuredText markup introduced to support Python documentation and how it
should be used, as well as the Sphinx build system.

.. _reStructuredText: http://docutils.sourceforge.net/rst.html
.. _docutils: http://docutils.sourceforge.net/
.. _Sphinx: http://sphinx.pocoo.org/

If you're interested in contributing to GetFEM's documentation, there's no
need to write reStructuredText if you're not so inclined; plain text
contributions are more than welcome as well. The main documentations are in the directory ``doc/sphinx/source`` of the project. A part of the documentation is automatic and comes from the sources of the project. This is in particular the case for the documentations of the interface commands which are located in the ``interface/src/gf_*.cc`` files.

It is highly recommending to document each created C++ class and exported function both in sources (for Oxygen documentation) and in the user documentation.

.. toctree::

   style.rst
   rest.rst
   markup.rst
   fromlatex.rst
