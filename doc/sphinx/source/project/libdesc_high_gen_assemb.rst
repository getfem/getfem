.. $Id: libdesc.rst 4917 2015-03-27 14:37:59Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _dp-libdesc_high_assemb:



The high-level generic assembly module in |gf|
-----------------------------------------------


Description
^^^^^^^^^^^

The high level generic assembly language of |gf| is a key module which allows to describe weak formulation of partial differential equation problems.

Files
^^^^^

.. csv-table::
   :header: "File(s)", "Description"
   :widths: 8, 15

   :file:`getfem_generic_assembly.h` and :file:`getfem_generic_assembly.cc`, "In order not to export all the internal of the generic assembly, all is implemented in a single file"

A few implementation details
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The assembly string is transformed in an assembly tree by a set of function in :file:`src/getfem\_generic\_assembly.cc`. The process has 4 steps:

 - Lexical analysis with the procedure ``ga_get_token``.

 - Syntax analysis and transformation into a syntax tree by ``ga_read_string``.

 - Semantic analysis, simplification (pre-computation) of constant expressions and enrichment of the tree.

 - Symbolic (automatic) differentiation.

 - Compilation in a sequence of instructions with optimization not to evaluate several time the same expression.

 - Execution of the sequence of instructions and assembly.

These steps are performed only once at the beginning of the assembly. The final tree is evaluated on each Gauss point of each element after compilation.




Predefined functions
^^^^^^^^^^^^^^^^^^^^

Some predefined scalar functions are available in |gf| generic assembly langage in order to describe a weak formulation (or also to make basic algebraic computations). This is limited to scalar functions of one or two arguments. Due to the automatic differentiation used to obtain the tangent system of described problems, the derivative each function have to be available. The principle retained is the following: For each predefined function is available:
  - A C++ function which computes the value given the argument(s).
  - The support of the function in the first each argument in term of a
    (possibly infinite) interval (this is for simplification of expressions).
  - The string corresponding of the derivative in terms of already known
    functions

A new predefined function is easy to add. See init_predefined_functions() in file :file:`src/getfem_generic_assembly.cc`. + describe how to give the derivative ...



State
^^^^^
Stable.

Perspectives
^^^^^^^^^^^^

- Is a certain extension to complex data possible ?

- More simplifications

- Integation of a tool allowing to compute inter-elements quantities (for a posteriori estimate for instance).
