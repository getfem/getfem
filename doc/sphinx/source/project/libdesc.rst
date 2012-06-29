.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _dp-libdesc:

Description of the different parts of the library
=================================================

gmm library
-----------

Description
^^^^^^^^^^^

|gmm| is a linear algebra library which was originally designed to make an
interface between the need in linear algebra procedures of |gf| and existing free
linear algebra libraries (MTL, Superlu, Blas, Lapack originally). It rapidly
evolves to an independent self-consistent library with its own vector and matrix
types. It is now used as a base linear algebra library by several other projects
(projet `KDE <http://websvn.kde.org/trunk/kdesupport/gmm>`_, for instance).

However, it preserves the characteristic to be a potential interface for more
specific packages. Any vector or matrix type having the minimum of compatibility
can be used by generic algorithms of |gmm| writing a ``linalg_traits`` structure.

A |gmm| standalone version is distributed since release 1.5 of |gf|. It is
however developed inside the |gf| project even though since release 3.0 it is
completely independent of any |gf| file.

In addition to the linear algebra procedures, it furnishes also the following
utilities to |gf|.

* Fix some eventual compatibility problems in :file:`gmm_std.h`.

* Error, warning and trace management in :file:`gmm_except.h`.

* Some extended math definitions in :file:`gmm_def.h`.

State
^^^^^

For the moment, |gmm| cover the needs of |gf| concerning the basic linear algebra
procedures.

Perspectives
^^^^^^^^^^^^

There is potentatialy several points to be improved in |gmm| (partial
introduction of expression template for some base types of matrix and vectors,
reflection on the way to represent in a more coherent manner sparse sub-vectors
and sub-matrices, introduction of C++ concepts, etc.). However, since |gmm|
globally cover the needs of |gf| and since there exists some other project like
`Glas <http://glas.sourceforge.net/>`_ to build a reference C++ library for
linear algebra, a global reflection seems not necessary for the moment. This part
is considered to be stabilized.

The current vocation of |gmm| is to continue to collect generic algorithms and
interfaces to some other packages in order to cover new needs of the whole
project. The library is now frequently used as a separate package and has also
the vocation to collect the contribution of any person who propose some
improvements, new algorithms or new interfaces.


MESH module
-----------

Description
^^^^^^^^^^^

This part of the library has the role to store and manage the meshes, i.e. a
collection of elements (real elements) connected to each other by some of their
faces. For that, it develops concepts of elements, elements of reference,
structure of meshes, collection of nodes, geometric transformations, subpart of
the boundary or subzone of the mesh.

There is no really effective meshing capabilities available for the moment in
|gf|. The meshes of complex objects must be imported from existing meshers such
as `Gmsh`_ or `GiD`_. Some importing functions of meshes have been written and
can be easily extended for other formats.

The object which represents a mesh declared in the file :file:`getfem_mesh.h` and
which is used as a basis for handling of the meshes in |gf| manages also the
possibility for the structures depending on a mesh (see MESHFEM and MESHIM
modules) to react to the evolution of the mesh (addition or removal of elements,
etc.).

State
^^^^^

The main C++ header files are

* :file:`bgeot_convex_structure.h`

  Describes the structure of an element disregarding the coordinates of its 
  vertices.

* :file:`bgeot_mesh_structure.h`

  Describes the structure of a mesh disregarding the coordinates of the nodes.

* :file:`bgeot_node_tab.h`

  A node container allowing the fast search of a node.

* :file:`bgeot_convex.h`

  Describes an element with its vertices.

* :file:`bgeot_convex_ref.h`

  Describe reference elements.

* :file:`bgeot_mesh.h`

  Describes a mesh with the collection of node (but without the description of 
  geometric transformations).

* :file:`bgeot_geometric_trans.h`

  Describes geometric transformations.

* :file:`bgeot_geotrans_inv.h`

  A tool to invert geometric transformations.

* :file:`getfem_mesh.h`

  Fully describes a mesh (with the geometric transformations, subparts of the 
  mesh, support for parallelization). Includes the Bank algorithm to refine a 
  mesh.

* :file:`getfem_mesher.h`

  An attempt to develop a mesher. To be use with care.

A prototype of mesher is in the files :file:`getfem_mesher.h` and
:file:`getfem_mesher.cc` which makes it possible to mesh geometries defined by
some level sets. However, the continuation of the development of this mesher is
not planned for the moment because the project |gf| has vocation to focus on the
finite element methods themselves.

Perspectives
^^^^^^^^^^^^

For the moment, the module is split into two parts which lie into two different
namespaces. Of course, It would be more coherent to gather the module in only one
namespace (``getfem``).

.. note::

   The file :file:`bgeot_mesh.h` could be renamed :file:`getfem_basic_mesh.h`.

A possible work to do on this part would be to examine the manner of storing the
meshes and possibly to make a bibliographical study on the manner of storing a
mesh (for instance see [remacle2002]_). It would be necessary to supplement
documentation and to examine also the management of the events and the way in
which the structures which depend on the mesh react to these events.


FEM module
----------

Description
^^^^^^^^^^^

The FEM module is the part of |gf| which describes the finite elements at the
element level and the degrees of freedom. Finite element methods can be of
different types. They could be scalar or vectorial, polynomial, piecewise
polynomial or non-polynomial, equivalent via the geometric transformation or not.
Moreover, the description of the degrees of freedom have to be such that it is
possible to gather the compatible degrees of freedom between two neighbor
elements in a generic way (for instance connecting a Lagrange 2D element to
another Lagrange 1D element).

State
^^^^^

The main files of the module are

* :file:`getfem_fem.h`

  Abstract definition of a finite element and a degree of freedom. Interface for 
  the exported functions of :file:`getfem_fem.cc` and 
  :file:`getfem_fem_composite.cc`.

* :file:`getfem_fem.cc`

  Definition of the polynomial finite elements and interface to get the 
  descriptor on these elements (function ``pfem fem_descriptor(std::string 
  name)``).

* :file:`getfem_fem_composite.cc`

  Definition of the piecewise polynomial finite elements.

The two files :file:`getfem_fem.cc` and :file:`getfem_fem_composite.cc` mainly
contains all the finite element description for basic elements. A exhaustive list
of the defined finite elements is given in :ref:`ud-appendixa`.

Some other files define some specific finite element such as
:file:`getfem_fem_level_set.h` which is a complex construction which allows to
"cut" a existing element by one or several level sets (see the LEVELSET module).

The manner to describe the degrees of freedom globally satisfies the needing
(connecting dof from an element to another in a generic way) but is a little bit
obscure and too much complicated.

Conversely, the way to represent non-equivalent elements with the supplementary
matrix ``M`` has proven its efficiency on several elements (Hermites elements,
Argyris, etc.).

Perspectives
^^^^^^^^^^^^

The principal dissatisfaction of this module is that description of the degrees
of freedom is not completely satisfactory. It is the principal reason why one
documentation on how to build an element from A to Z was not made for the moment
because description of the degrees of freedom was conceived to be temporary. An
effort of design is thus to be provided to completely stabilize this module
mainly thus with regard to the description of degrees of freedom but also perhaps
the description of finite elements which could be partially externalized in a
similar way to the cubature methods , at least for the simplest finite elements
(equivalent and polynomial finite elements).


CUBATURE module
---------------

Description
^^^^^^^^^^^

The CUBATURE module gives access to the numerical integration methods on
reference elements. In fact it does not only contain some cubature formulas
because it also give access to some exact integration methods. However, the exact
integration methods are only usable for polynomial element and affine geometric
transformations. This explain why exact integration methods are not widely used.
The description of cubature formulas is done either directly in the file
:file:`getfem_integration.h` or via a description file in the directory
``cubature`` of |gf|. The addition of new cubature formulas is then very simple,
it suffices to reference the element on which it is defined and the list of Gauss
points in a file and add it to this directory. Additionally, In order to
integrate terms defined on a boundary of a domain, the description should also
contains the reference to a method of same order on each face of the element.

State 
^^^^^

This module meets the present needs for the project and is considered as
stabilized. The list of available cubature formulas is given in 
:ref:`ud-appendixb`.

Perspectives 
^^^^^^^^^^^^

No change needed for the moment. An effort could be done on the documentation to
describe completely how to add a new cubature formula (format off descritption
files).


MESHFEM module
--------------

to be done

Description
^^^^^^^^^^^

State
^^^^^

Perspectives
^^^^^^^^^^^^

Parallelisation of dof numbering to be done. An optimal (an simple) algorithm
exits.


LEVELSET module
^^^^^^^^^^^^^^^

to be done

Description
^^^^^^^^^^^

State
^^^^^

Perspectives
^^^^^^^^^^^^


MESHIM module
-------------

to be done

Description
^^^^^^^^^^^

State
^^^^^

Perspectives
^^^^^^^^^^^^


INTEGELEM module
----------------

to be done

Description
^^^^^^^^^^^

State
^^^^^

Perspectives
^^^^^^^^^^^^


ASSEMBLE module
---------------

to be done

Description
^^^^^^^^^^^

State
^^^^^

Perspectives
^^^^^^^^^^^^


BRICK module
------------

to be done

Description
^^^^^^^^^^^

State
^^^^^

Perspectives
^^^^^^^^^^^^


Events management
-----------------

Description
^^^^^^^^^^^

The ``mesh``, |mf|, |mim| and |mo| description are linkedtogether in the sense
that there is some dependencies between them. For instance, when an element is
suppressed to a mesh, the |mf| object has to react.

State
^^^^^

The main tool to deal with simple dependence of object is in
:file:`getfem_context.h`. An object ``context_dependencies`` is defined there. In
order to deal with the dependencies of an object, the object
``context_dependencies`` needs to be a parent class of this object. It adds the
following methods to the object:

.. cfunction:: add_dependency(ct)

   Add an object (which has to have ``context_dependencies`` as a parent class)
   to the list of objects from which the current object depend.

.. cfunction:: touch()

   Indicates to the dependent objects that something has change in the object.

.. cfunction:: context_check()

   Check if the object has to be updated. if it is the case it makes first a
   check to the dependency list and call the update function of the object. (the
   update function of the dependencies are called before the update function of
   the current object).

.. cfunction:: context_valid()

   Says if the object has still a valid context, i.e. if the object in the
   dependency list still exist.

Moreover, the object has to define a method::

 ``void update_from_context(void) const``

which is called after a ``context_check()`` if the context has changed.

An additional system is present in the object |m|. Each individual element has a
version number in order for the objects |mf| and |mim| to detect which element
has changed between two calls.

Perspectives
^^^^^^^^^^^^

Some object do not manage satisfactorily events. This is the case for instance of
|mls|, |mfls|, |pmf|, etc.

This is clear that the event management still have to be tested and improved to
have a fully reactive system.


Python, Scilab and Matlab interfaces
------------------------------------

A simplified interface of |gf| is provided, so that it is possible to use getfem
in other languages.

Description
^^^^^^^^^^^
 
All sources are located in the :file:`interface/src` directory. The interface is
composed of one large library ``getfemint`` (which stands for getfem
interaction), which acts as a layer above the |gf| library, and is used by
the python, matlab and scilab interfaces.

This interface is not something that is generated automatically from c++ sources
(as that could be the case with tools such as swig). It is something that has
been designed as a simplified and consistent interface to getfem. Adding a new
language should be quite easy (assuming the language provides some structures for
dense arrays manipulations).

State
^^^^^

Here is a list of the various files, with a short description:

* :file:`getfem_interface.cc`.

  This is the bridge between the script language and the getfem interface. The 
  function getfem_interface_main is exported as an ``extern "C"`` function, so 
  this is a sort of c++ barrier between the script language and the getfem 
  interface (exporting only a C interface avoids many compilation problems).

* :file:`matlab/gfm_mex.c`.

  The matlab interface. The only thing it knows about getfem is in 
  :file:`getfem_interface.h`.

* :file:`python/getfem_python.c`.

  The python interface. The only thing it knows about getfem is in 
  :file:`getfem_interface.h`.

* :file:`gfi_array.h`, :file:`gfi_array.c`.

  Both :file:`gfm_mex.c` and :file:`getfem_python.c` need a simple convention on 
  how to send and receive arrays, and object handles, from 
  ``getfem_interface_main()``. This file provide such functionnality.

* :file:`getfemint_object.h`.

  Not all getfem objects are exported, only a selected subset, mostly |m|, |mim|, 
  |mf|, |sl|, |br|, etc. They are all wrapped in a common interface, which is 
  ``getfemint::getfem_object``.

* :file:`getfemint_mesh.h`, :file:`getfemint_mesh_fem.h`, etc.

  All the wrapped |gf| objects. Some of them are quite complicated 
  (getfemint_gsparse which export some kind of mutable sparse matrix that can 
  switch between different storage types, and real of complex elements).

* :file:`gf_workspace.cc`, :file:`gf_delete.cc`.

  Memory management for getfem objects. There is a layer in 
  ``getfemint::getfem_object`` which handles the dependency between for example a 
  ``getfemint_mesh`` and a ``getfemint_mesh_fem``. It makes sure that no object 
  will be destroyed while there is still another getfem_object using it. The goal 
  is to make sure that under no circumstances the user is able to crash getfem 
  (and the host program, matlab, scilab or python) by passing incorrect argument to the 
  getfem interface.

  It also provides a kind of workspace stack, which was designed to simplify 
  handling and cleaning of many getfem objects in matlab (since matlab does not 
  have "object destructors").

* :file:`getfemint.h`, :file:`getfemint.cc`.

  Define the ``mexarg_in``, ``mexarg_out`` classes, which are used to parse the 
  list of input and output arguments to the getfem interface functions. The name 
  is not adequate anymore since any reference to "mex" has been moved into 
  :file:`gfm_mex.c`.

* :file:`gf_mesh.cc`, :file:`gf_mesh_get.cc`, :file:`gf_mesh_set.cc`,
  :file:`gf_fem.cc`, etc.

  All the functions exported be the getfem interfaces, sorted by object type 
  (``gf_mesh*``, ``gf_mesh_fem*``, ``gf_fem*``), and then organized as one for 
  the object construction (``gf_mesh``), one for the object modification 
  (``gf_mesh_set``), and one for the object inquiry (``gf_mesh_get``). Each of 
  these files contain one main function, that receives a ``mexargs_in`` and 
  ``mexargs_out`` stack of arguments. It parses then, and usually interprets the 
  first argument as the name of a subfunction (``gf_mesh_get('nbpts')`` in 
  matlab, or ``Mesh.nbpts()`` in python).

* :file:`matlab/gfm_rpx_mexint.c`.

  An alternative to :file:`gfm_mex.c` which is used when the 
  ``--enable-matlab-rpc`` is passed to the ``./configure`` script. The main use 
  for that is debugging the interface, since in that case, the matlab interface 
  communicates via sockets with a "getfem_server" program, so it is possible to 
  debug that server program, and identify memory leaks or anything else without 
  having to mess with matlab (it is pain to debug).

* :file:`python/getfem.py`.

  The python interface is available as a ":file:`getfem.py`" file which is
  produced during compilation by the python script
  ":file:`bin/extract_doc.py`".



Objects, methods and functions of the interface 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The main concepts manipulated by the interface are a limited number of objects
(Fem, Mesh, MeshFem, Model ...), the associated methods and some functions defined on these objects.

A special effort has been done to facilitate the addition of new objects, methods and functions to the interface without doing it separetaly for each partsupported script language (Python, Scilab, Matlab).


All the information needed to build the interface for the different objects, methods and functions is contained in the files `interface/src/gf*.cc`. A python script (`bin/extract_doc`) produces all the necessary files from the information it takes there. In particular, it produces the python file getfem.py, the matlab m-files for the different functions and objects (including subdirectories) and it also produces the automatic documentations.

To make all the things work automatically, a certain number of rules have to be respected:


* An object have to be defined by three files on the interface

  - :file:`gf_objectname.cc` : contains the constructors of the object

  - :file:`gf_objectname_get.cc` : contains the methods which only get some information about the object (if any).

  - :file:`gf_objectname_set.cc` : contains the methods which transform the object (if any).

* A list of function is defined by only one file :file:`gf_commandname.cc`
  it contains a list of sub-comands.


* For each file, the main commentary on the list of functions or methods is delimited by the tags '/*@GFDOC' and '@*/'. For a file corresponding to the constructors of an object, the commentary should correspond to the description of the object.


* Each non trivial file gf_*.cc contains a macro allowing to define the
  methods of the object or the sub-commands. In particular, this system
  allows to have a efficient search of the called method/function.
  This macro allows to declare
  a new method/function with the following syntax::

   /*@GET val = ('method-name', params, ...)
      Documention of the method/function.
   @*/
   sub_command
   ("method-name", 0, 0, 0, 1,
     ...
     body of the method/function
     ...
   );

  The first three line are a c++ commentary which describes the call of the
  method/function with a special syntax and also gives a description of the
  method/function which will be included in the documentations. The first
  line of this commentary is important since it will be analyzed to produce
  the right interface for Python, Matlab and Scilab.

  The syntax for the description of the call of a method/function is the
  following: After ``/*@`` a special keyword should be present. It is either
  ``INIT``, ``GET``, ``SET``, ``RDATTR`` or ``FUNC``. The keyword
  ``INIT`` means that
  this is the description of a constructor of an object. ``RDATTR`` is for
  a short method allowing to get an attribut of an object. ``GET`` is for a
  method of an object which does not modify it. ``SET`` is for a method which
  modifies an object and ``FUNC`` is for the sub-command of a function list.

  If the method/function returns a value, then a name for the return value
  is given (which is arbitrary) followed by ``=``.

  The parameters of the method/function are described. For a method, the
  object itself is not mentionned. The first parameter should be the method
  or sub-command name between single quotes (a speical case is when
  this name begins with a dot; this means that it corresponds to a
  method/function where the command name is not required).

  The other parameters, if any, should be declared with a type. Predefined
  types are the following:

        - ``@CELL``   : a cell array,
        - ``@imat``   : matrix of integers,
        - ``@ivec``   : vector of integers,
        - ``@cvec``   : vector of complex values,
        - ``@dcvec``  : vector of complex values,
        - ``@dvec``   : vector of real values,
        - ``@vec``    : vector of real or complex values,
        - ``@dmat``   : matrix of real values,
        - ``@mat``    : matrix of real or complex values,
        - ``@str``    : a string,
        - ``@int``    : an integer,
        - ``@bool``   : a boolean,
        - ``@real``   : a real value,
        - ``@scalar`` : a real or complex value,
        - ``@list``   : a list.

  Moreover, ``@tobj`` refers to an object defined by the interface.
  For instance, ou can refer to ``@tmesh``, ``@tmesh_fem``, ``@tfem``, etc.
  There are some authorized abreviations:

        - ``@tmf``  for  ``@tmesh_fem``
        - ``@tbrick``  for  ``@tmdbrick``
        - ``@tstate``  for  ``@tmdstate``
        - ``@tgt``  for  ``@tgeotrans``
        - ``@tgf``  for  ``@tglobal_function``
	- ``@tmo``  for  ``@tmesher_object``
        - ``@tmls``  for  ``@tmesh_levelset``
	- ``@tmim``  for  ``@tmesh_im``
        - ``@tls``  for  ``@tlevelset``
        - ``@tsl``  for  ``@tslice``
        - ``@tsp``  for  ``@tspmat``
        - ``@tpre``  for  ``@tprecond``


  Three dots at the end of the parameter list (``...``) mean that 
  additional parameters are possible. Optional parameters can be described
  with brackets. For instance ``/*@SET v = ('name'[, @int i])``. But
  be carreful how it is interpreted by the :file:`extract_doc` script
  to build the python interface.

  The second to fifth parameters of the macro correspond respectively to
  the minimum number of input arguments, the maximum one, the minimum
  number of output arguments and the maximum number of output arguments. It
  is dynamically verified.

  Additional parameters for the function lists ....

  For unknown reasons, the body of the function cannot contain multiple
  declarations such as ``int a, b;`` (c++ believes that it is an additional
  parameter of the macro).

.. _reStructuredText: http://docutils.sourceforge.net/rst.html

* The parts of documentation included in the c++ commentaries should be in 
  `reStructuredText`_ format. In particular, math formulas can be included
  with \:math\:\`f(x) = 3x^2+2x+4\` or with::
  
    .. math::
 
      f(x) = 3x^2+2x+4

  It is possible to refer to another method or function of the interface
  with the syntax ``INIT::OBJNAME('method-name', ...)``,
  ``GET::OBJNAME('method-name', ...)``, ``SET::OBJNAME('method-name', ...)``,
  ``FUNC::FUNCNAME('subcommand-name', ...)``. This will be replaced with
  the right syntax depending on the language (Matlab, Scilab or Python).

* Still in the documentations, parts for a specific language can be added by
  ``@MATLAB{specific part ...}``, ``@SCILAB{specific part ...}`` and
  ``@PYTHON{specific part ...}``.
  If a method/sub-command is specific to an interface, it can be added,
  for instance for Matlab,
  replacing `GET` by `MATLABGET`, `FUNC` by `MATLABFUNC`, etc.
  If a specific code is needed for this additional function, it can be added
  with the tags ``/*@MATLABEXT``, ``/*@SCILABEXT``, ``/*@PYTHONEXT``. See
  for instance the file :file:`gf_mesh_fem_get.cc`.

* For Python and the Matlab object, if a `SET` method has the same name as
  a `GET` method, the `SET` method is prefixed by `set_`.







Adding a new function or object method to the getfem interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If one want to add a new function ``gf_mesh_get(m, "foobar", .)``, then the
main file to modify is :file:`gf_mesh_get.cc`. Remember to check every argument
passed to the function in order to make sure that the user cannot crash scilab, matlab or python when using that function. Use the macro defined in :file:`gf_mesh_get.cc` to add your function.

Do not forget to add documentation for that function: in :file:`gf_mesh_get.cc`,
this is the documentation that appears in the matlab/scilab/python help files (that is when on
type "``help gf_mesh_get``" at the matlab prompt), and in the getfem_python
autogenerated documentation.

IMPORTANT. Note that the array indices start at 0 in Python and 1 in Matlab and Scilab. A specific function::

   config::base_index()

whose value is 0 in python and 1 in Matlab and Scilab has to be used to exchange indices and array of indices. Take care not to make the correction twice. Some Array of indices are automatically shifted.

Adding a new object to the getfem interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to add a new object to the interface, you have to build the new corresponding sources :file:`gf_obj.cc`, :file:`gf_obj_get.cc` and :file:`gf_obj_set.cc`. Of course you can take the existing ones as a model.

A structure name `getfemint_object_name` has to be defined (see getfemint_mesh.h for instance).
Moreover, for the management of the object, you have to declare the class in :file:`getfemint.cc` and :file:`getfemint.h` and add the methods `is_object()`, `to_const_object()`, `to_object()` and `to_getfemint_object()`. You have to set its ``class_id`` in :file:`gfi_array.h` (with respect to the alphabetic order of its name).

You have also to add the call of the interface function in :file:`getfem_interface.cc` and modifiy the file :file:`bin/extract_doc` and run the configure file.

The methods ``get('char')`` and ``get('display')`` should be defined for each object. The first one should give a string allowing the object to be saved in a file and the second one is to give some information about the object. Additionnaly, a constructor from a string is necessary to load the object from a file.

For the Scilab interface the file :file:`sci_gateway/c/builder_gateway_c.sce.in` has to be modified and the files in the directory :file:`macros/overload`.

Perspectives
^^^^^^^^^^^^
The interface grows in conjunction with |gf|. The objective is to interface the maximum of the |gf| functionalities.

