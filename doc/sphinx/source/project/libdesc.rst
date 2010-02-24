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
interaction), which is acts as a layer above the |gf| library, and is used by
both the python and matlab interfaces.

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
  (and the host program, matlab or python) by passing incorrect argument to the 
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

* :file:`python/getfem.base.py`.

  The python interface is available as a ":file:`getfem.py`" file which is built 
  during compilation. Its source file is :file:`getfem.base.py`, it contains just 
  the list of classes, and for each class the names of the member functions.



Objects, methods and functions of the interface 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


The main concepts manipulated by the interface are a limited number of objects
(Fem, Mesh, MeshFem, Model .), the associated methods and some functions defined on these objects.

A special effort has been done to facilitate the addition of new objects, methods and functions to the interface without doing it separetaly for each part (Python, Scilab, Matlab).

The following rules have to be respected :

* An object have to be defined by three files on the interface

  - :file:`gf_objectname.cc` : contains the constructors of the object

  - :file:`gf_objectname_get.cc` : contains the methods which only get some information about the object.

  - :file:`gf_objectname_set.cc` : contains the methods which transform the object.

* A list of function is defined by only one file :file:`gf_commandname.cc`



How to document an object, a method or a function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In each file, the commentary between tags /*@GFDOC and @*/ is added to the documentation to the object or command. Specific documentation for each sub-command should be placed between the tags /*@SET('sub-command-name', .) and @*/ for 
for the set methods, /*@GET('sub-command-name', .) and @*/ for the get methods and /*@FUNC('sub-command-name', .) and @*/ for the functions. (/*@RDATTR a changer en GET ?)

On fait reference a une autre fonction par
OBJ:INIT, OBJ:GET, OBJ:SET ou ::FUNCTION  

--> Attention, il faut definir le format des fonctions pour qu'il soit dechiffrable et le transmettre a l'interface python par example
+ type des parametres .
La "declaration de fonction" doit tenir sur une seule ligne a l'exclusion de tout autre chose

A l'interieur, globalement, c'est la syntaxe rst qui doit predominer. Les formules mathematiques sont mise sous la forme   :math:`f(x)`



En python, si une commande set a le meme nom qu'une commande get, il doit etre accole un 'set_' au nom de la commande set.
(en matlab ?)


On peut ajouter des tags specifiques pour ajouter a l'interface scilab/matlab/python des parties de doc ou de commande pour regler les cas penibles (eval, MLABEXT, .)  -> les tags specifiques sont @MATLAB{} @PYTHON{} et @SCILAB{}

On peut faire reference a une autre commande dans la doc en mettant
nom_objet:get (ou set) ou ::nom_commande en majuscule .
suivant le langae destination, la syntaxe doit etre adaptee pour mettre la bonne commande (voir ce qui se passe en matlab)

Liste des types d'arguments pre-definis :
        d = string.replace(d, '@CELL',   '')
        d = string.replace(d, '@imat',   'imat')
        d = string.replace(d, '@ivec',   'ivec')
        d = string.replace(d, '@cvec',   'vec')
        d = string.replace(d, '@dcvec',  'vec')
        d = string.replace(d, '@dvec',   'vec')
        d = string.replace(d, '@vec',    'vec')
        d = string.replace(d, '@dmat',   'mat')
        d = string.replace(d, '@mat',    'mat')
        d = string.replace(d, '@str',    'str')
        d = string.replace(d, '@int',    'int')
        d = string.replace(d, '@bool',   'bool')
        d = string.replace(d, '@real',   'real')
        d = string.replace(d, '@scalar', 'scalar')
        d = string.replace(d, '@list',   'list')
- @tobj : ou obj est un objet de l'interface.
- @ivec : vecteur d'entiers
- @dvec : vecteur de doubles
- @cdvec : vecteur de complexe double doubles (a verif)
- @dmat : matrice de doubles
- @str  : chaine de caracteres
- @int  : entier

      # Authorized abbreviations
        d = string.replace(d, '@tmf',    'mesh_fem')
        d = string.replace(d, '@tbrick', 'mdbrick')
        d = string.replace(d, '@tstate', 'mdstate')
        d = string.replace(d, '@tgt',    'geotrans')
        d = string.replace(d, '@tgf',    'global_function')
        d = string.replace(d, '@tmls',   'mesh_levelset')
        d = string.replace(d, '@tls',    'levelset')
        d = string.replace(d, '@tsl',    'slice')
        d = string.replace(d, '@tsp',    'spmat')
        d = string.replace(d, '@tpre',   'precond')

Supprimer les <Par> (ou les gerer ?)

transformer les MESHFEM en MESH_FEM, les MEHSLEVELSET aussi

ATTENTION aux virgules dans les declaration de variables : les macros n'apprecient pas.

Attention au format de declaration des arguments. Python, ne comprends que des choses simples sinon il capitule (*args). Les trucs simples sont f(i, j), f(i[, j]). Les f({i|j}) ou f(i[, j[, k]]) sont transformes en *args.

Adding a new function or object method to the getfem interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If one want to add a new function ``gf_mesh_get(m, "foobar", .)``, then the
main file to modify is :file:`gf_mesh_get.cc`. Remember to check every argument
passed to the function in order to make sure that the user cannot crash matlab or python when using that function.

Do not forget to add documentation for that function: in :file:`gf_mesh_get.cc`,
this is the documentation that appears in the matlab/scilab/python help files (that is when on
type "``help gf_mesh_get``" at the matlab prompt), and in the getfem_python
autogenerated documentation. In order to have "foobar" as a member function of
the python |py_m| class, it is necessary to add it in the ``getfem.base.py``
file. It is also necessary to add documentation in the
``interface/doc/getfemmatlab.tex``, which was at the beginning the only
documentation available. It is still very matlab oriented, and a little bit
redundant with the documentation embedded in :file:`gf_mesh_get.cc`.


+ fonction specifiques a une seule interface 
LANGAGEEXT et LANGAGEFUNC (ou GET, SET, INIT)
examples : gf_compute, gf_mesh_fem_get .

+ pour un example de code (attention a etre compatible tout langage) le :: puis le decalage comme en rst standard.

+ regle python pour deux methodes get et set ayant un nom identique ...->  set_

Adding a new class to the getfem interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ajout dans pseudo_funcions dans 'Makefile.am' de src


 + ajout de la fonction dans getfem_interface.cc 
(du style void gf_mesh_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out); )

 + get('char')    -> sauvegarde d'un objet
 + get('display') -> affichage d'un objet
 + constructeur 'from string'

REMAINING-TODOLIST : completer l'approche Objet de Matlab en ajoutant automatiquement les sous-commandes dans le fichier subsref.m (d'autant plus qu'avec l'analyse pour Python, on devrait avoir tout le necessaire). Il y avait un embrion non automatique pour gfMesh, gfMeshFem, gfMeshIm, gfSlice, gfSpmat (voir les fichiers dans anciennes versions). Avec la version actuelle, les get sont a peu pres instancies, il manque les set.

Perspectives
^^^^^^^^^^^^


