.. $Id: libdesc.rst 4917 2015-03-27 14:37:59Z renard $

.. include:: ../replaces.txt

.. highlight:: c++

.. _dp-libdesc_interface:


Interface with scripting languages (Python, Scilab and Matlab)
--------------------------------------------------------------

A simplified (but rather complete) interface of |gf| is provided, so that it is possible to use getfem in some script languages.

Description
^^^^^^^^^^^

All sources are located in the :file:`interface/src` directory. The interface is
composed of one large library ``getfemint`` (which stands for getfem
interaction), which acts as a layer above the |gf| library, and is used by
the python, matlab and scilab interfaces.

This interface is not something that is generated automatically from c++ sources
(as that could be the case with tools such as swig). It is something that has
been designed as a simplified and consistent interface to getfem. Adding a new
language should be quite easy (assuming the language provides some structures
for dense arrays manipulations).

Files
^^^^^

All the files in the directory :file:`interface\src`. A short description of main files:

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

* :file:`getfemint_gsparse.h`, :file:`getfemint_gprecond.h`,
  :file:`getfemint_gmumps.h`, etc.

  Files specific to an interfaced object if needed.
  (getfemint_gsparse which export some kind of mutable sparse matrix that can
  switch between different storage types, and real of complex elements).

* :file:`gf_workspace.cc`, :file:`gf_delete.cc`.

  Memory management for getfem objects. There is a layer which handles the
  dependency between for example a ``mesh`` and a ``mesh_fem``.
  It makes sure that no object
  will be destroyed while there is still another getfem_object using it.
  The goal
  is to make sure that under no circumstances the user is able to crash getfem
  (and the host program, matlab, scilab or python) by passing incorrect
  argument to the getfem interface.

  It also provides a kind of workspace stack, which was designed to simplify
  handling and cleaning of many getfem objects in matlab (since matlab does not
  have "object destructors").

* :file:`getfemint.h`, :file:`getfemint.cc`.

  Define the ``mexarg_in``, ``mexarg_out`` classes, which are used to parse the
  list of input and output arguments to the getfem interface functions.
  The name  is not adequate anymore since any reference to "mex"
  has been moved into
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

A special effort has been done to facilitate the addition of new objects, methods and functions to the interface without doing it separately for each part supported script language (Python, Scilab, Octave, Matlab).


All the information needed to build the interface for the different objects, methods and functions is contained in the files `interface/src/gf*.cc`. A python script (`bin/extract_doc`) produces all the necessary files from the information it takes there. In particular, it produces the python file getfem.py, the matlab m-files for the different functions and objects (including subdirectories) and it also produces the automatic documentations.

To make all the things work automatically, a certain number of rules have to be respected:


* An object have to be defined by three files on the interface

  - :file:`gf_objectname.cc` : contains the constructors of the object

  - :file:`gf_objectname_get.cc` : contains the methods which only get some information about the object (if any).

  - :file:`gf_objectname_set.cc` : contains the methods which transform the object (if any).

* A list of function is defined by only one file :file:`gf_commandname.cc`.
  It contains a list of sub-commands.


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

  The first three lines are a C++ comment which describes the call of the
  method/function with a special syntax and also gives a description of the
  method/function which will be included in the documentations. The first
  line of this comment is important since it will be analyzed to produce
  the right interface for Python, Octave, Matlab and Scilab.

  The syntax for the description of the call of a method/function is the
  following: After ``/*@`` a special keyword should be present. It is either
  ``INIT``, ``GET``, ``SET``, ``RDATTR`` or ``FUNC``. The keyword
  ``INIT`` means that
  this is the description of a constructor of an object. ``RDATTR`` is for
  a short method allowing to get an attribute of an object. ``GET`` is for a
  method of an object which does not modify it. ``SET`` is for a method which
  modifies an object and ``FUNC`` is for the sub-command of a function list.

  If the method/function returns a value, then a name for the return value
  is given (which is arbitrary) followed by ``=``.

  The parameters of the method/function are described. For a method, the
  object itself is not mentionned. The first parameter should be the method
  or sub-command name between single quotes (a special case is when
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
  For instance, you can refer to ``@tmesh``, ``@tmesh_fem``, ``@tfem``, etc.
  There are some authorized abbreviations:

        - ``@tcs``  for  ``@tcont_struct``
        - ``@tmf``  for  ``@tmesh_fem``
        - ``@tgt``  for  ``@tgeotrans``
        - ``@tgf``  for  ``@tglobal_function``
        - ``@tmo``  for  ``@tmesher_object``
        - ``@tmls`` for  ``@tmesh_levelset``
        - ``@tmim`` for  ``@tmesh_im``
        - ``@tls``  for  ``@tlevelset``
        - ``@tsl``  for  ``@tslice``
        - ``@tsp``  for  ``@tspmat``
        - ``@tpre`` for  ``@tprecond``
        - ``@tmct`` for  ``@tmumps_context``


  Three dots at the end of the parameter list (``...``) mean that
  additional parameters are possible. Optional parameters can be described
  with brackets. For instance ``/*@SET v = ('name'[, @int i])``. But
  be careful how it is interpreted by the :file:`extract_doc` script
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
  the right syntax depending on the language (Octave, Matlab, Scilab or Python).

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
passed to the function in order to make sure that the user cannot crash Scilab, Octave, Matlab or Python when using that function. Use the macro defined in :file:`gf_mesh_get.cc` to add your function.

Do not forget to add documentation for that function: in :file:`gf_mesh_get.cc`,
this is the documentation that appears in the Octave/Matlab/Scilab/Python help files (that is when on
type "``help gf_mesh_get``" at the matlab prompt), and in the getfem_python
autogenerated documentation.

IMPORTANT. Note that the array indices start at 0 in Python and 1 in Octave, Matlab and Scilab. A specific function::

   config::base_index()

whose value is 0 in python and 1 in Octave, Matlab and Scilab has to be used to exchange indices and array of indices. Take care not to make the correction twice. Some Array of indices are automatically shifted.

Adding a new object to the getfem interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to add a new object to the interface, you have to build the new corresponding sources :file:`gf_obj.cc`, :file:`gf_obj_get.cc` and :file:`gf_obj_set.cc`. Of course you can take the existing ones as a model.

For the management of the object, you have to declare the class at the begining of :file:`getfemint.h` (with respect to the alphabetic order), and declare three functions::

  bool is_"name"_object(const mexarg_in &p);
  id_type store_"name"_object(const std::shared_ptr<object_class> &shp);
  object_class *to_"name"_object(const mexarg_in &p);

where "name" is the name of the object in the interface and ``object_class`` is the class name in getfem (for instance  ``getfem::mesh`` for the mesh object). Alternatively, for the object that are manipulated by a shared pointer in |gf|, the third function can return a shared pointer.

IMPORTANT: In order to be interfaced, a |gf| object has to derive from ``dal::static_stored_object``. However, if it is not the case, a wrapper class can be defined such as the one for ``bgeot::base_poly`` (see the end of :file:`getfemint.h`).

The previous three functions have to be implemented at the end of :file:`getfemint.cc`.It is possible to use one of the two macros defined in :file:`getfemint.cc`. The first macro is for a standard object and the second one for an object which is manipulated in |gf| with a shared pointer.

You have also to complete functions ``name_of_getfemint_class_id`` and ``class_id_of_object`` at the end of :file:`getfemint.cc`.


You have to add the call of the interface function in :file:`getfem_interface.cc` and modifiy the file :file:`bin/extract_doc` and run the configure file.

The methods ``get('char')`` and ``get('display')`` should be defined for each object. The first one should give a string allowing the object to be saved in a file and the second one is to give some information about the object. Additionaly, a constructor from a string is necessary to load the object from a file.

For the Scilab interface the file :file:`sci_gateway/c/builder_gateway_c.sce.in` has to be modified and the files in the directory :file:`macros/overload`.


State
^^^^^


Perspectives
^^^^^^^^^^^^
The interface grows in conjunction with |gf|. The main |gf| functionalities are interfaced.

