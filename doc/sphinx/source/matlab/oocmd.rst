.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: matlab

.. _mlab-oocmd:

|gfm| OO-commands
=================

The toolbox comes with a set of |Mlab| objects `mathworks-oo`_, (look at the
:file:`@gf*` sub-directories in the toolbox directory). These object are no more
than the getfem object handles, which are flagged by |mlab| as objects.

In order to use these objects, you have to call their constructors: ``gfMesh``,
``gfMeshFem``, ``gfGeoTrans``, ``gfFem``, ``gfInteg``.  These constructor just
call the corresponding |gfm| function (i.e.  ``gf_mesh``, ``gf_mesh_fem``, ...),
and convert the structure returned by these function into a |mlab| object. There
is also a ``gfObject`` function which converts any getfem handle into the
corresponding |mlab| object.

With such object, the most interesting feature is that you do not have to call
the "long" functions names ``gf_mesh_fem_get(obj,...)``,
``gf_slice_set(obj,...)`` etc., instead you just call the shorter
``get(obj,...)`` or ``set(obj,...)`` whatever the type of ``obj`` is.

A small number of "pseudo-properties" are also defined on these objects, for
example if ``m`` is a ``gfMesh`` object, you can use directly ``m.nbpts`` instead
of ``get(m, 'nbpts')``.

As an example::

  % classical creation of a mesh object
  >> m=gf_mesh('load', 'many_element.mesh_fem')
  m =
       id: 2
      cid: 0
  % conversion to a matlab object. the display function is overloaded for gfMesh.
  >> mm=gfMesh(m)
  gfMesh object ID=2 [11544 bytes], dim=3, nbpts=40, nbcvs=7
  % direct creation of a gfMesh object. Arguments are the same than those of gf_mesh
  >> m=gfMesh('load', 'many_element.mesh_fem')
  gfMesh object ID=3 [11544 bytes], dim=3, nbpts=40, nbcvs=7
  % get(m, 'pid_from_cvid') is redirected to gf_mesh_get(m,'pid from cvid')
  >> get(m, 'pid_from_cvid', 3)
  ans =
       8     9    11    15    17    16    18    10    12
  % m.nbpts is directly translated into gf_mesh_get(m,'nbpts')
  >> m.nbpts
  ans =
      40

  >> mf=gfMeshFem('load','many_element.mesh_fem')
  gfMeshFem object: ID=5 [1600 bytes], qdim=1, nbdof=99,
    linked gfMesh object: dim=3, nbpts=40, nbcvs=7
  >> mf.mesh
  gfMesh object ID=4 [11544 bytes], dim=3, nbpts=40, nbcvs=7
  % accessing the linked mesh object
  >> mf.mesh.nbpts
  ans =
      40
  >> get(mf.mesh, 'pid_from_cvid', 3)
  ans =
       8     9    11    15    17    16    18    10    12

  >> mf.nbdof
  ans =
      99

  % access to fem of convex 1
  >> mf.fem(2)
  gfFem object ID=0 dim=2, target_dim=1, nbdof=9,[EQUIV, POLY, LAGR], est.degree=4
   -> FEM_QK(2,2)
  >> mf.mesh.geotrans(1)
  gfGeoTrans object ID= 0 dim=2, nbpts= 6 : GT_PK(2,2)

Although this interface seems more convenient, you must be aware that this always
induce a call to a mex-file, and additional |mlab| code::

  >> tic; j=0; for i=1:1000, j=j+mf.nbdof; end; toc
  elapsed_time =
      0.6060
  >> tic; j=0; for i=1:1000, j=j+gf_mesh_fem_get(mf,'nbdof'); end; toc
  elapsed_time =
      0.1698
  >> tic; j=0;n=mf.nbdof;  for i=1:1000, j=j+n; end; toc
  elapsed_time =
      0.0088

Hence you should always try to store data in |mlab| arrays instead of
repetitively calling the getfem functions.

Avalaible object types are :envvar:`gfCvStruct`, :envvar:`gfGeoTrans`,
:envvar:`gfEltm`, :envvar:`gfInteg`, :envvar:`gfFem`, :envvar:`gfMesh`,
:envvar:`gfMeshFem`, :envvar:`gfMeshIm`, :envvar:`gfMdBrick`,
:envvar:`gfMdState`, :envvar:`gfModel`, :envvar:`gfSpmat`, :envvar:`gfPrecond`,
and :envvar:`gfSlice`.
