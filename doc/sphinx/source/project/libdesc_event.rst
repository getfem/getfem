.. $Id: libdesc.rst 4917 2015-03-27 14:37:59Z renard $

.. include:: ../replaces.txt

.. highlightlang:: none

.. _dp-libdesc_event:


Events management
-----------------

Description
^^^^^^^^^^^

The ``mesh``, |mf|, |mim| and |mo| description are linked together in the sense
that there is some dependencies between them. For instance, when an element is
suppressed to a mesh, the |mf| object has to react.

Files
^^^^^

.. csv-table::
   :header: "File(s)", "Description"
   :widths: 8, 15

   :file:`getfem_context.h` and :file:`getfem_context.cc`, "Define a class `context_dependencies` from which all object has to derive in order to manage events."


State
^^^^^

The main tool to deal with simple dependence of object is in
:file:`getfem_context.h`. An object ``context_dependencies`` is defined there. In
order to deal with the dependencies of an object, the object
``context_dependencies`` needs to be a parent class of this object. It adds the
following methods to the object:

.. c:function:: add_dependency(ct)

   Add an object (which has to have ``context_dependencies`` as a parent class)
   to the list of objects from which the current object depend.

.. c:function:: touch()

   Indicates to the dependent objects that something has change in the object.

.. c:function:: context_check()

   Check if the object has to be updated. if it is the case it makes first a
   check to the dependency list and call the update function of the object. (the
   update function of the dependencies are called before the update function of
   the current object).

.. c:function:: context_valid()

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

The event management of some objects should be analysed with care. This is the case for instance of |mls|, |mfls|, |pmf|, etc.

The event management still have to be improved to be a fully reactive system.

