.. $Id: libdesc.rst 4917 2015-03-27 14:37:59Z renard $

.. include:: ../replaces.txt

.. highlight:: c++

.. _dp-libdesc_dal:

Dal library
-----------



Description
^^^^^^^^^^^

In the very begining of |gf| (the first files was written in 1995) the S.T.L. was
not available and the containers defined in the ``dal`` namespace was used
everywhere. Now, in |gf|, the S.T.L. containers are mainly used. The remaining
uses of ``dal`` containers are eather historical or due to the specificities of
these containers. It is however clear that this is not the aim of the |gf|
project to develop new container concept. So, the use of the ``dal`` containers
has to be as much as possible reduced.

Furthermore, ``dal`` contains a certain number of basic algorithms to deal with static stored objects (description of finite element methods, intermediary structures for auxiliary computations ...).

Files
^^^^^

.. csv-table::
   :header: "File(s)", "Description"
   :widths: 8, 15

   :file:`dal_config.h`, Mainly load |gmm| header files
   :file:`dal_basic.h`, "A variable size array container, dal::dynamic_array<T>."
   :file:`dal_bit_vector.h` and :file:`dal_bit_vector.cc`,  "A improved bit vector container based on dal::dynamic_array<T>."
   :file:`dal_tas.h`, "A heap container based on dal::dynamic_array<T>."
   :file:`dal_tree_sorted.h`, "A balanced tree stored array based on dal::dynamic_array<T>."
   :file:`dal_static_stored_objects.h` and :file:`dal_static_stored_objects.cc`, "Allows to store some objects and dependencies between some objects. Used to store many things in |gf| (finite element methods, integration methods, pre-computations, ...)."
   :file:`dal_naming_system.h`, "A generic object to associate a name to a method descriptor and store the method descriptor. Used for finite element methods, integration methods and geometric transformations. Uses dal::static_stored_object."
   :file:`dal_shared_ptr.h`,  A simplified version of boost::shared_ptr.
   :file:`dal_singleton.h` and :file:`dal_singleton.cc`, "A simple singleton implementation which has been made thread safe for OpenMP (singletons are replicated n each thread)."
   :file:`dal_backtrace.h` and :file:`dal_backtrace.cc`, "For debugging, dump glibc backtrace."



State
^^^^^

Stable, not evolving too much.


Perspectives
^^^^^^^^^^^^

No plan.
