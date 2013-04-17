.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-contact-friction-large:

Large sliding/large deformation contact with friction bricks
------------------------------------------------------------

These bricks present some algorithms for contact and friction in the large sliding/large deformation framework. Of course, their computational cost is greatly higher than small sliding-small deformation bricks.

The multi-contact frame object
++++++++++++++++++++++++++++++

A |gf| object is dedicated to the computation of effective contact surfaces which is shared by all the bricks. This object stores the different potential contact surfaces. On most of methods, potential contact surface are classified into two categories: master and slave surface (see  :ref:`figure<ud-fig-masterslave>`).

.. _ud-fig-masterslave:

.. figure:: images/getfemusermodelmasterslave.png
   :align: center
   :scale: 60

The slave surface is the "contactor" and the master one the "target". Rigid obstacle are also considered. They are always master surfaces.  The basic rule is that the contact is considered between a slave surface and a master one. However, the multi-contact frame object and the |gf| bricks allow multi-contact situations, including contact between two master surfaces, self-contact of a master surface and an arbitrary number of slave and master surfaces. 

Basically, in order to detect the contact pairs, Gauss points or f.e.m. nodes of slave surfaces are projected on master surfaces (see  :ref:`figure<ud-fig-masterslave>`). If self-contact is considered, Gauss points or f.e.m. nodes of master surface are also projected on master surfaces.

The use of multi-contact frame object
*************************************

A multi-contact frame object is initialized as follows::

  multi_contact_frame mcf(size_type N, scalar_type release_distance,
                          int fem_nodes_mode = 0, bool use_delaunay = true,
                          bool ref_conf = false, bool self_contact = true,
                          scalar_type cut_angle = 0.3);

where `N` is the space dimension (typically, 2 or 3), `release_distance` is the limit distance beyond which two points are not considered in potential contact (should be typically comparable to element sizes). There is several optional parameters. If `fem_node_mode=0` (default value), then contact is considered on Gauss points, `fem_node_mode=1` then contact is considered on Gauss points for slave surfaces and on f.e.m. nodes for master surfaces (in that case, the f.e.m. should be of Lagrange type) and `fem_node_mode=2` then contact is considered on f.e.m. nodes for both slave and master surfaces. if `use_delaunay` is true (default value), then contact detection is done calling `Qhull <http://www.qhull.org>`_ package to perform a Delaunay triangulation on potential contact points. Otherwise, contact detection is performed by conputing some influences boxes of the element of master surfaces. If `ref_conf` is true, the contact detection is made on the reference configuration (without taking into account a displacement) CAUTION: not fully implemented for the moment.  If `self_contact` is true, the contact detection is also made between master surfaces and for a master surface with itself. The parameter `cut_angle` is an angle in radian wich is used for the simplification of unit normal cones in the case of f.e.m. node contact : if a contact cone has an angle less than `cut_angle` it is reduced to a mean unit normal to simplify the contact detection.

Once a multi-contact frame is build, one adds slave or master surfaces, or rigid obstacles. Note that rigid obstacles are defined by a level-set expression which is evaluated by the `MuParser <http://muparser.beltoforion.de/>`_ package. The methods of multi-contact frame object adding a contact boundary are::


  size_type add_obstacle(const std::string &obs);

  size_type add_slave_boundary(const getfem::mesh_im &mim,
                               const getfem::mesh_fem &mfu,
                               const model_real_plain_vector &U,
                               size_type region);

  size_type add_master_boundary(const getfem::mesh_im &mim,
                                const getfem::mesh_fem &mfu,
                                const model_real_plain_vector &U,
                                size_type region);

where `obs` is a string containing the expression of the level-set function which should be a signed distance to the obstacle (the coordinates are (`x`, `y`) in 2D, (`x`, `y`, `z`) in 3D and , (`x`, `y`, `z`, `w`) in 4D). `region` is the boundary number.


The contact pair detection algorithm
************************************

A contact pair is formed by a point of a slave (or master in case of self-contact) surface and a projected point on the nearest master surface (or rigid obstacle). The Algorithm used is summerized in :ref:`figure<ud-fig-algodetect>`

.. _ud-fig-algodetect:

.. figure:: images/getfemusermodeldetectcontact.png
   :align: center
   :scale: 100


It is impossible to distinguish between valid and invalid contact situations without a global topological criterion (such as in [Pantz2008]_), a fortiori for self-contact detection. However, this kind of criterion can be very costly to implement. Thus, one generally implements some simple heuristic criteria which cannot cover all the possible cases. We present such a set of criteria here. They are of course perfectible and subject to change. First, in :ref:`figure<ud-fig-invalidcontact>` one can see a certain number of situations of valid or invalid contact that criteria have to distinguish.


.. _ud-fig-invalidcontact:

.. figure:: images/getfemusermodelfalsecontact1.png
   :align: center
   :scale: 90


.. figure:: images/getfemusermodelfalsecontact2.png
   :align: center
   :scale: 90

Some details on the algorithm:

  - **Computation of influence boxes.** The influence box of an element is just
    an offset to its bounding box at a distance equal to the release distance.
    If this strategy is used, the release distance should not be too large
    compared to the element size. Otherwise, a point would correspond to a
    a large number of influence box which can considerably slow down the search
    of contact pairs. The influence boxes are stored in a region tree object
    in order to find the boxes containing a point with an algorithm having
    a mean complexity in :math:`O(log(N))`.
  
  - **What is a potential contact pair.** A potential contact pair is a pair
    slave point - master element face which will be investigated.
    The projection of the slave point on the master surface will be done
    and criteria will be applied.
 
  - **Projection algorithm.** The projection of the slave point onto a
    master element face is done by a parametrization of the surface on the
    reference element via the geometric transformation and the displasement
    field. During the projection, no constraint is applied to remain inside
    the element face, which means that the element face is prolongated
    analytically. The projection is performed by minimizing the distance
    between the slave point and the projected one using the parametrization
    and a BFGS algorithm.

The list of criteria:

  - **Criterion 1: the unit normal cone/vector should be compatible, and the
    two points do not share the same element.**
    Two unit normal vector are compatible if their scalar product are
    non-positive. In case of f.e.m. node contact, since a fem node is shared
    generally by several elements, a normal cone constituted of the unit normal
    vectors of each element is considered. Two normal cones are compatible if
    at least one pair of unit normal vector have their scalar product
    non-positive. In order to simplify the computation, a normal cone is
    reduced to a mean normal vector if the solid angle of the normal cone is
    less than `cut_angle` a parameter of the multi-contact frame object.
    This criterion allows to treat cases (B) and (K1).

  - **Criterion 2: the contact pair is eliminated when the search of the
    projection point do not converge.**
    When the BFGS algorithm used to compute the projection of the slave
    point on the master element surface fails to converge, the pair is
    not considered. A warning is generated.
    
  - **Criterion 3 : the projected point should be inside the element.**
    The slave point is projected on the surface of the master element by
    a BFGS algorithm without the constraint to remain inside the face
    (which means that the face is prolongated). If the orthogonal
    projection is outside the face, the pair is not considered. This
    is the present state, however, to treat case (J3) an aditional
    treatment will have to be considered (projection on the face with
    the constraint to remain inside it and test of the normal cone at
    this point)
    This criterion allows to treat cases (F2), (K2), (M1) and (M2).

  - **Criterion 4 : the release distance is applied.**
    If the distance between the slave point and its projection on the master
    surface is greater than the release distance, the contact pair is not
    considered. This can treat cases (C), (E), (F1), (G), (H) if the release
    distance is adapted and the deformation not too important.

  - **Criterion 5 : comparison with rigid obstacles.**
    If the signed distance between the slave point and its projection on
    the master surface is greater than the one with a rigid obstacle
    (considering that the release distance is also first applied to rigid
    obstacle) then the contact pair is not considered.

  - **Criterion 6 : for self-contact only : apply a test on
    unit normals in reference configuration.**
    In case of self contact, a contact pair is eliminated when the slave point
    and the master element belong to the same mesh and if the slave point is
    behind the master surface (with respect to its unit outward normal vector)
    and not four times farther than the release distance.
    This can treat cases (A), (C), (D), (H).

  - **Criterion 7 : smallest signed distance on contact pairs.**
    Between the retained contact pairs (or rigid obstacle) the one
    corresponding to the smallest signed distance is retained.




The available bricks
++++++++++++++++++++

Sorry, for the moment no brick is fully working. Coming soon ...