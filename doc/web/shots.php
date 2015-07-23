<?php $thisPage="Screenshots"; include("header.inc") ?>
  <div id="content">
  <h1>GetFEM++ in action..</h1>

  <a name="genmesh"></a><h3 class="sshot">Generic mesh handling</h3>
  <p>
    The first images illustrate the general mesh handling of getfem. The <a
      href="strange.mesh_fem">mesh description</a> is hand-made, and involves many
    different element types and convex types, as you can see (the mesh, and a random
    function interpolated on the mesh):
  </p>
  <p align="center">
    <a href="images/strangemesh.png"><img src="images/strangemesh_small.png" alt="strange mesh" border="none"></a>
    <a href="images/strangernd.png"><img src="images/strangernd_small.png" alt="strange mesh with random data" border="none"></a>
  </p>
  <p>
    The mesh is 3D. There is a quadrangle, a curved quadrangle/triangle, a kind of
    curved prism and hexahedron, and a very curved (geometrical transformation of
    degree 3) quadrangle.</p>

  <a name="linelast"></a><h3 class="sshot">Linear elasticity</h3>
  <p>
    A tripod is fixed on the ground
    and loaded with a vertical force on its top. The mesh was generated with <a
      href="http://gid.cimne.upc.es/">GiD</a>, using quadratic (i.e. curved)
    tetrahedrons. The solution is computed on a P2 FEM (i.e. <em>P2 isoparametric FEM</em>).
    Below is the Von Mises stress, represented on the deformed tripod. The source
    code of this example uses the matlab interface, and can be found <a
      href="demo_tripod.html">here</a>.
  </p>
  <p align="center"><a href="images/tripodvonmiseswithmesh.png"><img
	src="images/tripodvonmiseswithmesh_small.png" alt="tripod" border="none"></a>
  </p>
  <p>
    If you want to see what is inside the tripod, download the following animation (mpeg-4 movie, 6MB, 45secs) <a href="http://download.gna.org/getfem/misc/tripod_slice.avi">tripod_slice.avi</a>
  </p>

  <a name="stokes2d"></a><h3 class="sshot">Stokes equation</h3>
  <p>
    An incompressible viscous fluid
    flows in a 2D tube. The mesh is made of curved triangles, and the solution is
    computed on a mixed P2+/P1 FEM (P2 with a cubic bubble for the velocity field,
    and discontinuous P1 for the pressure field). The source code is <a
      href="demo_stokes_tube2D.html">here</a>. 
  </p>
  <p align="center"><a href="images/tube.png"><img src="images/tube_small.png" alt="2D tube" border="none"></a>
  </p>

  <a name="stokes3d"></a>
  <p>The next example is still the Stokes problem, inside a 3D cylindrical
    tank. The picture show the norm of the fluid velocity, with some streamlines.
  </p>

  <p align="center"><a href="images/cuve.png"><img
	src="images/cuve_3D_streamlines_small.png" alt="3D tank" border="none"></a></p>

  <a name="helmholtz"></a><h3 class="sshot">Helmholtz equation</h3>
  <p>This is a basic 2D
    scattering example. An incoming plane wave is scaterred by a perfectly
    reflective circular obstacle. The mesh is made of only 25 quadrangles
    whose geometric transformations are polynomials of degree
    6. Computations are done with a P10 FEM, hence it is possible to have
    2 wavelength per element ! (with a P1 fem, the rule is at least 6
    elements per wavelength). The source is <a
      href="demo_wave2D.html">here</a>.
  </p>
  <p align="center">
    <img src="images/helm_mesh_k7_P10_gt6.png" alt="helmholtz mesh" border="none">
    <img src="images/helm_k7_P10_gt6.png" border="none" alt="the real part of the scaterred field">
  </p>

  <a name="paolo"></a><h3 class="sshot">Eigenmodes of a structure (thanks to Paolo Bertolo)</h3>
  <p align="center"><img src="images/modestructure_paolo_small.png" border="none"
      alt="eigenmode of a vibrating structure"></p>
  You can look at a small movie showing the 24 first modes of the structure: <a
    href="http://download.gna.org/getfem/misc/oggetto_modes.mpeg">(mpeg1, 4MB)</a> or <a
    href="http://download.gna.org/getfem/misc/oggetto_modes.avi">(mpeg4, 8MB)</a>.

  <a name="donut"></a>
  <h3 class="sshot">Contact with friction problem (Houari Khenous)</h3>

  <p>
    This example shows the deformation of a tire under its own weight. The tire is meshed with one layer of
    regular hexahedric cells (384 cells), whose geometric transformation is of order
    2, and a Q2 FEM. This picture shows the Von Mises criterion on the deformed
    tire.
  </p>
  <p align="center">
    <img src="images/pneu_Q2_vonmises_small.png" border="none"
      alt="contact problem">
  </p>
  <p> An animation of a (soft) elastic disk is also available (mpeg-4 movie, 4MB, 12secs) <a href="http://download.gna.org/getfem/misc/disk_in_contact.avi">(mpeg1, 4MB)</a> (mpeg-4 movie, 1MB, 12secs) <a href="http://download.gna.org/getfem/misc/disk_in_contact.avi">(mpeg1, 1MB)</a> (A newmark scheme adapted for the unilateral contact condition) 
  </p>

  <a name="xfem"></a>
  <h3 class="sshot">Xfem cracks in a beam</h3>

  <p>
    Here we used XFem to handle cracks in a beam. XFem is an enrichment of the classical finite element space (a P2 FEM was used for this example) with
  </p>
  <ul>
    <li>A discontinuous function. Thanks to this function, the crack path does not have to follow the original mesh. Note how the crack cross elements on the mesh below.</li>
    <li>Four singular functions, which form a basis for asymptotical solution to the linear elasticity problem near the crack tips.
  </ul>
  <p align="center">
    <img src="images/xfembeammesh.png" title="The original mesh, with the 1D meshes of the cracks" alt="xfem mesh of a cracked beam">
  </p>
  <p align="center">
    <img src="images/xfembeam.png" title="The Tresca criterion on the deformed beam" alt="Tresca criterion on a cracked beam">
  </p>

  
  <h3 class="sshot">a 3D crack, made via level-set</h3>
  <p>
    In this example, the mesh was a simple cartesian mesh <tt>20x20x1</tt>, and the crack geometry was defined implicitely via a levelset.
  </p>
  <p align="center">
    <img src="images/fissure_3d_de_traviole.png" title="a 3D crack" alt="a 3D crack">
  </p>
  
  <a name="nonlinelast"></a>
  <h3 class="sshot">Large strain</h3>
  <p>
    In this example, a bar is twisted. Each step is solved with a Newton method. The material law is a "Ciarlet Geymonat" one. A P2 FEM is used. The source code for this example can be found in the <tt>tests/nonlinear_elastostatic.C</tt> file of getfem++ package. This picture was made with OpenDX.
  <p align="center">
    <img src="images/torsion034.png" title="Torsion of a rubber bar" alt="">
  </p>
  <p>
    A short animation is also available: (mpeg-4 movie, 3MB) <a href="http://download.gna.org/getfem/misc/torsion.avi">torsion.avi</a>.
  </p>

 <a name="shapeoptimization"></a>
  <h3 class="sshot">Shape and topological optimization</h3>
  <p>
    This images were obtained with the script interface/tests/matlab/demo_structural_optimization.m (Alassane SY and Yves Renard). It represents a shape optimization of a structure submitted to a vertical load at the right and clambed at the left. A (Xfem like) fictitious domain approach is used together with both a shape gradient and a topological gradient.
  <p align="center">
    <img src="images/shape1.png" title="Shape optimization, remaining surface 1.039 / 2" alt="">
    <img src="images/shape2.png" title="Shape optimization, remaining surface 0.954 / 2" alt="">
  </p>
  The first image corresponds to an initial structure with pre-existing holes. For the second one the holes are initiated by the topological optimization. The two following images correspond to a 3D case.
    <p align="center">
    <img src="images/shape3.png" title="3D shape optimization" alt="" height="100%">
    <img src="images/shape4.png" title="3D shape optimization" alt="" height="100%">
  </p>
  </div>
<?php include("footer.inc") ?>

