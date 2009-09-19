<?php $thisPage="News"; include("header.inc") ?>
  <div id="content">			  
  <h1>What's new ?</h1>

  <div class="gfnews">
    <h2>2009/09/19 Getfem++-4.0.0. released</h2>
     <p>
      This is a major update to Getfem++. The main changes is the
      introduction of a new model 
      bricks system. The old system is kept and compatibility with 3.x
      releases is globally ensured. However some functionalities are
      deprecated.

      The main changes are:
      </p>
      <ul>
        <li>
	The mesh_fem object has undergone significant changes. Now it is possible to perform linear combination of degrees of freedom in order to describe some special finite element spaces. The main application is to obtain a finite element space reduced on a boundary or a curve. But it can be used also to prescibe directly some matching condition. The main change in the use of the mesh_fem object is the introduction of "basic" and "reduced" dofs. See the documentation.
      </li>
      <li>
        A new algorithm gmm_range_basis allows to select a basis between the columns of a matrix. It has been specially designed to select a basis of the trace on a boundary of a finite element space. 
      </li>
      <li>
	The partial_mesh_fem object has been completely changed. It is now a lighter object which is intensively used in the new model bricks to obtain finite element spaces on a boundary. 
      </li>
      <li>
	Introduction of the new model brick system. The bricks are more simple to build and it is now really designed to the representation of coupled/multiphysics models. A generic manner to deals with time dependent models from static models is also introduced.
      </li>
      <li>
        Python interface uses Numpy instead of Numarray.
      </li>
    </ul>

All the old bricks have not been rewritten into new bricks. This will be done
gradually in the near future. A Scilab interface is close to be finished and should be included in the future release. 

    </p>
    <h2>2008/09/09 Getfem++-3.1. minor version</h2>
     <p>
      A certain number of small bug fixed in Getfem++ and Gmm++.
      Clarification of copyrights. 
    </p>
    <h2>2007/07/12 Getfem++-3.0.1. minor update</h2>
     <p>
      Two bugs were fixed: a memory leakage problem
      and a bad identification of some dofs.
    </p>
    <h2>2007/06/27 Getfem++-3.0 released</h2>
      <p>Getfem++ 3.0 is now available !</p>
      <p>Not so many changes, but some of them are incompatible with getfem 2.0:</p>
      <ul>
         <li>The Getfem and Gmm header files have been moved into their respective subdirectories. So, as a consequence, the include directives have to be updated:

<p><tt>#include "gmm_xxx.h"</tt> should be replaced with <tt>#include "gmm/gmm_xxx.h"</tt></p>
<p><tt>#include "getfem_xxx.h"</tt> should be replaced with <tt>#include "getfem/getfem_xxx.h"</tt></p>

         <li>The Getfem interface (python and matlab) is now included in the Getfem tar.gz file, in the '<tt>interface</tt>' subdirectory. They can be enabled with the '<tt>--enable-python</tt>' or '<tt>--enable-matlab</tt>' switch of the <tt>configure</tt> script</li>
         <li>Some C1 composite elements have been added (triangles and quadrilaterals)</li>
         <li>Levelset support has been improved</li>
      </ul>
      The full list of changes is available in the <a href="http://svn.gna.org/viewcvs/getfem/trunk/getfem%2B%2B/ChangeLog?rev=2640&amp;view=auto">ChangeLog</a>.
    <h2>2006/11/10 Getfem++-2.0.2, minor update</h2>
    <p>
      The GMSH mesh import has been fixed.
    </p>
    <h2>2006/04/06 Getfem++-2.0.1, minor update</h2>
    <p>
      Two bugs were fixed which could be toggled in particular conditions with nonlinear terms.
    </p>
    <h2>2006/03/20 Getfem++-2.0, Gmm++-2.0 and Getfem-Interface 2.0 released</h2>
    <p>
      This is a major update to Getfem++, which make some backward-incompatible changes:
    </p>
    <ul>
      <li>
	the old <tt>mesh_fem</tt> has been split into two disjoint
	classes: <tt>mesh_fem</tt> which handles all that is related
	to FEM, and <tt>mesh_im</tt> which handles the integration
	methods on a mesh.
      </li>
      <li>
	the old <tt>getfem::getfem_mesh</tt> class has been renamed to <tt>getfem::mesh</tt>
      </li>
      <li>
	the "boundaries" which were attached to a <tt>mesh_fem</tt> in
	previous versions, are now attached to a <tt>mesh</tt>, and
	they are now called "regions" (because they can stored
	boundaries, and also sets of convexes).
      </li>
      <li>
	the model bricks have been reworked -- especially the Dirichlet conditions.
      </li>
    </ul>
    <p>
    Some news features have been introduced in this release:
    </p>
    <ul>
      <li>
	introduction of level-set objects. Integration methods can be cut
	with respect to these level-set and discontinuous elements
	across the level-set are provided.
      </li>
      <li>
	parallelization of the assembly.
      </li>
      <li>
	interface to MUMPS.
      </li>
      <li>
	many news elements, Hermite and vectorial elements are now fully supported: 1D, 2D and 3D hermite, Argyris triangle, HCT triangle, RT0 and Nedelec elements are now available.
      </li>
      <li>
	automatic mesh refinement.
      </li>
    </ul>

    <p>
      Major changes for the matlab and python interface: they follow
      the changes that occured in Getfem. An interface to the Getfem++
      model bricks has been added.
    </p>

    <p>
      Next releases of Getfem++ will try to maintain backward compatability with this release.
    </p>

    <h2>2005/01/05 Getfem++ 1.7, Gmm++ 1.7 and Getfem-Interface 1.7 released</h2>  
    <p>
      An important number of improvements have been done on Getfem++ 1.7. Note that the next release will be Getfem 2.0, some of its changes won't maintain backward compatibility with getfem++-1.7. 
    </p>
    <ul>
      <li>
	Introduction of the "model brick" system, which provides a general framework for the solution of common PDEs. Each brick is dedicated to a specific task (i.e. "handle Dirichlet conditions", "assembly of the Stokes Problem", "solve a linear system", etc.). These bricks are then connected to each other. Examples of use can be found in the "tests/" directory of Getfem++.
      </li>
      <li>
	New models : Small strain plasticity, <a href="torsion034.png">large strain elasticity</a>,
	contact and friction conditions, linearized plates,
	incompressibility in small and large strain elasticity.
      </li>
      <li>
	Simplifications and optimizations in elementary computations.
      </li>
      <li>
	A direct sparse solver (<a href="http://crd.lbl.gov/~xiaoye/SuperLU/">SuperLU 3.0</a>) is available "out of the box"
      </li>
      <li>
	Ability to export results to <a href="http://www.vtk.org">VTK</a> and <a href="http://www.opendx.org">OpenDX</a>.
      </li>
    </ul>
    <p>
      Major changes in Gmm++ 1.7:
    </p>
    <ul>
      <li>New preconditionner ILUTP.</li>
      <li>A BFGS algorithm has been developped.</li>
      <li>gmm++ now handles (valid) operations mixing complex and scalars.</li>
      <li>gmm::real_part(V) and gmm::imag_part(V) gives a possibly writable reference on the real and imaginary part of a complex vector or matrix.</li>
      <li>the SuperLU interface has been updated for SuperLU 3.0.</li>
    </ul>
    <p>
      getfem-matlab has been renamed "getfem-interface" since it now provides an interface for Matlab and <a href="http://www.python.org">Python</a> (with the <a href="http://www.stsci.edu/resources/software_hardware/numarray">Numarray</a> package). Note that, while it is <a href="getfem_python_reference.html">documented</a> and working, the python interface is still considered a "work in progress". You have to enable it explicitly with "./configure --enable-python". An example of use can be found <a href="demo_tripod.py.html">here</a>. An interface to the gmm++ sparse matrices and solvers is also provided.
    </p>
  </div>

  <div class="gfnews">
    <h2>2004/01/23 Getfem++ 1.6 and Gmm++ 1.6 released</h2>  
    <p>
      Getfem++ 1.6 is mostly a bugfix and performance improvements
      release.
    </p>
    <ul>
      <li>
	Some new integration methods were added (high order methods for
	triangles such as &quot;IM_TRIANGLE(19)&quot; from <it>P. Solin, K. Segeth
      and I. Dolezel: &quot;Higher-Order &amp; Finite Element Methods&quot;, Chapman
      & Hall/CRC Press, 2003</it>).
  </li>
    <li>
      Performance of interpolation and geometric transformation
      inversion was much improved.
    </li>
    <li>
      Support for emc2 meshes
    </li>
  </ul>
    <p>
      The Gmm++ library has been much improved version 1.6 and version 1.5. We have especially focused on its robustness.
    </p>
    <ul>
      <li>Many bugs were fixed, especially for complex matrices.</li>
      <li>QR algorithms were introduced for dense matrices.</li>
      <li>A LAPACK/ATLAS interface is available.</li>
      <li>SuperLU 2.0 interface.</li>
      <li>Small simplification in <code>linalg_traits</code> structure.</li>
      <li>Generic resize procedures for vector and matrices were introduced.</li>
      <li>It is possible to use a column or row matrix view of a vector with <code>gmm::row_vector</code> and <code>gmm::col_vector</code>.</li>
      <li>Generic <code>gmm::reshape</code> and <code>gmm::conjugated</code> functions.</li>
      <li>Intensive tests with random type of matrices and vectors.</li>
    </ul>
  </div>
  <div class="gfnews">    
    <h2>2003/07/25 Getfem++ 1.5 and Gmm++ 1.5 released</h2>

    First standalone release of Gmm++, which now includes some preconditioners and harwell-boeing/matrix-market data file support. 
    It is now possible to use high precision computations of elementary integrals with the (optional) QD library.
    Quadrature data has been moved into data files in the <tt>cubature/</tt> directory.
    Initial support for XFem. 
    Mesh slices in Getfem++ and getfem-matlab. The Matlab interface was merged into a single giant mex-file.
  </div>
  
  <div class="gfnews">
    <h2>2003/03/03 Getfem++ 1.4 released</h2>
    
    The Matlab interface is now fully working and documented. Huge speed
    improvement on elementary computations. New generic assembly
    procedures. Introduction of Gmm++.
  </div>

  <div class="gfnews">
    <h2>2002/09/24 Getfem++ 1.3 released</h2>
    Introduction of hierarchical and composite FEMs and integration methods.
  </div>
  
  <div class="gfnews">
    <h2>2002/08/21 Getfem++ 1.2 released</h2>
    
    Introduction of the Hermite element (not fully working). Support for
    non-tau-equivalent elements. Introduction of a consistent naming
    system for FEMs, geometric transformations and integration methods.
  </div>
  
  <div class="gfnews">
    <h2>2002/07/18 Getfem++ 1.1 released</h2> Many improvements.
    Introduction of the Matlab interface.
  </div>
  <div class="gfnews">
    <h2>2002/06/28 Getfem++ 1.0 released</h2> First public release.
  </div>
  </div>
<?php include("footer.inc") ?>

