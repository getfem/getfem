<?php $thisPage="Getfem++ HomePage"; include("header.inc") ?>
<META NAME="Keywords" CONTENT=" finite element library, finite element package, finite element software, finite elements">
  <div id="content">
    <div style="text-align:center;"><img src="images/logo_getfem.png" alt="the Getfem++ logo"></div>
    <h1>What is Getfem++</h1>

    <p>The Getfem++ project focuses on the development of a generic C++
       finite element library which aims to offer the widest range of
       finite element methods and elementary matrix computations for
       the approximation of linear or non-linear problems, possibly
       in hybrid form and possibly coupled. The dimension of the problem
       is arbitrary and may be a parameter of the problem. Getfem++ offers
       a description of models in the form of bricks whose objective is to
       enable reusability of the approximations made. The system of bricks,
       now mature, is used to assemble components such as standard models
       (elasticity in small and large deformations, Helmholtz problem,
       scalar elliptic problem ...) to components representing the
       boundary conditions (Neumann, Dirichlet, Fourier-Robin, contact,
       friction...), also to components representing constraints
       (incompressibility, removing rigid motions ...) and to coupling
       components for coupled models.</p>


    <p>Two strong points of Getfem++ are structural mechanics (in particular contact mechanics) and taking into account discontinuities by fictitious domain methods of XFEM type (eg cracking).</p>

    <p>It is proposed three interfaces (with Scilab, Matlab and Python) that allow to use of the main features of the software without the need of C++ programming and allowing graphical post-processing.</p>

    <p>Getfem++ offers a complete separation between the integration methods (exact or approximated), geometric transformations (linear or not) and finite element methods of arbitrary degree. The library can help to write more integrated finite element codes in relieving the basic technical calculations.</p>

    <p>Examples of families of finite elements available are: Pk on simplices of arbitrary degree and dimension, Qk on parallelepipeds, P1, P2 with bubble functions,
 Hermite elements, Argyris element, HCT and FVS, elements with hierarchical basis (for multigrid methods for instance), discontinuous Pk and Qk, XFEM methods, vector elements (RT0, Nédélec) ...</p>

    <p>The addition of a new finite element method is relatively easy. A description on the reference element must be provided (in most cases it is the description of the basic functions and nothing more). Extensions are provided to describe Hermite elements, piecewise polynomial or non-polynomial elements, vector elements and XFEM.</p>

    <p>The library also includes the usual tools for finite elements such as assembly procedures for classical PDEs, interpolation methods, the calculation of norms, mesh operations (including automatic refinement), management of boundary conditions, post-treatment with a tool to make arbitrary cuts ...</p>

    <p>
      Getfem++ can be used to build very general finite elements
      codes, where the finite elements, integration methods, dimension
      of the meshes, are just some parameters that can be changed very
      easily, thus allowing a large spectrum of
      experimentations. Several examples are provided (see the <a
      href="shots.html">screenshot section</a>).
    </p>


    <p>
      Getfem++ has no meshing capabilities (apart regular meshes and a small attempt),
      hence it is necessary to import meshes. Imports formats
      currently known by getfem are <a
      href="http://gid.cimne.upc.es/">GiD</a> , <a
      href="http://www.geuz.org/gmsh/">GmSH</a> and <a
      href="http://pauillac.inria.fr/cdrom/www/emc2/eng.htm">emc2</a>
      mesh files.  However, given a mesh, it is possible to refine it automatically.
    </p>

    <h1>Gmm++</h1>

    <p>
      Getfem++ includes a <a href="gmm_intro.html" title="Generic Matrix Methods">generic matrix template</a> library inspired by <a
      href="http://www.osl.iu.edu/research/mtl/" title="Matrix Template Library">MTL</a> and <a
      href="http://www.osl.iu.edu/research/itl/" title="Iterative Template Library">ITL</a>.  
    </p>
      
    <h1>Matlab interface</h1>

    <p>
      A <a href="http://www.mathworks.com">Matlab</a>&reg; interface
      to this library is also provided. Hence it is possible to use
      Getfem++ without any knowledge of C++. Moreover, this interface
      provides some post-processing functions which use matlab
      graphics. A great effort has been made to offer pictures that
      represent precisely the finite element solution (i.e. preserve
      discontinuities across elements, preserve polynomial order of
      f.e.m,..), and tools are provided for slice views with respect
      to a plane/half space/cylindar/sphere, and streamlines.  The
      pictures from the <a href="shots.html">screenshot</a> section were
      all generated via the Matlab interface.
    </p>

    <h1>Python interface</h1>

      A <a href="http://www.python.org">python</a> interface is also
      available, it is very similar to the Matlab one. Coming soon:
      the graphical post_processing unctions with tvtk.

    <h1>Scilab interface</h1>

      A <a href="http://www.scilab.org">scilab</a> interface has also
      more recently been developped. It is also very similar to the
      Matlab one and in particular offers the same graphical capabilities.

    <h1>Awards</h1>
      
      Getfem++ has been awarded by the second price at the
      <a href="http://www.tropheesdulibre2007.org/"> "Trophées du Libre 2007"</a> in the category of scientific softwares.

    <h1>Licence</h1>

    <p>
      Getfem++ is freely distributed under the terms of the
      <a href="http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html">
      Gnu Lesser General Public License, either version 2.1 of the license or any later version</a>.
    </p>
    <!-- Mise en commentaire de la suite
    <div class="footer">  
      <p>This project is supported by the <a
					     href="http://www-gmm.insa-toulouse.fr/">Department of
	  Mathematics</a> of the INSA Toulouse and by the <a
							     href="http://mip.ups-tlse.fr">MIP laboratory</a>.</p>

      <p>Authors: <a href="mailto:Yves.Renard@insa-lyon.fr">Yves
	  Renard</a>, <a
			 href="mailto:Julien.Pommier@insa-toulouse.fr">Julien
	  Pommier</a>.</p>

      <p>Many thanks to Pierre Saramito, Xavier Montagutelli and
	Thomas Montfort for their help. <br><br> Many thanks to Philippe
	    Guillaume, Alain Huard, Faker Ben Belgacem, Michel
	    Fourni&eacute; and Abderrahmanne Bendali for their
	    collaboration. </p>
    </div>
     -->
<a href="http://gna.org/projects/getfem"><img
  src="images/hostedbygna.png" alt="gna" border="none"></a>
 </div>
<?php include("footer.inc") ?>

