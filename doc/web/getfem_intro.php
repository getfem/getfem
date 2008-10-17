<?php $thisPage="Getfem++ HomePage"; include("header.inc") ?>
<META NAME="Keywords" CONTENT=" finite element library, finite element package, finite element software, finite elements">
  <div id="content">
    <div style="text-align:center;"><img src="images/logo_getfem.png" alt="the Getfem++ logo"></div>
    <h1>What is Getfem++</h1>
    <p>The Getfem++ project focuses on the development of a generic and
      efficient C++ library for finite element methods. The
      goal is to provide a library allowing the computation of any elementary
      matrix (even for mixed finite element methods) on the largest class of
      methods and elements, and for arbitrary dimension
      (i.e. not only 2D and 3D problems). 
    </p>
    <p>
      It offers a complete separation between integration methods
      (exact or approximated), <span class="popup"
      title="transformation from the reference element to the real element">geometric transformations</span> (linear or not) and
      finite element methods of arbitrary degrees. It can really
      relieve a more integrated finite element code of technical
      difficulties of elementary computations.
    </p>
    <p>
      Examples of available finite element method are : Pk on simplices in
      arbitrary degrees and dimensions, Qk on parallelepipeds, P1, P2 with
				       bubble functions, Hermite elements, Argyris element, elements with hierarchic basis
				       (for multigrid methods for instance), discontinuous Pk or Qk, XFem, vectorial elements (RT0, Nedelec) ...
    </p>
    <p>
      The addition of a new finite element method is
      relatively easy. Its description on the 
      reference element must be provided (in most of
      the cases, this is the description of the basis functions, and
      nothing more).  Extensions are provided for Hermite elements,
      piecewise polynomial, non-polynomial, vectorial elements and XFem.
    </p>
    <p>
      The library also includes the usual tools for finite
      elements such as assembly procedures for classical PDEs,
				       interpolation methods, computation of norms, mesh operations (including automatic refinement), boundary
      conditions, post-processing tools such as extraction of slices from a mesh ... 
    </p>

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
      the graphical fpost_processing unctions with tvtk.

    <h1>Awards</h1>
      
      Getfem++ has been awarded by the second price at the
      <a href="http://www.tropheesdulibre.org/"> "Trophées du Libre 2007"</a> int the category of scientific softwares.

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

