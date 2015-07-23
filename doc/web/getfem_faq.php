<?php $thisPage="Faq"; include("header.inc") ?>
  <div id="content">
    <h1>Faq</h1>

    <div class="faqq">
      What are the main differences between GetFEM++ and Deal II <a
      href="http://gaia.iwr.uni-heidelberg.de/~deal">http://www.dealii.org</a>
    </div>

    <div class="faqr">
      <p>
	Of course, every single package have its own
	specificities and advantages. I think, the logic is slightly
	different in the sense that the main goal of GetFEM++ is to be
	able to handle virtually any FEM, in any number of dimensions.
      </p>
      <ul>
	<li> The Deal.II library is restricted by design to lines/quadrangles/hexahedrons. It provides tool for mesh generation, parallelization and mesh refinements, has adaptivity as the fundamental principle of the library, and 
	  is working on hp-methods (1 <= p <= 4).
	</li>
	<li>
	  <p>
	  GetFEM++ provide a large set of pre-programmed methods. It
	  is possible to use GetFEM++ without knowing the details of
	  implementation of finite element methods since they are
	  described with a character string like "PK(N,K)" and the
	  MATLAB interface hides all the c++ internals for people who
	  don't want to deal with C++.</p>
	  <p>
	    GetFEM++ is more
	  flexible, since it provides separate basic descriptions of
	  Finite Element methods, geometric description and
	  integration methods. This means that you can either choose
	  any pre-programmed fem with any geometric transformation
	  (linear, quadratic ...) and any integration method defined
	  on the same geometric element or define your own methods. If
	  you define properly a new finite element method on the
	  reference element you will be able to use it with any
	  geometric transformation.</p>
	  <p>GetFEM++ can handle FEMs of
	  any dimension, and this dimension is not fixed at
	  compile-time (it is not a template parameter). 
</p>
<p>On the other
	  side GetFEM++ is not parallelized and does not have mesh
	  generations tools.
	  </p>
	</li>
      </ul>	
    </div>

    <div class="faqq">
      The 3D graphics from getfem-matlab are ugly, there are many artifacts.
    </div>
    <div class="faqr">
      <p>
	You should disable OpenGL rendering. It won't slow the drawing,
	but these artifacts ( inconsistant orientation of faces) will disappear.
      </p>
      <p>
	<tt>
	  set(gcf,'Renderer','zbuffer');
	</tt>
      </p><p>
      If you want to completely disable OpenGL rendering in Matlab, you can put <tt>opengl neverselect</tt>
      in your <tt>~/matlab/startup.m</tt>.
      </p>
    </div>

    <div class="faqq">
      Matlab crashes very frequently when I use the getfem-matlab toolbox
    </div>
    <div class="faqr">
      <p>
	Unfortunatly, linking a big c++ library with matlab via
	mex-files has proven to be quite unstable. There are many issues
	with dynamic libraries, exceptions, dynamic casting etc.
      </p>
      
      Starting with getfem 1.5, two options are available for getfem-matlab:
      <ul>
	<li> A giant-mex C++ file (default), containing
	  everything. It works in most of the cases, but not all (for
	  example icc won't build a correct mex-file with matlab
	  6.5)</li> 
	  
	<li> A very small C mex-file, which communicates with an
	  external process (the getfem_sever
	  executable). Communications between matlab and the
	  getfem_server use RPC (Remote Procedure Calls). The
	  advantage is that getfem and matlab process are completly
	  separated (they could even run on different machines). Hence
	  it is much easier to pin-point problems in getfem or matlab,
	  and to debug them.
      </ul>      
    </div>
    <div class="faqq">
      When the getfem-matlab interface does work: it says '<i>libgcc_s.so.1: version `GCC_3.4' not found"</i>'
    </div>
    <div class="faqr">
      <p>
        The fix for that problem, using LD_PRELOAD,  is explained <a href="https://mail.gna.org/public/getfem-users/2007-03/msg00014.html">here</a>.
      </p>
  </div>
<?php include("footer.inc") ?>

