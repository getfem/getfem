<?php $thisPage="Getfem-python interface"; include("header.inc") ?>
  <div id="content">
  <h1>Getfem-python interface</h1>

  <h2>Introduction</h2>

  As of version 1.7, getfem++ provides an interface to the <a
    href="http://www.python.org">Python</a> scripting language. Python is
  a nice, cross-platform, and free language. With the addition of the <a
    href="http://www.stsci.edu/resources/software_hardware/numarray">numarray</a>
  package, python provides a basic subset of Matlab functionalities
  (i.e. dense arrays). The <a href="http://public.kitware.com/VTK/">VTK</a> toolkit may provide visualization tools
  via its python interface (or via <a href="http://mayavi.sourceforge.net/">mayavi</a>), and data files for <a href="http://www.opendx.org/">openDX</a>
  may be exported. The sparse matrix routines are provided by the getfem
  interface.

  <h2>Building the python interface</h2>

  Use <tt>./configure --enable-python=yes</tt> when the
  getfem-interface is built. It requires the python developpement
  files (<tt>python.h</tt> etc.) to be available, and also the
  numarray package to be installed (its installation is
  straightforward if it not provided by your linux distribution).

  <h2>The getfem module</h2> 

  The python interface is available via a python module getfem.py. In
  order to use the interface you have to load it with

  <pre>
import getfem;
m=getfem.Mesh('cartesian', range(0, 3), range(0,3))
  </pre>
  or 
  <pre>
from getfem import *;
m=Mesh('cartesian', range(0, 3), range(0,3))
  </pre>

  <p>
    If the <tt>getfem.py</tt> (and the internal getfem_.so) module is not installed in a standard location for python, you may have to set the <tt>PYTHONPATH</tt> environnement variable to its location.
  </p>

  <p>
    A nice command-line python shell is <a href="http://ipython.scipy.org/">ipython</a>.
  </p>

  <h2>getfem-python Classes</h2>
  The general organization of the python-interface is the following:
  <ul>

    <li> Each class from the matlab interface has a corresponding class in
      the python interface: the gfMesh class becomes the getfem.Mesh class
      in python, the gfSlice becomes the getfem.Slice etc.

    <li> Each get and set method of the matlab interface has been
      translated into a method of the corresponding class in the python
      interface. For example

      <code> gf_mesh_get(m, 'outer faces'); gf_mesh_get(m, 'pts'); </code>
      becomes 
      <code>m.outer_faces(); m.pts();</code>

      Some methods have been renamed when there was ambiguity, for example 
      <code>gf_mesh_set(m, 'pts', P)</code> is <code>getfem.Mesh.set_pts(P)</code>

    <li> 
      
      The other getfem-matlab function function have a very simple
      mapping to their python equivalent:
  <table style="margin:1em;">
    <tr>
      <td width="40%">
	<tt>gf_compute(mf, U, 'foo',...)</tt>
      </td>
      <td width="40%">
	<tt>getfem.compute_foo(mf, U)</tt> or <tt>getfem.compute('foo',...)</tt>
      </td>
    </tr>
    <tr>
      <td>
	<tt>gf_asm('foobar',...)</tt>
      </td>
      <td>
	<tt>getfem.asm_foobar(...)</tt> or <tt>getfem.asm('foobar',...)</tt>
      </td>
    </tr>
    <tr>
      <td>
	<tt>gf_linsolve('gmres',...)</tt>
      </td>
      <td>
	<tt>getfem.linsolve_gmres(...)</tt>
	or
	<tt>getfem.linsolve('gmres',...)</tt>
      </td>
    </tr>
  </table>
</ul> 

  <h2>memory management</h2>

  <p>A nice advantage over the Matlab interface is that you do not have
    to explicitely delete objects that are not used any more, this is done
    automagically. You can however inspect the content of the getfem
    workspace with the function <tt>getfem.memstats()</tt>.
  </p>
  
  <h2>Documentation</h2>

  The getfem.py module is largely documented. This documentation has
  been extracted into the <a
    href="getfem_python_reference.html">getfem-python reference</a>. The
  getfem-matlab user guide may also be used, as 95% of its content
  translates quite directly into python (with the exception of the
  plotting functions, which are specific to matlab).


  <h2>Examples</h2>
  <ul>
    <li>
      <tt><a href="demo_tripod.py.html">tests/python/demo_tripod.py</a></tt> : this is the python equivalent of the matlab demo_tripod. There is also a <tt>demo_tripod_alt.py</tt> which does not use the model bricks.
    </li>
    <li>
      <tt><a href="demo_plate.py.html">tests/python/demo_plate.py</a></tt> : an example of use of the linear plate model bricks.
    </li>
  </ul>
</div>
<?php include("footer.inc") ?>

