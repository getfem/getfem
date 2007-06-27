<?php $thisPage="Download"; include("header.inc") ?>
  <div id="content">
    <h1> Download Getfem++</h1>
    <p>
      Getfem++ is freely distributed under the terms of the <a
      href="http://www.gnu.org/copyleft/lesser.html">Gnu Lesser
      General Public License</a>.
    </p>

    <p>
      The latest <span class="embg">stable</span> release of getfem++:
      <ul>
	<li> core getfem++ library: <a href="http://download.gna.org/getfem/stable/getfem++-2.0.2.tar.gz">getfem++-2.0.2.tar.gz</a> (includes gmm++)</li>
<li> Matlab and Python interface: <a href="http://download.gna.org/getfem/stable/getfem-interface-2.0.tar.gz">getfem-interface-2.0.tar.gz</a> (updated 2006/03/29 to avoid a <a href="http://gcc.gnu.org/bugzilla/show_bug.cgi?id=16849">bug</a> in gcc-3.3/3.4, <b>updated 2006/04/03</b> for matlab-R13 compat of tests programs and addition of missing gfPrecond.m file)</li>
	<li> gmm++ standalone: <a href="http://download.gna.org/getfem/stable/gmm-2.02.tar.gz">gmm-2.02.tar.gz</a></li>
      </ul>
    </p>
    <p>
    The latest <span class="embg">unstable</span> releases (cvs snapshot) can be found <a href="http://download.gna.org/getfem/unstable/">here</a>.
    </p>
    <p>
    For older releases, look <a href="http://download.gna.org/getfem/stable/">here</a>.
    </p>
    <p>
      Building a portable c++ library is not an easy task. We try to
      build it with many combinations of OS and compilers.
      The last stable version has been tested on the following configurations:
      <ul>
        <li>Linux/x86 with g++ 3.x, g++ 4.x</li>
        <li>Intel C++ Compiler 8.0</li>
        <li>Linux/Itanium with g++</li>
        <li>MacOS X Tiger (with the python and matlab interface)</li>
        <li>Windows with <a href="http://www.mingw.org/">MinGW</a> and <a href="http://www.mingw.org/msys.shtml">MSys</a> (getfem++ only -- see <a href="http://download.gna.org/getfem/misc/getfem-matlab-2.0_R14_win32.README.txt">specific notes</a> for the matlab interface)</li> 
      </ul>
    <p>


  Some "not-so-easy" pre-built binaries are also available:
  <ul>
    <li>
      (2006/04/03) binary for the matlab-interface, for matlab-R13 on linux/i386 only (crashes with matlab-R14): <a
	href="http://download.gna.org/getfem/misc/getfem-matlab-2.0_R13_i386.bin.tar.gz">getfem-matlab-2.0_R13_i386.bin.tar.gz</a>
      (and some <a
	href="http://download.gna.org/getfem/misc/getfem-matlab-2.0_R13_i386.README.txt">notes</a>
      on how it was built).
    </li>
    <li> T(2006/04/18) binary for the matlab-interface for matlab-R14 on windows:
      <a href="http://download.gna.org/getfem/misc/getfem-matlab-2.0_R14_win32.zip">getfem-matlab-2.0_R14_win32.zip</a> (and some <a href="http://download.gna.org/getfem/misc/getfem-matlab-2.0_R14_win32.README.txt">notes</a>).</li>
  </ul>
  
    </p>
  </div>
<?php include("footer.inc") ?>
