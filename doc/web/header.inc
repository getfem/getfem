<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">
<html>
  <head>
    <title><?php echo $thisPage?></title>
    <link rel="stylesheet" type="text/css" href="style.css" media="screen" title="Normal">
  </head>
  <body>
    <div id="menubloc">
      <div id="menu">
	<div style="margin-bottom:2em;"><img src="images/logo_getfem_small.png" alt="the getfem logo"></div>
	  <ul>      
	    <li<?php if ($thisPage=="News") echo " id=\"currentpage\""; ?>>
	      <a href="news.html">What's New</a>
	    </li>
	  </ul>
	  <h3>GetFEM++</h3>
	  <ul>
	    <li<?php if ($thisPage=="GetFEM++ HomePage") echo " id=\"currentpage\""; ?>>
	      <a href="getfem_intro.html">Introduction</a>
	    </li>
	    <li<?php if ($thisPage=="Screenshots") echo " id=\"currentpage\""; ?>>
	      <a href="shots.html">Screenshots</a>
	      <?php if ($thisPage=="Screenshots") : ?>
	      <ul class="menusub">
		<li><a href="#genmesh">Generic mesh</a></li>
		<li><a href="#linelast">Linear elasticity</a></li>
		<li><a href="#stokes2d">Stokes equation (2D)</a></li>
		<li><a href="#stokes3d">Stokes equation (3D), mesh slicing</a></li>
		<li><a href="#helmholtz">Scattering &amp; high order FEM</a></li>
		<li><a href="#paolo">Eigenmodes of a structure</a></li>
		<li><a href="#donut">Contact/friction donut</a></li>
		<li><a href="#xfem">XFem crack</a></li>
		<li><a href="#nonlinelast">Non-linear elasticity</a></li>
	      </ul>
	      <?php endif ?>
	    </li>
	    <li<?php if ($thisPage=="Support") echo " id=\"currentpage\""; ?>>
	      <a href="support.html">Support/FAQ</a>
	      <?php if ($thisPage=="Support") : ?>
	      <ul class="menusub">
		<li><a href="getfem_faq.html">FAQ</a></li>
	      </ul>
	      <?php endif ?>
	    </li>
	    <li<?php if ($thisPage=="Roadmap") echo " id=\"currentpage\""; ?>>
	      <a href="roadmap.html">Future</a>
	    </li>
	  </ul>
	  <h3>Gmm++</h3>
	  <ul>
	    <li<?php if ($thisPage=="Gmm++ Home Page") echo " id=\"currentpage\""; ?>>
	      <a href="gmm_intro.html">Introduction</a>
	    </li>
	    <li<?php if ($thisPage=="Gmm++ Support/FAQ") echo " id=\"currentpage\""; ?>>
	      <a href="gmm_faq.html">Support/Faq</a>
	    </li>
	  </ul>
	  <h3></h3>
	  <ul>
	    <li<?php if ($thisPage=="Download") echo " id=\"currentpage\""; ?>>
	      <a href="download.html">Download</a>
	      <?php if ($thisPage=="Download") : ?>
	      <ul class="menusub">
		<li><a href="http://download.gna.org/getfem/stable/">Stable</a></li>
		<li><a href="http://download.gna.org/getfem/unstable/">Unstable</a></li>
	      </ul>
	      <?php endif ?>
	    </li>
	    <li<?php if ($thisPage=="Documentation") echo " id=\"currentpage\""; ?>>
	      <a href="doc.html">Documentation</a>
	      <?php if ($thisPage=="Documentation") : ?>
	      <ul class="menusub">
		<li><a href="http://download.gna.org/getfem/doc/getfemuser/getfemuser.html">GetFEM++ user guide</a></li>
		<li><a href="http://download.gna.org/getfem/doc/getfem_matlab/gfm.html">Matlab interface</a></li>
		<li><a href="http://download.gna.org/getfem/doc/getfem_python_reference.html">Python interface</a></li>
		<li><a href="http://download.gna.org/getfem/doc/gmmuser/gmmuser.html">Gmm++ user guide</a></li>
	      </ul>
	      <?php endif ?>
	    </li>
            <li<?php if ($thisPage=="Links") echo " id=\"currentpage\""; ?>>
	      <a href="links.html">Links</a>
	    </li>
	  </ul>
      </div>

      <div class="foot">
     	<a href="http://www.insa-lyon.fr"><img src="images/logoINSAt.png" border="none" alt="INSA logo"></a>
	<a href="http://math.univ-lyon1.fr"><img src="images/icj.png" border="none" alt="Lamcos logo"></a>
	<a href="http://lamcos.insa-lyon.fr"><img src="images/logolamcos.png" border="none" alt="Lamcos logo"></a>
      </div>
    </div>

