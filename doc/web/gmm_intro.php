<?php $thisPage="Gmm++ HomePage"; include("header.inc") ?>
  <div id="content">
    <h1>What is GMM++</h1>
    <div id="biglogo"><img src="images/logo_gmm.png" alt="the Gmm++ logo"></div>
    <p>
      <abbr title="Generic Matrix Methods">GMM++</abbr> is a
      generic C++ template library for sparse, dense and skyline
      matrices. It is built as a set of generic algorithms (mult, add,
      copy, sub-matrices, dense and sparse solvers ...) for any
      interfaced vector type or matrix type. It can be view as a glue
      library allowing cooperation between several vector and matrix
      types. However, basic sparse, dense and skyline matrix/vector types are built
      in Gmm++, hence it can be used as a standalone linear algebra
      library.

      Interfacing a vector or matrix type means writing "traits" objects called
      "<code>linalg_traits</code>", which describe their properties. The library
      offers predefined dense, sparse and skyline matrix types.
    </p>
    <p>
      The goal is to create a general, adaptable and easy to use
      framework of pre-defined methods for matrix computation. When a
      vector or a matrix type has been interfaced (i.e. its
      <code>linalg_traits</code> has been filled), all generic algorithms works on
      it. However, it is always possible (and easy) to specialize some
      generic algorithms for efficiency reason. Major generic
      algorithms are
    </p>
      <ul>
	<li> A set of miscellaneous generic commands (clear, clean,
	scalar product, scale, norms, ...)</li> 
	
	<li> Vector-Vector addition with the possibility to mix
	formats (sparse, dense, skyline)</li>
	
	<li> Matrix-Vector mult for any format.</li>

	<li> Matrix-Matrix mult with the possibility to mix formats
	(sparse, dense, skyline, row major, column major, ...)</li>

	<li> Generic linear solvers (<abbr
	title="Conjugate Gradient">cg</abbr>, <abbr
	title="bi-Conjugated Gradient">bicgstag</abbr>, <abbr
	title="Quasi-Minimal Residual Method">qmr</abbr>, <abbr title="Generalized Minimum Residual Method">gmres</abbr> ...) with preconditioners for sparse matrices
	(<abbr title="Incomplete LU factorization with fill-in and threshold">ILUT</abbr>, <abbr title="Incomplete LU factorization with fill-in, threshold and column pivoting">ILUTP</abbr>, <abbr title="Incomplete LDLT factorization">ILDLT</abbr>, ...). Some of them are imported form <a href="http://www.osl.iu.edu/research/itl/" title="Iterative Template Library">ITL</a> (eventually corrected and optimized), some of them are new. </li>



	<li> Reference to sub-matrices (with sub-interval, sub-slice
	or sub-index) for any sparse dense or skyline matrix for read
	or write operations.</li>
	
	<li> LU and QR factorizations for dense matrices.</li>
	
	<li> Eigenvalues computation for dense matrices.</li>
      </ul>
    <p>
      The structure of GMM++ is largely inspired from <a
      href="http://www.osl.iu.edu/research/mtl/" title="Matrix Template Library">MTL</a>. The major
      differences are : simpler use, built as an interface for existing
      matrix types, sub-matrices for any matrix types. The efficiency
      is comparable (see <a href="http://grh.mur.at/misc/sparselib_benchmark/">
      http://grh.mur.at/misc/sparselib_benchmark/</a> for instance).
    </p>

    <p>
      NOTE : For performance reason, an interface with <a
      href="http://www.netlib.org/lapack/">LAPACK</a> or <a
      href="http://math-atlas.sourceforge.net/" title="Automatically Tuned Linear Algebra Software">ATLAS</a> is provided
      for dense matrices. See the <a href="http://download.gna.org/getfem/doc/gmmuser/gmmuser.html">documentation</a> (if you make some
      benchmarks, do not forget to use optimization compiler options,
      at least -O3 and you should disable checks with
      -dNDEBUG).
    </p>
    
    <p>
      A small interface to <a
      href="http://crd.lbl.gov/~xiaoye/SuperLU/">SuperLU 3.0</a>
      (sparse matrix direct solver) is also proposed for sparse
      matrices.
    </p>

    <p>
      GMM++ has been tested with
      <a href="http://www.cs.berkeley.edu/~yozo">QD</a>
      an
      efficient library for double double and quadruple double
      precision. See on the documentation how to link QD. This means
      that GMM++ should work with any reasonable arbitrary precision
      floating point library.
    </p>

    <h1>Licence</h1>     
    GMM++ is freely distributed under the terms of the <a
    href="http://www.gnu.org/copyleft/lesser.html">Gnu Lesser General
    Public License</a>.

    <h1>Contribute to GMM++</h1>
    <p>
      GMM++ offers a framework to develop efficient methods for linear algebra. This library is open-source and will remain open-source. Here are some examples of possible extensions:
    </p>
      <ul>
	<li>Specialize some algorithms to optimize them for particular matrix implementation.</li>
	<li>New solvers and preconditioners.</li>
	<li>Eigenvalues computation for sparse matrices. </li>
	<li> ...</li>
      </ul>

    <h1>Random test procedures</h1>
    <p>
      A problem with generic programming is to be sure that every
      configuration has been fully tested. This is why there is now
      a random generator of
      tests. This means that a number of test procedures will
      be called with random parameters, i.e. random type of vector,
      sub-vector, matrix or sub-matrix types, with random base type
      (float, double, long double, std::complex&lt;float&gt;,
      std::complex&lt;double&gt;, dd_real ...) and random size and filling, testing all
      the possibilities of mixing formats in operations such as mult,
      add ...
    </p>
    
    <p> 
      You are encouraged to test them, runing a "make
      check" on the distribution of GMM++ and sending us a bug report
      if it fails. We will also appreciate if you send us new test
      procedures.
    </p>
  </div>
<?php include("footer.inc") ?>

