<?php $thisPage="Gmm++ Support/FAQ"; include("header.inc") ?>
  <div id="content">
    <h1>Bug reports</h1>
    <p>
      If you find any bug or misbehaviour of Gmm++, please send a mail
      with a full report to <a
      href="mailto:Yves.Renard@insa-lyon.fr">Yves.Renard@insa-lyon.fr</a>.
    </p>

    <h1>Known bugs</h1>
    <p>
    <ul>
<li> The computation of ILDLTT preconditioner is very very slow. Fixed in version gmm-1.7-20040408.tar.gz and later.
    </ul>
    </p>

    <h1>Gmm++ Faq</h1>
    <div class="faqq">
      GMM++ seems to crash frequently, is it bugged ?
    </div>
    <div class="faqr">
      <p>
	Remember that you have to CATCH ERRORS in your main
	procedure. GMM++ uses "throw" when an error occurs such as
	uncompatible dimensions. See the documentation.
      </p>
    </div>
    <div class="faqq">
      In a function such as mult(A, B, C), is the dimensions of the
      result (C) adapted to the dimensions of the input parameters ?
    </div>
    <div class="faqr"> 
      No, You have to declare C with the right
      dimensions. This is the case for all the algorithms in GMM++. The
      reason is that C could be a reference such as a sub matrix, a sub
      vector or a line of a matrix for instance. With last version 1.6
      of GMM++ you can use resize(C, ..) or reshape(C, ..) to change the
      size of a vector or the dimensions of a matrix, but of course only
      if it is not a reference.
    </div>
    <div class="faqq">
      How do I release the memory used by a vector or matrix ? <tt>resize(0)</tt> does not work.
    </div>
    <div class="faqr"> 
      If the object cannot be destroyed, the usual way is to swap its content with the content
      of a short-lived empty object:
      <pre>
     std::vector&lt;double&gt; v(100000000);
     ...
     { 
       dense_vector&lt;double&gt; w; 
       v.swap(w); // or swap(v,w)
     } // w is destroyed, and v.capacity() == 0
      </pre>
    </div>

    <div class="faqq">
      How does Gmm++ compare with other c++ linear algebra libraries (MTL, Pooma, uBlas, ..) ?
    </div>
    <div class="faqr"> 
      <p>
	The main difference is that Gmm++ primary aim is not to be a
	standalone linear algebra library, but is more aimed at
	interoperability between several linear algebra packages. It
	started as a glue code between the several vector and matrix
	classes used in getfem++.
      </p>
      <p>
	Its code size is kept as small as possible, and no
	attempt has been made (for now) to parallelize it.
      </p>
    </div>
    <div class="faqq"> 
      Why GMM++ defines add, mult, scale procedures instead of
      overloaded operator +, *, - ...
    </div>
    <div class="faqr"> 
      This is a big discussion. The choice in GMM++ is to be able to
      have reasonably optimized operations in all mixed cases
      (operations mixing sparse, skyline and dense matrices and
      vectors), to accept various format (for instance for spares
      matrices) and finally to be able to interface already existing
      matrix and vector types. This seems not to be possible with
      overloaded operator (since it is not possible to overload
      operator = outside of a class), and in our opinion this does not
      offer a big advantage (but a big complexity !).
    </div>

    <div class="faqq"> 
      How can I interface a CSC or CSR matrix coming form a Fortran or C
      code ?
    </div>
    <div class="faqr">
      An interface exists in the file "gmm_interface.h". The usage is
      <pre>
      gmm::csc_matrix_ref&lt;PT1, PT2, PT3, shift&gt;  M1(pr, ir, jc, nrows, ncols)
      gmm::csr_matrix_ref&lt;PT1, PT2, PT3, shift&gt;  M2(pr, ir, jc, nrows, ncols)
      </pre>
      where PT1 is the type of pointer to the data (double * for instance),
      PT1 and PT2 the types of pointers to the indices (int * for instance)
      and shift is 1 for matrices coming from Fortran codes and 0 for the ones
      coming from C or MATLAB codes. This is a read_only reference.
      If you want to modify your matrix you have to copy it first in a
      writable matrix such as gmm::col_matrix&lt;gmm::rsvector&lt;T&gt; &gt;.
    </div>

  </div>
<?php include("footer.inc") ?>

