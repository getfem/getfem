.. $Id: gasm_high.rst -1   $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: asm, generic assembly

.. _ud-gasm-high:

Compute arbitrary elementary matrices - high-level generic assembly procedures
==============================================================================

This is a work in progress for |gf| 5.0. Not fully working yet.




This section presents the second version of generic assembly which has been implemented in |gf| 5.0. It is a high-level generic assembly in the sense that the language used to describe the assembly to be done is quite close to the weak formulation of boundary value problems of partial differential equations. It mainly has been developped to circumvent the difficulties with the older generic assembly (see :ref:`ud-gasm-low`) for which nonlinear term are quite difficult to take into account. Conversly, an automatic differentiation algorithm is used with this version to simplify the writing of new nonlinear terms.

Differences in execution time between high and low level generic assembly
-------------------------------------------------------------------------
For basic linear assembly terms, the high and low level generic assembly procedures have approximately the same efficiency in term of computational time. Both have been thoroughly optimized. On the one hand, the fact that the high-level generic assembly incorporate a compilation in basic optimized instructions and operates simplifications makes that it can be really faster especially on complex terms. On the other hand, the fact that the low-level generic assembly incorporates a mechanism to precompute on the reference element the linear term for elements with a linear transformation makes that it can be faster on simple linear terms in that situation. But even in that case, the high level generic assembly is sometime faster. Of course, a possibility would be to incorporate the ability to precompute on the reference element the linear term for linear transformations in the high level generic assembly. However, it would be rather complicated due to the high genericity of the language. 


The assembly language
---------------------

A specific language has been developped to describe the weak formulation of boundary value problems. It aims to be close to the structure of a standard weak formulation. The components are the following:

  - A certain number of predefined scalar functions (sin(t), cos(t), pow(t,u), sqrt(t), sqr(t), Heaviside(t) to be described ...). A scalar function can be appile to scalar or vector/matrix/tensor expression. It applies componentwise. For functions having two arguments (pow(t,u), min(t,u) ...) if two non-scalar arguments are passed, the dimension have to be the same. For instance "max([1;2],[0;3])" will return "[0;3]".

  - The constant pi and some special functions: meshdim(u) (for u a fem variable, the dimension of the corresponding mesh), qdim(u) (for u a mesh variable, the dimension of the vector field corresponding to the variable, 1 for a scalar field), Id(n) (indentity nxn matrix).

  - Variable names: A list of variables should be given. The variables are described on a finite element method or can be a simple vector of unknows. For instance "u", "v", "p", "pressure", "electric_field" are valid variable names. A variable name should begin by a letter (case sensitive) or un underscore followed by a letter, a number or an undescore. The name should not begin by "Test\_", "Grad\_" or "Hess\_". The variable name should not correspond to a predefined function (sin, cos, acos ...). In case of conflict, the name is understand as a function name.

  - Test functions corresponding to the variables. It is identified by the prefix "Test\_" followed by the variable name. For instance  "Test_u", "Test_v", "Test_p", "Test_pressure", "Test_electric_field".
  
  - The gradient of a variable or of test functions are identified by "Grad\_" followed by the variable name or by "Test\_" followed itself by the variable name. This is available for fem variables only. For instance "Grad_u", "Grad_v", "Grad_p", "Grad_pressure", "Grad_electric_field" and "Grad_Test_u", "Grad_Test_v", "Grad_Test_p", "Grad_Test_pressure", "Grad_Test_electric_field". The gradient is either a vector for scalar variables or a matrix for vector field variables. In the latter case, the first index corresponds to the vector field dimension and the second one to the index of the partial derivative.

  - The Hessian of a variable or of test functions are identified by "Hess\_" followed by the variable name or by "Test\_" followed itself by the variable name. This is available for fem variables only. For instance "Hess_u", "Hess_v", "Hess_p", "Hess_pressure", "Hess_electric_field" and "Hess_Test_u", "Hess_Test_v", "Hess_Test_p", "Hess_Test_pressure", "Hess_Test_electric_field". The Hessian is either a matrix for scalar variables or a third order tensor for vector field variables. In the latter case, the first index corresponds to the vector field dimension and the two remaining to the indices of partial derivatives.

  - Constant names: A list of constants could be given. The rule are the same as for the variables but no test function can be associated to constants.

  - Constant expressions: numbers, C_PI, mesh_dim,  ...

  - A certain number of operators: "+", "-", "*", "/", ":", ".", ".*", "./", "@", "'".
    
    - "+" and "-" are the standard addition and substraction of scalar, vector, matrix or tensors.

    - "*" stands for the scalar, matrix-vector, matrix-matrix or (fourth order tensor)-matrix multiplication.

    - "/" stands for the division by a scalar.

    - "." stands for the scalar product of vectors, or more generally to the reduction of a tensor with respect to the last index with a vector. Note that "*" and "." are equivalent for matrix-vector multiplication.

    - ":" stands for the the Frobenius product of matrices or more generally to the reduction of a tensor with respect to the two last indices with a matrix. Note that "*" and ":" are equivalent for (fourth order tensor)-matrix multiplication.

    - ".*" stands for the multiplication of two vectors/matrix/tensor componentwise.

    - "./" stands for the division of two vectors/matrix/tensor componentwise.

    - "@" stands for the tensor product.

    - "'" stands for the transpose of a matrix or line view of a vector.

  - Parentheses can be used to gather expressions. For instance "(1+2)*4" or "(u+v)*Test_u" are correct. 

  - Constant matrices. The syntax is close to Matlab one but the white space cannot be used to separate components (instead only the comma can be used). For instance "[1,2;3,4]" denotes a 2x2 matrix. Note that the components can be the result of an expression. For instance "[1+2,2;3,4+a]" is correct if "a" is a declared constant or variable.

  - Constant vectors. A constant vector is a constant matrix with only one column. For instance  "[1;2;3;4]" is a constant vector of size four.

  - Constant fourth order tensors. It is also possible to give a constant fourth order tensor with a syntax close to the one for the matrices. Additionnal separators for the two supplementary dimensions are ',,' and ';;'. For instance "[1,1;1,2,,1,1;1,2;;1,1;1,2,,1,1;1,2]" is a 2x2x2x2 valid tensor. Note that constant fourth order tensors can also be obtained by the tensorial product of two constant matrices. Note also that constant third order tensors are not supported.



  - The access to a component of a vector/matrix/tensor can be done by following a term by a left parenthesis, the list of components and a right parenthesis. For instance "[1,1,2](3)" is correct and will return "2". Note that indices are assumed to begin by 1 for the compatibility with matlab (even in C++ and with the python interface). The expressions "[1,1;2,3](2,2)" and "Grad_u(2,2)" are also correct provided that "u" is a vector valued declared variable. Note that the components can be the result of a constant computation. For instance "[1,1;2,3](1+1,a)" is correct provided that "a" is a declared constant but not if it is declared as a variable. A colon can replace the value of an index in a matlab like syntax for instance to access to a line or a column of a matrix. "[1,1;2,3](1,:)" denotes the first line of the matrix "[1,1;2,3]". It can also be used for a fourth order tensor.

  - Trace operator (can be applied to test functions).

  - Print command : for debugging, print the tensor and pass it unchanged. "Grad_u.Print(Grad_Test_u)" will have the same effect as "Grad_u.Grad_Test_u" but printing the tensor "Grad_Test_u" for each Gauss point of each element. Note that constant terms are printed only once at the begining of the assembly. Note also that the expression could be derived so that the derivative of the term may be printed instead of the term itself.

  - A certain number of predefined nonlinear operator (tr, norm, Idmat, det, target_dim(u) ...). Cannot be applied to test functions.


Classical examples ("u" the variable name, "a" a coefficient):

  - Simple Laplace operator for a scalar field "Grad_u.Grad_Test_u"
  - Simple componentwize Laplace operator for a vector field "Grad_u:Grad_Test_u"
  - Laplace operator with a scalar coefficient for a scalar field "a*Grad_u.Grad_Test_u"
  - Laplace operator with a matrix coefficient for a vector field "(a*Grad_u).Grad_Test_u" or "([2,1;1,4]*Grad_u).Grad_Test_u" for a constant coefficient in 2D.
  - Linear isotropic elasticity "(lambda*Trace(Grad_u)*Id(mesh_dim(u)) + mu*(Grad_u+Grad_u')):Grad_Test_u" ou "lambda*Trace(Grad_u)*Trace(Grad_Test_u) + mu*(Grad_u + Grad_u'):Grad_Test_u"
  

The assembly string is transformed in an assembly tree by a set of function in :file:`src/getfem_generic_assembly.cc`. The process has 4 steps:

 - Lexical analysis with the procedure `ga_get_token`.

 - Syntax analysis and transformation into a syntax tree by `ga_read_string`.

 - Semantic analysis, simplification (pre-computation) of constant expressions and enrichment of the tree.

 - Symbolic (automatic) differentiation.

 - Compilation in a sequence of instructions with optimisation not to evaluate several time the same expression.

 - Execution of the sequence of instructions and assembly.

These steps are perfromed only once at the begining of the assembly. The final tree is evaluated on each Gauss point of each element ...



Main error messages of the generic assembly:
...
