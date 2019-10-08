% FUNCTION [...] = gf_model_get(model M, [operation [, args]])
%
%   Get information from a model object.
%   
%
%   * b = gf_model_get(model M, 'is_complex')
%   Return 0 is the model is real, 1 if it is complex.
%
%   * T = gf_model_get(model M, 'nbdof')
%   Return the total number of degrees of freedom of the model.
%
%   * dt = gf_model_get(model M, 'get time step')
%   Gives the value of the time step.
%
%   * t = gf_model_get(model M, 'get time')
%   Give the value of the data `t` corresponding to the current time.
%   
%
%   * T = gf_model_get(model M, 'tangent_matrix')
%   Return the tangent matrix stored in the model .
%
%   * gf_model_get(model M, 'rhs')
%   Return the right hand side of the tangent problem.
%
%   * gf_model_get(model M, 'brick term rhs', int ind_brick[, int ind_term, int sym, int ind_iter])
%   Gives the access to the part of the right hand side of a term
%   of a particular nonlinear brick. Does not account of the eventual
%   time dispatcher. An assembly of the rhs has to be done first.
%   `ind_brick` is the brick index. `ind_term` is the index of the
%   term inside the brick (default value : 1).
%   `sym` is to access to the second right hand side of for symmetric
%   terms acting on two different variables (default is 0).
%   `ind_iter` is the iteration number when time dispatchers are
%   used (default is 1).
%   
%
%   * z = gf_model_get(model M, 'memsize')
%   Return a rough approximation of the amount of memory (in bytes) used by
%   the model.
%
%   * gf_model_get(model M, 'variable list')
%   print to the output the list of variables and constants of the model.
%
%   * gf_model_get(model M, 'brick list')
%   print to the output the list of bricks of the model.
%
%   * gf_model_get(model M, 'list residuals')
%   print to the output the residuals corresponding to all terms
%   included in the model.
%
%   * V = gf_model_get(model M, 'variable', string name)
%   Gives the value of a variable or data.
%
%   * V = gf_model_get(model M, 'interpolation', string expr, {mesh_fem mf | mesh_imd mimd | vec pts,  mesh m}[, int region[, int extrapolation[, int rg_source]]])
%   Interpolate a certain expression with respect to the mesh_fem `mf`
%   or the mesh_im_data `mimd` or the set of points `pts` on mesh `m`.
%   The expression has to be valid according to the high-level generic
%   assembly language possibly including references to the variables
%   and data of the model.
%   
%   The options `extrapolation` and `rg_source` are specific to
%   interpolations with respect to a set of points `pts`.
%
%   * V = gf_model_get(model M, 'local_projection', mesh_im mim, string expr, mesh_fem mf[, int region])
%   Make an elementwise L2 projection of an expression with respect
%   to the mesh_fem `mf`. This mesh_fem has to be
%   a discontinuous one.
%   The expression has to be valid according to the high-level generic
%   assembly language possibly including references to the variables
%   and data of the model.
%
%   * mf = gf_model_get(model M, 'mesh fem of variable', string name)
%   Gives access to the `mesh_fem` of a variable or data.
%
%   * name = gf_model_get(model M, 'mult varname Dirichlet', int ind_brick)
%   Gives the name of the multiplier variable for a Dirichlet brick.
%   If the brick is not a Dirichlet condition with multiplier brick,
%   this function has an undefined behavior
%
%   * I = gf_model_get(model M, 'interval of variable', string varname)
%   Gives the interval of the variable `varname` in the linear system of
%   the model.
%
%   * V = gf_model_get(model M, 'from variables')
%   Return the vector of all the degrees of freedom of the model consisting
%   of the concatenation of the variables of the model (useful
%   to solve your problem with you own solver).
%
%   * gf_model_get(model M, 'assembly'[, string option])
%   Assembly of the tangent system taking into account the terms
%   from all bricks. `option`, if specified, should be 'build_all',
%   'build_rhs', 'build_matrix'.
%   The default is to build the whole
%   tangent linear system (matrix and rhs). This function is useful
%   to solve your problem with you own solver.
%
%   * {nbit, converged} = gf_model_get(model M, 'solve'[, ...])
%   Run the standard getfem solver.
%   
%   Note that you should be able to use your own solver if you want
%   (it is possible to obtain the tangent matrix and its right hand
%   side with the gf_model_get(model M, 'tangent matrix') etc.).
%   
%   Various options can be specified:
%   
%   - 'noisy' or 'very_noisy'
%   the solver will display some information showing the progress
%   (residual values etc.).
%   - 'max_iter', int NIT
%   set the maximum iterations numbers.
%   - 'max_res', @float RES
%   set the target residual value.
%   - 'diverged_res', @float RES
%   set the threshold value of the residual beyond which the iterative
%   method is considered to diverge (default is 1e200).
%   - 'lsolver', string SOLVER_NAME
%   select explicitely the solver used for the linear systems (the
%   default value is 'auto', which lets getfem choose itself).
%   Possible values are 'superlu', 'mumps' (if supported),
%   'cg/ildlt', 'gmres/ilu' and 'gmres/ilut'.
%   - 'lsearch', string LINE_SEARCH_NAME
%   select explicitely the line search method used for the linear systems (the
%   default value is 'default').
%   Possible values are 'simplest', 'systematic', 'quadratic' or 'basic'.
%   
%   Return the number of iterations, if an iterative method is used.
%   
%   Note that it is possible to disable some variables
%   (see gf_model_set(model M, 'disable variable') ) in order to
%   solve the problem only with respect to a subset of variables (the
%   disabled variables are then considered as data) for instance to
%   replace the global Newton strategy with a fixed point one.
%   
%   
%
%   * gf_model_get(model M, 'test tangent matrix'[, scalar EPS[, int NB[, scalar scale]]])
%   Test the consistency of the tangent matrix in some random positions
%   and random directions (useful to test newly created bricks).
%   `EPS` is the value of the small parameter for the finite difference
%   computation of the derivative is the random direction (default is 1E-6).
%   `NN` is the number of tests (default is 100). `scale` is a parameter
%   for the random position (default is 1, 0 is an acceptable value) around
%   the current position.
%   Each dof of the random position is chosen in the range
%   [current-scale, current+scale].
%   
%
%   * gf_model_get(model M, 'test tangent matrix term', string varname1, string varname2[, scalar EPS[, int NB[, scalar scale]]])
%   Test the consistency of a part of the tangent matrix in some
%   random positions and random directions
%   (useful to test newly created bricks).
%   The increment is only made on variable `varname2` and tested on the
%   part of the residual corresponding to `varname1`. This means that
%   only the term (`varname1`, `varname2`) of the tangent matrix is tested.
%   `EPS` is the value of the small parameter for the finite difference
%   computation of the derivative is the random direction (default is 1E-6).
%   `NN` is the number of tests (default is 100). `scale` is a parameter
%   for the random position (default is 1, 0 is an acceptable value)
%   around the current position.
%   Each dof of the random position is chosen in the range
%   [current-scale, current+scale].
%   
%
%   * expr = gf_model_get(model M, 'Neumann term', string varname, int region)
%   Gives the assembly string corresponding to the Neumann term of
%   the fem variable `varname` on `region`. It is deduced from the
%   assembly string declared by the model bricks.
%   `region` should be the index of a boundary region
%   on the mesh where `varname` is defined. Care to call this function
%   only after all the volumic bricks have been declared.
%   Complains, if a brick
%   omit to declare an assembly string.
%
%   * V = gf_model_get(model M, 'compute isotropic linearized Von Mises or Tresca', string varname, string dataname_lambda, string dataname_mu, mesh_fem mf_vm[, string version])
%   Compute the Von-Mises stress or the Tresca stress of a field (only
%   valid for isotropic linearized elasticity in 3D). `version` should
%   be  'Von_Mises' or 'Tresca' ('Von_Mises' is the default).
%   Parametrized by Lame coefficients.
%   
%
%   * V = gf_model_get(model M, 'compute isotropic linearized Von Mises pstrain', string varname, string data_E, string data_nu, mesh_fem mf_vm)
%   Compute the Von-Mises stress  of a displacement field for isotropic
%   linearized elasticity in 3D or in 2D with plane strain assumption.
%   Parametrized by Young modulus and Poisson ratio.
%   
%
%   * V = gf_model_get(model M, 'compute isotropic linearized Von Mises pstress', string varname, string data_E, string data_nu, mesh_fem mf_vm)
%   Compute the Von-Mises stress  of a displacement field for isotropic
%   linearized elasticity in 3D or in 2D with plane stress assumption.
%   Parametrized by Young modulus and Poisson ratio.
%   
%
%   * V = gf_model_get(model M, 'compute Von Mises or Tresca', string varname, string lawname, string dataname, mesh_fem mf_vm[, string version])
%   Compute on `mf_vm` the Von-Mises stress or the Tresca stress of a field
%   for nonlinear elasticity in 3D. `lawname` is the constitutive law which
%   could be 'SaintVenant Kirchhoff', 'Mooney Rivlin', 'neo Hookean' or
%   'Ciarlet Geymonat'.
%   `dataname` is a vector of parameters for the constitutive law. Its length
%   depends on the law. It could be a short vector of constant values or a
%   vector field described on a finite element method for variable coefficients.
%   `version` should be  'Von_Mises' or 'Tresca' ('Von_Mises' is the default).
%   
%
%   * V = gf_model_get(model M, 'compute finite strain elasticity Von Mises',  string lawname, string varname, string params, mesh_fem mf_vm[, int region])
%   Compute on `mf_vm` the Von-Mises stress of a field `varname`
%   for nonlinear elasticity in 3D. `lawname` is the constitutive law which
%   should be a valid name. `params` are the parameters law. It could be
%   a short vector of constant values or may depend on data or variables
%   of the model.
%   Uses the high-level generic assembly.
%   
%
%   * V = gf_model_get(model M, 'compute second Piola Kirchhoff tensor', string varname, string lawname, string dataname, mesh_fem mf_sigma)
%   Compute on `mf_sigma` the second Piola Kirchhoff stress tensor of a field
%   for nonlinear elasticity in 3D. `lawname` is the constitutive law which
%   could be 'SaintVenant Kirchhoff', 'Mooney Rivlin', 'neo Hookean' or
%   'Ciarlet Geymonat'.
%   `dataname` is a vector of parameters for the constitutive law. Its length
%   depends on the law. It could be a short vector of constant values or a
%   vector field described on a finite element method for variable
%   coefficients.
%   
%
%   * gf_model_get(model M, 'elastoplasticity next iter', mesh_im mim, string varname, string previous_dep_name, string projname, string datalambda, string datamu, string datathreshold, string datasigma)
%   Used with the old (obsolete) elastoplasticity brick to pass from an
%   iteration to the next one.
%   Compute and save the stress constraints sigma for the next iterations.
%   'mim' is the integration method to use for the computation.
%   'varname' is the main variable of the problem.
%   'previous_dep_name' represents the displacement at the previous time step.
%   'projname' is the type of projection to use. For the moment it could only be 'Von Mises' or 'VM'.
%   'datalambda' and 'datamu' are the Lame coefficients of the material.
%   'datasigma' is a vector which will contain the new stress constraints values.
%
%   * gf_model_get(model M, 'small strain elastoplasticity next iter', mesh_im mim,  string lawname, string unknowns_type [, string varnames, ...] [, string params, ...] [, string theta = '1' [, string dt = 'timestep']] [, int region = -1])
%   Function that allows to pass from a time step to another for the
%   small strain plastic brick. The parameters have to be exactly the
%   same than the one of `add_small_strain_elastoplasticity_brick`,
%   so see the documentation of this function for the explanations.
%   Basically, this brick computes the plastic strain
%   and the plastic multiplier and stores them for the next step.
%   Additionaly, it copies the computed displacement to the data
%   that stores the displacement of the previous time step (typically
%   'u' to 'Previous_u'). It has to be called before any use of
%   `compute_small_strain_elastoplasticity_Von_Mises`.
%   
%
%   * V = gf_model_get(model M, 'small strain elastoplasticity Von Mises', mesh_im mim, mesh_fem mf_vm, string lawname, string unknowns_type [, string varnames, ...] [, string params, ...] [, string theta = '1' [, string dt = 'timestep']] [, int region])
%   This function computes the Von Mises stress field with respect to
%   a small strain elastoplasticity term, approximated on `mf_vm`,
%   and stores the result into `VM`.  All other parameters have to be
%   exactly the same as for `add_small_strain_elastoplasticity_brick`.
%   Remember that `small_strain_elastoplasticity_next_iter` has to be called
%   before any call of this function.
%   
%
%   * V = gf_model_get(model M, 'compute elastoplasticity Von Mises or Tresca', string datasigma, mesh_fem mf_vm[, string version])
%   Compute on `mf_vm` the Von-Mises or the Tresca stress of a field for plasticity and return it into the vector V.
%   `datasigma` is a vector which contains the stress constraints values supported by the mesh.
%   `version` should be  'Von_Mises' or 'Tresca' ('Von_Mises' is the default).
%
%   * V = gf_model_get(model M, 'compute plastic part', mesh_im mim, mesh_fem mf_pl, string varname, string previous_dep_name, string projname, string datalambda, string datamu, string datathreshold, string datasigma)
%   Compute on `mf_pl` the plastic part and return it into the vector V.
%   `datasigma` is a vector which contains the stress constraints values supported by the mesh.
%
%   * gf_model_get(model M, 'finite strain elastoplasticity next iter', mesh_im mim, string lawname, string unknowns_type, [, string varnames, ...] [, string params, ...] [, int region = -1])
%   Function that allows to pass from a time step to another for the
%   finite strain plastic brick. The parameters have to be exactly the
%   same than the one of `add_finite_strain_elastoplasticity_brick`,
%   so see the documentation of this function for the explanations.
%   Basically, this brick computes the plastic strain
%   and the plastic multiplier and stores them for the next step.
%   For the Simo-Miehe law which is currently the only one implemented,
%   this function updates the state variables defined in the last two
%   entries of `varnames`, and resets the plastic multiplier field given
%   as the second entry of `varnames`.
%   
%
%   * V = gf_model_get(model M, 'compute finite strain elastoplasticity Von Mises', mesh_im mim, mesh_fem mf_vm, string lawname, string unknowns_type, [, string varnames, ...] [, string params, ...] [, int region = -1])
%   Compute on `mf_vm` the Von-Mises or the Tresca stress of a field for plasticity and return it into the vector V.
%   The first input parameters ar as in the function 'finite strain elastoplasticity next iter'.
%   
%
%   * V = gf_model_get(model M, 'sliding data group name of large sliding contact brick', int indbrick)
%   Gives the name of the group of variables corresponding to the
%   sliding data for an existing large sliding contact brick.
%
%   * V = gf_model_get(model M, 'displacement group name of large sliding contact brick', int indbrick)
%   Gives the name of the group of variables corresponding to the
%   sliding data for an existing large sliding contact brick.
%
%   * V = gf_model_get(model M, 'transformation name of large sliding contact brick', int indbrick)
%   Gives the name of the group of variables corresponding to the
%   sliding data for an existing large sliding contact brick.
%
%   * V = gf_model_get(model M, 'sliding data group name of Nitsche large sliding contact brick', int indbrick)
%   Gives the name of the group of variables corresponding to the
%   sliding data for an existing large sliding contact brick.
%
%   * V = gf_model_get(model M, 'displacement group name of Nitsche large sliding contact brick', int indbrick)
%   Gives the name of the group of variables corresponding to the
%   sliding data for an existing large sliding contact brick.
%
%   * V = gf_model_get(model M, 'transformation name of Nitsche large sliding contact brick', int indbrick)
%   Gives the name of the group of variables corresponding to the
%   sliding data for an existing large sliding contact brick.
%
%   * M = gf_model_get(model M, 'matrix term', int ind_brick, int ind_term)
%   Gives the matrix term ind_term of the brick ind_brick if it exists
%   
%
%   * s = gf_model_get(model M, 'char')
%   Output a (unique) string representation of the model.
%   
%   This can be used to perform comparisons between two
%   different model objects.
%   This function is to be completed.
%   
%
%   * gf_model_get(model M, 'display')
%   displays a short summary for a model object.
%
%
function [varargout]=gf_model_get(varargin)
  if (nargout),
    [varargout{1:nargout}]=gf_matlab('model_get', varargin{:});
  else
    gf_matlab('model_get', varargin{:});
    if (exist('ans', 'var') == 1), varargout{1}=ans; end;
  end;
% autogenerated mfile;
