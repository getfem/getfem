% FUNCTION [...] = gf_mesh_fem_get(mesh_fem MF, [operation [, args]])
%
%   General function for inquiry about mesh_fem objects.
%   
%
%   * n = gf_mesh_fem_get(mesh_fem MF, 'nbdof')
%   Return the number of degrees of freedom (dof) of the mesh_fem.
%
%   * n = gf_mesh_fem_get(mesh_fem MF, 'nb basic dof')
%   Return the number of basic degrees of freedom (dof) of the mesh_fem.
%
%   * DOF = gf_mesh_fem_get(mesh_fem MF, 'dof from cv',mat CVids)
%   Deprecated function. Use gf_mesh_fem_get(mesh_fem MF, 'basic dof from cv') instead.
%
%   * DOF = gf_mesh_fem_get(mesh_fem MF, 'basic dof from cv',mat CVids)
%   Return the dof of the convexes listed in `CVids`.
%   
%   WARNING: the Degree of Freedom might be returned in ANY order, do
%   not use this function in your assembly routines. Use 'basic dof from cvid'
%   instead, if you want to be able to map a convex number with its
%   associated degrees of freedom.
%   
%   One can also get the list of basic dof on a set on convex faces, by
%   indicating on the second row of `CVids` the faces numbers (with
%   respect to the convex number on the first row).
%
%   * {DOFs, IDx} = gf_mesh_fem_get(mesh_fem MF, 'dof from cvid'[, mat CVids])
%   Deprecated function. Use gf_mesh_fem_get(mesh_fem MF, 'basic dof from cvid') instead.
%   
%
%   * {DOFs, IDx} = gf_mesh_fem_get(mesh_fem MF, 'basic dof from cvid'[, mat CVids])
%   Return the degrees of freedom attached to each convex of the mesh.
%   
%   If `CVids` is omitted, all the convexes will be considered (equivalent
%   to `CVids = 1 ... gf_mesh_get(mesh M, 'max cvid')`).
%   
%   `IDx` is a row vector, `length(IDx) = length(CVids)+1`.
%   `DOFs` is a row vector containing the concatenated list
%   of dof of each convex in `CVids`. Each entry of `IDx` is the position
%   of the corresponding convex point list in `DOFs`. Hence, for example,
%   the list of points of the second convex is DOFs(IDx(2):IDx(3)-1).
%   
%   If `CVids` contains convex #id which do not exist in the mesh, their
%   point list will be empty.
%
%   * gf_mesh_fem_get(mesh_fem MF, 'non conformal dof'[, mat CVids])
%   Deprecated function. Use gf_mesh_fem_get(mesh_fem MF, 'non conformal basic dof') instead.
%   
%
%   * gf_mesh_fem_get(mesh_fem MF, 'non conformal basic dof'[, mat CVids])
%   Return partially linked degrees of freedom.
%   
%   Return the basic dof located on the border of a convex and which belong
%   to only one convex, except the ones which are located on the border
%   of the mesh.  For example, if the convex 'a' and 'b' share a common
%   face, 'a' has a P1 FEM, and 'b' has a P2 FEM, then the basic dof on the
%   middle of the face will be returned by this function (this can be
%   useful when searching the interfaces between classical FEM and
%   hierarchical FEM).
%
%   * gf_mesh_fem_get(mesh_fem MF, 'qdim')
%   Return the dimension Q of the field interpolated by the mesh_fem.
%   
%   By default, Q=1 (scalar field). This has an impact on the dof numbering.
%
%   * {FEMs, CV2F} = gf_mesh_fem_get(mesh_fem MF, 'fem'[, mat CVids])
%   Return a list of FEM used by the mesh_fem.
%   
%   `FEMs` is an array of all fem objects found in the convexes
%   given in `CVids`. If `CV2F` was supplied as an output argument,
%   it contains, for each convex listed in `CVids`, the index of its
%   correspounding FEM in `FEMs`.
%   
%   Convexes which are not part of the mesh, or convexes which do not
%   have any FEM have their correspounding entry in `CV2F` set to -1.
%   
%   Example::
%   
%   cvid=gf_mesh_get(mf,'cvid');
%   [f,c2f]=gf_mesh_fem_get(mf, 'fem');
%   for i=1:size(f), sf{i}=gf_fem_get('char',f(i)); end;
%   for i=1:size(c2f),
%   disp(sprintf('the fem of convex %d is %s',...
%   cvid(i),sf{i}));
%   end;
%   
%
%   * CVs = gf_mesh_fem_get(mesh_fem MF, 'convex_index')
%   Return the list of convexes who have a FEM.
%
%   * bB = gf_mesh_fem_get(mesh_fem MF, 'is_lagrangian'[, mat CVids])
%   Test if the mesh_fem is Lagrangian.
%   
%   Lagrangian means that each base function Phi[i] is such that
%   Phi[i](P[j]) = delta(i,j), where P[j] is the dof location of
%   the jth base function, and delta(i,j) = 1 if i==j, else 0.
%   
%   If `CVids` is omitted, it returns 1 if all convexes in the mesh
%   are Lagrangian. If `CVids` is used, it returns the convex indices
%   (with respect to `CVids`) which are Lagrangian.
%
%   * bB = gf_mesh_fem_get(mesh_fem MF, 'is_equivalent'[, mat CVids])
%   Test if the mesh_fem is equivalent.
%   
%   See gf_mesh_fem_get(mesh_fem MF, 'is_lagrangian')
%
%   * bB = gf_mesh_fem_get(mesh_fem MF, 'is_polynomial'[, mat CVids])
%   Test if all base functions are polynomials.
%   
%   See gf_mesh_fem_get(mesh_fem MF, 'is_lagrangian')
%
%   * bB = gf_mesh_fem_get(mesh_fem MF, 'is_reduced')
%   Return 1 if the optional reduction matrix is applied to the dofs.
%
%   * bB = gf_mesh_fem_get(mesh_fem MF, 'reduction matrix')
%   Return the optional reduction matrix.
%
%   * bB = gf_mesh_fem_get(mesh_fem MF, 'extension matrix')
%   Return the optional extension matrix.
%
%   * Vr = gf_mesh_fem_get(mesh_fem MF, 'reduce vector', vec V)
%   Multiply the provided vector V with the extension matrix of the mesh_fem.
%
%   * Ve = gf_mesh_fem_get(mesh_fem MF, 'extend vector', vec V)
%   Multiply the provided vector V with the reduction matrix of the mesh_fem.
%
%   * DOFs = gf_mesh_fem_get(mesh_fem MF, 'basic dof on region',mat Rs)
%   Return the list of basic dof (before the optional reduction) lying on one
%   of the mesh regions listed in `Rs`.
%   
%   More precisely, this function returns the basic dof whose support is
%   non-null on one of regions whose #ids are listed in `Rs` (note
%   that for boundary regions, some dof nodes may not lie exactly
%   on the boundary, for example the dof of Pk(n,0) lies on the center
%   of the convex, but the base function in not null on the convex
%   border).
%
%   * DOFs = gf_mesh_fem_get(mesh_fem MF, 'dof on region',mat Rs)
%   Return the list of dof (after the optional reduction) lying on one
%   of the mesh regions listed in `Rs`.
%   
%   More precisely, this function returns the basic dof whose support is
%   non-null on one of regions whose #ids are listed in `Rs` (note
%   that for boundary regions, some dof nodes may not lie exactly
%   on the boundary, for example the dof of Pk(n,0) lies on the center
%   of the convex, but the base function in not null on the convex
%   border).
%   
%   For a reduced mesh_fem
%   a dof is lying on a region if its potential corresponding shape
%   function is nonzero on this region. The extension matrix is used
%   to make the correspondance between basic and reduced dofs.
%
%   * DOFpts = gf_mesh_fem_get(mesh_fem MF, 'dof nodes'[, mat DOFids])
%   Deprecated function. Use gf_mesh_fem_get(mesh_fem MF, 'basic dof nodes') instead.
%
%   * DOFpts = gf_mesh_fem_get(mesh_fem MF, 'basic dof nodes'[, mat DOFids])
%   Get location of basic degrees of freedom.
%   
%   Return the list of interpolation points for the specified
%   dof #IDs in `DOFids` (if `DOFids` is omitted, all basic dof are
%   considered).
%
%   * DOFP = gf_mesh_fem_get(mesh_fem MF, 'dof partition')
%   Get the 'dof_partition' array.
%   
%   Return the array which associates an integer (the partition number)
%   to each convex of the mesh_fem. By default, it is an all-zero array.
%   The degrees of freedom of each convex of the mesh_fem are connected
%   only to the dof of neighbouring convexes which have the same
%   partition number, hence it is possible to create partially
%   discontinuous mesh_fem very easily.
%
%   * gf_mesh_fem_get(mesh_fem MF, 'save',string filename[, string opt])
%   Save a mesh_fem in a text file (and optionaly its linked mesh object
%   if `opt` is the string 'with_mesh').
%
%   * gf_mesh_fem_get(mesh_fem MF, 'char'[, string opt])
%   Output a string description of the mesh_fem.
%   
%   By default, it does not include the description of the linked mesh
%   object, except if `opt` is 'with_mesh'.
%
%   * gf_mesh_fem_get(mesh_fem MF, 'display')
%   displays a short summary for a mesh_fem object.
%
%   * m = gf_mesh_fem_get(mesh_fem MF, 'linked mesh')
%   Return a reference to the mesh object linked to `mf`.
%
%   * m = gf_mesh_fem_get(mesh_fem MF, 'mesh')
%   Return a reference to the mesh object linked to `mf`.
%   (identical to gf_mesh_get(mesh M, 'linked mesh'))
%
%   * gf_mesh_fem_get(mesh_fem MF, 'export to vtk',string filename, ... ['ascii'], U, 'name'...)
%   Export a mesh_fem and some fields to a vtk file.
%   
%   The FEM and geometric transformations will be mapped to order 1
%   or 2 isoparametric Pk (or Qk) FEMs (as VTK does not handle higher
%   order elements). If you need to represent high-order FEMs or
%   high-order geometric transformations, you should consider
%   gf_slice_get(slice S, 'export to vtk').
%
%   * gf_mesh_fem_get(mesh_fem MF, 'export to dx',string filename, ...['as', string mesh_name][,'edges']['serie',string serie_name][,'ascii'][,'append'], U, 'name'...)
%   Export a mesh_fem and some fields to an OpenDX file.
%   
%   This function will fail if the mesh_fem mixes different convex types
%   (i.e. quads and triangles), or if OpenDX does not handle a specific
%   element type (i.e. prism connections are not known by OpenDX).
%   
%   The FEM will be mapped to order 1 Pk (or Qk) FEMs. If you need to
%   represent high-order FEMs or high-order geometric transformations,
%   you should consider gf_slice_get(slice S, 'export to dx').
%
%   * gf_mesh_fem_get(mesh_fem MF, 'export to pos',string filename[, string name][[,mesh_fem mf1], mat U1, string nameU1[[,mesh_fem mf2], mat U2, string nameU2,...]])
%   Export a mesh_fem and some fields to a pos file.
%   
%   The FEM and geometric transformations will be mapped to order 1
%   isoparametric Pk (or Qk) FEMs (as GMSH does not handle higher
%   order elements).
%
%   * gf_mesh_fem_get(mesh_fem MF, 'dof_from_im',mesh_im mim[, int p])
%   Return a selection of dof who contribute significantly to the
%   mass-matrix that would be computed with `mf` and the integration
%   method `mim`.
%   
%   `p` represents the dimension on what the integration method
%   operates (default `p = mesh dimension`).
%   
%   IMPORTANT: you still have to set a valid integration method on
%   the convexes which are not crosses by the levelset!
%
%   * U = gf_mesh_fem_get(mesh_fem MF, 'interpolate_convex_data',mat Ucv)
%   
%   Interpolate data given on each convex of the mesh to the mesh_fem dof.
%   The mesh_fem has to be lagrangian, and should be discontinuous (typically
%   a FEM_PK(N,0) or FEM_QK(N,0) should be used).
%   
%   The last dimension of the input vector Ucv should have
%   gf_mesh_get(mesh M, 'max cvid') elements.
%   
%   Example of use: gf_mesh_fem_get(mesh_fem MF, 'interpolate_convex_data', gf_mesh_get(mesh M, 'quality'))
%
%   * z = gf_mesh_fem_get(mesh_fem MF, 'memsize')
%   Return the amount of memory (in bytes) used by the mesh_fem object.
%   
%   The result does not take into account the linked mesh object.
%
%   * gf_mesh_fem_get(mesh_fem MF, 'has_linked_mesh_levelset')
%   Is a mesh_fem_level_set or not.
%
%   * gf_mesh_fem_get(mesh_fem MF, 'linked_mesh_levelset')
%   if it is a mesh_fem_level_set gives the linked mesh_level_set.
%
%   * U = gf_mesh_fem_get(mesh_fem MF, 'eval', expr [, DOFLST])
%   
%   Call gf_mesh_fem_get_eval. This function interpolates an expression on a
%   lagrangian mesh_fem (for all dof except if DOFLST is specified).
%   The expression can be a
%   numeric constant, or a cell array containing numeric constants, string
%   expressions or function handles. For example::
%   
%   U1=gf_mesh_fem_get(mf,'eval',1)
%   U2=gf_mesh_fem_get(mf,'eval',[1;0]) % output has two rows
%   U3=gf_mesh_fem_get(mf,'eval',[1 0]) % output has one row, only valid if qdim(mf)==2
%   U4=gf_mesh_fem_get(mf,'eval',{'x';'y.*z';4;@myfunctionofxyz})
%   
%   
%
%
function [varargout]=gf_mesh_fem_get(varargin)

  if (nargin>=2 & strcmpi(varargin{2},'eval')),
    [varargout{1:nargout}]=gf_mesh_fem_get_eval(varargin{[1 3:nargin]}); return;
  end;
  
  if (nargout),
    [varargout{1:nargout}]=gf_matlab('mesh_fem_get', varargin{:});
  else
    gf_matlab('mesh_fem_get', varargin{:});
    if (exist('ans', 'var') == 1), varargout{1}=ans; end;
  end;
% autogenerated mfile;
