% FUNCTION [...] = gf_compute(mesh_fem MF, vec U, [operation [, args]])
%
%   
%   Various computations involving the solution U to a finite element problem.
%   
%
%   * n = gf_compute(mesh_fem MF, vec U, 'L2 norm', mesh_im mim[, mat CVids])
%   Compute the L2 norm of the (real or complex) field `U`.
%   
%   If `CVids` is given, the norm will be computed only on the listed
%   elements.
%
%   * n = gf_compute(mesh_fem MF, vec U, 'L2 dist', mesh_im mim, mesh_fem mf2, vec U2[, mat CVids])
%   Compute the L2 distance between `U` and `U2`.
%   
%   If `CVids` is given, the norm will be computed only on the listed
%   elements.
%
%   * n = gf_compute(mesh_fem MF, vec U, 'H1 semi norm', mesh_im mim[, mat CVids])
%   Compute the L2 norm of grad(`U`).
%   
%   If `CVids` is given, the norm will be computed only on the listed
%   elements.
%
%   * n = gf_compute(mesh_fem MF, vec U, 'H1 semi dist', mesh_im mim, mesh_fem mf2, vec U2[, mat CVids])
%   Compute the semi H1 distance between `U` and `U2`.
%   
%   If `CVids` is given, the norm will be computed only on the listed
%   elements.
%
%   * n = gf_compute(mesh_fem MF, vec U, 'H1 norm', mesh_im mim[, mat CVids])
%   Compute the H1 norm of `U`.
%   
%   If `CVids` is given, the norm will be computed only on the listed
%   elements.
%
%   * n = gf_compute(mesh_fem MF, vec U, 'H2 semi norm', mesh_im mim[, mat CVids])
%   Compute the L2 norm of D^2(`U`).
%   
%   If `CVids` is given, the norm will be computed only on the listed
%   elements.
%
%   * n = gf_compute(mesh_fem MF, vec U, 'H2 norm', mesh_im mim[, mat CVids])
%   Compute the H2 norm of `U`.
%   
%   If `CVids` is given, the norm will be computed only on the listed
%   elements.
%
%   * DU = gf_compute(mesh_fem MF, vec U, 'gradient', mesh_fem mf_du)
%   Compute the gradient of the field `U` defined on mesh_fem `mf_du`.
%   
%   The gradient is interpolated on the mesh_fem `mf_du`, and returned in
%   `DU`. For example, if `U` is defined on a P2 mesh_fem, `DU` should be
%   evaluated on a P1-discontinuous mesh_fem. `mf` and `mf_du` should
%   share the same mesh.
%   
%   `U` may have any number of dimensions (i.e. this function is not
%   restricted to the gradient of scalar fields, but may also be used
%   for tensor fields). However the last dimension of `U` has to be
%   equal to the number of dof of `mf`. For example, if `U` is a
%   [3x3xNmf] array (where Nmf is the number of dof of `mf`), `DU` will
%   be a [Nx3x3[xQ]xNmf_du] array, where N is the dimension of the mesh,
%   Nmf_du is the number of dof of `mf_du`, and the optional Q dimension
%   is inserted if `Qdim_mf != Qdim_mf_du`, where Qdim_mf is the Qdim of
%   `mf` and Qdim_mf_du is the Qdim of `mf_du`.
%
%   * HU = gf_compute(mesh_fem MF, vec U, 'hessian', mesh_fem mf_h)
%   Compute the hessian of the field `U` defined on mesh_fem `mf_h`.
%   
%   See also gf_compute('gradient', mesh_fem mf_du).
%
%   * UP = gf_compute(mesh_fem MF, vec U, 'eval on triangulated surface', int Nrefine, [vec CVLIST])
%   [OBSOLETE FUNCTION! will be removed in a future release]
%   Utility function designed for 2D triangular meshes : returns a list
%   of triangles coordinates with interpolated U values. This can be
%   used for the accurate visualization of data defined on a
%   discontinous high order element. On output, the six first rows of UP
%   contains the triangle coordinates, and the others rows contain the
%   interpolated values of U (one for each triangle vertex) CVLIST may
%   indicate the list of convex number that should be consider, if not
%   used then all the mesh convexes will be used. U should be a row
%   vector.
%   
%
%   * Ui = gf_compute(mesh_fem MF, vec U, 'interpolate on', {mesh_fem mfi | slice sli | vec pts})
%   Interpolate a field on another mesh_fem or a slice or a list of points.
%   
%   - Interpolation on another mesh_fem `mfi`:
%   `mfi` has to be Lagrangian. If `mf` and `mfi` share the same
%   mesh object, the interpolation will be much faster.
%   - Interpolation on a slice `sli`:
%   this is similar to interpolation on a refined P1-discontinuous
%   mesh, but it is much faster. This can also be used with
%   gf_slice('points') to obtain field values at a given set of
%   points.
%   - Interpolation on a set of points `pts`
%   
%   See also gf_asm('interpolation matrix')
%   
%
%   * Ue = gf_compute(mesh_fem MF, vec U, 'extrapolate on', mesh_fem mfe)
%   Extrapolate a field on another mesh_fem.
%   
%   If the mesh of `mfe` is stricly included in the mesh of `mf`, this
%   function does stricly the same job as gf_compute('interpolate_on').
%   However, if the mesh of `mfe` is not exactly included in `mf`
%   (imagine interpolation between a curved refined mesh and a coarse
%   mesh), then values which are outside `mf` will be
%   extrapolated.
%   
%   See also gf_asm('extrapolation matrix')
%
%   * E = gf_compute(mesh_fem MF, vec U, 'error estimate', mesh_im mim)
%   Compute an a posteriori error estimate.
%   
%   Currently there is only one which is available: for each convex,
%   the jump of the normal derivative is integrated on its faces.
%
%   * E = gf_compute(mesh_fem MF, vec U, 'error estimate nitsche', mesh_im mim, int GAMMAC, int GAMMAN, scalar lambda_, scalar mu_, scalar gamma0, scalar f_coeff, scalar vertical_force)
%   Compute an a posteriori error estimate in the case of Nitsche method.
%   
%   Currently there is only one which is available: for each convex,
%   the jump of the normal derivative is integrated on its faces.
%
%   * gf_compute(mesh_fem MF, vec U, 'convect', mesh_fem mf_v, vec V, scalar dt, int nt[, string option[, vec per_min, vec per_max]])
%   Compute a convection of `U` with regards to a steady state velocity
%   field `V` with a Characteristic-Galerkin method. The result is returned
%   in-place in `U`.
%   This method is restricted to pure Lagrange fems for U. `mf_v` should
%   represent a continuous finite element method. `dt` is the integration time
%   and `nt` is the number of integration step on the caracteristics. `option`
%   is an option for the part of the boundary where there is a re-entrant
%   convection.
%   `option = 'extrapolation'` for an extrapolation on the nearest element,
%   `option = 'unchanged'` for a constant value on that boundary or
%   `option = 'periodicity'` for a peridiodic boundary. For this latter option
%   the two vectors per_min, per_max has to be given and represent the limits
%   of the periodic domain (on components where per_max[k] < per_min[k]
%   no operation is done).
%   This method is rather dissipative, but stable.
%   
%
%   * [U2[,MF2,[,X[,Y[,Z]]]]] = gf_compute(mesh_fem MF, vec U, 'interpolate on Q1 grid', {'regular h', hxyz | 'regular N', Nxyz | X[,Y[,Z]]})
%   
%   Creates a cartesian Q1 mesh fem and interpolates U on it. The
%   returned field U2 is organized in a matrix such that in can be drawn
%   via the MATLAB command 'pcolor'. The first dimension is the Qdim of
%   MF (i.e.  1 if U is a scalar field)
%   
%   example (mf_u is a 2D mesh_fem):
%   >> Uq=gf_compute(mf_u, U, 'interpolate on Q1 grid', 'regular h', [.05, .05]);
%   >> pcolor(squeeze(Uq(1,:,:)));
%   
%
%
function [varargout]=gf_compute(varargin)

  if (nargin>=3 & strcmpi(varargin{3}, 'interpolate on Q1 grid')),
    [varargout{1:nargout}]=gf_compute_Q1grid_interp(varargin{[1 2 4:nargin]});
    return;
  end;
  
  if (nargout),
    [varargout{1:nargout}]=gf_matlab('compute', varargin{:});
  else
    gf_matlab('compute', varargin{:});
    if (exist('ans', 'var') == 1), varargout{1}=ans; end;
  end;
% autogenerated mfile;
