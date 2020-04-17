% FUNCTION [...] = gf_slice_get(slice S, [operation [, args]])
%
%   General function for querying information about slice objects.
%   
%
%   * d = gf_slice_get(slice S, 'dim')
%   Return the dimension of the slice (2 for a 2D mesh, etc..).
%
%   * a = gf_slice_get(slice S, 'area')
%   Return the area of the slice.
%
%   * CVids = gf_slice_get(slice S, 'cvs')
%   Return the list of convexes of the original mesh contained in the slice.
%
%   * n = gf_slice_get(slice S, 'nbpts')
%   Return the number of points in the slice.
%
%   * ns = gf_slice_get(slice S, 'nbsplxs'[, int dim])
%   Return the number of simplexes in the slice.
%   
%   Since the slice may contain points (simplexes of dim 0), segments
%   (simplexes of dimension 1), triangles etc., the result is a vector
%   of size gf_slice_get(slice S, 'dim')+1, except if the optional argument `dim`
%   is used.
%
%   * P = gf_slice_get(slice S, 'pts')
%   Return the list of point coordinates.
%
%   * {S, CV2S} = gf_slice_get(slice S, 'splxs',int dim)
%   Return the list of simplexes of dimension `dim`.
%   
%   On output, S has 'dim+1' rows, each column contains the point
%   numbers of a simplex.  The vector `CV2S` can be used to find the
%   list of simplexes for any convex stored in the slice. For example
%   'S(:,CV2S(4):CV2S(5)-1)'
%   gives the list of simplexes for the fourth convex.
%
%   * {P, E1, E2} = gf_slice_get(slice S, 'edges')
%   Return the edges of the linked mesh contained in the slice.
%   
%   `P` contains the list of all edge vertices, `E1` contains
%   the indices of each mesh edge in `P`, and `E2` contains the
%   indices of each "edges" which is on the border of the slice.
%   This function is useless except for post-processing purposes.
%
%   * Usl = gf_slice_get(slice S, 'interpolate_convex_data', mat Ucv)
%   Interpolate data given on each convex of the mesh to the slice nodes.
%   
%   The input array `Ucv` may have any number of dimensions, but its
%   last dimension should be equal to gf_mesh_get(mesh M, 'max cvid').
%   
%   Example of use: gf_slice_get(slice S, 'interpolate_convex_data', gf_mesh_get(mesh M, 'quality')).
%
%   * m = gf_slice_get(slice S, 'linked mesh')
%   Return the mesh on which the slice was taken.
%
%   * m = gf_slice_get(slice S, 'mesh')
%   Return the mesh on which the slice was taken
%   (identical to 'linked mesh')
%
%   * z = gf_slice_get(slice S, 'memsize')
%   Return the amount of memory (in bytes) used by the slice object.
%
%   * gf_slice_get(slice S, 'export to vtk', string filename, ...)
%   Export a slice to VTK.
%   
%   Following the `filename`, you may use any of the following options:
%   
%   - if 'ascii' is not used, the file will contain binary data
%   (non portable, but fast).
%   - if 'edges' is used, the edges of the original mesh will be
%   written instead of the slice content.
%   
%   More than one dataset may be written, just list them. Each dataset
%   consists of either:
%   
%   - a field interpolated on the slice (scalar, vector or tensor),
%   followed by an optional name.
%   - a mesh_fem and a field, followed by an optional name.
%   
%   Examples:
%   
%   - gf_slice_get(slice S, 'export to vtk', 'test.vtk', Usl, 'first_dataset', mf,
%   U2, 'second_dataset')
%   - gf_slice_get(slice S, 'export to vtk', 'test.vtk', 'ascii', mf,U2)
%   - gf_slice_get(slice S, 'export to vtk', 'test.vtk', 'edges', 'ascii', Uslice)
%
%   * gf_slice_get(slice S, 'export to pov', string filename)
%   Export a the triangles of the slice to POV-RAY.
%
%   * gf_slice_get(slice S, 'export to dx', string filename, ...)
%   Export a slice to OpenDX.
%   
%   Following the `filename`, you may use any of the following
%   options:
%   
%   - if 'ascii' is not used, the file will contain binary data
%   (non portable, but fast).
%   - if 'edges' is used, the edges of the original mesh will be
%   written instead of the slice content.
%   - if 'append' is used, the opendx file will not be overwritten,
%   and the new data will be added at the end of the file.
%   
%   More than one dataset may be written, just list them. Each dataset
%   consists of either:
%   
%   - a field interpolated on the slice (scalar, vector or tensor),
%   followed by an optional name.
%   - a mesh_fem and a field, followed by an optional name.
%
%   * gf_slice_get(slice S, 'export to pos', string filename[, string name][[,mesh_fem mf1], mat U1, string nameU1[[,mesh_fem mf1], mat U2, string nameU2,...])
%   Export a slice to Gmsh.
%   
%   More than one dataset may be written, just list them.
%   Each dataset consists of either:
%   
%   - a field interpolated on the slice (scalar, vector or tensor).
%   - a mesh_fem and a field.
%
%   * s = gf_slice_get(slice S, 'char')
%   Output a (unique) string representation of the slice.
%   
%   This can be used to perform comparisons between two
%   different slice objects.
%   This function is to be completed.
%   
%
%   * gf_slice_get(slice S, 'display')
%   displays a short summary for a slice object.
%
%
function [varargout]=gf_slice_get(varargin)
  if (nargout),
    [varargout{1:nargout}]=gf_octave('slice_get', varargin{:});
  else
    gf_octave('slice_get', varargin{:});
    if (exist('ans', 'var') == 1), varargout{1}=ans; end;
  end;
% autogenerated mfile;
