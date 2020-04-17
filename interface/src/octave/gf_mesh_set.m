% FUNCTION [...] = gf_mesh_set(mesh M, [operation [, args]])
%
%   General function for modification of a mesh object.
%   
%
%   * PIDs = gf_mesh_set(mesh M, 'pts', mat PTS)
%   Replace the coordinates of the mesh points with those given in `PTS`.
%
%   * PIDs = gf_mesh_set(mesh M, 'add point', mat PTS)
%   Insert new points in the mesh and return their #ids.
%   
%   `PTS` should be an ``nxm`` matrix , where ``n`` is the mesh
%   dimension, and ``m`` is the number of points that will be
%   added to the mesh. On output, `PIDs` contains the point #ids
%   of these new points.
%   
%   Remark: if some points are already part of the mesh (with a small
%   tolerance of approximately ``1e-8``), they won't be inserted again,
%   and `PIDs` will contain the previously assigned #ids of these
%   points.
%
%   * gf_mesh_set(mesh M, 'del point', ivec PIDs)
%   Removes one or more points from the mesh.
%   
%   `PIDs` should contain the point #ids, such as the one returned by
%   the 'add point' command.
%
%   * CVIDs = gf_mesh_set(mesh M, 'add convex', geotrans GT, mat PTS)
%   Add a new convex into the mesh.
%   
%   The convex structure (triangle, prism,...) is given by `GT`
%   (obtained with gf_geotrans('...')), and its points are given by
%   the columns of `PTS`. On return, `CVIDs` contains the convex #ids.
%   `PTS` might be a 3-dimensional array in order to insert more than
%   one convex (or a two dimensional array correctly shaped according
%   to Fortran ordering).
%
%   * gf_mesh_set(mesh M, 'del convex', mat CVIDs)
%   Remove one or more convexes from the mesh.
%   
%   `CVIDs` should contain the convexes #ids, such as the ones
%   returned by the 'add convex' command.
%
%   * gf_mesh_set(mesh M, 'del convex of dim', ivec DIMs)
%   Remove all convexes of dimension listed in `DIMs`.
%   
%   For example; ``gf_mesh_set(mesh M, 'del convex of dim', [1,2])`` remove
%   all line segments, triangles and quadrangles.
%
%   * gf_mesh_set(mesh M, 'translate', vec V)
%   Translates each point of the mesh from `V`.
%
%   * gf_mesh_set(mesh M, 'transform', mat T)
%   Applies the matrix `T` to each point of the mesh.
%   
%   Note that `T` is not required to be a ``NxN`` matrix (with
%   ``N = gf_mesh_get(mesh M, 'dim')``). Hence it is possible to transform
%   a 2D mesh into a 3D one (and reciprocally).
%
%   * gf_mesh_set(mesh M, 'boundary', int rnum, mat CVFIDs)
%   DEPRECATED FUNCTION. Use 'region' instead.
%
%   * gf_mesh_set(mesh M, 'region', int rnum, mat CVFIDs)
%   Assigns the region number `rnum` to the set of convexes or/and convex
%   faces provided in the matrix `CVFIDs`.
%   
%   The first row of `CVFIDs` contains convex #ids, and the second row
%   contains a face number in the convex (or 0
%   for the whole convex (regions are usually used to store a list of
%   convex faces, but you may also use them to store a list of convexes).
%   
%   If a vector is provided (or a one row matrix) the region will represent
%   the corresponding set of convex.
%
%   * gf_mesh_set(mesh M, 'extend region', int rnum, mat CVFIDs)
%   Extends the region identified by the region number `rnum` to include
%   the set of convexes or/and convex faces provided in the matrix
%   `CVFIDs`, see also ``gf_mesh_set(mesh M, 'set region)``.
%
%   * gf_mesh_set(mesh M, 'region intersect', int r1, int r2)
%   Replace the region number `r1` with its intersection with region number `r2`.
%
%   * gf_mesh_set(mesh M, 'region merge', int r1, int r2)
%   Merge region number `r2` into region number `r1`.
%
%   * gf_mesh_set(mesh M, 'region subtract', int r1, int r2)
%   Replace the region number `r1` with its difference with region
%   number `r2`.
%
%   * gf_mesh_set(mesh M, 'delete boundary', int rnum, mat CVFIDs)
%   DEPRECATED FUNCTION. Use 'delete region' instead.
%
%   * gf_mesh_set(mesh M, 'delete region', ivec RIDs)
%   Remove the regions whose #ids are listed in `RIDs`
%
%   * gf_mesh_set(mesh M, 'merge', mesh m2[, scalar  tol])
%   Merge with the mesh `m2`.
%   
%   Overlapping points, within a tolerance radius `tol`, will not be
%   duplicated. If `m2` is a mesh_fem object, its linked mesh will be used.
%
%   * gf_mesh_set(mesh M, 'optimize structure'[, int with_renumbering])
%   Reset point and convex numbering.
%   
%   After optimisation, the points (resp. convexes) will
%   be consecutively numbered from 1 to gf_mesh_get(mesh M, 'max pid')
%   (resp. gf_mesh_get(mesh M, 'max cvid')).
%
%   * gf_mesh_set(mesh M, 'refine'[, ivec CVIDs])
%   Use a Bank strategy for mesh refinement.
%   
%   If `CVIDs` is not given, the whole mesh is refined. Note
%   that the regions, and the finite element methods and
%   integration methods of the mesh_fem and mesh_im objects linked
%   to this mesh will be automagically refined.
%
%
function [varargout]=gf_mesh_set(varargin)
  if (nargout),
    [varargout{1:nargout}]=gf_octave('mesh_set', varargin{:});
  else
    gf_octave('mesh_set', varargin{:});
    if (exist('ans', 'var') == 1), varargout{1}=ans; end;
  end;
% autogenerated mfile;
