.. include:: ../replaces.txt

.. highlightlang:: matlab

.. _scilab-plotcmdref:

Draw Command reference
======================


gf_colormap
-------------------------------------------

**Synopsis**

::

  c=gf_colormap(name)


**Description :**

  return a colormap, or change the current colormap.
  name can be: 'tripod', 'chouette', 'froid', 'tank'
  or 'earth'.


gf_plot
-------------------------------------------

**Synopsis**

::

  [hsurf, hcontour, hquiver, hmesh, hdefmesh]=gf_plot(mesh_fem mf, U, ...)

  The options are specified as pairs of "option name"/"option value"

  'zplot',{'off'|'on'}       : values of ``U`` are mapped on the $z$-axis (only possible when qdim=1, mdim=2).
  'norm', {'off'|'on'}       : if qdim >= 2, color-plot the norm of the field
  'dir',[]	              : or the scalar product of the field with 'dir' (can be a vector, or 'x', 'y' etc..)
  'refine',8		      : nb of refinments for curved edges and surface plots
  'interpolated',{'off'|'on'}: if triangular patch are interpolated
  'pcolor',{'on'|'off'}      : if the field is scalar, a color plot of its values is plotted
  'quiver',{'on'|'off'}      : if the field is vector, represent arrows 	
  'quiver_density',50        : density of arrows in quiver plot
  'quiver_scale',1           : scaling of arrows (0=>no scaling)
  'mesh',{'off'|'on'}	      : show the mesh ?
  'meshopts',{cell(0)}	         : cell array of options passed to gf_plot_slice for the mesh
  'deformed_mesh', {'off'|'on'} : shows the deformed mesh (only when qdim == mdim)
  'deformed_meshopts', {cell(0)}: cell array of options passed to gf_plot_slice for the deformed mesh
  'deformation',[]	      : plots on the deformed object (only when qdim == mdim)
  'deformation_mf',[]        : plots on the deformed object (only when qdim == mdim)
  'deformation_scale','10%'  : indicate the amplitude of the deformation. Can be a percentage of the mesh width if given as a string, or an absolute value if given as a number
  'cvlst',[]		      : list of convexes to plot (empty=>all convexes)
  'title',[]                 : set the title
  'contour',[]               : list of contour values
  'disp_options', {'off'|'on'} : shows the option or not.



**Description :**


  The function expects ``U`` to be a row vector. If ``U`` is a scalar
  field, then ``gf_plot(mf,U)`` will fill the mesh with colors
  representing the values of ``U``. If ``U`` is a vector field, then
  the default behavior of ``gf_plot`` is to draw vectors representing
  the values of ``U``.

  On output, this function returns the handles to the various
  graphical objects created: ``hmesh`` is the handles to the mesh
  lines, ``hbound`` is the handles to the edges of the boundaries, ``hfill``
  is the handle of the patch objects of faces, ``hvert`` (resp
  ``hconv``, ``hdof``) is the handles of the vertices (resp. convexes,
  dof) labels.

  For example, plotting a scalar field on the border of a 3D mesh can be done with ::

    % load the 'strange.mesh_fem' (found in the getfem_matlab/tests directory)
    mf=gf_mesh_fem('load', 'strange.mesh_fem')
    U=rand(1, gf_mesh_fem_get(mf, 'nbdof')); # random field that will be drawn
    gf_plot(mf, U, 'refine', 25, 'cvlst', gf_mesh_get(mf,'outer faces'), 'mesh','on');




gf_plot_1D
-------------------------------------------

**Synopsis**

::

  gf_plot_1D(mesh_fem mf, U, ...)

  Available options are specified as pairs of "option name"/"option value"

  'style', 'bo-'       : line style and dof marker style (same syntax as in the Scilab command 'plot');
  'color', ''          : override line color (by a given color name);
  'dof_color', ''      : override color of dof markers;
  'width', 2           : line width.


**Description :**


  This function plots a 1D finite element field.


gf_plot_mesh
-------------------------------------------

**Synopsis**

::

  gf_plot_mesh(m, ...)

  'vertices', {'off'|'on'}    : displays also vertices numbers.
  'convexes', {'off'|'on'}    : displays also convexes numbers.
  'dof',{'off'|'on'}          : displays also finite element nodes. In that case, ``m`` should be a ``mesh_fem`` identifier.
  'regions',BLST              : displays the boundaries listed in BLST.
  'cvlst',CVLST               : display only the listed convexes. If CVLST has two rows, display only the faces listed in the second row.
  'edges', {'on' | 'off'}     : display edges ?
  'faces', {'off'|'on'}       : fills each 2D-face of the mesh
  'curved', {'off'|'on'}      : displays curved edges
  'refine',N                  : refine curved edges and filled faces N times
  'deformation', Udef         : optionnal deformation applied to the mesh (M must be a mesh_fem object)
  'edges_color',[.6 .6 1]     : RGB values for the color of edges
  'edges_width',1             : width of edges
  'faces_color',[.75 .75 .75]): RGB values for the color of faces
  'quality',{ 'off' | 'on' }  : Display the quality of the mesh.


**Description :**

  This function is used to display a mesh.

  Example ::

    % the mesh is in the tests directory of the distribution
    m=gf_mesh('import','gid','donut_with_quadratic_tetra_314_elements.msh');
    gf_plot_mesh(m,'refine',15,'cvlst',gf_mesh_get(m,'outer faces'),'faces','on',\ldots, 'faces_color',[1. .9 .2],'curved','on','edges_width',2);
    camlight % turn on the light!



gf_plot_slice
-------------------------------------------

**Synopsis**

::

  gf_plot_slice(sl, ...)

  The options are specified as pairs of "option name"/"option value"


  data    []          : data to be plotted (one value per slice node)
  convex_data    []   : data to be plotted (one value per mesh convex)
  mesh, ['auto']      : 'on' -> show the mesh (faces of edges), 'off' -> ignore mesh
  mesh_edges, ['on']  : show mesh edges ?
  mesh_edges_color, [0.60 0.60 1] : color of mesh edges
  mesh_edges_width, [0.70] : width of mesh edges
  mesh_slice_edges, ['on'] : show edges of the slice ?
  mesh_slice_edges_color, [0.70 0 0] : color of slice edges
  mesh_slice_edges_width, [0.50] : width of slice edges
  mesh_faces, ['off'] : 'on' -> fill mesh faces (otherwise they are transparent)
  mesh_faces_color, [0.75 0.75 0.75]
  pcolor, ['on']      : if the field is scalar, a color plot of its values is plotted
  quiver, ['on']      : if the field is vector, represent arrows
  quiver_density, 50  : density of arrows in quiver plot
  quiver_scale, 1     : density of arrows in quiver plot
  tube, ['on']        : use tube plot for 'filar' (1D) parts of the slice
  tube_color, ['red'] : color of tubes (ignored if 'data' is not empty and 'pcolor' is on)
  tube_radius, ['0.5%'] : tube radius; you can use a constant, or a percentage (of the mesh size) or a vector of nodal values
  showoptions, ['on'] : display the list of options

  the 'data' and 'convex_data' are mutually exclusive.


**Description :**

  This function can be used to plot mesh slices. It is also used by the ``gf_plot_mesh`` and ``gf_plot`` functions.


  Example : consider that you have a 3D mesh_fem ``mf`` and a vector field ``U`` defined on this mesh_fem, solution of the Stokes problem in a tank (see the demo ``demo_stokes_3D_tank_draw.m`` in the tests directory). ::

    figure;
    % slice the mesh with two half spaces, and take the boundary of the resulting quarter-cylinder
    sl=gf_slice(\{'boundary',\{'intersection',\{'planar',+1,[0;0;0],[0;1;0]\},\ldots
                                              \{'planar',+1,[0;0;0],[1;0;0]\}\}\},m,6);
    Usl=gf_compute(pde.mf_u,U,'interpolate on', sl);  % interpolate the solution on the slice
    % show the norm of the displacement on this slice
    gf_plot_slice(sl,'mesh','on','data',sqrt(sum(Usl.^2,1)),'mesh_slice_edges','off');

    % another slice: now we take the lower part of the mesh
    sl=gf_slice(\{'boundary',\{'intersection',\{'planar',+1,[0;0;6],[0;0;-1]\},\ldots
                                            \{'planar',+1,[0;0;0],[0;1;0]\}\}\},m,6);
    Usl=gf_compute(pde.mf_u,U,'interpolate on', sl);
    hold on;
    gf_plot_slice(sl,'mesh','on','data',sqrt(sum(Usl.^2,1)),'mesh_slice_edges','off');

    % this slice contains the transparent mesh faces displayed on the picture
    sl2=gf_slice(\{'boundary',\{'planar',+1,[0;0;0],[0;1;0]\}\},\ldots
                m,6,setdiff(all_faces',TOPfaces','rows')');
    gf_plot_slice(sl2,'mesh_faces','off','mesh','on','pcolor','off');

    % last step is to plot the streamlines
    hh=[1 5 9 12.5 16 19.5]; % vertical position of the different starting points of the streamlines
    H=[zeros(2,numel(hh));hh];

    % compute the streamlines
    tsl=gf_slice('streamlines',pde.mf_u,U,H);
    Utsl=gf_compute(pde.mf_u,U,'interpolate on', tsl);

    % render them with "tube plot"
    [a,h]=gf_plot_slice(tsl,'mesh','off','tube_radius',.2,'tube_color','white');
    hold off;
    % use a nice colormap
    caxis([0 .7]);
    c=[0 0 1; 0 .5 1; 0 1 .5; 0 1 0; .5 1 0; 1 .5 0; 1 .4 0; 1 0 0; 1 .2 0; 1 .4 0; 1 .6 0; 1 .8 0];
    colormap(c);
