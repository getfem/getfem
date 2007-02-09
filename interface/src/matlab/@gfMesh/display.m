function display(m)
% gfMesh/display displays a short summary for a gfMesh object
% see gfMesh/char for an exhaustive description of the object.

  for i=1:numel(m)
    disp(sprintf('gfMesh object ID=%u [%d bytes], dim=%d, nbpts=%d, nbcvs=%d',double(m(i).id),...
                 gf_mesh_get(m(i),'memsize'), gf_mesh_get(m(i),'dim'), ...
                 gf_mesh_get(m(i),'nbpts'), gf_mesh_get(m(i),'nbcvs')));
  end;
