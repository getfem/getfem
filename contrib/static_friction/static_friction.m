% -*- matlab -*- (enables emacs matlab mode)
% addpath ~/source++/getfem_toolbox
% addpath ~/source++/getfem++/contrib/static_friction/

gf_workspace('clear all');
mesh = gf_mesh('load', 'static_friction.meshfem');
mf = gf_mesh_fem('load', 'static_friction.meshfem', mesh);
mf_vm = gf_mesh_fem('load', 'static_friction.meshfem_vm', mesh);
U = load('static_friction.U')';
VM = load('static_friction.VM')';
% VM = max(1, VM);
% VM = log(VM);
% VM = max(1E6, VM);
% VM = min(1E8, VM);
% for i = 1:size(VM, 2),
% if (VM(i) > 1E5)
%   VM(i) = 1E5;
% end;
% end;
% clear VM2; VM2(1:2:2*size(VM, 2)) = VM; VM2(2*size(VM, 2)) = 0;
gf_plot(mf_vm, VM, 'norm', 'on', 'refine', 1, 'deformation', U, 'deformation_mf', mf, 'deformed_mesh','on', 'deformation_scale', 1.0);
%caxis([1E3 2e8]);
% axis([-0.2 1.2 -0.2 1.2]);
colorbar;
