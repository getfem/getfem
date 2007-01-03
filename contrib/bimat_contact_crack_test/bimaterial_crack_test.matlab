
gf_workspace('clear all');
mf = gf_mesh_fem('load', 'bimaterial_crack.meshfem'); 
U = load('bimaterial_crack.U')';
gf_plot(mf, U, 'refine', 2, 'deformation', U, 'deformed_mesh','on', 'deformation_scale', '100%');