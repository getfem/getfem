mesh=gf_mesh('load','xfem_stab_unilat_contact_friction.meshfem');    
mf=gf_mesh_fem('load','xfem_stab_unilat_contact_friction.meshfem', mesh);
U=load('xfem_stab_unilat_contact_friction.U');   
gf_plot(mf,U','mesh','off','norm','on','deformed_mesh','on','deformation_scale',1, 'deformation_mf', mf, 'deformation', U')
%caxis([0 0.009])