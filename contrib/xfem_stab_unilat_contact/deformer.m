mesh=gf_mesh('load','xfem_stab_unilat_contact.meshfem');    
mf=gf_mesh_fem('load','xfem_stab_unilat_contact.meshfem', mesh);
U=load('xfem_stab_unilat_contact.U');   
gf_plot(mf,U','mesh','off','norm','on','deformed_mesh','on','deformation_scale',1, 'deformation_mf', mf, 'deformation', U')
