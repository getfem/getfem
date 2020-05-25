// Scilab GetFEM interface
//
// Copyright (C) 2011 Mariama Ndiaye, Yves Renard, Yann Collette.
//
// This file is a part of GetFEM
//
// GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//  Transient Navier-Stokes equation on a driven cavity with two numerical
//  scheme : a projection and a semi-implicit scheme. Without using the
//  bricks.
//
//  This program is used to check that matlab-getfem is working. This is
//  also a good example of use of GetFEM.
//

lines(0);
stacksize('max');

if getos()=='Windows' then
  // Under Windows, all the trace messages are available in the dos console
  // Under Linuxs, all the trace messages are redirected to the Scilab console
  consolebox('on');
end
gf_util('trace level',3);
gf_util('warning level',3);

gf_workspace('clear all');

NX = 20;      // Space resolution
Dt = 0.02;    // Time step
T  = 10;      // Time interval
nu = 0.005;   // Viscosity
v  = 30;      // Driven velocity
scheme = 1;   // 1 : Projection scheme.
              // 2 : Semi-implicit scheme
rho = 2;      // density
g = 9.81;     // gravity constant

//m = gf_mesh('cartesian', [0:1/NX:1],[0:1/NX:1]); 
m = gf_mesh('triangles grid',[0:1/NX:1],[0:1/NX:1]);
border = gf_mesh_get(m,'outer faces');
// mark it as boundary #1
gf_mesh_set(m, 'boundary', 1, border);
drawlater; clf();
gf_plot_mesh(m, 'regions', [1]); // the boundary edges appears in red
drawnow;

// create mesh_fem objects
mf_u = gfMeshFem(m,2);  // For the velocity.
mf_p = gfMeshFem(m,1);  // For the pressure.
mf_f = gfMeshFem(m,1);  // For the external forces.

//mim = gf_mesh_im(m,  gf_integ('IM_GAUSS_PARALLELEPIPED(2, 2)'));
mim = gf_mesh_im(m, gf_integ('IM_HCT_COMPOSITE(IM_TRIANGLE(13))'));

// assign the fems to all convexes of the mesh_fems
if (scheme == 1) then
  gf_mesh_fem_set(mf_u,'fem',gf_fem('FEM_PK(2,1)'));
  gf_mesh_fem_set(mf_p,'fem',gf_fem('FEM_PK(2,1)'));
  gf_mesh_fem_set(mf_f,'fem',gf_fem('FEM_PK(2,1)'));
elseif (scheme == 2) then
  gf_mesh_fem_set(mf_u,'fem',gf_fem('FEM_PK(2,2)'));
  gf_mesh_fem_set(mf_p,'fem',gf_fem('FEM_PK(2,1)'));
  gf_mesh_fem_set(mf_f,'fem',gf_fem('FEM_PK(2,1)'));
end

// Assembly

Fd = [gf_mesh_fem_get_eval(mf_f, list(0)); gf_mesh_fem_get_eval(mf_f, list(rho*g))];
FD = Fd;
U  = zeros(gf_mesh_fem_get(mf_u, 'nbdof'), 1); // initial condition

M  = gf_asm('mass matrix',mim,mf_u) / Dt;
K  = nu*gf_asm('volumic','M(#1,#1)+=comp(vGrad(#1).vGrad(#1))(:,i,j,:,i,j)', mim, mf_u);
F  = gf_asm('volumic source',mim,mf_u,mf_f,Fd); //#1 methode d'elmt fini 1, vBase vecteur de base de methode d'EF 1, vGrad grad vect
Kp = gf_asm('volumic','M(#1,#1)+=comp(Grad(#1).Grad(#1))(:,i,:,i)', mim, mf_p);
D  = gf_asm('volumic','M(#1,#2)+=comp(Base(#1).vGrad(#2))(:,:,i,i)', mim, mf_p, mf_u) / Dt;
B  = gf_asm('volumic','M(#1,#2)+=comp(vBase(#1).Grad(#2))(:,i,:,i)', mim, mf_u, mf_p);

// for the vorticity computation
Mo  = gf_asm('mass matrix', mim, mf_f);
MVo = gf_asm('volumic','t=comp(Base(#1).vGrad(#2));M(#1,#2)+=t(:,:,1,2)-t(:,:,2,1)', mim, mf_f, mf_u);


UBOUND = gf_mesh_fem_get(mf_u, 'dof on region', 1);// supplies the list of dofs on the boundary
UNODES = gf_mesh_fem_get(mf_u, 'basic dof nodes');
Kp(1, :) = 0; // In order to fix the pressure on a node for scheme 1.
Kp(1, 1) = 1;
Ndofu = size(D,2);          // Dof number for the velocity
Ndofp = size(D,1);          // Dof number for the pressure

h = get("current_figure");
h.color_map = jetcolormap(255);

for t=0:Dt:T
  printf('Time step = %f / %f\n', t, T);
  
  if (scheme == 1) then
    
    C = gf_asm('volumic','a=data(#1);M(#1,#1)+=comp(vBase(#1).vGrad(#1).vBase(#1))(i,j,:,k,j,:,k).a(i)', mim,mf_u, U);
    A = M + K + C;
    L = F + M * U;
  
    for i=UBOUND     // Boundary conditions
        A(i, :) = 0; A(i,i) = 1; L(i) = 0;
        if (modulo(double(i), 2) == 1) then
            node = UNODES(:, i);
            if (abs(node(2)-1) < 1e-10 & abs(node(1)-0.5) < 0.499) then
               L(i) = v * node(1) * (1-node(1));
            end
        end
    end

    U1_2 = A \ L;
 
    L2 = -D * U1_2;
    L2(1) = 0;
    P =  Kp \ L2 ;
    U = M \ (M * U1_2 - B * P);
  
  elseif (scheme == 2) then
      
      C = gf_asm('volumic','a=data(#1);M(#1,#1)+=comp(vBase(#1).vGrad(#1).vBase(#1))(i,j,:,k,j,:,k).a(i)', mim,mf_u, U);
      C = C+gf_asm('volumic','a=data(#1);M(#1,#1)+=comp(vBase(#1).vGrad(#1).vBase(#1))(:,i,j,k,k,:,i).a(j)', mim,mf_u, U)/2;
      
      A = [M+K+C, (-Dt*D)'; -Dt*D, zeros(Ndofp)];
      L = F + M * U;
      
      for i=UBOUND     // Boundary conditions
        A(i, :) = 0; A(i,i) = 1; L(i) = 0;
        if (modulo(double(i), 2) == 1) then
            node = UNODES(:, i);
            if (abs(node(2)-1) < 1e-10 & abs(node(1)-0.5) < 0.499) then
               L(i) = v * node(1) * (1-node(1));
            end
        end
      end
      
      A(Ndofu+1, :) = 0;   // In order to fix the pressure on a node.
      A(Ndofu+1, Ndofu+1) = 1;
      
      UP = A \ [L; zeros(Ndofp,1)];
      U  = UP(1:Ndofu);
      P  = UP(Ndofu+1:Ndofu+Ndofp);
      
  end

  Vo = Mo \ (MVo * U); // Vorticity projected on mf_f. 
  
  drawlater; clf();
  gf_plot(mf_p,P','refine',1);
  //a = get("current_axes"); a.data_bounds = [0 1 0 1];
  //hold on;
  gf_plot(mf_f,Vo','refine',1,'contour',[-40,-20,-10,10,20,40,80], 'pcolor', 'off');
  gf_plot(mf_u, U','mesh','off', 'quiver_density', 15, 'quiver_scale',0.3);
  //hold off;
  colorbar(min(P'),max(P'));
  title(sprintf('Quiver plot of U, with color plot of the pressure and vorticity contour lines, t=%g', t));
  drawnow;
  sleep(1000);
   
end
