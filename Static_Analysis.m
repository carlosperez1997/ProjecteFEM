
run('input_data');

input.CN = CN_global(input.xpoints, input.T);

%Plot(input.xpoints, input.T)

% FORCES
F_gravity = Gravity_force(input.xpoints, input.T, input.A, input.rho, input.g);
F = assembly_F(input.xpoints, input.T, input.CN, input.Fext);
F_grav = assembly_F(input.xpoints, input.T, input.CN, F_gravity);

% Força termica ja expressada en coord. globals
F_termic = Termical_force(input.xpoints, input.T, input.CN, input.E, input.A, input.Delta_T, input.alpha);

% Centro de masas 
[CG, M] = CentroMasas(input.xpoints, input.rho, input.A, input.T, input.mass)
I_CM = Inercia_CM (input.xpoints, input.rho, input.A, input.T, input.mass, CG);

% MATRIZ KG y MG
[Kg, Mg] = assembly_KG_MG(input.xpoints, input.T, input.E, input.A, input.CN, input.rho);

F = F + F_grav; 

% Solver
[U, R, vl, vr] = solver (input.xpoints, input.T, Kg, F, input.fixnodes);

% Sigma
output.sx = compute_sigma(input.xpoints, input.T, input.CN, input.E, U);

U_vector = Convert2vector(input.xpoints, input.T, U);

Final_structure = input.xpoints+U_vector;

X = input.xpoints;
Y = Final_structure;
T = input.T;

Ndim = size(X,2);
Nnodes = size(X,1);

% Conectivities for bars and cables
%Tbar   = T(Tmat==2,:); %Tcable = T(Tmat==1,:);


for i = 1:Ndim
    Xbar{i} = reshape(X(T,i),size(T))';
    Xbar2{i} = reshape(Y(T,i),size(T))';
    %Xcable{i} = reshape(x(Tcable,i),size(Tcable))';
end

%Xbar = reshape(input.xpoints(input.T,:),size(input.T))';

figure('visible','off','color','w','Name','Parachute');

% Initial plot
p_bar = patch(Xbar{:},zeros(size(Xbar{1})),'edgecolor',[0.5,0.5,0.5],'linewidth',2);
hold on;
p_cable = patch(Xbar2{:},zeros(size(Xbar2{1})),'edgecolor','red','linewidth',0.5);
%p_bar_ref = patch(Xbar{:},zeros(size(Xbar{1})),'edgecolor',[0.7,0.7,0.7],'linewidth',2);
%p_cable_ref = patch(Xcable{:},zeros(size(Xcable{1})),'edgecolor',[0.8,0.8,0.8],'linewidth',0.5);
view(30,25);
axis equal;
colormap jet; 
cbar = colorbar;
set(cbar,'visible','off');
set(gcf,'visible','on');

sum(R)