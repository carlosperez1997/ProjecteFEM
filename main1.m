%--------------------------------------------------------------------------
% MAIN CODE
%--------------------------------------------------------------------------
% Authors: Paula Sorolla and Carlos Pérez
% Date: 2/4/2020
%--------------------------------------------------------------------------

clear all;
close all;
clc

run('input_data');

input.CN = CN_global(input.xpoints, input.T);

Plot(input.xpoints, input.T)

F_gravity = Gravity_force(input.xpoints, input.T, input.A, input.rho, input.g);

% Centro de masas 
[CG, M] = CentroMasas(input.xpoints, input.rho, input.A, input.T, input.mass);

I_CM = Inercia_CM (input.xpoints, input.rho, input.A, input.T, input.mass, CG);

[Kg, Mg] = assembly_KG_MG(input.xpoints, input.T, input.E, input.A, input.CN, input.rho);

F = assembly_F(input.xpoints, input.T, input.CN, input.Fext);
F_grav = assembly_F(input.xpoints, input.T, input.CN, F_gravity);
[U, R, vl, vr] = solver (input.xpoints, input.T, Kg, F, input.fixnodes);
%[U, R, vl, vr] = solveSys(input.xpoints, input.T, Kg, F, input.fixnodes);

%elementcount(input.T, 155)

output.sx = compute_sigma(input.xpoints, input.T, input.CN, input.E, U);



    Kll = Kg(vl,vl);
    Mll = Mg(vl,vl);
    [phi,lambda] = eigs(Kll,Mll);
    lambda = sqrt(diag(lambda));
    freq = sqrt(lambda)/2/pi;
    
    % Direct method
    X_dir = direct(Kg, Mg, F, vl, lambda);
    % Modal method (2 and 1 eigenvalues)
    X_mod = modal(Kg, Mg, F, vl, phi, lambda);


disp('Finish');



