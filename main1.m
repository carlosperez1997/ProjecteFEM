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

CN = CN_global(input.xpoints, input.T);

Plot(input.xpoints, input.T)

coord_motor = input.xpoints(33,:);

% Centro de masas
[CG, M] = CentroMasas(input.xpoints, input.rho, input.A, input.T, coord_motor, Masa_motor);

I_CM = Inercia_CM (input.xpoints, input.rho, input.A, input.T, coord_motor, Masa_motor, CG);

[Kg, Mg] = assembly_KG_MG(input.xpoints, input.T, input.E, input.A, CN, input.rho);

F = assembly_F(input.xpoints, input.T, CN, input.Fext);
%[U, R, vl, vr] = solveSys(input.xpoints, input.T, Kg, F, input.fixnodes);

%elementcount(input.T, 155)

NdofsxNode = 3;
Ndofs = 3*size(input.xpoints,1);
[U, R, vl, vr] = solveSys(NdofsxNode, Ndofs,input.fixnodes,Kg,F);
%[u ,R, vl, vr] = solveSys(NdofsXnode,Ndofs,fixNod,KG,f)
output.sx = compute_sigma(input.xpoints, input.T, CN, input.E, U);



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



