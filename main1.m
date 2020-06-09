%--------------------------------------------------------------------------
% MAIN CODE
%--------------------------------------------------------------------------
% Authors: Paula Sorolla and Carlos Pérez
% Date: 2/4/2020
%--------------------------------------------------------------------------

clear all;
close all;
clc

%% STATIC

Static_Analysis


% INERTIA RELIEF

F_grav2 = F_grav;
F_ext2 = assembly_F(input.xpoints, input.T, input.CN, input.Fext);
F = F_ext2 + F_termic + F_grav2;

acceleration = F / M;

F_inercial = acceleration * M;




%% DYNAMIC
    
    Kll = Kg(vl,vl);
    Mll = Mg(vl,vl);
    [phi,lambda] = eigs(Kll,Mll);
    lambda = sqrt(diag(lambda));
    freq = sqrt(lambda)/2/pi;
    
    
    
    %% Direct method
    %X_dir = direct(Kg, Mg, F, vl, lambda);
    %% Modal method (2 and 1 eigenvalues)
    %X_mod = modal(Kg, Mg, F, vl, phi, lambda);

%%DAVID
% Eigenmodes computation

    Mnn = Mg(vl, vl);
    Knn = Kg(vl, vl);

    [phi, w_n] = eig((Mnn)^(-1)*Knn); % Structure's eigenmodes

    % Frequency response
    
    masa_motor = 65;    
    disp = 1*10^-3;
    
    F_dyn = [
        33, 1, 1;
        33, 2, 1;
        33, 3, 1;
    ];

    F_dynamic = assembly_F(input.xpoints, input.T, input.CN, F_dyn);
        
    F = F_dynamic * disp;
    
    %% DIRECT FREQUENCY RESPONSE
    tic
    w_excs = 0.01:0.05:30; % Excitation frequencies range
    %f = zeros(168,1);
    a_0 = 10^-6;
    
    for i = 1:length(w_excs)
        
        w_exc = w_excs(i);
        F_dir = - w_exc^2 * masa_motor * F;

            Knn1 = Knn * (1+0.02*i);
        
        X_DF = dirfr(Knn1, Mnn, vl, w_exc, F_dir);
        
        % Node Piloto
        X = X_DF(36);
        
        acc (i) = - w_exc^2 * X;
        acc_log(i) = 20*log(abs(acc(i))/a_0);
        
    end
    
    figure;
    plot(w_excs,acc_log);
    
    % for i=1:length(X_DF)
%     
    % vecnorm1_squared = sqrt(X_DF(1,:).^2 + X_DF(2,:).^2);
% 
% end



    vecnorm1_squared = zeros(1,length(X_DF));
    vecnorm2_squared = zeros(length(phi),length(X_DF));

    for i=1:size(X_DF,1)
    
        vecnorm1_squared = vecnorm1_squared + X_DF(i,:).^2;

    end

vecnorm1 = sqrt(vecnorm1_squared);

% MODAL FREQUENCY RESPONSE

for i = 1:length(w_excs)
    for n=1:floor(length(phi)) % Only the last one considered

        w_exc = w_excs(i);
        F_dir = - w_exc^2 * masa_motor * F;
        
        freq = w_exc / (2*pi);
        
        if freq < 10
            Bnn = Knn*(0.01*i);
            Knn1 = Knn*(1+0.01*i);
        else
            Bnn = Knn*(0.05*i);
            Knn1 = Knn*(1+0.05*i);
        end

        X_MF = modfr2(phi, Mnn, Bnn, Knn1, w_exc, n, vl,F_dir);

        X =  X_MF(36);
        acc2 (i) = - w_exc^2 * X;
        acc_log2(i) = 20*log(abs(acc2(i))/a_0);
        %for i=1:size(X_MF,1)

        %    vecnorm2_squared(n,:) = vecnorm2_squared(n,:) + X_MF(i,:).^2;

        %end

    end
end

figure;
    plot(w_excs,acc_log); hold on;
    plot(w_excs,acc_log2);

vecnorm2 = sqrt(vecnorm2_squared);

toc


freq_plot(w_exc,vecnorm1,vecnorm2,phi)



disp('Finish');



