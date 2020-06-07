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

%% DYNAMIC

    Kll = Kg(vl,vl);
    Mll = Mg(vl,vl);
    [phi,lambda] = eigs(Kll,Mll);
    lambda = sqrt(diag(lambda));
    freq = sqrt(lambda)/2/pi;
    
    % Direct method
    X_dir = direct(Kg, Mg, F, vl, lambda);
    % Modal method (2 and 1 eigenvalues)
    X_mod = modal(Kg, Mg, F, vl, phi, lambda);


% Eigenmodes computation

Mnn = Mg(vl, vl);
Knn = Kg(vl, vl);

[phi, w_n] = eig((Mnn)^(-1)*Knn); % Structure's eigenmodes

% Frequency response

% DIRECT FREQUENCY RESPONSE
tic
w_exc = 0:1:1500; % Excitation frequencies range
f = zeros(168,1);
X_DF = dirfr(Knn, Mnn, vl, w_exc, f);
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

for n=1:floor(length(phi)) % Only the last one considered
   
X_MF = modfr(phi, Mnn, Knn, w_exc, n, vl,f);

for i=1:size(X_MF,1)
    
vecnorm2_squared(n,:) = vecnorm2_squared(n,:) + X_MF(i,:).^2;

end

end

vecnorm2 = sqrt(vecnorm2_squared);

toc


freq_plot(w_exc,vecnorm1,vecnorm2,phi)



disp('Finish');



