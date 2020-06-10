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

%Static_Analysis


%% DYNAMIC
    
% Apartado 6 

%     Kll = Kg(vl,vl);
%     Mll = Mg(vl,vl);
%     [phi1,lambda1] = eigs(Kll,Mll);
%     lambda1 = sqrt(diag(lambda1));
%     freq1 = sqrt(lambda1)/2/pi;
% 
    [phi, lambda] = eig((Mg)^(-1)*Kg);
    %[phi,lambda] = eigs(Kg,Mg);
    
    lambda = sqrt(diag(lambda));
    freq = sqrt(lambda)/2/pi;
    
    % 150; 151; 152; 153; 154; 155; 158; 159
    sol.a1 = phi(:,152) * (1/max(phi(:,152))); 
    %Aircraft_plot(input.xpoints,input.T,input.le,sol.a1);
    %freq(152,1)

    % Modal 
    
    M_modal = phi'*Mg*phi;
    K_modal = phi'*Kg*phi;
    
    [phi, lambda] = eig((M_modal)^(-1)*K_modal);
    
    H_modal = - (lambda').^2 * M_modal + K_modal;
    
    %U_vector = Static_plot(input.xpoints, input.T, input.le, phi(:,159), input.Fext, type);
    
%% Direct method

    Mnn = Mg(vl, vl); Knn = Kg(vl, vl);
    
    [phi, w_n] = eig((Mnn)^(-1)*Knn); % Structure's eigenmodes

    % Frequency response
    
    masa_motor = 65; disp = 1*10^-3;
    
    F_dyn = [ 
        33, 1, 1; %node motor: 33
        33, 2, 1*i;
        33, 3, -1*i;
    ];

    F_dynamic = assembly_F(input.xpoints, input.T, input.CN, F_dyn);
        
    F = F_dynamic * disp;
    
    %% DIRECT FREQUENCY RESPONSE
    tic
    w_excs = 0.01:0.05:30; % Excitation frequencies range
    a_0 = 10^-6; 
    ur = zeros(size(input.fixnodes,1),1);
    
    for i = 1:length(w_excs)
        
        w_exc = w_excs(i);
        F_dir = - w_exc^2 * masa_motor * F;

        Knn1 = Knn * (1+0.02*i);
        X_DF = dirfr(Knn1, Mnn, vl, w_exc, F_dir);
        
        % Displacements : Node Pilot 36
        
        U(vl) = X_DF;
        U(vr) = ur;
        
        U = Convert2vector(input.xpoints, input.T, U);
        X = U(36,:);
        
        acc.x (i) = - w_exc^2 * X(1);
        acc.y (i) = - w_exc^2 * X(2);
        acc.z (i) = - w_exc^2 * X(3);
        acc_log.x(i) = 20*log(abs(real(acc.x(i)))/a_0);
        acc_log.y(i) = 20*log(abs(real(acc.y(i)))/a_0);
        acc_log.z(i) = 20*log(abs(real(acc.z(i)))/a_0);
        
    end
    
    %figure;
    %plot(w_excs,acc_log.x); hold on;
    %plot(w_excs,acc_log.y); 
    %plot(w_excs,acc_log.z);
    
    
    
%% GOOD DIRECT FREQUENCY ANALYSIS:

%Mass Matrix Calculation: Model: Consistent
x = input.xpoints;
T = input.T;

Ndim = size(x,1);           % Number of dimensions (DOFs for each node)
Nnodes = size(x,2);         % Number of nodes
NnodesXelement = size(T,1); % Number of nodes for each element
Nelements = size(T,2);      % Number of elements
Ndofs = Nnodes*Ndim ;       % Total number of degrees of freedom

Count=1;
eta=0.02; %Percentage of Structural Damping

% GLOBAL STIFFNESS MATRIX PLUS STRUCTURAL DAMPING:

% Aquí has de cridar la definida a Static, però de moment, anem programant
% amb això

KG = Knn;
%KG=rand(Ndofs, Ndofs);
KG=KG*(1+1i*eta);

MTOT = Mnn;

% EXTERNAL FORCE MATRIX:

Off=1e-3;
Meng=65;

%Fext is Centrifugue Force: F=m*r*(w^2)

        %Node   %DOF        %Magnitude      
Fext=[  33,      1,         Meng*Off;
        33,      2,         1i*Meng*Off;
        33,      3,         Meng*Off;
    ]';

Forces=zeros(Ndofs, 1);
Nforces=size(Fext, 2);
for e=1:Nforces
    Forces(Fext(1, e)*Fext(2, e), 1)=Fext(3, e);
end

Forces = Forces(vl);

[Psi, Lambda]=eig(KG, MTOT);
%wn=sqrt(diag(Lambda));
%fn=wn/(2*pi);

tic
w_sweeps=0.01:0.01:30;
a_0=10e-6;

for i=1:length(w_sweeps)
    
    w_sweep=w_sweeps(i);
    F_syst=-Forces*w_sweep^2;
    
    H=KG-MTOT*w_sweep^2;
    X_DF=H\F_syst;
    %X=X_DF(45);
    
    U(vl) = X_DF;
    U(vr) = ur;
    
    U = Convert2vector(input.xpoints, input.T, U);
    X = U(36,:);

    Disp(i,:) = real(X);
    ACC(i,:)=w_sweep^2*X;
    ACC_LOG(i,:)=20*log10(abs(real(ACC(i,:)))/a_0);
    
end

figure; hold on;
plot(w_sweeps, ACC_LOG(:,1));
plot(w_sweeps, ACC_LOG(:,2));
plot(w_sweeps, ACC_LOG(:,3));


figure; hold on;
plot(w_sweeps, Disp(:,1));
plot(w_sweeps, Disp(:,2));
plot(w_sweeps, Disp(:,3));


% MODAL FREQUENCY RESPONSE
% lambda = sqrt(diag(lambda));
%  
% for i=1:length(w_sweeps)
%         
%         freq = w_sweeps(i)/(2*pi);
%         w_exc = w_sweeps(i);
%         
%         F_dir = - w_exc^2 * masa_motor * F;
%         
%         if freq < 10
%             chi = 0.01;
%         else
%             chi = 0.05;
%         end
%         
%         for j = 1:size(phi,1)
%         
%             for k=1:size(phi,2)
%                 Cnn(j,k) = 2 * Mnn(j,k) * wn(j) *chi;
%             end
%         
%         end
% 
%         X_MF = modfr2(phi, Mnn, Cnn, Knn, chi, lambda, w_exc, vl, F_dir);
%         
%         U(vl) = X_MF;
%         U(vr) = ur;
%         
%         U = Convert2vector(input.xpoints, input.T, U);
%         X = U(36,:);
%         
%         Disp2(i,:) = abs(X);
%         ACC2(i,:)=w_sweep^2*X;
%         ACC_LOG2(i,:)=20*log10(abs(real(ACC2(i,:)))/a_0);
%         
% end
% 
% figure; hold on;
% plot(w_sweeps, ACC_LOG2(:,1));
% plot(w_sweeps, ACC_LOG2(:,2));
% plot(w_sweeps, ACC_LOG2(:,3));
% 
% 
% figure; hold on;
% plot(w_sweeps, Disp2(:,1));
% plot(w_sweeps, Disp2(:,2));
% plot(w_sweeps, Disp2(:,3));
%     A = 1;
    
%% MARGA

wn = sqrt(Lambda);

tic
w_sweeps=0.01:0.01:30;
a_0=1e-6;

for i=1:length(w_sweeps)
    
    w_sweep=w_sweeps(i);
    F_syst=-Forces*w_sweep^2;
    Fmodal=phi'*F_syst;
    
    if (w_sweep/(2*pi))<10
        modal_damping(i)=0.01;
    else
        modal_damping(i)=0.05;
    end 
    
    Hmod=inv(wn^2 - w_sweep^2 + 2*1i*w_sweep*wn*modal_damping(i) );
    Q=Hmod*Fmodal;
    X_M=phi*Q;
    
    U(vl) = X_M;
    U(vr) = ur;
    
    U = Convert2vector(input.xpoints, input.T, U);
    X = U(36,:);

    Disp2(i,:) = real(X);
    ACC2(i,:)=w_sweep^2*X;
    ACC_LOG2(i,:)=20*log10(abs(real(ACC2(i,:)))/a_0);
    
    %X=X_M(54); %Node under pilot
    %XdisM(i)=abs(X);
    %YdisM(i)=angle(X);
    
    %accM(i)=w_sweep^2*X;
    %acc_logM(i)=20*log10(abs(real(accM(i)/a_0)));
    
end

figure (5)
plot(w_sweeps, ACC_LOG2(:,1)); hold on;
plot(w_sweeps, ACC_LOG(:,1));

figure(6)
plot(w_sweeps, XdisM);

figure(7)
plot(w_sweeps, YdisM);

disp('Finish');


%% NO TOCAR QUE VA!

% wn = sqrt(Lambda);
% 
% tic
% w_sweeps=0.01:0.05:30;
% a_0=1e-6;
% 
% 
% for i=1:length(w_sweeps)
%     
%     w_sweep=w_sweeps(i);
%     F_syst=-Forces*w_sweep^2;
%     Fmodal=phi'*F_syst;
%     
%     if (w_sweep/(2*pi))<10
%         modal_damping(i)=0.01;
%     else
%         modal_damping(i)=0.05;
%     end 
%     
%     Hmod=inv(wn^2 - w_sweep^2 + 2*1i*w_sweep*wn*modal_damping(i) );
%     Q=Hmod*Fmodal;
%     X_M=phi*Q;
%     
%     X=X_M(54); %Node under pilot
%     XdisM(i)=abs(X);
%     YdisM(i)=angle(X);
%     
%     accM(i)=w_sweep^2*X;
%     acc_logM(i)=20*log10(abs(real(accM(i)/a_0)));
%     
% end

