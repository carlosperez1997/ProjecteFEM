function Plot(X,T)

Ndim = size(X,2);
Nnodes = size(X,1);


% Conectivities for bars and cables
%Tbar   = T(Tmat==2,:);
%Tcable = T(Tmat==1,:);


for i = 1:Ndim
    Xbar{i} = reshape(X(T,i),size(T))';
    %Xcable{i} = reshape(x(Tcable,i),size(Tcable))';
end

%Xbar = reshape(input.xpoints(input.T,:),size(input.T))';

figure('visible','off','color','w','Name','Parachute');

% Initial plot
p_bar = patch(Xbar{:},zeros(size(Xbar{1})),'edgecolor',[0.5,0.5,0.5],'linewidth',2);
%p_cable = patch(Xcable{:},zeros(size(Xcable{1})),'edgecolor',[0.0,0.0,0.0],'linewidth',0.5);
%p_bar_ref = patch(Xbar{:},zeros(size(Xbar{1})),'edgecolor',[0.7,0.7,0.7],'linewidth',2);
%p_cable_ref = patch(Xcable{:},zeros(size(Xcable{1})),'edgecolor',[0.8,0.8,0.8],'linewidth',0.5);
view(30,25);
axis equal;
colormap jet; 
cbar = colorbar;
set(cbar,'visible','off');
set(gcf,'visible','on');