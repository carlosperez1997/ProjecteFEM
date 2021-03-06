function F_termic = Termical_force(x, T, CN, Es, A, Delta_T, alpha)
    % Dimensions
    Ndim = size(x,2); % Number of dimensions (DOFs for each node)
    Nnodes = size(x,1); % Number of nodes
    NnodesXelement = size(T,2); % Number of nodes for each element
    Nelements = size(T,1); % Number of elements
    Ndofs = Nnodes*Ndim; % Total number of degrees of freedom

    % Computation of element stiffness matrices
    F_termic = zeros(Ndofs,1);
    F_e_global = zeros(Nelements,NnodesXelement*Ndim);
    
    for e=1:Nelements
        Area=A(e);  E=Es(e); delta_T = Delta_T(e);
        x1=x(T(e,1),1);  x2=x(T(e,2),1);
        y1=x(T(e,1),2);  y2=x(T(e,2),2);
        z1=x(T(e,1),3);  z2=x(T(e,2),3);
        
        le=sqrt((x2-x1)^2+ (y2-y1)^2 + (z2-z1)^2);
        
        Re = 1/le * [ x2-x1, y2-y1, z2-z1, 0, 0, 0 ;...  
            0, 0, 0, x2-x1, y2-y1, z2-z1];
 
        vt = [-1; 0; 0; 1; 0; 0];
        
        F_e_global(e,:) = Area*E*delta_T*alpha* Re'* [-1 1; 1 -1]*Re*vt;        
        
    end
    
for e=1:Nelements
        F_termic(CN(e,:)) = F_termic(CN(e,:)) + F_e_global(e,:)';  
end