function [KG, MG] = assembly_KG_MG(x, T, Es, A, CN, rhos)

    % Dimensions
    Ndim = size(x,2); % Number of dimensions (DOFs for each node)
    Nnodes = size(x,1); % Number of nodes
    NnodesXelement = size(T,2); % Number of nodes for each element
    Nelements = size(T,1); % Number of elements
    Ndofs = Nnodes*Ndim; % Total number of degrees of freedom

    % Computation of element stiffness matrices
    Kel = zeros(NnodesXelement*Ndim,NnodesXelement*Ndim,Nelements);
    for e=1:Nelements
        Area=A(e);  E=Es(e); rho = rhos(e);
        x1=x(T(e,1),1);  x2=x(T(e,2),1);
        y1=x(T(e,1),2);  y2=x(T(e,2),2);
        z1=x(T(e,1),3);  z2=x(T(e,2),3);
        %le=sqrt((x2-x1)^2+(y2-y1)^2);
        le=sqrt((x2-x1)^2+(y2-y1)^2 + (z2-z1)^2);
        %s=(y2-y1)/le;  %c=(x2-x1)/le;
        %Ke=Area*E/le*[c^2, c*s, -c^2, -c*s; c*s, s^2, -c*s, -s^2; -c^2, -c*s, c^2, c*s; -c*s, -s^2, c*s, s^2];
        Re = 1/le * [ x2-x1, y2-y1, z2-z1, 0, 0, 0 ;...
            0, 0, 0, x2-x1, y2-y1, z2-z1];
        K_e_prima = Area * E / le * [1 -1 ; -1 1];
        Ke = Re'*K_e_prima*Re;
        
        ML = rho*le*Area*0.5*eye(6);
        
        for r=1:NnodesXelement*Ndim
            for s=1:NnodesXelement*Ndim
                Kel(r,s,e)=Ke(r,s);
                Mel(r,s,e)=ML(r,s);
            end
        end
    end

    % Global Stiffness matrix assembly
    KG = zeros(Ndofs,Ndofs);
    MG = zeros(Ndofs,Ndofs);
    for e=1:Nelements
        for i=1:NnodesXelement*Ndim
            I=CN(e,i);
            for j=1:NnodesXelement*Ndim
                J=CN(e,j);
                KG(I,J)=KG(I,J)+Kel(i,j,e);
                MG(I,J)=MG(I,J)+Mel(i,j,e);
            end        
        end
    end
end

    
    
   
    
