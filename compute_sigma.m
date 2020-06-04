function [sig] = compute_sigma (x, T, CN, E, U)

    % Dimensions
    Ndim = size(x,2); % Number of dimensions (DOFs for each node)
    Nnodes = size(x,1); % Number of nodes
    NnodesXelement = size(T,2); % Number of nodes for each element
    Nelements = size(T,1); % Number of elements
    Ndofs = Nnodes*Ndim; % Total number of degrees of freedom

    % Compute stress for each element
    sig = zeros(Nelements,1);
    de=zeros(NnodesXelement*Ndim,1);
    for e=1:Nelements
        x1e=x(T(e,1),1);
        x2e=x(T(e,2),1);
        y1e=x(T(e,1),2);
        y2e=x(T(e,2),2);
        z1e=x(T(e,1),3);
        z2e=x(T(e,2),3);
        
        le=sqrt((x2e-x1e)^2+(y2e-y1e)^2+(z2e-z1e)^2);
        Re=1/le*[ x2e-x1e, y2e-y1e, z2e-z1e, 0,0,0;
           0, 0, 0,x2e-x1e, y2e-y1e, z2e-z1e];
    
        %obtain displacements
        for r=1:NnodesXelement*Ndim
               p=CN(e,r);
               de(r)=U(p);
        end
 
        dprim=Re*de;
        %calculate Strain and Stress
        Epsilon=(1/le)*[-1, 1]*dprim;
        sig(e,1)= E(e) * Epsilon;
    end
    
    