function F_gravity = Gravity_force (x, T, A, rhos, g)

    % Dimensions
    Ndim = size(x,2); % Number of dimensions (DOFs for each node)
    Nnodes = size(x,1); % Number of nodes
    NnodesXelement = size(T,2); % Number of nodes for each element
    Nelements = size(T,1); % Number of elements
    Ndofs = Nnodes*Ndim; % Total number of degrees of freedom

    % Computation of element stiffness matrices
    F_gravity = zeros(Nnodes,3);
    
    i = 1;
    
    for e=1:Nelements
        Area=A(e); rho = rhos(e);
        x1=x(T(e,1),1);  x2=x(T(e,2),1);
        y1=x(T(e,1),2);  y2=x(T(e,2),2);
        z1=x(T(e,1),3);  z2=x(T(e,2),3);
        %le=sqrt((x2-x1)^2+(y2-y1)^2);
        le=sqrt((x2-x1)^2+ (y2-y1)^2 + (z2-z1)^2);

        peso = Area*rho*le*g;
        
        F_gravity (T(e,1),:) = [T(e,1), 3, -peso/2];
        F_gravity (T(e,2),:) = [T(e,2), 3, -peso/2];        
        
    end