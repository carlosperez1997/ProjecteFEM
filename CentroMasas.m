function CG = CentroMasas(x, rhos, A, T)

CG = zeros(3,1);

    Ndim = size(x,2); % Number of dimensions (DOFs for each node)
    Nnodes = size(x,1); % Number of nodes
    NnodesXelement = size(T,2); % Number of nodes for each element
    Nelements = size(T,1); % Number of elements
    Ndofs = Nnodes*Ndim; % Total number of degrees of freedom
    
for e = 1:Nelements

   Area=A(e);  rho = rhos(e);
        x1=x(T(e,1),1); x2=x(T(e,2),1);
        y1=x(T(e,1),2); y2=x(T(e,2),2);
        z1=x(T(e,1),3); z2=x(T(e,2),3);
        %le=sqrt((x2-x1)^2+(y2-y1)^2);
   le=sqrt((x2-x1)^2+(y2-y1)^2 + (z2-z1)^2);
        
   peso(e) = rho*le*Area;
        
   x_medio(e) = (x1 + x2 )/ 2;
   y_medio(e) = (y1 + y2 )/ 2;
   z_medio(e) = (z1 + z2 )/ 2;
        
end

num_x = 0; num_y = 0; num_z = 0;
for e = 1:Nelements
    num_x = num_x + x_medio(e) * peso(e);
    num_y = num_y + y_medio(e) * peso(e);
    num_z = num_z + z_medio(e) * peso(e);
end

CG(1) = num_x / sum(peso)
CG(2) = num_y / sum(peso)
CG(3) = num_z / sum(peso)

CG(1) = 0;
