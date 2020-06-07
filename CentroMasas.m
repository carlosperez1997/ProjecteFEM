function [CG, M] = CentroMasas(x, rhos, A, T, punctual_mass)

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

if isempty(punctual_mass) 
    M = sum(peso);
else
    
    for j = 1:size(punctual_mass,1)
        num_x = num_x + x(punctual_mass(j,1),1)*punctual_mass(j,2);
        num_y = num_y + x(punctual_mass(j,1),2)*punctual_mass(j,2);
        num_z = num_z + x(punctual_mass(j,1),3)*punctual_mass(j,2);
    end
    
    M = sum(peso) + sum(punctual_mass(:,2));
end

CG(1) = num_x / M;
CG(2) = num_y / M;
CG(3) = num_z / M;



%CG(1) = 0;
