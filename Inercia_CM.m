function [I_CM] = Inercia_CM(x, rhos, A, T, punctual_mass, CG)

I_CM = zeros(3);

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

% Inercias

num_xx = 0; num_xy = 0; num_xz = 0;
num_yx = 0; num_yy = 0; num_yz = 0;
num_zx = 0; num_zy = 0; num_zz = 0;

for e = 1:Nelements
    
    x_dist_cm(e) = x_medio(e) - CG(1);
    y_dist_cm(e) = y_medio(e) - CG(2);
    z_dist_cm(e) = z_medio(e) - CG(3);
    
    num_xx = num_xx + (y_dist_cm(e)^2 + z_dist_cm(e)^2 )* peso(e);
    num_xy = num_xy - x_dist_cm(e) * y_dist_cm(e) * peso(e);
    num_xz = num_xz - x_dist_cm(e) * z_dist_cm(e) * peso(e);
    
    num_yx = num_yx - y_dist_cm(e) * x_dist_cm(e) * peso(e);
    num_yy = num_yy + (x_dist_cm(e)^2 + z_dist_cm(e)^2) * peso(e);
    num_yz = num_yz + y_dist_cm(e) * z_dist_cm(e) * peso(e);
    
    num_zx = num_zx - z_dist_cm(e) * x_dist_cm(e) * peso(e);
    num_zy = num_zy - z_dist_cm(e) * y_dist_cm(e) * peso(e);
    num_zz = num_zz + (x_dist_cm(e)^2 + y_dist_cm(e)^2) * peso(e);
    
end



for j = 1:size(punctual_mass,1)
    
    x_dist_cm = x(punctual_mass(j,1),1) - CG(1);
    y_dist_cm = x(punctual_mass(j,1),2) - CG(2);
    z_dist_cm = x(punctual_mass(j,1),3) - CG(3);
    masa = punctual_mass(j,2);
    
    num_xx = num_xx + (y_dist_cm^2 + z_dist_cm^2) * masa; 
    num_xy = num_xy - (x_dist_cm * y_dist_cm) * masa; 
    num_xz = num_xz - (x_dist_cm * z_dist_cm) * masa; 

    num_yx = num_yx - (y_dist_cm * x_dist_cm) * masa; 
    num_yy = num_yy + (x_dist_cm^2 + z_dist_cm^2) * masa; 
    num_yz = num_yz - (y_dist_cm * z_dist_cm) * masa; 

    num_zx = num_zx - (z_dist_cm * x_dist_cm) * masa; 
    num_zy = num_zy - (z_dist_cm * y_dist_cm) * masa; 
    num_zz = num_zz + (x_dist_cm^2 + y_dist_cm^2) * masa; 
end

I_CM = [num_xx, num_xy, num_xz; ...
    num_yx, num_yy, num_yz; ...
    num_zx, num_zy, num_zz ];
