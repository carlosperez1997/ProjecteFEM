function CN = CN_global(x, T)

% Dimensions
Ndim = size(x,2); % Number of dimensions (DOFs for each node)
Nnodes = size(x,1); % Number of nodes
NnodesXelement = size(T,2); % Number of nodes for each element
Nelements = size(T,1); % Number of elements
Ndofs = Nnodes*Ndim; % Total number of degrees of freedom

CN = zeros(Nelements, NnodesXelement*Ndim);

for i = 1:Nelements
    for j = 1:NnodesXelement
        for k = 1:Ndim
            J = Ndim*(j-1)+k;
            I = T(i,j);
            CN(i,J) = Ndim*(I-1)+k;
        end
    end
end
