function X_vector = Convert2vector(x, T, X)

    % Dimensions
    Ndim = size(x,2); % Number of dimensions (DOFs for each node)
    Nnodes = size(x,1); % Number of nodes
    NnodesXelement = size(T,2); % Number of nodes for each element
    Nelements = size(T,1); % Number of elements
    Ndofs = Nnodes*Ndim; % Total number of degrees of freedom
    
    X_vector = zeros(Nnodes, Ndim);
    for i = 1:Nnodes
        for j = 1:Ndim
            I = Ndim*(i-1)+j;
            X_vector(i,j) = X(I);
        end
    end