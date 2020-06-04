function [F] = assembly_F (x, T, CN, Fext)
    
    % Dimensions
    Ndim = size(x,2); % Number of dimensions (DOFs for each node)
    Nnodes = size(x,1); % Number of nodes
    NnodesXelement = size(T,2); % Number of nodes for each element
    Nelements = size(T,1); % Number of elements
    Ndofs = Nnodes*Ndim; % Total number of degrees of freedom
    
    % Global force vector assembly
    Nforces = size(Fext,1); % Number of imposed forces
    F = zeros(Ndofs,1);

    for k=1:Nforces
        I=Fext(k,1);
        i=Fext(k,2);
        p=Ndim*(I-1)+i;
        F(p,1)=F(p,1)+Fext(k,3);
    end
