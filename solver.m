function [U, R] = solver (x, T, KG, F, fixnodes);

    % Dimensions
    Ndim = size(x,2); % Number of dimensions (DOFs for each node)
    Nnodes = size(x,1); % Number of nodes
    NnodesXelement = size(T,2); % Number of nodes for each element
    Nelements = size(T,1); % Number of elements
    Ndofs = Nnodes*Ndim; % Total number of degrees of freedom
    
    % Restricted degrees of freedom
    NfixNod = size(fixnodes,1); % Number of fixed degrees of freedom
    vr = zeros(1,NfixNod);
    ur = zeros(NfixNod,1);
    
    for k=1:NfixNod
        A = fixnodes(k,1);
        i = fixnodes(k,2);
        vr(1,k) = Ndim*(A-1)+i;
        ur(k,1) = fixnodes(k,3);
    end

    % Free degrees of freedom
    vl = setdiff(1:Ndofs,vr);


    % Reduced system matrices
    Kll = KG(vl,vl);
    Klr = KG(vl,vr);
    Krl = KG(vr,vl);
    Krr = KG(vr,vr);
    fl = F(vl,1);
    fr = F(vr,1);

    % Solve linear system
    ul = Kll\(fl - Klr*ur);
    Rr = Krr*ur + Krl*ul - fr;

    % Assembly of global displacements
    U(vl) = ul;
    U(vr) = ur;

    % Assembly of global reactions
    R(vl) = 0;
    R(vr) = Rr;