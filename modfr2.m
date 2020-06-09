function  X_MF  = modfr2( phi, Mnn, Bnn, Knn, w_exc, n, vl,f)
%--------------------------------------------------------------------------
% This function calculates the displacements of the structure using the
% modal frequency response.
%--------------------------------------------------------------------------
   
    x = zeros(size(phi,1), length(w_exc));
    F=f(vl);

    m = phi(:,1:n)'*Mnn*phi(:,1:n);
    b = phi(:,1:n)'*Bnn*phi(:,1:n);
    k = phi(:,1:n)'*Knn*phi(:,1:n);

    force = phi(:,1:n)'*F;
     
    q = - w_exc^2*m + k + i*b*w_exc;
    e = q\force;
    X_MF(:) = phi(:,1:n)*e;
    
end