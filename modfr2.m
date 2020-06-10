function  X_MF  = modfr2( phi, Mnn, Cnn, Knn, chi, wn, w_exc, vl,f)
%--------------------------------------------------------------------------
% This function calculates the displacements of the structure using the
% modal frequency response.
%--------------------------------------------------------------------------
    
    for i = 1:size(phi,1)
        x = zeros(size(phi,1), length(w_exc));
        F=f(vl);
    
        m = phi(:,i)'*Mnn*phi(:,i);
        b = phi(:,i)'*Cnn*phi(:,i);
        k = phi(:,i)'*Knn*phi(:,i);

        force = phi(:,i)'*F;
    
        q = - w_exc^2*m + k + i*b*w_exc;
        e = force/q;
        
        x(i) = force/q;
        
    end
    
    X_MF(:) = phi'*x;
    
    a = 1;
end