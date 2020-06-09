function  X_MF  = modfr( phi, Mnn, Knn, w_exc, n, vl,f)
%--------------------------------------------------------------------------
% This function calculates the displacements of the structure using the
% modal frequency response.
%--------------------------------------------------------------------------
   
x = zeros(size(phi,1), length(w_exc));
F=f(vl);

m = phi(:,1:n)'*Mnn*phi(:,1:n);
k = phi(:,1:n)'*Knn*phi(:,1:n);
force = phi(:,1:n)'*F;
fprintf('Iteraci�n modal n�mero %d \n',n);
for i = 1:length(w_exc)
     
    q = k - w_exc(i)^2*m;
    e = q\force;
    X_MF(:,i) = phi(:,1:n)*e;
end
end