function X_DF = dirfr( Knn, Mnn, vl, w_exc, f)
%--------------------------------------------------------------------------
% This function calculates the displacements of the structure using the
% direct frequency response.
%--------------------------------------------------------------------------
    X_DF = zeros(length(Knn),length(w_exc));
 
for i = 1:length(w_exc)
  Q = Knn-Mnn*(w_exc(i))^2;
  X_DF(:,i) = Q\f(vl); % Displacements
  
end

end