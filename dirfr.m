function X_DF = dirfr( Knn, Mnn, vl, w_exc, f)
%--------------------------------------------------------------------------
% This function calculates the displacements of the structure using the
% direct frequency response.
%--------------------------------------------------------------------------
    X_DF = zeros(length(Knn),1);

  Q = Knn-Mnn*w_exc^2;
  X_DF = Q\f(vl); % Displacements
  
end