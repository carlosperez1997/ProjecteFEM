%--------------------------------------------------------------------------
% Authors: Paula Sorolla and Carlos Pérez
% Date: 2nd April 2020
%--------------------------------------------------------------------------
% MODAL FUNCTION:
% This funciton applies the modal method response
%--------------------------------------------------------------------------

function X_mod_xi = modal(KG, MG, F, vl, phi, lambda)
    
    Kll = KG(vl,vl);
    Mll = MG(vl,vl);
    Fl = F(vl,1);
    
    omega = linspace(0,5,sqrt(max(lambda))*1.5);
%     X_mod_xi = zeros(length(vl),length(omega));
    
    % Modal method for lowest eig vector
    Mxi = phi(:,2)' * Mll * phi(:,2);
    Kxi = phi(:,2)' * Kll * phi(:,2);
    Fxi = phi(:,2)' * Fl;
    xi = zeros(1,length(omega));
    for ind = 1:length(omega)
        Q = -(omega(ind)^2 * Mxi) + Kxi;
        xi(ind) = Fxi / Q;
    end
    X_mod_xi = phi(:,2)*xi;

    
    %% Plot the response for the Modal Method
    figure;
    semilogy(omega,abs(X_mod_xi(2,:)./Fl(2)),'lineWidth',1);
    hold on;
    y1=get(gca,'ylim');
%     semilogy([lambda(1) lambda(1)],y1,'k--');
    semilogy([lambda(2) lambda(2)],y1,'k--');
    xlabel('\omega [rad/s]');
    ylabel('X/F');
    title('Modal method response for lowest eigenvector');
end