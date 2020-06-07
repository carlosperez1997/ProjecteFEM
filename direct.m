%--------------------------------------------------------------------------
% Authors: Paula Sorolla and Carlos Pérez
% Date: 2nd April 2020
%--------------------------------------------------------------------------
% DIRECT FUNCTION:
% This funciton applies the direct method response
%--------------------------------------------------------------------------

function X_mod = direct(KG, MG, F, vl, lambda)
    
    Kll = KG(vl,vl);
    Mll = MG(vl,vl);
    Fl = F(vl,1);

    omega = linspace(0,1,100);
    X_mod = zeros(length(vl),length(omega));
    
    for ind = 1:length(omega)
        Q = -(Mll*omega(ind)^2) + Kll;
        X_mod(:,ind) = Q\Fl;
    end
    
    
    %% Plot the response for the Direct Method
    figure;
    semilogy(omega,abs(X_mod(2,:)./Fl(2)),'lineWidth',1);
    hold on;
    y1=get(gca,'ylim');
    semilogy([lambda(1) lambda(1)],y1,'k--');
    semilogy([lambda(2) lambda(2)],y1,'k--');
%     grid on;
    xlabel('\omega [rad/s]');
    ylabel('X/F');
    title('Direct method response');
end