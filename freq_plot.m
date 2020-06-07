function [] = freq_plot(w_exc,vecnorm1,vecnorm2,phi)
% FREQUENCY PLOTS
    
    %Direct Frequency respose 
    figure
    semilogy(w_exc/(2*pi), vecnorm1,'k');
    xlabel("Freq [Hz]", 'interpreter', 'latex')
    ylabel("|X| [m]", 'interpreter', 'latex')
    legend("DFR", 'interpreter', 'latex')
    xlim([0 max(w_exc/(2*pi))])
    
%     % Modal frequency response using the first eigenvalue
%     figure
%     semilogy(w_exc/(2*pi), vecnorm2(1,:),'k');
%     xlabel("Freq [Hz]", 'interpreter', 'latex')
%     ylabel("|X| [m]", 'interpreter', 'latex')
%     legend("MFR - 1st", 'interpreter', 'latex')
%     xlim([0 max(w_exc/(2*pi))])
    
    % Modal frequency response considering all eigenvalues
    figure
    semilogy(w_exc/(2*pi), vecnorm2,'k');
    xlabel("Freq [Hz]", 'interpreter', 'latex')
    ylabel("|X| [m]", 'interpreter', 'latex')
    legend("MFR - All Eigenvalues", 'interpreter', 'latex')
    xlim([0 max(w_exc/(2*pi))])
    
    % Methods comparison
    figure
    semilogy(w_exc/(2*pi), vecnorm1,'k');
    hold on
    semilogy(w_exc/(2*pi), vecnorm2(length(phi),:),'r--');
    xlabel("Freq [Hz]", 'interpreter', 'latex')
    ylabel("|X| [m]", 'interpreter', 'latex')
    legend("DFR","MFR", 'interpreter', 'latex')
    hold off
end