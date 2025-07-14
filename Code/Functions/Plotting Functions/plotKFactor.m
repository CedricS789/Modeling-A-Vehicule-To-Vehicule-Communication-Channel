function plotKFactor(distances, K_factor_dB)
% Plots the Rician K-factor in dB versus distance.
%
% The K-factor is the ratio of the power in the dominant, line-of-sight (LOS)
% component to the power in the NLOS components.
% A high K-factor indicates a channel dominated by the LOS path.
%
% INPUTS:
%   distances   - Vector of distances for the x-axis.
%   K_factor_dB - Vector of Rician K-factors in dB.

    figure('Name', 'Rician K-Factor vs. Distance', 'NumberTitle', 'off');
    
    semilogx(distances, K_factor_dB, 'k-', 'LineWidth', 2);
    
    grid on;
    title('Rician K-Factor vs. Distance', 'FontSize', 16);
    xlabel('Distance, d (m)', 'FontSize', 14);
    ylabel('Rician K-Factor (dB)', 'FontSize', 14);
    xlim([min(distances), max(distances)]);

end
