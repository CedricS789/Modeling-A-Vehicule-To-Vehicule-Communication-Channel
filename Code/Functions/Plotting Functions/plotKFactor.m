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

    figure('Name', 'Rician K-Factor vs. Distance', 'NumberTitle', 'off', 'Position', [300 300 800 600]);
    
    semilogx(distances, K_factor_dB, 'k-', 'LineWidth', 2);
    
    grid on;
    title('Rician K-Factor vs. Distance (5 Reflections)', 'FontSize', 16);
    xlabel('Distance, d (m)', 'FontSize', 12);
    ylabel('Rician K-Factor (dB)', 'FontSize', 12);
    xlim([min(distances), max(distances)]);
end
