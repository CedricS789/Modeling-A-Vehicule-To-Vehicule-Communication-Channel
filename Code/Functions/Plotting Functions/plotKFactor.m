function plotKFactor(distances, k_factor_dB)
% plotKFactor - Plots the Rician K-factor in dB versus distance.
%
% The K-factor is the ratio of the power in the dominant, line-of-sight (LOS)
% component to the power in the scattered, non-line-of-sight (NLOS) components.
% A high K-factor indicates a channel dominated by the LOS path.
%
% INPUTS:
%   distances   - Vector of distances for the x-axis.
%   k_factor_dB - Vector of Rician K-factors in dB.

    figure('Name', 'Rician K-Factor vs. Distance', 'NumberTitle', 'off', 'Position', [300 300 800 600]);
    
    semilogx(distances, k_factor_dB, 'k-', 'LineWidth', 2);
    
    grid on;
    title('Rician K-Factor vs. Separation Distance');
    xlabel('Separation Distance, d (m)');
    ylabel('Rician K-Factor (dB)');
    xlim([min(distances), max(distances)]);
end
