function plotKFactorVsDistance(distances, K_factors_dB)
% PLOTKFACTORvsDISTANCE - Plots the Rician K-factor in dB vs. distance.
%
% INPUTS:
%   distances    - Vector of distances for the x-axis.
%   K_factors_dB - Vector of Rician K-factors in dB.

    figure('Name', 'Rician K-Factor vs. Distance', 'NumberTitle', 'off', 'Position', [150 150 800 600]);
    semilogx(distances, K_factors_dB, 'k-', 'LineWidth', 2);
    grid on;
    title('Rician K-Factor vs. Distance');
    xlabel('Separation Distance d (m)');
    ylabel('Rician K-Factor (dB)');
    xlim([1, 1000]);
end
