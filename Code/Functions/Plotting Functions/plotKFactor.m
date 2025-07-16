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
    grid minor;
    
    title('$K = f(d)$', 'FontSize', 20, 'Interpreter', 'latex');
    xlabel('$d$ (m)', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$K$-factor (dB)', 'FontSize', 18, 'Interpreter', 'latex');
    % xlim([min(distances), max(distances)]);
    
    axis tight;
end
