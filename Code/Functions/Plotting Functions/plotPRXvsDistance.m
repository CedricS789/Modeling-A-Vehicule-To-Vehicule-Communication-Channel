function plotPRXvsDistance(distances, total_RX_power_dBm, los_power_dBm)
% Plots received power vs. distance.
%
% This function plots the total received power from all multipath components
% and compares it against the power from the Line-of-Sight (LOS) path alone,
% which follows the Friis transmission formula.
%
% INPUTS:
%   distances          - Vector of distances for the x-axis.
%   total_RX_power_dBm - Vector of total received power in dBm (multipath).
%   los_power_dBm      - Vector of LOS-only received power in dBm (Friis).

    figure('Name', 'Received Power vs. Distance', 'NumberTitle', 'off');
    
    % Use a semilogarithmic scale for the x-axis
    semilogx(distances, total_RX_power_dBm, 'b-', 'LineWidth', 2, 'DisplayName', 'Total Received Power');
    hold on;
    
    % Plot the theoretical LOS power as a dashed red line
    semilogx(distances, los_power_dBm, 'r--', 'LineWidth', 2, 'DisplayName', 'Friis Formula');
    
    grid on;
    title('$P_{RX} = f(d)$', 'FontSize', 20, 'Interpreter', 'latex');
    xlabel('$d$ (m)', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$P_{RX}$ (dBm)', 'FontSize', 18, 'Interpreter', 'latex');
    legend('show', 'Location', 'northeast', 'FontSize', 18, 'Interpreter', 'latex');
    axis tight;
    hold off;
end