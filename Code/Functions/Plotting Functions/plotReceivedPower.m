function plotReceivedPower(distances, total_RX_power_dBm, los_power_dBm)
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

    figure('Name', 'Received Power vs. Distance', 'NumberTitle', 'off', 'Position', [150 150 800 600]);
    
    % Use a semilogarithmic scale for the x-axis
    semilogx(distances, total_RX_power_dBm, 'b-', 'LineWidth', 2, 'DisplayName', 'Total Received Power (Multipath)');
    hold on;
    
    % Plot the theoretical LOS power as a dashed red line
    semilogx(distances, los_power_dBm, 'r--', 'LineWidth', 2, 'DisplayName', 'Friis Formula');
    
    grid on;
    title('Received Power vs. Distance ($K=5$ Reflections)', 'FontSize', 16, 'Interpreter', 'latex');
    xlabel('Distance, $d$ (m)', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('Received Power, $P_{RX}$ (dBm)', 'FontSize', 12, 'Interpreter', 'latex');
    legend('show', 'Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex');
    xlim([min(distances), max(distances)]);
    ylim([-110, -20]);
    hold off;
end