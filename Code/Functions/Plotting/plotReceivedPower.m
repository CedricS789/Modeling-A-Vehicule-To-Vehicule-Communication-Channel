function plotReceivedPower(distances, total_rx_power_dBm, los_power_dBm)
% plotReceivedPower - Plots received power vs. distance.
%
% This function plots the total received power from all multipath components
% and compares it against the power from the Line-of-Sight (LOS) path alone,
% which follows the Friis transmission formula.
%
% INPUTS:
%   distances          - Vector of distances for the x-axis.
%   total_rx_power_dBm - Vector of total received power in dBm (multipath).
%   los_power_dBm      - Vector of LOS-only received power in dBm (Friis).

    figure('Name', 'Received Power vs. Distance', 'NumberTitle', 'off', 'Position', [150 150 800 600]);
    
    % Plot the total received power from the full channel model
    semilogx(distances, total_rx_power_dBm, 'b-', 'LineWidth', 2, 'DisplayName', 'Total Received Power (Multipath)');
    hold on;
    
    % Plot the LOS power, which represents the Friis formula
    semilogx(distances, los_power_dBm, 'r--', 'LineWidth', 2, 'DisplayName', 'LOS Power (Friis Formula)');
    
    grid on;
    title('Received Power vs. Separation Distance');
    xlabel('Separation Distance, d (m)');
    ylabel('Received Power, P_{RX} (dBm)');
    legend('show', 'Location', 'northeast');
    xlim([min(distances), max(distances)]);
    ylim([-110, -20]); % Adjust ylim for better visualization
    hold off;
end
