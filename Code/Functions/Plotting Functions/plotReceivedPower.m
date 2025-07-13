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

    figure('Name', 'Received Power vs. Distance', 'NumberTitle', 'off', 'Position', [150 150 800 600]); % Create a new, dedicated window for the plot with a descriptive title and size.
    
    % Use a semilogarithmic scale for the x-axis (distance) to clearly visualize
    % power trends over several orders of magnitude of distance.
    semilogx(distances, total_rx_power_dBm, 'b-', 'LineWidth', 2, 'DisplayName', 'Total Received Power (Multipath)');
    hold on;                                                                                          % Retain the current plot to overlay the next data series.
    
    % Plot the theoretical LOS power as a dashed red line to serve as a clear
    % baseline for comparison against the more complex multipath model.
    semilogx(distances, los_power_dBm, 'r--', 'LineWidth', 2, 'DisplayName', 'LOS Power (Friis Formula)');
    
    grid on;                                                                                          % Add a grid to make it easier to read power values at specific distances.
    title('Received Power vs. Separation Distance');                                                  % Set a clear and descriptive title for the plot.
    xlabel('Separation Distance, d (m)');                                                             % Label the x-axis to indicate the physical quantity and its units.
    ylabel('Received Power, P_{RX} (dBm)');                                                           % Label the y-axis to indicate the physical quantity and its units.
    legend('show', 'Location', 'northeast');                                                          % Display the legend in the top-right corner to identify each line.
    xlim([min(distances), max(distances)]);                                                           % Ensure the plot's x-axis spans the exact range of simulated distances.
    ylim([-110, -20]);                                                                                % Set a fixed y-axis range to focus on the most relevant power levels.
    hold off;                                                                                         % Release the plot hold, a good practice for subsequent plotting commands.
end