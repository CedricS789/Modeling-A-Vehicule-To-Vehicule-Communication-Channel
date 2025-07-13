function plotCoverageHeatmap(ax, rx_x_coordinates, rx_y_coordinates, Prx_heatmap_dBm, tx_pos, walls, sens_dBm)
% PLOTCOVERAGEHEATMAP - Generates a 2D heatmap of the received power on a given axes.
%
% INPUTS:
%   ax              - Axes handle to plot on. If empty, a new figure is created.
%   rx_x_coordinates- Vector of x-coordinates for the heatmap grid.
%   rx_y_coordinates- Vector of y-coordinates for the heatmap grid.
%   Prx_heatmap_dBm - Matrix of received power values in dBm.
%   tx_pos          - 1x2 position of the transmitter.
%   walls           - Struct array defining wall geometry.
%   sens_dBm        - Receiver sensitivity in dBm. Points below this will not be shown.

    if isempty(ax)                                                                        % Allows the function to either draw on an existing plot or create its own.
        figure('Name', 'Received Power Heatmap', 'NumberTitle', 'off', 'Position', [300 300 1000 400]);
        ax = gca;                                                                         % Get the handle to the current axes to direct all plotting commands.
    end
    
    h_img = imagesc(ax, rx_x_coordinates, rx_y_coordinates, Prx_heatmap_dBm);              % Display the power data as a scaled color image on the specified axes.
    set(ax, 'YDir', 'normal');                                                             % Ensure the y-axis direction is bottom-up, matching standard Cartesian coordinates.
    
    % Create a logical mask to make areas with power below the receiver's
    % sensitivity level transparent, effectively showing coverage gaps.
    set(h_img, 'AlphaData', (Prx_heatmap_dBm >= sens_dBm));
    set(ax, 'Color', 'w');                                                                % Set the background for transparent areas to white for better visibility.
    
    hold(ax, 'on');                                                                       % Retain the current plot so that geometry can be overlaid.
    
    % Draw the physical environment on top of the heatmap.
    for i = 1:length(walls)
        plot(ax, walls(i).coordinates(:,1), walls(i).coordinates(:,2), 'k-', 'LineWidth', 2); % Plot each wall as a thick black line for clarity.
    end
    
    % Mark the location of the signal source.
    plot(ax, tx_pos(1), tx_pos(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'); % Use a distinct red circle marker for the transmitter.
    
    % Add descriptive labels to make the plot informative.
    title(ax, 'Received Power Heatmap');
    xlabel(ax, 'X Coordinate (m)');
    ylabel(ax, 'Y Coordinate (m)');
    c = colorbar(ax);                                                                     % Add a color bar to serve as a legend for the power levels.
    ylabel(c, 'Received Power (dBm)');
    caxis(ax, [-90 -40]);                                                                 % Set fixed color limits for consistent visualization across different plots.
    
    % Configure the axes for a clean, true-to-scale presentation.
    axis(ax, 'equal', 'tight');                                                            % 'equal' prevents geometric distortion; 'tight' removes extra white space.
    xlim(ax, [min(rx_x_coordinates), max(rx_x_coordinates)]);                              % Fit the plot's x-limits to the simulation grid.
    ylim(ax, [min(rx_y_coordinates), max(rx_y_coordinates)]);                              % Fit the plot's y-limits to the simulation grid.
    colormap(ax, 'jet');                                                                   % Apply a perceptually uniform colormap for intuitive data interpretation.
    
    hold(ax, 'off');                                                                      % Release the plot hold as good practice for subsequent plotting commands.
end