function plotCoverageHeatmap(ax, rx_x_coords, rx_y_coords, Prx_heatmap_dBm, tx_pos, walls, sens_dBm)
% PLOTCOVERAGEHEATMAP - Generates a 2D heatmap of the received power on a given axes.
%
% INPUTS:
%   ax              - Axes handle to plot on. If empty, a new figure is created.
%   rx_x_coords     - Vector of x-coordinates for the heatmap grid.
%   rx_y_coords     - Vector of y-coordinates for the heatmap grid.
%   Prx_heatmap_dBm - Matrix of received power values in dBm.
%   tx_pos          - 1x2 position of the transmitter.
%   walls           - Struct array defining wall geometry.
%   sens_dBm        - Receiver sensitivity in dBm. Points below this will not be shown.

    % If no axes handle is provided, create a new figure.
    if isempty(ax)
        figure('Name', 'Received Power Heatmap', 'NumberTitle', 'off', 'Position', [300 300 1000 400]);
        ax = gca;
    end
    
    % Plot the heatmap image on the specified axes
    h_img = imagesc(ax, rx_x_coords, rx_y_coords, Prx_heatmap_dBm);
    set(ax, 'YDir', 'normal');
    
    % Apply transparency for values below sensitivity
    set(h_img, 'AlphaData', (Prx_heatmap_dBm >= sens_dBm));
    set(ax, 'Color', 'w'); 
    
    hold(ax, 'on');
    
    % Plot walls
    for i = 1:length(walls)
        plot(ax, walls(i).coords(:,1), walls(i).coords(:,2), 'k-', 'LineWidth', 2);
    end
    
    % Plot transmitter
    plot(ax, tx_pos(1), tx_pos(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
    
    % Set titles, labels, and colorbar for the specified axes
    title(ax, 'Received Power Heatmap');
    xlabel(ax, 'X Coordinate (m)');
    ylabel(ax, 'Y Coordinate (m)');
    c = colorbar(ax);
    ylabel(c, 'Received Power (dBm)');
    caxis(ax, [-90 -40]); % Adjust color axis limits as needed
    
    % Configure axes properties
    axis(ax, 'equal', 'tight');
    xlim(ax, [min(rx_x_coords), max(rx_x_coords)]);
    ylim(ax, [min(rx_y_coords), max(rx_y_coords)]);
    colormap(ax, 'jet');
    
    hold(ax, 'off');
end