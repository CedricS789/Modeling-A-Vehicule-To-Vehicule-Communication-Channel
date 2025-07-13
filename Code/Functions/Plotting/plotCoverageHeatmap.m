function plotCoverageHeatmap(rx_x_coords, rx_y_coords, rx_power_heatmap_dBm, tx_pos, walls)
% plotCoverageHeatmap - Generates a 2D heatmap of the received power.
%
% INPUTS:
%   rx_x_coords          - Vector of x-coordinates for the heatmap grid.
%   rx_y_coords          - Vector of y-coordinates for the heatmap grid.
%   rx_power_heatmap_dBm - Matrix of received power values in dBm.
%   tx_pos               - 1x2 position of the transmitter.
%   walls                - Struct array defining wall geometry.

    figure('Name', '2D Received Power Coverage Map', 'NumberTitle', 'off', 'Position', [250 250 1000 450]);
    
    % Use imagesc to display the power matrix as an image.
    imagesc(rx_x_coords, rx_y_coords, rx_power_heatmap_dBm);
    
    % By default, the y-axis is inverted in imagesc. 'normal' fixes this.
    set(gca, 'YDir', 'normal');
    hold on;
    
    % Overlay the walls on the heatmap.
    for i = 1:length(walls)
        plot(walls(i).coords(:,1), walls(i).coords(:,2), 'w-', 'LineWidth', 3);
    end
    
    % Overlay the transmitter position.
    plot(tx_pos(1), tx_pos(2), 'w^', 'MarkerSize', 12, 'MarkerFaceColor', '#D95319', 'DisplayName', 'Transmitter');
    
    title('Received Power Coverage Map');
    xlabel('X Coordinate (m)');
    ylabel('Y Coordinate (m)');
    
    % Add and label the color bar.
    c = colorbar;
    ylabel(c, 'Received Power (dBm)');
    
    % Set the color axis limits to a reasonable range for dBm values.
    caxis([-90 -30]);
    
    axis equal tight; % Use tight, equal aspect ratio for the axes.
    xlim([min(rx_x_coords), max(rx_x_coords)]);
    ylim([min(rx_y_coords), max(rx_y_coords)]);
    
    % Use a perceptually uniform colormap.
    colormap('jet');
    
    hold off;
end
