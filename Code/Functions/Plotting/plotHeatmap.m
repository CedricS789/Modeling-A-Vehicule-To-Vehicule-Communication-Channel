function plotHeatmap(rx_x_coords, rx_y_coords, Prx_heatmap_dBm, tx_pos, walls)
% PLOTHEATMAP - Generates a 2D heatmap of the received power.
%
% INPUTS:
%   rx_x_coords     - Vector of x-coordinates for the heatmap grid.
%   rx_y_coords     - Vector of y-coordinates for the heatmap grid.
%   Prx_heatmap_dBm - Matrix of received power values in dBm.
%   tx_pos          - 1x2 position of the transmitter.
%   walls           - Struct array defining wall geometry.

    figure('Name', 'Received Power Heatmap', 'NumberTitle', 'off', 'Position', [300 300 1000 400]);
    
    imagesc(rx_x_coords, rx_y_coords, Prx_heatmap_dBm);
    set(gca, 'YDir', 'normal');
    hold on;
    
    for i = 1:length(walls)
        plot(walls(i).coords(:,1), walls(i).coords(:,2), 'w-', 'LineWidth', 3);
    end
    
    plot(tx_pos(1), tx_pos(2), 'w^', 'MarkerSize', 12, 'MarkerFaceColor', 'red');
    
    title('Received Power Heatmap in the Urban Canyon');
    xlabel('X Coordinate (m)');
    ylabel('Y Coordinate (m)');
    c = colorbar;
    ylabel(c, 'Received Power (dBm)');
    caxis([-90 -40]);
    axis equal tight;
    xlim([min(rx_x_coords), max(rx_x_coords)]);
    ylim([min(rx_y_coords), max(rx_y_coords)]);
    colormap('jet');
    hold off;
end
