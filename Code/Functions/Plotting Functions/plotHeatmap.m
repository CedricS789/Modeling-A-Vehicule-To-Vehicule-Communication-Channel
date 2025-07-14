function plotHeatmap(ax, RX_x_coordinates, RX_y_coordinates, PRX_dBm_domain, tx_pos, walls, sens_dBm)
% Generates a 2D heatmap of the received power on a given axes.
%
% INPUTS:
%   ax              - Axes handle to plot on. If empty, a new figure is created.
%   RX_x_coordinates- Vector of x-coordinates for the heatmap grid.
%   RX_y_coordinates- Vector of y-coordinates for the heatmap grid.
%   PRX_heatmap_dBm - Matrix of received power values in dBm.
%   tx_pos          - 1x2 position of the transmitter.
%   walls           - Struct array defining wall geometry.
%   sens_dBm        - Receiver sensitivity in dBm. Points below this will not be shown.

    if isempty(ax)
        figure('Name', 'Received Power Heatmap', 'NumberTitle', 'off', 'Position', [300 300 1000 400]);
        ax = gca;
    end
    
    h_img = imagesc(ax, RX_x_coordinates, RX_y_coordinates, PRX_dBm_domain);
    set(ax, 'YDir', 'normal');
    set(h_img, 'AlphaData', (PRX_dBm_domain >= sens_dBm));
    set(ax, 'Color', 'w');
    hold(ax, 'on');
    
    for i = 1:length(walls)
        plot(ax, walls(i).coordinates(:,1), walls(i).coordinates(:,2), 'k-', 'LineWidth', 2);
    end
    
    plot(ax, tx_pos(1), tx_pos(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', '$T_x$');
    
    title(ax, 'Received Power Heatmap', 'Interpreter', 'latex');
    xlabel(ax, '$x$ axis (m)', 'Interpreter', 'latex');
    ylabel(ax, '$y$ axis (m)', 'Interpreter', 'latex');
    c = colorbar(ax);
    ylabel(c, 'Received Power (dBm)', 'Interpreter', 'latex');
    caxis(ax, [-90 -40]);
    
    axis(ax, 'equal', 'tight');
    xlim(ax, [min(RX_x_coordinates), max(RX_x_coordinates)]);   
    ylim(ax, [min(RX_y_coordinates), max(RX_y_coordinates)]);
    
    % Use 256 color levels from the 'jet' colormap for a finer gradient
    colormap(ax, jet(256));
    
    hold(ax, 'off');
end