function plotRays(walls, transmitter_position, receiver_position, all_rays_data, max_reflection_order)
% plotRays - Visualizes the environment and all found ray rays.
%
% This function is dedicated to creating a 2D plot of the simulation results.
%
% INPUTS:
%   walls                - Struct array defining wall geometry.
%   transmitter_position - 1x2 vector [x, y] for the transmitter's position.
%   receiver_position    - 1x2 vector [x, y] for the receiver's position.
%   all_rays_data       - Cell array of structs with ray information.
%   max_reflection_order - Maximum number of reflections, used for the plot title.

    figure('Name', 'V2V Ray Tracing Results', 'NumberTitle', 'off');
    hold on; 
    grid on; 
    ax = gca;

    % --- 1. Plot Walls ---
    for i = 1:length(walls)
        plot(ax, walls(i).coords(:,1), walls(i).coords(:,2), 'k', 'LineWidth', 2.5, 'DisplayName', 'wall');
    end

    % --- 2. Plot Transmitter and Receiver ---
    plot(ax, transmitter_position(1), transmitter_position(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', '#0077BE', 'DisplayName', 'TX');
    plot(ax, receiver_position(1), receiver_position(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', '#D95319', 'DisplayName', 'RX');

    % --- 3. Plot All Found rays ---
    colors = lines(max_reflection_order); 
    
    for i = 1:length(all_rays_data)
        current_ray = all_rays_data{i};
        
        if strcmp(current_ray.type, 'LOS')
            plot(ax, current_ray.path(:,1), current_ray.path(:,2), 'g-', 'LineWidth', 1.5, 'DisplayName', 'LOS');
        else
            reflection_count = sscanf(current_ray.type, '%d-Refl');
            if ~isempty(reflection_count) && reflection_count <= size(colors, 1)
                plot_color = colors(reflection_count,:);
            else
                plot_color = 'm';
            end
            plot(ax, current_ray.path(:,1), current_ray.path(:,2), '-', 'Color', plot_color, 'LineWidth', 1, 'DisplayName', current_ray.type);
        end
    end
    
    hold off;
    
    % --- 4. Add Labels and a Clean Legend ---
    title(ax, sprintf('Ray Tracing with up to %d Reflections', max_reflection_order));
    xlabel(ax, 'X Coordinate (m)');
    ylabel(ax, 'Y Coordinate (m)');
    
    all_plot_handles = findobj(gcf,'Type','line');
    if ~isempty(all_plot_handles)
        [~, unique_indices] = unique(get(all_plot_handles, 'DisplayName'), 'stable');
        legend(all_plot_handles(unique_indices));
    end
end
