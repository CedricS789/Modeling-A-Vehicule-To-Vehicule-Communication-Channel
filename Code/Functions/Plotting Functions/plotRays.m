function plotRays(walls, transmitter_position, receiver_position, all_rays_data, M)
% Visualizes the environment and all found ray paths.
%
% This function is dedicated to creating a 2D plot of the simulation results.
%
% INPUTS:
%   walls                - Struct array defining wall geometry.
%   transmitter_position - 1x2 vector [x, y] for the transmitter's position.
%   receiver_position    - 1x2 vector [x, y] for the receiver's position.
%   all_rays_data        - Cell array of structs with ray information.
%   M                    - Maximum number of reflections, used for the plot title.

    figure('Name', 'V2V Ray Tracing Results', 'NumberTitle', 'off');
    hold on;
    grid on;
    ax = gca;

    % Plot Walls
    for i = 1:length(walls)
        plot(ax, walls(i).coordinates(:,1), walls(i).coordinates(:,2), 'k', 'LineWidth', 7, 'DisplayName', 'Wall');
    end

    % Plot Transmitter and Receiver
    plot(ax, transmitter_position(1), transmitter_position(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', '$T_x$');
    plot(ax, receiver_position(1), receiver_position(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', '#0077BE', 'DisplayName', '$R_x$');
    
    % Plot All Found Rays
    colors = lines(M);
    
    for i = 1:length(all_rays_data)
        current_ray = all_rays_data{i};
        
        if strcmp(current_ray.type, 'LOS')
            plot(ax, current_ray.coordinates(:,1), current_ray.coordinates(:,2), 'g-', 'LineWidth', 2, 'DisplayName', 'LOS');
        else
            reflection_count = sscanf(current_ray.type, '%d-Refl');
            if ~isempty(reflection_count) && reflection_count <= size(colors, 1)
                plot_color = colors(reflection_count,:);
            else
                plot_color = 'm';
            end
            plot(ax, current_ray.coordinates(:,1), current_ray.coordinates(:,2), '-', 'Color', plot_color, 'LineWidth', 1, 'DisplayName', current_ray.type); 
        end
    end
    
    hold off;
    
    title(ax, sprintf('Ray Tracing with up to $M = %d$ Reflections', M), 'FontSize', 16, 'Interpreter', 'latex');
    xlabel(ax, '$x$-axis (m)', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel(ax, '$y$-axis (m)', 'FontSize', 14, 'Interpreter', 'latex');
    
    all_plot_handles = findobj(gcf,'Type','line');
    if ~isempty(all_plot_handles)
        [~, unique_indices] = unique(get(all_plot_handles, 'DisplayName'), 'stable');
        legend(all_plot_handles(unique_indices), 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
    end
end