function plotRays(walls, transmitter_position, receiver_position, all_rays_data, K)
% plotRays - Visualizes the environment and all found ray paths.
%
% This function is dedicated to creating a 2D plot of the simulation results.
%
% INPUTS:
%   walls                - Struct array defining wall geometry.
%   transmitter_position - 1x2 vector [x, y] for the transmitter's position.
%   receiver_position    - 1x2 vector [x, y] for the receiver's position.
%   all_rays_data        - Cell array of structs with ray information.
%   K                    - Maximum number of reflections, used for the plot title.

    figure('Name', 'V2V Ray Tracing Results', 'NumberTitle', 'off');                      % Create a new figure window with a specific name.
    hold on;                                                                              % Ensure all subsequent plot commands draw on the same axes.
    grid on;                                                                              % Add a grid for better readability of coordinates.
    ax = gca;                                                                             % Get the current axes handle for direct manipulation.

    % --- Plot Walls ---
    for i = 1:length(walls)
        plot(ax, walls(i).coordinates(:,1), walls(i).coordinates(:,2), 'k', 'LineWidth', 2.5, 'DisplayName', 'wall'); % Draw walls as thick black lines for high contrast.
    end

    % --- Plot Transmitter and Receiver ---
    plot(ax, transmitter_position(1), transmitter_position(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', '#0077BE', 'DisplayName', 'TX'); % Mark the transmitter with a distinct blue-filled circle.
    plot(ax, receiver_position(1), receiver_position(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', '#D95319', 'DisplayName', 'RX'); % Mark the receiver with a distinct orange-filled circle.

    % --- Plot All Found Rays ---
    colors = lines(K);                                                                    % Generate a set of K distinct colors to represent each reflection order.
    
    for i = 1:length(all_rays_data)
        current_ray = all_rays_data{i};
        
        if strcmp(current_ray.type, 'LOS')
            plot(ax, current_ray.coordinates(:,1), current_ray.coordinates(:,2), 'g-', 'LineWidth', 1.5, 'DisplayName', 'LOS'); % Plot the Line-of-Sight path in green for special emphasis.
        else
            reflection_count = sscanf(current_ray.type, '%d-Refl');                       % Extract the integer number of reflections from the ray's type string.
            if ~isempty(reflection_count) && reflection_count <= size(colors, 1)
                plot_color = colors(reflection_count,:);                                  % Assign the pre-defined color for this reflection order.
            else
                plot_color = 'm';                                                         % Use a default color (magenta) for any unexpected case.
            end
            plot(ax, current_ray.coordinates(:,1), current_ray.coordinates(:,2), '-', 'Color', plot_color, 'LineWidth', 1, 'DisplayName', current_ray.type); % Draw the multi-segment reflected path.
        end
    end
    
    hold off;                                                                             % Release the plot hold, a good practice for figure manipulation.
    
    % --- Add Labels and a Clean Legend ---
    title(ax, sprintf('Ray Tracing with up to %d Reflections', K), 'FontSize', 16);                       % Set a descriptive title for the plot.
    xlabel(ax, 'x axis(m)', 'FontSize', 12);
    ylabel(ax, 'y axis(m)', 'FontSize', 12);
    
    % This section creates a clean legend by showing only one entry for each unique
    all_plot_handles = findobj(gcf,'Type','line');                                        % Find all line objects that have been plotted in the figure.
    if ~isempty(all_plot_handles)
        [~, unique_indices] = unique(get(all_plot_handles, 'DisplayName'), 'stable');     % Identify the first occurrence of each unique display name.
        legend(all_plot_handles(unique_indices));                                         % Generate the legend using only the handles for these unique items.
    end
end