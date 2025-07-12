function ray_tracing(w, d, k, tx_pos, rx_pos)
% Inputs:
%    w      - Half-width of the street (m)
%    d      - Distance between TX and RX (m)
%    tx_pos - Coordinates of the Transmitter as [x, y]
%    rx_pos - Coordinates of the Receiver as [x, y]
%    k      - Maximum number of reflections to trace


%% ---- Create the Plot ----
figure('Name', 'Recursive Image Method', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 800]);
hold on;
grid on;
ax = gca;

% Plot environment
plot(ax, [0, d], [w, w], 'k-', 'LineWidth', 1.5, 'DisplayName', 'Building Wall');
plot(ax, [0, d], [-w, -w], 'k-', 'LineWidth', 1.5, 'HandleVisibility','off');
plot(ax, tx_pos(1), tx_pos(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', '#0077BE', 'MarkerEdgeColor', '#0077BE', 'DisplayName', 'Transmitter (TX)');
plot(ax, rx_pos(1), rx_pos(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', '#D95319', 'MarkerEdgeColor', '#D95319', 'DisplayName', 'Receiver (RX)');
plot(ax, [tx_pos(1), rx_pos(1)], [tx_pos(2), rx_pos(2)], 'g-', 'LineWidth', 1, 'DisplayName', 'LOS Path');

%% ---- Start Recursive Ray Tracing ----
trace_recursive(ax, tx_pos, rx_pos, w, 1, k);

%% ---- Final Plot Formatting ----
hold off;
axis equal;
ylim(ax, [-w*2, w*2]);
xlim(ax, [-d/10, d + d/10]);
xlabel(ax, 'Distance Along Street (m)');
ylabel(ax, 'Distance Across Street (m)');
title(ax, sprintf('V2V Ray Tracing (d=%.0fm, k=%d reflections)', d, k));
legend(ax, 'show', 'Location', 'best');
set(ax, 'Color', [0.98 0.98 0.98]);
end

%% ---- Functions ----
function trace_recursive(ax, tx_pos, rx_pos, w, n, k)
    % Base case for the recursion: stop when n exceeds the max reflections
    if n > k
        return;
    end
    % Generate distinct colors for each reflection order
    colors = lines(k);
    
    % ---- Trace and plot paths for the current reflection order 'n' ----
    for path_type = 1:2 % Two path families: starting top or starting bottom
        
        % Generate the sequence of n image sources for the current path
        images = cell(n, 1);
        current_pos = tx_pos;
        for i = 1:n
            % Determine if the current reflection is off the top or bottom wall
            is_top_reflection = (path_type == 1 && mod(i, 2) == 1) || (path_type == 2 && mod(i, 2) == 0);
            if is_top_reflection
                current_pos = [current_pos(1), 2*w - current_pos(2)];
            else
                current_pos = [current_pos(1), -2*w - current_pos(2)];
            end
            images{i} = current_pos;
        end
        
        % Determine all reflection points by working backwards from the receiver
        path_points = zeros(n + 2, 2);
        path_points(1, :) = tx_pos;
        path_points(end, :) = rx_pos;
        next_point_in_path = rx_pos;
        is_path_valid = true;
        
        for i = n:-1:1
            image_source = images{i};
            is_top_reflection = (path_type == 1 && mod(i, 2) == 1) || (path_type == 2 && mod(i, 2) == 0);
            
            if is_top_reflection
                refl_point = intersect_line_y(image_source, next_point_in_path, w);
            else
                refl_point = intersect_line_y(image_source, next_point_in_path, -w);
            end
            
            % Validate that the reflection occurs physically between TX and RX
            if refl_point(1) <= tx_pos(1) || refl_point(1) >= rx_pos(1)
                is_path_valid = false;
                break;
            end
            
            path_points(i + 1, :) = refl_point;
            next_point_in_path = refl_point;
        end
        
        % Plot the full path if it's valid
        if is_path_valid
            % Check if a legend entry for this reflection order already exists
            h_plots = findobj(ax, 'Type', 'line');
            legend_str = sprintf('%d-Reflection', n);
            if ~any(strcmpi(get(h_plots, 'DisplayName'), legend_str))
                plot(ax, path_points(:,1), path_points(:,2), '-', 'Color', colors(n,:), 'LineWidth', 1, 'DisplayName', legend_str);
            else
                plot(ax, path_points(:,1), path_points(:,2), '-', 'Color', colors(n,:), 'LineWidth', 1, 'HandleVisibility','off');
            end
        end
    end
    
    % ---- Recursive call for the next order of reflections ----
    trace_recursive(ax, tx_pos, rx_pos, w, n + 1, k);
end

function p_intersect = intersect_line_y(p1, p2, y_line)
    % Finds the intersection of the line defined by points p1 and p2
    % with the horizontal line y = y_line.
    if p1(1) == p2(1) % Handle vertical line case
        x_intersect = p1(1);
    else
        m = (p2(2) - p1(2)) / (p2(1) - p1(1));
        x_intersect = p1(1) + (y_line - p1(2)) / m;
    end
    p_intersect = [x_intersect, y_line];
end
