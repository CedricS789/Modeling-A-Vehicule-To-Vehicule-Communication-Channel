function ray_tracing(w, d, k, tx_pos, RX_pos)
% Inputs:
%    w      - Half-width of the street
%    d      - Distance between TX and RX
%    tx_pos - Coordinates of the Transmitter as [x, y]
%    RX_pos - Coordinates of the Receiver as [x, y]
%    k      - Maximum number of reflections to trace


figure('Name', 'Recursive Image Method', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 800]);
hold on;
grid on;
ax = gca;

plot(ax, [0, d], [w, w], 'k-', 'LineWidth', 1.5, 'DisplayName', 'Building Wall');
plot(ax, [0, d], [-w, -w], 'k-', 'LineWidth', 1.5, 'HandleVisibility','off');
plot(ax, tx_pos(1), tx_pos(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', '#0077BE', 'MarkerEdgeColor', '#0077BE', 'DisplayName', 'Transmitter (TX)');
plot(ax, RX_pos(1), RX_pos(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', '#D95319', 'MarkerEdgeColor', '#D95319', 'DisplayName', 'Receiver (RX)');
plot(ax, [tx_pos(1), RX_pos(1)], [tx_pos(2), RX_pos(2)], 'g-', 'LineWidth', 1, 'DisplayName', 'LOS ray');

trace_recursive(ax, tx_pos, RX_pos, w, 1, k);

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

function trace_recursive(ax, tx_pos, RX_pos, w, n, k)
    if n > k
        return;
    end
    colors = lines(k);
    
    for ray_type = 1:2
        
        images = cell(n, 1);
        current_pos = tx_pos;
        for i = 1:n
            is_top_reflection = (ray_type == 1 && mod(i, 2) == 1) || (ray_type == 2 && mod(i, 2) == 0);
            if is_top_reflection
                current_pos = [current_pos(1), 2*w - current_pos(2)];
            else
                current_pos = [current_pos(1), -2*w - current_pos(2)];
            end
            images{i} = current_pos;
        end
        
        ray_points = zeros(n + 2, 2);
        ray_points(1, :) = tx_pos;
        ray_points(end, :) = RX_pos;
        next_point_in_ray = RX_pos;
        is_ray_valid = true;
        
        for i = n:-1:1
            image_source = images{i};
            is_top_reflection = (ray_type == 1 && mod(i, 2) == 1) || (ray_type == 2 && mod(i, 2) == 0);
            
            if is_top_reflection
                refl_point = intersect_line_y(image_source, next_point_in_ray, w);
            else
                refl_point = intersect_line_y(image_source, next_point_in_ray, -w);
            end
            
            if refl_point(1) <= tx_pos(1) || refl_point(1) >= RX_pos(1)
                is_ray_valid = false;
                break;
            end
            
            ray_points(i + 1, :) = refl_point;
            next_point_in_ray = refl_point;
        end
        
        if is_ray_valid
            h_plots = findobj(ax, 'Type', 'line');
            legend_str = sprintf('%d-Reflection', n);
            if ~any(strcmpi(get(h_plots, 'DisplayName'), legend_str))
                plot(ax, ray_points(:,1), ray_points(:,2), '-', 'Color', colors(n,:), 'LineWidth', 1, 'DisplayName', legend_str);
            else
                plot(ax, ray_points(:,1), ray_points(:,2), '-', 'Color', colors(n,:), 'LineWidth', 1, 'HandleVisibility','off');
            end
        end
    end
    
    trace_recursive(ax, tx_pos, RX_pos, w, n + 1, k);
end

function p_intersect = intersect_line_y(p1, p2, y_line)
    if p1(1) == p2(1)
        x_intersect = p1(1);
    else
        m = (p2(2) - p1(2)) / (p2(1) - p1(1));
        x_intersect = p1(1) + (y_line - p1(2)) / m;
    end
    p_intersect = [x_intersect, y_line];
end
