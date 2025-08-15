function plotRays(walls, tx_pos, rx_pos, all_rays_data, M)
    figure('Name', 'Ray Tracing Results', 'NumberTitle', 'off');
    hold on;
    grid on;
    ax = gca;

    for i = 1:length(walls)
        plot(ax, walls(i).coordinates(:,1), walls(i).coordinates(:,2), 'k', 'LineWidth', 5, 'DisplayName', 'Wall');
    end

    plot(ax, tx_pos(1), tx_pos(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', '$T_x$');
    plot(ax, rx_pos(1), rx_pos(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', '#0077BE', 'DisplayName', '$R_x$');
    
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
    
    title(ax, sprintf('Ray Tracing with up to $M = %d$ Reflections', M), 'FontSize', 20, 'Interpreter', 'latex');
    xlabel(ax, '$x$-axis (m)', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel(ax, '$y$-axis (m)', 'FontSize', 18, 'Interpreter', 'latex');

    w = 20;
    xlim(ax, [tx_pos(1)-5, rx_pos(1)+20]);
    ylim(ax, [-w/2-5, w/2+5]);

    all_plot_handles = findobj(gcf,'Type','line');
    if ~isempty(all_plot_handles)
        [~, unique_indices] = unique(get(all_plot_handles, 'DisplayName'), 'stable');
        legend(all_plot_handles(unique_indices), 'Interpreter', 'latex', 'FontSize', 18, 'Location', 'northeast');
    end
end