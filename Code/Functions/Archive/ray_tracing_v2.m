function [alphas, rays_data] = ray_tracing_v2(walls, k_max, tx_pos, RX_pos, params)
 rays_data = {};
    is_ray_clear = true;
    for i = 1:length(walls)
        intersection = findSegmentIntersection(tx_pos, RX_pos, walls(i).coordinates(1,:), walls(i).coordinates(2,:));
        if ~isempty(intersection)
            is_ray_clear = false; 
            break;
        end
    end
    
    if is_ray_clear
        los_ray.coordinates = [tx_pos; RX_pos];
        los_ray.type = 'LOS';
        los_ray.distance_total = norm(RX_pos - tx_pos);
        los_ray.gamma_prod = 1;
        los_ray.alpha_n = calculate_alpha(los_ray, params);
        rays_data{end+1} = los_ray;
    end

    for order = 1:k_max
        reflected_rays = findReflectionrays(tx_pos, RX_pos, walls, order, [], {tx_pos}, params);
        rays_data = [rays_data, reflected_rays];
    end
    
    num_rays = length(rays_data);
    alphas = zeros(1, num_rays, 'like', 1i);

    for i = 1:num_rays
        current_ray_struct = rays_data{i};
        alphas(i) = current_ray_struct.alpha_n;
    end

    plot_environment(walls, tx_pos, RX_pos, rays_data, k_max);

end



function found_rays = findReflectionrays(current_tx_pos, RX_pos, walls, reflections_remaining, wall_index_sequence, image_source_sequence, params)
    found_rays = {};

    if reflections_remaining == 0
        ray_data = tracePhysicalray(RX_pos, walls, wall_index_sequence, image_source_sequence, params);
        
        if ~isempty(ray_data)
            found_rays = {ray_data};
        end
        return;
    end

    last_wall_index = 0;
    if ~isempty(wall_index_sequence), last_wall_index = wall_index_sequence(end); end

    for i = 1:length(walls)
        if i == last_wall_index, continue; end 
        
        new_image_source = reflectPointAcrossWall(current_tx_pos, walls(i).coordinates);
        sub_rays = findReflectionrays(new_image_source, RX_pos, walls, reflections_remaining - 1, [wall_index_sequence, i], [image_source_sequence, {new_image_source}], params);
        
        found_rays = [found_rays, sub_rays];
    end
end

function ray_data = tracePhysicalray(RX_pos, walls, wall_index_sequence, image_source_sequence, params)
    ray_data = [];
    K = length(image_source_sequence) - 1;
    if K < 1, return; end

    ray_points = zeros(K + 2, 2);
    ray_points(1,:) = image_source_sequence{1};
    ray_points(end,:) = RX_pos;
    
    is_ray_valid = true;
    current_target = RX_pos;

    for i = K:-1:1
        image_for_current_segment = image_source_sequence{i+1};
        reflecting_wall_coordinates = walls(wall_index_sequence(i)).coordinates;

        reflection_point = findSegmentIntersection(image_for_current_segment, current_target, reflecting_wall_coordinates(1,:), reflecting_wall_coordinates(2,:));
        
        if isempty(reflection_point)
            is_ray_valid = false; break;
        end
        ray_points(i+1,:) = reflection_point;

        for j = 1:length(walls)
            is_start_wall = (j == wall_index_sequence(i));
            is_end_wall = (i < K && j == wall_index_sequence(i+1));

            if ~is_start_wall && ~is_end_wall
                if ~isempty(findSegmentIntersection(reflection_point, current_target, walls(j).coordinates(1,:), walls(j).coordinates(2,:)))
                    is_ray_valid = false; break;
                end
            end
        end
        if ~is_ray_valid, break; end
        
        current_target = reflection_point;
    end

    if is_ray_valid
        for j = 1:length(walls)
            if j ~= wall_index_sequence(1) && ~isempty(findSegmentIntersection(ray_points(1,:), ray_points(2,:), walls(j).coordinates(1,:), walls(j).coordinates(2,:)))
                is_ray_valid = false; break;
            end
        end
    end

    if is_ray_valid
        total_dist = 0;
        cumulative_gamma = 1;

        for i = 1:(K + 1)
            p1 = ray_points(i,:);
            p2 = ray_points(i+1,:);
            total_dist = total_dist + norm(p2 - p1);

            if i <= K
                incident_vec = p2 - p1;
                wall = walls(wall_index_sequence(i));
                wall_vec = wall.coordinates(2,:) - wall.coordinates(1,:);
                normal_vec = [wall_vec(2), -wall_vec(1)];
                
                cos_theta_i = abs(dot(incident_vec/norm(incident_vec), normal_vec/norm(normal_vec)));
                sin_theta_i_sq = 1 - cos_theta_i^2;
                eps_r = wall.eps_r;
                
                gamma_num = cos_theta_i - sqrt(eps_r - sin_theta_i_sq);
                gamma_den = cos_theta_i + sqrt(eps_r - sin_theta_i_sq);
                gamma_perp = gamma_num / gamma_den;
                
                cumulative_gamma = cumulative_gamma * gamma_perp;
            end
        end
        
        ray_data.coordinates = ray_points;
        ray_data.type = sprintf('%d-Refl', K);
        ray_data.distance_total = total_dist;
        ray_data.gamma_prod = cumulative_gamma;
        ray_data.alpha_n = calculate_alpha(ray_data, params);
    end
end

function alpha_n = calculate_alpha(ray_data, params)
    d_n = ray_data.distance_total;
    gamma_prod = ray_data.gamma_prod;
    fc = params.fc;
    c = params.c;
    Z_0 = params.Z_0;
    R_a = params.R_a;
    lambda = c / fc;
    
    tau_n = d_n / c;
    phase_shift = exp(-1j * 2 * pi * fc * tau_n);
    amplitude = (lambda * Z_0) / (4 * pi^2 * R_a * d_n);
    
    alpha_n = 1j * amplitude * phase_shift * gamma_prod;
end

function reflected_point = reflectPointAcrossWall(point, wall_coordinates)
    wall_p1 = wall_coordinates(1,:);
    wall_vector = wall_coordinates(2,:) - wall_p1;
    normal_vector = [wall_vector(2), -wall_vector(1)];
    reflected_point = point - 2 * dot(point - wall_p1, normal_vector) / dot(normal_vector, normal_vector) * normal_vector;
end

function intersection_point = findSegmentIntersection(seg1_p1, seg1_p2, seg2_p1, seg2_p2)
    v1 = seg1_p2 - seg1_p1;
    v2 = seg2_p2 - seg2_p1;
    cross_product = v1(1)*v2(2) - v1(2)*v2(1);
    if abs(cross_product) < 1e-10, intersection_point = []; return; end

    start_points_vector = seg2_p1 - seg1_p1;
    t = (start_points_vector(1)*v2(2) - start_points_vector(2)*v2(1)) / cross_product;
    u = (start_points_vector(1)*v1(2) - start_points_vector(2)*v1(1)) / cross_product;

    if t >= -1e-9 && t <= 1+1e-9 && u >= -1e-9 && u <= 1+1e-9
        intersection_point = seg1_p1 + t * v1;
    else
        intersection_point = [];
    end
end

function plot_environment(walls, tx_pos, RX_pos, rays_data, k_max)
    figure('Name', 'V2V Ray Tracing', 'NumberTitle', 'off');
    hold on; grid on; axis equal;
    ax = gca;

    for i = 1:length(walls)
        plot(ax, walls(i).coordinates(:,1), walls(i).coordinates(:,2), 'k', 'LineWidth', 2.5);
    end
    plot(ax, tx_pos(1), tx_pos(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', '#0077BE', 'DisplayName', 'TX');
    plot(ax, RX_pos(1), RX_pos(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', '#D95319', 'DisplayName', 'RX');

    colors = lines(k_max + 1); 
    for i = 1:length(rays_data)
        ray = rays_data{i};
        if strcmp(ray.type, 'LOS')
            plot(ax, ray.coordinates(:,1), ray.coordinates(:,2), 'g-', 'LineWidth', 1.5, 'DisplayName', 'LOS');
        else
            num_refl = sscanf(ray.type, '%d-Refl');
            plot(ax, ray.coordinates(:,1), ray.coordinates(:,2), '-', 'Color', colors(num_refl,:), 'LineWidth', 1, 'DisplayName', ray.type);
        end
    end
    
    hold off;
    title(ax, sprintf('Ray Tracing with up to %d Reflections', k_max));
    xlabel(ax, 'x axis(m)');
    ylabel(ax, 'y axis(m)');
    
    h = findobj(gcf,'Type','line');
    if ~isempty(h)
        [~, unique_indices] = unique(get(h, 'DisplayName'), 'stable');
        legend(h(unique_indices));
    end
end
