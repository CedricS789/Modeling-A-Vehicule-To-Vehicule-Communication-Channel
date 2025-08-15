function [all_alphas, all_rays_data] = runRayTracing(walls, M, tx_pos, RX_pos, params)
% Coordinates the ray-tracing process.
    
    all_rays_data = {};

    is_los_obstructed = false;
    for i = 1:length(walls)
        intersection_point = findSegmentIntersection(tx_pos, RX_pos, walls(i).coordinates(1,:), walls(i).coordinates(2,:));
        if ~isempty(intersection_point)
            is_los_obstructed = true;
            break;
        end
    end

    if ~is_los_obstructed
        los_ray.coordinates = [tx_pos; RX_pos];
        los_ray.type = 'LOS';
        los_ray.distance_total = norm(RX_pos - tx_pos);
        los_ray.tau_n =  los_ray.distance_total / params.c;
        los_ray.theta_n = 0;
        los_ray.gamma_tot_n = 1;
        los_ray.alpha_n = calculateAlpha_n(los_ray, params);
        all_rays_data{end+1} = los_ray;
    end

    for order = 1:M
        initial_source_pos = tx_pos;
        initial_image_sequence = {initial_source_pos};
        
        reflected_rays = findReflectedRaysRecursive(initial_source_pos, RX_pos, walls, ...
            order, [], initial_image_sequence, params);
        
        all_rays_data = [all_rays_data, reflected_rays];
    end

    if ~isempty(all_rays_data)
        all_alphas = cellfun(@(ray) ray.alpha_n, all_rays_data);
    else
        all_alphas = [];
    end
end