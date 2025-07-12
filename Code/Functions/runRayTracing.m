function [alphas, all_rays_data] = runRayTracing(walls, max_reflection_order, transmitter_position, receiver_position, simulation_parameters)
% runRayTracing - Main engine to find all propagation rays.
%
% This function serves as the high-level coordinator for the entire ray-tracing
% process. It systematically finds the Line-of-Sight (LOS) ray and all
% possible reflected rays up to the specified maximum reflection order.
%
% INPUTS:
%   walls                - Struct array defining wall geometry and properties.
%   max_reflection_order - Maximum number of reflections to consider.
%   transmitter_position - 1x2 vector [x, y] for the transmitter's position.
%   receiver_position    - 1x2 vector [x, y] for the receiver's position.
%   simulation_parameters- Struct with simulation parameters.
%
% OUTPUTS:
%   alphas         - 1xN complex vector of channel gains for the N found rays.
%   all_rays_data - 1xN cell array of structs, each with detailed ray info.

    all_rays_data = {};

    % --- 1. LINE-OF-SIGHT (LOS) ray CALCULATION ---
    is_los_ray_obstructed = false;
    for i = 1:length(walls)
        intersection_point = findSegmentIntersection(transmitter_position, receiver_position, walls(i).coords(1,:), walls(i).coords(2,:));
        if ~isempty(intersection_point)
            is_los_ray_obstructed = true;
            break;
        end
    end

    if ~is_los_ray_obstructed
        los_ray_data.coordinates = [transmitter_position; receiver_position];
        los_ray_data.type = 'LOS';
        los_ray_data.dist = norm(receiver_position - transmitter_position);
        los_ray_data.gamma_prod = 1;
        los_ray_data.gain = calculateAlpha_n(los_ray_data, simulation_parameters);
        all_rays_data{end+1} = los_ray_data;
    end

    % --- 2. REFLECTED rayS CALCULATION ---
    for order = 1:max_reflection_order
        reflected_rays = findReflectedRaysRecursive(transmitter_position, receiver_position, walls, order, [], {transmitter_position}, simulation_parameters);
        all_rays_data = [all_rays_data, reflected_rays];
    end

    % --- 3. FINAL OUTPUT PREPARATION ---
    num_rays_found = length(all_rays_data);
    alphas = zeros(1, num_rays_found, 'like', 1i);
    for i = 1:num_rays_found
        alphas(i) = all_rays_data{i}.alpha_n;
    end
end
