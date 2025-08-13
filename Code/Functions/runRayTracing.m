function [all_alphas, all_rays_data] = runRayTracing(walls, M, tx_pos, RX_pos, params)
% This function coordinates the ray-tracing process. It first checks for a
% Line-of-Sight path and then initiates a recursive search for all
% reflected paths up to the specified maximum reflection order.
%
% INPUTS:
%   walls                - Struct array defining wall geometry and properties.
%   M                    - Number of reflections to consider.
%   tx_pos               - [x, y] for the transmitter's position.
%   RX_pos               - [x, y] for the receiver's position.
%   params               - Struct with simulation parameters.
%
% OUTPUTS:
%   all_alphas           - 1xN complex vector of channel gains for the N found rays.
%   all_rays_data        - 1xN cell array of structs, each with detailed ray info.
    
    all_rays_data = {};

    % LOS PATH CALCULATION
    is_los_obstructed = false;
    for i = 1:length(walls)
        intersection_point = findSegmentIntersection(tx_pos, RX_pos, walls(i).coordinates(1,:), walls(i).coordinates(2,:));
        if ~isempty(intersection_point)
            is_los_obstructed = true;
            break;                                                              % If one wall obstructs, no need to check others.
        end
    end

    % If not obstructed, calculate the properties of the LOS ray.
    if ~is_los_obstructed
        los_ray.coordinates = [tx_pos; RX_pos];
        los_ray.type = 'LOS';
        los_ray.distance_total = norm(RX_pos - tx_pos);
        los_ray.tau_n =  los_ray.distance_total / params.c;
        los_ray.theta_n = 0;
        los_ray.gamma_tot_n = 1;                                                 % Reflection coefficient product is 1 for LOS.
        los_ray.alpha_n = calculateAlpha_n(los_ray, params);
        all_rays_data{end+1} = los_ray;
    end

    % Recursively find all valid reflected paths for each reflection order.
    for order = 1:M
        
        initial_source_pos = tx_pos;                                            % Start the recursion from the original transmitter position.
        initial_image_sequence = {initial_source_pos};
        
        reflected_rays = findReflectedRaysRecursive(initial_source_pos, RX_pos, walls, ...
            order, [], initial_image_sequence, params);
        
        % Append the valid rays found from this order to the main list.
        all_rays_data = [all_rays_data, reflected_rays];
    end

    % Extract the complex gains from the ray data structs.
    if ~isempty(all_rays_data)
        all_alphas = cellfun(@(ray) ray.alpha_n, all_rays_data);
    else
        all_alphas = [];
    end
end
