function [alphas, rays_data] = ray_tracing_v2(walls, k_max, tx_pos, rx_pos, params)
% This function implements the image method to find all propagation rays
% between a transmitter (TX) and receiver (RX) for a given 2D environment.
% It handles the Line-of-Sight (LOS) ray and reflected rays up to k_max reflections.
% For each valid ray, it calculates the complex channel gain coefficient 'alpha_n'
% based on ray length and Fresnel reflection coefficients.
%
% SYNTAX:
%   [alphas, rays_data] = ray_tracing_v2(walls, k_max, tx_pos, rx_pos, params)
%
% INPUTS:
%   walls        - Struct array defining wall geometry and properties.
%                  walls(i).coordinates: 2x2 matrix [x1 y1; x2 y2] for endpoints.
%                  walls(i).eps_r: Relative permittivity of the wall.
%   k_max        - Maximum number of reflections to consider (e.g., 3).
%   tx_pos       - 1x2 vector [x, y] for the transmitter's position.
%   rx_pos       - 1x2 vector [x, y] for the rx_pos's position.
%   params   - Struct with simulation parameters (.fc, .c, .Z_0, .R_a).
%
% OUTPUTS:
%   alphas       - 1xN complex vector of gain coefficients for N found rays.
%   rays_data   - 1xN cell array of structs, each with detailed ray info.
%

    % =================================================================
    % INITIALIZE
    % =================================================================
    % This cell array will store the detailed results for every valid
    % ray that is found.
    rays_data = {};
    
    % =================================================================
    % CALCULATE LINE-OF-SIGHT (LOS) ray
    % =================================================================
    % The LOS ray is the simplest case: a straight line from TX to RX.
    % We must first check if this ray is blocked by any walls.
    is_ray_clear = true;
    for i = 1:length(walls)
        % Check for intersection between the TX-RX line and each wall.
        intersection = findSegmentIntersection(tx_pos, rx_pos, walls(i).coordinates(1,:), walls(i).coordinates(2,:));
        if ~isempty(intersection)
            is_ray_clear = false; % ray is blocked.
            break;                 % No need to check other walls.
        end
    end
    
    % If the ray is not blocked, calculate its properties.
    if is_ray_clear
        los_ray.coordinates = [tx_pos; rx_pos]; % Geometric ray points
        los_ray.type = 'LOS';
        los_ray.total_distance = norm(rx_pos - tx_pos); % ray length
        los_ray.gamma_prod = 1; % Product of reflection coefficients is 1 (no reflections)
        
        % Calculate the complex gain 'alpha_n' for this ray.
        los_ray.alpha_n = calculate_alpha(los_ray, params);
        
        % Add the completed LOS ray data to our list of results.
        rays_data{end+1} = los_ray;
    end

    % =================================================================
    % FIND ALL REFLECTED rayS
    % =================================================================
    % We loop from 1 to the maximum number of reflections (k_max).
    % In each iteration, we find all rays with exactly 'order' reflections.
    for order = 1:k_max
        % This function starts the recursive search for rays with 'order' reflections.
        % It begins with the original transmitter position.
        reflected_rays = findReflectionrays(tx_pos, rx_pos, walls, order, [], {tx_pos}, params);
        
        % Add any valid rays found to our main list of results.
        rays_data = [rays_data, reflected_rays];
    end

    % =================================================================
    % PREPARE FINAL OUTPUTS AND PLOT
    % =================================================================
    % Extract just the complex gains from the 'rays_data' cell array 
    % into a simple vector for easy use. This uses a standard for-loop for clarity.
    
    num_rays = length(rays_data);
    alphas = zeros(1, num_rays, 'like', 1i); % Pre-allocate a complex array

    for i = 1:num_rays
        % For each ray found, get its struct from the cell array.
        current_ray_struct = rays_data{i};
        % Extract the 'gain' value and place it in the alphas vector.
        alphas(i) = current_ray_struct.alpha_n;
    end

    % Call the plotting function to visualize the environment and all found rays.
    plot_environment(walls, tx_pos, rx_pos, rays_data, k_max);

end

% =========================================================================
% HELPER FUNCTIONS 
% =========================================================================

%% --- Recursive ray Finding Engine ---
function found_rays = findReflectionrays(current_tx_pos, rx_pos, walls, reflections_remaining, wall_index_sequence, image_source_sequence, params)
    % This is the core of the image method. It recursively builds a tree of
    % all possible image source locations.
    found_rays = {};

    % BASE CASE: If no more reflections are needed, our image chain is complete.
    if reflections_remaining == 0
        % Now that we have a full chain of image sources, we can try to
        % trace the physical ray from the final image to the rx_pos.
        ray_data = tracePhysicalray(rx_pos, walls, wall_index_sequence, image_source_sequence, params);
        
        % If a valid physical ray was found, it is returned.
        if ~isempty(ray_data)
            found_rays = {ray_data};
        end
        return;
    end

    % RECURSIVE STEP: If we still need more reflections...
    last_wall_index = 0;
    if ~isempty(wall_index_sequence), last_wall_index = wall_index_sequence(end); end

    % Iterate through every wall in the environment to create the next image.
    for i = 1:length(walls)
        % A ray cannot reflect off the same wall twice in a row.
        if i == last_wall_index, continue; end 
        
        % Create the next image source by reflecting the current source across wall 'i'.
        new_image_source = reflectPointAcrossWall(current_tx_pos, walls(i).coordinates);

        % --- Make the recursive call ---
        % We now search for rays from this *new* image source, with one
        % less reflection required. We also update our lists of which walls
        % we've hit and the positions of the image sources.
        sub_rays = findReflectionrays(new_image_source, rx_pos, walls, reflections_remaining - 1, [wall_index_sequence, i], [image_source_sequence, {new_image_source}], params);
        
        % Add any rays found by the recursive call to our list.
        found_rays = [found_rays, sub_rays];
    end
end

%% --- Physical ray Tracer and Validator ---
function ray_data = tracePhysicalray(rx_pos, walls, wall_index_sequence, image_source_sequence, params)
    % Takes a complete image chain and determines if a valid, unobstructed
    % physical ray exists. If so, it calculates its properties.
    ray_data = [];
    K = length(image_source_sequence) - 1;
    if K < 1, return; end

    % This array will store the [TX, R1, R2, ..., RX] coordinates.
    ray_points = zeros(K + 2, 2);
    ray_points(1,:) = image_source_sequence{1}; % ray starts at the original TX.
    ray_points(end,:) = rx_pos;               % ray ends at the RX.
    
    is_ray_valid = true;
    current_target = rx_pos;

    % TRACE BACKWARDS: From RX to TX, find each reflection point.
    for i = K:-1:1
        image_for_current_segment = image_source_sequence{i+1};
        reflecting_wall_coordinates = walls(wall_index_sequence(i)).coordinates;

        % The reflection point is the intersection of that line and wall 'i'.
        reflection_point = findSegmentIntersection(image_for_current_segment, current_target, reflecting_wall_coordinates(1,:), reflecting_wall_coordinates(2,:));
        
        % VALIDATION 1: The intersection must be on the physical wall segment.
        if isempty(reflection_point)
            is_ray_valid = false; break;
        end
        ray_points(i+1,:) = reflection_point; % Store valid reflection point.

        % VALIDATION 2: The new ray segment must not be blocked by other walls.
        for j = 1:length(walls)
            % A segment is allowed to hit its start and end walls. We must
            % check for intersections with all OTHER walls.
            is_start_wall = (j == wall_index_sequence(i));
            
            % The end wall is the wall for the *next* reflection in the sequence.
            % This doesn't apply to the last segment (reflection -> rx_pos).
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

    % FINAL VALIDATION: Check the first segment (TX to first reflection point).
    if is_ray_valid
        for j = 1:length(walls)
            if j ~= wall_index_sequence(1) && ~isempty(findSegmentIntersection(ray_points(1,:), ray_points(2,:), walls(j).coordinates(1,:), walls(j).coordinates(2,:)))
                is_ray_valid = false; break;
            end
        end
    end

    % --- If ray is fully validated, calculate its physical properties ---
    if is_ray_valid
        total_dist = 0;
        cumulative_gamma = 1; % Start with a product of 1.

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
        ray_data.total_distance = total_dist;
        ray_data.gamma_prod = cumulative_gamma;
        ray_data.alpha_n = calculate_alpha(ray_data, params);
    end
end

%% --- Complex Gain (Alpha) Calculation ---
function alpha_n = calculate_alpha(ray_data, params)
    % Implements the formula for the complex gain from the project report.
    d_n = ray_data.total_distance;
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

%% --- Geometric Helper: Point Reflection ---
function reflected_point = reflectPointAcrossWall(point, wall_coordinates)
    wall_p1 = wall_coordinates(1,:);
    wall_vector = wall_coordinates(2,:) - wall_p1;
    normal_vector = [wall_vector(2), -wall_vector(1)];
    reflected_point = point - 2 * dot(point - wall_p1, normal_vector) / dot(normal_vector, normal_vector) * normal_vector;
end

%% --- Geometric Helper: Segment Intersection ---
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

%% --- Plotting Function ---
function plot_environment(walls, tx_pos, rx_pos, rays_data, k_max)
    figure('Name', 'V2V Ray Tracing', 'NumberTitle', 'off');
    hold on; grid on; axis equal;
    ax = gca;

    for i = 1:length(walls)
        plot(ax, walls(i).coordinates(:,1), walls(i).coordinates(:,2), 'k', 'LineWidth', 2.5);
    end
    plot(ax, tx_pos(1), tx_pos(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', '#0077BE', 'DisplayName', 'TX');
    plot(ax, rx_pos(1), rx_pos(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', '#D95319', 'DisplayName', 'RX');

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
    xlabel(ax, 'X Coordinate (m)');
    ylabel(ax, 'Y Coordinate (m)');
    
    h = findobj(gcf,'Type','line');
    if ~isempty(h)
        [~, unique_indices] = unique(get(h, 'DisplayName'), 'stable');
        legend(h(unique_indices));
    end
end
