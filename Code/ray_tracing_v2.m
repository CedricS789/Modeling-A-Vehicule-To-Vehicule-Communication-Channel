function [alphas, paths_data] = ray_tracing_v2(walls, k_max, tx_pos, rx_pos, sim_params)
% RAY_TRACING_V2 - Performs ray tracing for a V2V communication scenario.
%
% This function identifies all propagation paths between a transmitter (TX)
% and a receiver (RX) up to a maximum number of reflections (k_max). It
% calculates the complex gain coefficient 'alpha' for each path based on
% the physical properties of the channel.
%
% SYNTAX:
%   [alphas, paths_data] = ray_tracing_v2(walls, k_max, tx_pos, rx_pos, sim_params)
%
% INPUTS:
%   walls        - Struct array defining the geometry and properties of walls.
%                  Each element should have:
%                  walls(i).coords: A 2x2 matrix [x1 y1; x2 y2] for wall endpoints.
%                  walls(i).eps_r: The relative permittivity of the wall material.
%   k_max        - Maximum number of reflections to consider.
%   tx_pos       - 1x2 vector [x, y] of the transmitter's position.
%   rx_pos       - 1x2 vector [x, y] of the receiver's position.
%   sim_params   - Struct containing simulation parameters:
%                  .fc: Carrier frequency (Hz)
%                  .c: Speed of light (m/s)
%                  .Z_0: Impedance of free space (Ohms)
%                  .R_a: Antenna radiation resistance (Ohms)
%
% OUTPUTS:
%   alphas       - 1xN complex vector of gain coefficients for N found paths.
%   paths_data   - 1xN cell array of structs, each containing detailed info
%                  about a path (gain, geometric points, type, distance).
%
% Based on the project "Modeling A Vehicle-to-Vehicle Communication Channel"
% ELEC-H415, ULB/VUB.

    % --- Setup Environment ---
    paths_data = {};
    
    % --- Line-of-Sight (LOS) Path Calculation ---
    is_path_clear = true;
    for i = 1:length(walls)
        if ~isempty(findSegmentIntersection(tx_pos, rx_pos, walls(i).coords(1,:), walls(i).coords(2,:)))
            is_path_clear = false;
            break;
        end
    end
    
    if is_path_clear
        path_data.path = [tx_pos; rx_pos];
        path_data.type = 'LOS';
        path_data.dist = norm(rx_pos - tx_pos);
        path_data.gamma_prod = 1; % No reflections for LOS
        
        % Calculate complex gain for LOS path
        path_data.gain = calculate_alpha(path_data, sim_params);
        paths_data{end+1} = path_data;
    end

    % --- Reflected Paths Calculation ---
    % Initiate recursive search for reflected paths for each reflection order
    for order = 1:k_max
        reflected_paths = findReflectionPaths(tx_pos, rx_pos, walls, order, [], {tx_pos}, sim_params);
        paths_data = [paths_data, reflected_paths];
    end

    % --- Final Output Preparation ---
    if ~isempty(paths_data)
        alphas = cellfun(@(c) c.gain, paths_data);
    else
        alphas = [];
    end

    % --- Plotting ---
    plot_environment(walls, tx_pos, rx_pos, paths_data, k_max);

end

%% --- Recursive Path Finding ---
function found_paths = findReflectionPaths(current_source, receiver, walls, reflections_remaining, wall_index_sequence, image_source_sequence, sim_params)
    % Recursively finds all possible image source chains and traces their physical paths.
    found_paths = {};
    if reflections_remaining == 0
        % Base case: A full image chain is formed, now trace the physical path.
        path_data = tracePhysicalPath(receiver, walls, wall_index_sequence, image_source_sequence, sim_params);
        if ~isempty(path_data)
            found_paths = {path_data};
        end
        return;
    end

    last_wall_index = 0;
    if ~isempty(wall_index_sequence), last_wall_index = wall_index_sequence(end); end

    for i = 1:length(walls)
        if i == last_wall_index, continue; end % Avoid reflecting off the same wall twice in a row.
        
        new_image_source = reflectPointAcrossWall(current_source, walls(i).coords);

        % Recursive call with one less reflection remaining
        sub_paths = findReflectionPaths(new_image_source, receiver, walls, reflections_remaining - 1, [wall_index_sequence, i], [image_source_sequence, {new_image_source}], sim_params);
        found_paths = [found_paths, sub_paths];
    end
end

%% --- Physical Path Tracing and Validation ---
function path_data = tracePhysicalPath(receiver, walls, wall_index_sequence, image_source_sequence, sim_params)
    % From a complete image chain, this function traces the path backwards from
    % the receiver to the transmitter, validates it, and calculates its properties.
    path_data = [];
    num_reflections = length(image_source_sequence) - 1;
    if num_reflections < 1, return; end

    path_points = zeros(num_reflections + 2, 2);
    path_points(1,:) = image_source_sequence{1}; % Original TX
    path_points(end,:) = receiver;               % RX
    
    is_path_valid = true;
    current_target = receiver;

    % Trace backwards from RX to find all reflection points
    for i = num_reflections:-1:1
        image_for_current_segment = image_source_sequence{i+1};
        reflecting_wall_coords = walls(wall_index_sequence(i)).coords;

        reflection_point = findSegmentIntersection(image_for_current_segment, current_target, reflecting_wall_coords(1,:), reflecting_wall_coords(2,:));
        
        if isempty(reflection_point) % VALIDATION 1: Must intersect the physical wall segment.
            is_path_valid = false; break;
        end
        path_points(i+1,:) = reflection_point;

        % VALIDATION 2: Check for obstructions by other walls.
        for j = 1:length(walls)
            is_start_wall = (j == wall_index_sequence(i));
            is_end_wall = (i < num_reflections && j == wall_index_sequence(i+1));
            if ~is_start_wall && ~is_end_wall
                if ~isempty(findSegmentIntersection(reflection_point, current_target, walls(j).coords(1,:), walls(j).coords(2,:)))
                    is_path_valid = false; break;
                end
            end
        end
        if ~is_path_valid, break; end
        current_target = reflection_point;
    end

    % Final validation for the first segment (TX to first reflection)
    if is_path_valid
        for j = 1:length(walls)
            if j ~= wall_index_sequence(1) && ~isempty(findSegmentIntersection(path_points(1,:), path_points(2,:), walls(j).coords(1,:), walls(j).coords(2,:)))
                is_path_valid = false; break;
            end
        end
    end

    % If the path is valid, calculate all its physical properties.
    if is_path_valid
        % --- Calculate Geometry and Reflection Coefficients ---
        total_dist = 0;
        total_gamma = 1;

        for i = 1:(num_reflections + 1)
            p1 = path_points(i,:);
            p2 = path_points(i+1,:);
            total_dist = total_dist + norm(p2 - p1);

            % If this segment ends in a reflection, calculate the coefficient
            if i <= num_reflections
                incident_vec = p2 - p1;
                wall = walls(wall_index_sequence(i));
                wall_vec = wall.coords(2,:) - wall.coords(1,:);
                normal_vec = [wall_vec(2), -wall_vec(1)];
                
                % Angle of incidence w.r.t. the normal
                cos_theta_i = abs(dot(incident_vec/norm(incident_vec), normal_vec/norm(normal_vec)));
                sin_theta_i_sq = 1 - cos_theta_i^2;
                eps_r = wall.eps_r;
                
                % Fresnel reflection coefficient for perpendicular polarization
                gamma_num = cos_theta_i - sqrt(eps_r - sin_theta_i_sq);
                gamma_den = cos_theta_i + sqrt(eps_r - sin_theta_i_sq);
                gamma_perp = gamma_num / gamma_den;
                total_gamma = total_gamma * gamma_perp;
            end
        end
        
        % --- Package results ---
        path_data.path = path_points;
        path_data.type = sprintf('%d-Refl', num_reflections);
        path_data.dist = total_dist;
        path_data.gamma_prod = total_gamma;
        path_data.gain = calculate_alpha(path_data, sim_params);
    end
end

%% --- Complex Gain Calculation ---
function alpha = calculate_alpha(path_data, sim_params)
    % Calculates the complex gain 'alpha' for a given path.
    % Unpack parameters
    d_n = path_data.dist;
    gamma_prod = path_data.gamma_prod;
    fc = sim_params.fc;
    c = sim_params.c;
    Z_0 = sim_params.Z_0;
    R_a = sim_params.R_a;
    lambda = c / fc;
    
    % Propagation delay
    tau_n = d_n / c;
    
    % Complex gain formula from the report
    % alpha_n = (j * (lambda * Z_0) / (4*pi^2 * R_a * d_n) * exp(-j*2*pi*fc*tau_n)) * (Gamma_prod)
    phase_shift = exp(-1j * 2 * pi * fc * tau_n);
    amplitude = (lambda * Z_0) / (4 * pi^2 * R_a * d_n);
    
    alpha = 1j * amplitude * phase_shift * gamma_prod;
end

%% --- Geometric Helper Functions ---
function reflected_point = reflectPointAcrossWall(point, wall_coords)
    wall_p1 = wall_coords(1,:);
    wall_vector = wall_coords(2,:) - wall_p1;
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

    if t >= -1e-10 && t <= 1+1e-10 && u >= -1e-10 && u <= 1+1e-10
        intersection_point = seg1_p1 + t * v1;
    else
        intersection_point = [];
    end
end

%% --- Plotting Function ---
function plot_environment(walls, tx_pos, rx_pos, paths_data, k_max)
    % --- Setup the Plotting Environment ---
    figure('Name', 'V2V Ray Tracing', 'NumberTitle', 'off');
    hold on; grid on; axis equal;
    ax = gca;

    % --- Plot Geometry ---
    for i = 1:length(walls)
        plot(ax, walls(i).coords(:,1), walls(i).coords(:,2), 'k', 'LineWidth', 2);
    end
    plot(ax, tx_pos(1), tx_pos(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', '#0077BE', 'DisplayName', 'TX');
    plot(ax, rx_pos(1), rx_pos(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', '#D95319', 'DisplayName', 'RX');

    % --- Plot Paths ---
    colors = lines(k_max + 1); % +1 for LOS color
    for i = 1:length(paths_data)
        path = paths_data{i};
        if strcmp(path.type, 'LOS')
            plot(ax, path.path(:,1), path.path(:,2), 'g-', 'LineWidth', 1.5, 'DisplayName', 'LOS');
        else
            num_refl = sscanf(path.type, '%d-Refl');
            plot(ax, path.path(:,1), path.path(:,2), '-', 'Color', colors(num_refl,:), 'LineWidth', 1, 'DisplayName', path.type);
        end
    end
    
    % --- Finalize Plot ---
    hold off;
    title(ax, sprintf('Ray Tracing with up to %d Reflections', k_max));
    xlabel(ax, 'X Coordinate (m)');
    ylabel(ax, 'Y Coordinate (m)');
    
    % Clean up legend
    h = findobj(gcf,'Type','line');
    if ~isempty(h)
        legend_handles = findobj(gca, 'Type', 'Line', '-and', '-not', {'DisplayName', ''});
        [~, unique_indices] = unique(get(legend_handles, 'DisplayName'), 'stable');
        if ~isempty(unique_indices)
            legend(legend_handles(unique_indices));
        end
    end
end
