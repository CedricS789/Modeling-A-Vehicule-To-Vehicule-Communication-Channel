%% V2V DOPPLER & TIME VARIANCE ANALYSIS
clear; close all; clc;

%% ========================================================================
%                       MAIN SCRIPT
% =========================================================================

%% PARAMETERS
fprintf('Defining simulation parameters\n');
params.fc = 5.9e9;
params.c = 3e8;
params.lambda = params.c / params.fc;
params.Z0 = 377;
params.Ra = 73.1;
params.PTX = 0.1;

M = 10;
w = 20;
L = 1000;
eps_r = 4;

d_fixed = 100;
v_kmh = 50;
v_ms = v_kmh * 1000 / 3600;
params.v = v_ms;

f_m_max = v_ms / params.lambda;
fs = 5 * f_m_max;
t_duration = 0.1;
t = 0:1/fs:t_duration;

TX_pos = [0, 0];
RX_pos_initial = [d_fixed, 0];
walls(1).coordinates = [[0, w/2];  [L, w/2]];  walls(1).eps_r = eps_r;
walls(2).coordinates = [[0, -w/2]; [L, -w/2]];  walls(2).eps_r = eps_r;

%% 1. PERFORM STATIC RAY TRACING
fprintf('Performing ray tracing for d = %.1f m and M = %d reflections\n', d_fixed, M);
[~, all_rays] = runRayTracing(walls, M, TX_pos, RX_pos_initial, params);

if isempty(all_rays)
    error('No rays found for the specified geometry. Cannot proceed with Doppler analysis.');
end

fprintf('   - Found %d valid MPCs.\n', length(all_rays));

%% 2. CALCULATE DOPPLER SHIFTS
fprintf('Calculating Doppler shifts for each MPC...\n');
doppler_results = calculateDoppler(all_rays, params);
f_m = doppler_results.f_m;
fprintf('   - Maximum Doppler Shift (fm): %.2f Hz\n', f_m);

%% 3. COMPUTE TIME-VARYING CHANNEL
fprintf('Computing time-variant channel gain h_NB(t)...\n');
h_nb_t = zeros(size(t));
for i = 1:length(t)
    time_phasors = exp(-1j * 2 * pi * doppler_results.f_D_n * t(i));
    h_nb_t(i) = sum(doppler_results.alphas_t0 .* time_phasors);
end

Tc = 1 / (2 * f_m);
fprintf('   - Calculated Coherence Time (Tc): %.4f ms\n', Tc * 1000);


%% 4. PLOT RESULTS
fprintf('Plotting Doppler analysis results...\n');
plotDopplerResults(t, h_nb_t, doppler_results, Tc);


%% ========================================================================
%                       LOCAL FUNCTIONS
% =========================================================================

function [all_alphas, all_rays_data] = runRayTracing(walls, M, tx_pos, RX_pos, params)
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

function valid_rays_found = findReflectedRaysRecursive(current_tx_pos, RX_pos, walls, num_reflections_remaining, processed_wall_indices, all_images_pos, params)
    valid_rays_found = {};

    if num_reflections_remaining == 0
        ray_coordinates = validateRayPath(RX_pos, walls, processed_wall_indices, all_images_pos); 
        if ~isempty(ray_coordinates)
            ray_data = calculatePhysicalProperties(ray_coordinates, processed_wall_indices, walls, params);
            valid_rays_found = {ray_data};
        end
        return;
    end

    last_wall_index = 0;
    if ~isempty(processed_wall_indices)
        last_wall_index = processed_wall_indices(end);
    end
    
    for i = 1:length(walls)
        if i == last_wall_index
            continue;
        end
        
        new_image_source = findSymmetricAcrossLine(current_tx_pos, walls(i).coordinates);
        
        rays_from_branch = findReflectedRaysRecursive(new_image_source, RX_pos, walls, ...
            num_reflections_remaining - 1, [processed_wall_indices, i], ...
            [all_images_pos, {new_image_source}], params);
        
        valid_rays_found = [valid_rays_found, rays_from_branch];
    end
end

function final_ray_coordinates = validateRayPath(RX_pos, walls, processed_wall_indices, all_images_pos)
    final_ray_coordinates = [];
    M = length(all_images_pos) - 1;
    if M < 1
        return;
    end
    potential_ray_coordinates = zeros(M + 2, 2);
    original_tx_pos = all_images_pos{1};
    potential_ray_coordinates(1,:) = original_tx_pos;
    potential_ray_coordinates(end,:) = RX_pos;
    is_path_valid = true;
    current_target_point = RX_pos;

    for k = M:-1:1
        image_source_for_segment = all_images_pos{k+1};
        reflecting_wall_index = processed_wall_indices(k);
        reflecting_wall_coordinates = walls(reflecting_wall_index).coordinates;

        reflection_point = findSegmentIntersection(...
            image_source_for_segment, current_target_point, ...
            reflecting_wall_coordinates(1,:), reflecting_wall_coordinates(2,:));
        if isempty(reflection_point)
            is_path_valid = false;
            break;
        end
        potential_ray_coordinates(k+1,:) = reflection_point;

        for j = 1:length(walls)
            is_this_segments_wall = (j == reflecting_wall_index);
            is_next_segments_wall = (k < M && j == processed_wall_indices(k+1));
            if ~is_this_segments_wall && ~is_next_segments_wall
                if ~isempty(findSegmentIntersection(reflection_point, current_target_point, ...
                        walls(j).coordinates(1,:), walls(j).coordinates(2,:)))
                    is_path_valid = false;
                    break;
                end
            end
        end
        if ~is_path_valid
            break;
        end
        current_target_point = reflection_point;
    end

    if is_path_valid
        first_reflection_wall_idx = processed_wall_indices(1);
        tx_pos = potential_ray_coordinates(1,:);
        first_reflection_point = potential_ray_coordinates(2,:);
        for j = 1:length(walls)
            if j ~= first_reflection_wall_idx
                if ~isempty(findSegmentIntersection(tx_pos, first_reflection_point, ...
                        walls(j).coordinates(1,:), walls(j).coordinates(2,:)))
                    is_path_valid = false;
                    break;
                end
            end
        end
    end

    if is_path_valid
        final_ray_coordinates = potential_ray_coordinates;
    end
end

function ray_data = calculatePhysicalProperties(ray_coordinates, hit_walls_indices, walls, params)
    c = params.c;
    M = size(ray_coordinates, 1) - 2;
    d_tot = 0;
    gamma_tot_n = 1;

    for m = 1:(M + 1)
        start_segment = ray_coordinates(m,:);
        end_segment = ray_coordinates(m+1,:);
        
        d_tot = d_tot + norm(end_segment - start_segment);
        
        is_reflection_segment = (m <= M);
        if is_reflection_segment
            incident_vector = end_segment - start_segment;
            reflecting_wall = walls(hit_walls_indices(m));
            
            wall_vector = reflecting_wall.coordinates(2,:) - reflecting_wall.coordinates(1,:);
            normal_vector = [wall_vector(2), -wall_vector(1)];
            
            cos_theta_n = abs(dot(incident_vector, normal_vector)) / (norm(incident_vector) * norm(normal_vector));
            sin_theta_n_sq = 1 - cos_theta_n^2;
            epsilon_r = reflecting_wall.eps_r;
            
            gamma_num = cos_theta_n - sqrt(epsilon_r - sin_theta_n_sq);
            gamma_den = cos_theta_n + sqrt(epsilon_r - sin_theta_n_sq);
            gamma_m = gamma_num / gamma_den;
            
            gamma_tot_n = gamma_tot_n * gamma_m;
        end
    end
    
    ray_data.coordinates = ray_coordinates;
    ray_data.type = sprintf('%d-Refl', M);
    ray_data.distance_total = d_tot;
    ray_data.tau_n = d_tot / c;
    ray_data.gamma_tot_n = gamma_tot_n;
    ray_data.alpha_n = calculateAlpha_n(ray_data, params);
end

function alpha_n = calculateAlpha_n(ray_data, params)
    d_n = ray_data.distance_total;
    gamma_tot_n = ray_data.gamma_tot_n;
    fc = params.fc;
    c = params.c;
    Z0 = params.Z0;
    Ra = params.Ra;
    lambda = params.lambda;
    
    tau_n = d_n / c;
    phase_shift = exp(-1j * 2 * pi * fc * tau_n);
    amplitude = (lambda * Z0) / (4 * pi^2 * Ra * d_n);
    alpha_n = 1j * amplitude * phase_shift * gamma_tot_n;
end

function doppler_results = calculateDoppler(all_rays, params)
    num_rays = length(all_rays);
    alphas_t0 = zeros(1, num_rays);
    theta_n_rad = zeros(1, num_rays);
    
    motion_vector = [1, 0]; 
    
    for i = 1:num_rays
        ray = all_rays{i};
        alphas_t0(i) = ray.alpha_n;
        
        arrival_vector = ray.coordinates(end,:) - ray.coordinates(end-1,:);
        cos_theta = dot(arrival_vector, motion_vector) / (norm(arrival_vector) * norm(motion_vector));
        cos_theta = max(min(cos_theta, 1), -1);
        
        theta_n_rad(i) = acos(cos_theta);
    end
    
    v = params.v;
    lambda = params.lambda;
    f_m = v / lambda;
    f_D_n = f_m * cos(theta_n_rad);
    
    doppler_results.alphas_t0 = alphas_t0;
    doppler_results.theta_n_rad = theta_n_rad;
    doppler_results.f_D_n = f_D_n;
    doppler_results.f_m = f_m;
end

function symmetric_point = findSymmetricAcrossLine(point, line_coordinates)
    line_point_1 = line_coordinates(1,:);
    line_direction_vector = line_coordinates(2,:) - line_point_1;
    normal_vector = [line_direction_vector(2), -line_direction_vector(1)];
    vector_to_point = point - line_point_1;
    scale = 2 * dot(vector_to_point, normal_vector) / dot(normal_vector, normal_vector);
    symmetric_point = point - scale * normal_vector;
end

function intersection_point = findSegmentIntersection(p1, p2, p3, p4)
    v1 = p2 - p1;
    v2 = p4 - p3;
    v1_cross_v2 = v1(1)*v2(2) - v1(2)*v2(1);

    if abs(v1_cross_v2) < 1e-10
        intersection_point = [];
        return;
    end

    p3_minus_p1 = p3 - p1;
    t = (p3_minus_p1(1)*v2(2) - p3_minus_p1(2)*v2(1)) / v1_cross_v2;
    u = (p3_minus_p1(1)*v1(2) - p3_minus_p1(2)*v1(1)) / v1_cross_v2;

    if (t >= -1e-9 && t <= 1+1e-9) && (u >= -1e-9 && u <= 1+1e-9)
        intersection_point = p1 + t * v1;
    else
        intersection_point = [];
    end
end

function plotDopplerResults(t, h_nb_t, doppler_results, Tc)
    figure('Name', 'Time Fading due to Doppler', 'NumberTitle', 'off');
    plot(t * 1000, 20*log10(abs(h_nb_t)), 'b-', 'LineWidth', 1.5);
    grid on; grid minor;
    title('Time Fading of Channel Magnitude $|h_{NB}(t)|$', 'Interpreter', 'latex', 'FontSize', 16);
    xlabel('Time (ms)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$|h_{NB}(t)|$ (dB)', 'Interpreter', 'latex', 'FontSize', 14);
    xlim([0, t(end)*1000]);
    legend(sprintf('Coherence Time $T_c \\approx %.2f$ ms', Tc*1000), 'Interpreter', 'latex', 'Location', 'southwest');
    
    figure('Name', 'Power Doppler Spectrum', 'NumberTitle', 'off');
    powers_dB = 10*log10(abs(doppler_results.alphas_t0).^2);
    stem(doppler_results.f_D_n, powers_dB, 'filled', 'LineWidth', 1.5, 'MarkerSize', 8);
    
    grid on; grid minor;
    title('Power Doppler Spectrum $P(f_D)$', 'Interpreter', 'latex', 'FontSize', 16);
    xlabel('Doppler Shift $f_D$ (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Power (dB)', 'Interpreter', 'latex', 'FontSize', 14);
    xlim([-doppler_results.f_m*1.1, doppler_results.f_m*1.1]);
    
    hold on;
    xline(doppler_results.f_m, '--r', '$f_m$', 'Interpreter', 'latex', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'FontSize', 14);
    xline(-doppler_results.f_m, '--r', '$-f_m$', 'Interpreter', 'latex', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'FontSize', 14);
    hold off;
end
