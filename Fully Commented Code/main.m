%% V2V NARROWBAND CHANNEL SIMULATION
clear; close all; clc;

addpath('Functions');
addpath('Functions/Plotting Functions');

%% PARAMATERES
fprintf('STEP 1: Defining simulation parameters...\n');

% --- Physical & Channel Parameters ---
params.fc = 5.9e9;              % Carrier frequency in Hz
params.c = 3e8;                 % Speed of light in m/s
params.Z_0 = 377;               % Impedance of free space in Ohms
params.R_a = 73.1;              % Radiation resistance for a half-wave dipole
params.PTX = 0.1;              % Transmit power in Watts
params.PRX_sens_dBm = -70;     % Receiver sensitivity in dBm

% Derived parameters
params.PTX_dBm = 10 * log10(params.PTX * 1000);
params.lambda = params.c / params.fc;

% --- Ray-Tracing Configuration ---
M = 1;                          % Maximum number of reflections to consider

% --- Environment Geometry ---
w = 20;
L = 1e10;                 % Length of wall in meters
eps_r = 4;                  % Relative permittivity building walls

walls(1).coordinates = [[0, w/2]; [L, w/2]];    walls(1).eps_r = eps_r;
walls(2).coordinates = [[0, -w/2]; [L, -w/2]];  walls(2).eps_r = eps_r;
fprintf('   ...Configuration complete.\n\n');


%% STEP 2 - LOS-ONLY ANALYSIS
fprintf('STEP 2: Performing LOS-only analysis for d = %.1f ...\n', d);
d = 150;
tx_pos = [0, 0];
RX_pos = [d, 0];

% Run ray tracing for LOS path only (M = 0)
[~, LOS_rays_data] = runRayTracing(walls, 0, tx_pos, RX_pos, params);

if ~isempty(LOS_rays_data)
    LOS_ray = LOS_rays_data{1};
    LOS_delay = LOS_ray.distance_total / params.c;
    h_nb_LOS = LOS_ray.alpha_n;
    
    % Calculate received power using the channel transfer function
    PRX_LOS = params.PTX * abs(h_nb_LOS)^2;
    PRX_LOS_dBm = 10 * log10(PRX_LOS * 1000);
    
    fprintf('   - For d = %.1f m:\n', d);
    fprintf('     - \tau_LOS = %.3e s\n', LOS_delay);
    fprintf('     - Narrowband Gain |h_NB|: %.3e, Angle: %.2f deg\n', abs(h_nb_LOS), rad2deg(angle(h_nb_LOS)));
    fprintf('     - Received Power (LOS): %.2f dBm (should match Friis formula)\n', PRX_LOS_dBm);
else
    fprintf('   - No LOS path found for d = %.1f m.\n', d);
end
fprintf('   ...Step 2 analysis complete.\n\n');


%% STEP 3 - FULL MULTIPATH CHANNEL ANALYSIS
fprintf('STEP 3: Performing full multipath channel analysis...\n');

% --- Analyze MPCs for one distance ---
fprintf('   - Analyzing MPCs for d = %.1f m...\n', d);
[all_alphas, all_rays] = runRayTracing(walls, M, tx_pos, RX_pos, params);

% Plot the ray paths for the distance
fprintf('   - Plotting ray-tracing visualization...\n');
plotRays(walls, tx_pos, RX_pos, all_rays, M);
xlim([tx_pos(1)-20, RX_pos(1)+20]);
ylim([-w, w]);

% Display properties of each found ray
for i = 1:length(all_rays)
    ray = all_rays{i};
    fprintf('Ray %2d: Type = %-7s          d_%d = %7.2f m          |alpha_%d| = %.4e         arg(alpha_%2d) = %7.2f° \n', i, ray.type, i, ray.distance_total, i, abs(ray.alpha_n), i, rad2deg(angle(ray.alpha_n)));
end

% --- Calculate total received power from all paths ---
h_nb_total = sum(all_alphas);
PRX_total = params.PTX * abs(h_nb_total)^2;
PRX_total_dBm = 10 * log10(PRX_total * 1000);
fprintf('\n   - Total Narrowband Gain h_NB at %.1f m: %.3e\n', d, h_nb_total);
fprintf('   - Total Received Power PRX at %.1f m: %.2f dBm\n\n', d, PRX_total_dBm);


% --- Simulate over a range of distances ---
fprintf('   - Simulating over distances ...\n');
distances_domain = logspace(0, log10(L), 5e10); % points from 1m to 10^...
num_distances = length(distances_domain);
PRX_LOS_dBm_vs_dist = zeros(1, num_distances);
PRX_total_dBm_vs_dist = zeros(1, num_distances);
K_factor_dB_vs_dist = zeros(1, num_distances);

sim_waitbar = waitbar(0, 'Running Simulation vs. Distance...');
for i = 1:num_distances
    d = distances_domain(i);
    current_tx_pos = [0, 0];
    current_RX_pos = [d, 0];
    
    [alphas, rays] = runRayTracing(walls, M, current_tx_pos, current_RX_pos, params);
    
    if isempty(alphas)
        PRX_total_dBm_vs_dist(i) = -Inf;
        PRX_LOS_dBm_vs_dist(i) = -Inf;
        K_factor_dB_vs_dist(i) = -Inf;
        continue;
    end
    
    % Calculate total power
    h_nb = sum(alphas);
    PRX_total = params.PTX * abs(h_nb)^2;
    PRX_total_dBm_vs_dist(i) = 10 * log10(PRX_total * 1000);
    
    % Separate LOS and NLOS power for K-factor calculation
    alpha_LOS = 0;
    P_NLOS = 0;
    for j = 1:length(rays)
        if strcmp(rays{j}.type, 'LOS') % If type is LOS
            alpha_LOS = rays{j}.alpha_n;
        else
            P_NLOS = P_NLOS + params.PTX * abs(rays{j}.alpha_n)^2;
        end
    end
    
    P_LOS = params.PTX * abs(alpha_LOS)^2;
    PRX_LOS_dBm_vs_dist(i) = 10 * log10(P_LOS * 1000);
    
    if P_NLOS > 0
        K_factor = P_LOS / P_NLOS;
        K_factor_dB_vs_dist(i) = 10 * log10(K_factor);
    else
        K_factor_dB_vs_dist(i) = Inf; % Pure LOS case
    end
    
    waitbar(i/num_distances, sim_waitbar);
end
close(sim_waitbar);
fprintf('   ...Data computation complete.\n\n');

% --- Plotting Results ---
fprintf('   - Plotting results...\n');
plotPRXvsDistance(distances_domain, PRX_total_dBm_vs_dist, PRX_LOS_dBm_vs_dist);
plotKFactor(distances_domain, K_factor_dB_vs_dist);
% [n_pl, sigma_L, PL_d0] = plotPathLoss(distances, PRX_total_dBm_vs_dist, params);
% fprintf('     - Calculated Path Loss Exponent n = %.2f\n', n_pl);
% fprintf('     - Shadowing Standard Deviation sigma_L = %.2f dB\n', sigma_L);
% plotCellRangeAnalysis(distances, n_pl, PL_d0, sigma_L, params);
% fprintf('   ...Step 3 analysis complete.\n\n');



%% HEATMAP
fprintf('STEP 4: Generating 2D coverage heatmap...\n');

% --- Heatmap Control Parameters ---
tx_pos     = [0-10, 0];   % Transmitter position [x, y] in meters
x_start    = 0-20;        % Starting x-coordinate of the heatmap
x_end      = 200;       % Ending x-coordinate of the heatmap
y_start    = -w/2-10;      % Starting y-coordinate
y_end      = w/2+10;       % Ending y-coordinate
resolution = 0.5;          % Meters per pixel.

% --- Calculate simulation grid from control parameters ---
fprintf('   - Generating %.1fm x %.1fm grid with %.1fm resolution...\n', ...
        (x_end - x_start), (y_end - y_start), resolution);
        
num_x_points = round((x_end - x_start) / resolution) + 1;
num_y_points = round((y_end - y_start) / resolution) + 1;

RX_x_coordinates = linspace(x_start, x_end, num_x_points);
RX_y_coordinates = linspace(y_start, y_end, num_y_points);

% --- Run Simulation for each point on the grid ---
PRX_inst_dBm = zeros(length(RX_y_coordinates), length(RX_x_coordinates));
PRX_local_dBm = zeros(length(RX_y_coordinates), length(RX_x_coordinates));
% PRX_global_dBm = zeros(length(RX_y_coordinates), length(RX_x_coordinates));

heatmap_waitbar = waitbar(0, 'Generating Heatmap Data...');
total_points = length(PRX_inst_dBm);
point_count = 0;
% d0 = 1; % Reference distance for path loss model is 1m

for i = 1:length(RX_x_coordinates)
    for j = 1:length(RX_y_coordinates)
        RX_pos = [RX_x_coordinates(i), RX_y_coordinates(j)];
        dist_from_tx = norm(RX_pos - tx_pos);
        
        % Avoid calculating at the transmitter's exact location
        if dist_from_tx < resolution
            PRX_inst_dBm(j, i) = NaN;
            PRX_local_dBm(j, i) = NaN;
            % PRX_global_dBm(j, i) = NaN;
        else
            % --- Global Average Power <<PRX>> ---
            % path_loss_dB = PL_d0 + 10 * n_pl * log10(dist_from_tx / d0);
            % PRX_global_dBm(j, i) = params.PTX_dBm - path_loss_dB;
            
            % --- Ray Tracing for Instantaneous and Local Average Power ---
            [alphas, ~] = runRayTracing(walls, M, tx_pos, RX_pos, params);
            if ~isempty(alphas)

                % --- Instantaneous Power PRX ---
                h_nb = sum(alphas);
                PRX_inst = params.PTX * abs(h_nb)^2;
                PRX_inst_dBm(j, i) = 10 * log10(PRX_inst * 1000);
                
                % --- Local Average Power <PRX> ---
                power_of_each_path = params.PTX * abs(alphas).^2;
                PRX_local_avgs = sum(power_of_each_path);                      %TODO: This is not what the teacher wants I think, will resolve later
                PRX_local_dBm(j, i) = 10 * log10(PRX_local_avgs * 1000);
            else
                PRX_inst_dBm(j, i) = NaN;
                PRX_local_dBm(j, i) = NaN;
            end
        end
        point_count = point_count + 1;
        waitbar(point_count/total_points, heatmap_waitbar);
    end
end
close(heatmap_waitbar);

% --- Plotting Results ---
fprintf('   - Plotting heatmaps...\n');
figure('Name', 'Power Heatmap Comparison', 'NumberTitle', 'off', 'Position', [50, 200, 1800, 400]);

% Plot 1: Instantaneous Power
ax1 = subplot(2, 1, 1);
plotCoverageHeatmap(ax1, RX_x_coordinates, RX_y_coordinates, PRX_inst_dBm, tx_pos, walls, params.PRX_sens_dBm);
title('Instantaneous Power (P_{RX})');

% Plot 2: Local Average Power
ax2 = subplot(2, 1, 2);
plotCoverageHeatmap(ax2, RX_x_coordinates, RX_y_coordinates, PRX_local_dBm, tx_pos, walls, params.PRX_sens_dBm);
title('Local Average Power (<P_{RX}>)');

% Plot 3: Global Average Power
% ax3 = subplot(3, 1, 3);
% plotCoverageHeatmap(ax3, RX_x_coordinates, RX_y_coordinates, PRX_global_dBm, tx_pos, walls, params.PRX_sens_dBm);
% title('Global Average Power (<<P_{RX}>>)');


fprintf('   ...Step 4 analysis complete.\n\n');