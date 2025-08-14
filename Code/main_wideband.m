%% V2V NARROWBAND CHANNEL SIMULATION - Parallelism Implementation
% Parallel Processing toolbox required
clear; close all; clc;
addpath('Functions');
addpath('Functions/Plotting Functions');



%% START A PARLLEL POOL
if isempty(gcp('nocreate'))
    parpool(); 
end



%% PARAMATERES
fprintf('Defining simulation parameters\n');
params.fc = 5.9e9;
params.c = 3e8;
params.Z0 = 377;
params.Ra = 73;
params.PTX = 0.1;
params.BRF = 100e6;
params.Ltaps = 80;
params.resolution = 1 / params.BRF;
params.PRX_sens_dBm = -70;
params.PTX_dBm = 10 * log10(params.PTX * 1000);
params.Gain = params.Z0 / (pi * params.Ra);
params.G_dBi = 10 * log10(params.Gain);
params.lambda = params.c / params.fc;

PTX = params.PTX;
sens_dBm = params.PRX_sens_dBm;
lambda = params.lambda;

M = 3;                      % Maximum number of reflections to consider
w = 20;
L = 1000;                % Length of wall in meters
eps_r = 4;                  % Relative permittivity building walls

TX_pos = [0, 0];
walls(1).coordinates = [[0, w/2];  [L, w/2]];  walls(1).eps_r = eps_r;
walls(2).coordinates = [[0, -w/2]; [L, -w/2]];  walls(2).eps_r = eps_r;

d_fixed = 100;
RX_pos = [d_fixed, 0];



%% WIDEBAND LOS ANALYSIS - SAME AS NARROWBAND
% Run ray tracing for LOS path : M = 0
[~, all_rays] = runRayTracing(walls, 0, TX_pos, RX_pos, params);
fprintf('\nPerforming LOS analysis for d = %.1fm \n', d_fixed);

if ~isempty(all_rays)
    LOS_ray = all_rays{1};
    LOS_delay = LOS_ray.tau_n;
    h_nb_LOS = LOS_ray.alpha_n;
    
    % Calculate received power using the channel transfer function
    PRX_LOS = PTX * abs(h_nb_LOS)^2;
    PRX_LOS_dBm = 10 * log10(PRX_LOS * 1000);
    
    fprintf('   - For d = %.1fm:\n', d_fixed);
    fprintf('      * tau_LOS = %.3es\n', LOS_delay);
    fprintf('      * Narrowband Gain |h_NB_LOS| = %.3e, arg(h_NB_LOS) = %.2f°\n', abs(h_nb_LOS), rad2deg(angle(h_nb_LOS)));
    fprintf('      * PRX_LOS: %.2fdBm\n', PRX_LOS_dBm);
else
    fprintf('   - No LOS path found for d = %.1fm.\n', d_fixed);
end


%% Full Channel - NARROWBAND
fprintf('\nPerforming full multipath channel analysis\n');

% Analyze MPCs for one distance
fprintf('   - Analyzing MPCs for d = %.1f m, M = %.0f reflections\n', d_fixed, M);
[all_alphas, all_rays] = runRayTracing(walls, M, TX_pos, RX_pos, params);

% Plot the ray paths for the distance
fprintf('   - Plotting ray-tracing visualization\n');
plotRays(walls, TX_pos, RX_pos, all_rays, M);
    
% Display properties of each found ray
for i = 1:length(all_rays)
    ray = all_rays{i};
    fprintf('      * Ray %2d: Type = %-7s          d_%.2d = %1.2f m          tau_%.2d = %5.2f ns          theta_%.2d = %5.2f°          gamma_tot_%.2d = %9.2e          |alpha_%.2d| = %8.4e         arg(alpha_%.2d) = %8.2f° \n', ... 
        i, ray.type, i, ray.distance_total, i, ray.tau_n*1e9, i, ray.theta_n, i, ray.gamma_tot_n, i, abs(ray.alpha_n), i, rad2deg(angle(ray.alpha_n)));
end

% Calculate total received power from all paths
h_nb_total = sum(all_alphas);
PRX_total = PTX * abs(h_nb_total)^2;
PRX_total_dBm = 10 * log10(PRX_total * 1000);
fprintf('\n   - Total Narrowband Gain at %.1fm: |h_NB| = %.3e\n', d_fixed, abs(h_nb_total));
fprintf('   - Total Received Power at %.1fm: PRX = %.2fdBm\n', d_fixed, PRX_total_dBm);


plotWidebandChannel(all_rays, params);  
