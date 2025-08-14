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
params.Ra = 73.1;
params.PTX = 0.1;
params.BRF = 100e6;
params.Ltaps = 90;
params.resolution = 1 / params.BRF;
params.PRX_sens_dBm = -70;
params.PTX_dBm = 10 * log10(params.PTX * 1000);
params.Gain = params.Z0 / (pi * params.Ra);
params.G_dBi = 10 * log10(params.Gain);
params.lambda = params.c / params.fc;

PTX = params.PTX;
sens_dBm = params.PRX_sens_dBm;
lambda = params.lambda;

M = 10;                      % Maximum number of reflections to consider
w = 20;
L = 1e6;                % Length of wall in meters
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

if ~isempty(all_rays)
    plotLOSChannel(LOS_ray, params);  % h(τ) = α_LOS δ(τ − τ_LOS)
end
