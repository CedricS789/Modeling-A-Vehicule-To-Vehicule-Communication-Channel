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
params.PRX_sens_dBm = -1000;
params.PTX_dBm = 10 * log10(params.PTX * 1000);
params.Gain = params.Z0 / (pi * params.Ra);
params.G_dBi = 10 * log10(params.Gain);
params.lambda = params.c / params.fc;

PTX = params.PTX;
sens_dBm = params.PRX_sens_dBm;
lambda = params.lambda;

M = 10;                      % Maximum number of reflections to consider
w = 20;
L = 10;                % Length of wall in meters
eps_r = 4;                  % Relative permittivity building walls

TX_pos = [0, 0];
walls(1).coordinates = [[0, w/2];  [L, w/2]];  walls(1).eps_r = eps_r;
walls(2).coordinates = [[0, -w/2]; [L, -w/2]];  walls(2).eps_r = eps_r;

d_fixed = 100;
RX_pos = [d_fixed, 0];


%% NARROWBAND LOS ANALYSIS

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



%% Full Channel - NARROWBAND PRX vs Distance and K-factor
% Simulate over a range of distances
fprintf('   - Simulating over distances\n');
distance_domain = logspace(0, log10(L), 50000); % points from 1m to 10^
PRX_LOS_dBm_domain = zeros(1, length(distance_domain));
PRX_total_dBm_domain = zeros(1, length(distance_domain));
K_factor_dB_domain = zeros(1, length(distance_domain));

queue_p_v_d = parallel.pool.DataQueue;
waitbar_p_v_d = waitbar(0, 'P_{RX} vs. Distance');
set(waitbar_p_v_d, 'UserData', 0);
afterEach(queue_p_v_d, @(~) updateWaitbar(waitbar_p_v_d, length(distance_domain), 'P_{RX} vs. Distance'));

tic;
parfor i = 1:length(distance_domain)
    d_i = distance_domain(i);
    current_tx_pos = [0, 0];
    current_RX_pos = [d_i, 0];
    [all_alphas, all_rays] = runRayTracing(walls, M, current_tx_pos, current_RX_pos, params);
    if isempty(all_alphas)
        PRX_total_dBm_domain(i) = -Inf;
        PRX_LOS_dBm_domain(i) = -Inf;
        K_factor_dB_domain(i) = -Inf;
        send(queue_p_v_d, 1);
        continue;
    end

    % Calculate total power
    h_nb = sum(all_alphas);
    PRX_total = PTX * abs(h_nb)^2;
    PRX_total_dBm_domain(i) = 10 * log10(PRX_total * 1000);
    
    % Use logical indexing to separate LOS and NLOS components
    isLOS = cellfun(@(r) strcmp(r.type, 'LOS'), all_rays);
 
    alpha_LOS = all_alphas(isLOS);

    if isempty(alpha_LOS)
        alpha_LOS = 0; 
    end % Ensure alpha_LOS is not empty if no LOS path exists
    
    P_LOS = PTX * sum(abs(alpha_LOS).^2); 
    P_NLOS = sum(PTX * abs(all_alphas(~isLOS)).^2);
    

    PRX_LOS_dBm_domain(i) = 10 * log10(P_LOS * 1000);
    if P_NLOS > 0
        K_factor = P_LOS / P_NLOS;
        K_factor_dB_domain(i) = 10 * log10(K_factor);
    else
        K_factor_dB_domain(i) = Inf; % Pure LOS case
    end
    send(queue_p_v_d, 1);
end
fprintf('      * PRX vs Distance simulation took %.2f seconds.\n', toc);

% Cleanup waitbar
if exist('waitbar_p_v_d', 'var') && ishandle(waitbar_p_v_d)
    close(waitbar_p_v_d);
end

fprintf('      * Data computation complete.\n');

% Plotting Results
fprintf('   - Plotting results\n');
plotPRXvsDistance(distance_domain, PRX_total_dBm_domain, PRX_LOS_dBm_domain);
plotKFactor(distance_domain, K_factor_dB_domain);



%% NARROWBAND PRX vs PTX
fprintf('\nAnalyzing Received Power vs. Transmitted Power\n');
d_fixed = 1000;

% Set a fixed distance between Tx and Rx
TX_pos = [0, 0];
RX_pos = [d_fixed, 0];
fprintf('   - Simulating for a fixed distance d = %.1fm\n', d_fixed);

% Calculate the channel transfer function just once for this fixed geometry.
[all_alphas, all_rays] = runRayTracing(walls, M, TX_pos, RX_pos, params);

% Calculate the total channel transfer function from all paths
h_nb_total = sum(all_alphas);

% Find the LOS-only channel transfer function using logical indexing
isLOS = cellfun(@(r) strcmp(r.type, 'LOS'), all_rays);
alpha_LOS = all_alphas(isLOS);
if isempty(alpha_LOS)
    alpha_LOS = 0; % Handle case where no LOS path is found
end

PTX_dBm_domain = -70:1:10; 
PRX_total_dBm_domain = zeros(1, length(PTX_dBm_domain));
PRX_LOS_dBm_domain = zeros(1, length(PTX_dBm_domain));

% Progress bar setup
queue_pr_v_pt = parallel.pool.DataQueue;
waitbar_pr_v_pt = waitbar(0, 'P_{RX} vs. P_{TX}');
set(waitbar_pr_v_pt, 'UserData', 0);
afterEach(queue_pr_v_pt, @(~) updateWaitbar(waitbar_pr_v_pt, length(PTX_dBm_domain), 'P_{RX} vs. P_{TX}'));

tic;
% Loop over the transmit power domain and calculate Rx power for each
parfor i = 1:length(PTX_dBm_domain)
    current_PTX_dbm = PTX_dBm_domain(i);
    current_PTX = 10^((current_PTX_dbm - 30) / 10);

    % Calculate total received power
    PRX_total = current_PTX * abs(h_nb_total)^2;
    PRX_total_dBm_domain(i) = 10 * log10(PRX_total * 1000);

    % Calculate LOS-only received power
    PRX_los = current_PTX * abs(alpha_LOS)^2;
    PRX_LOS_dBm_domain(i) = 10 * log10(PRX_los * 1000);
    send(queue_pr_v_pt, 1);
end
fprintf('      * PRX vs PTX simulation took %.2f seconds.\n', toc);

% Cleanup waitbar
if exist('waitbar_pr_v_pt', 'var') && ishandle(waitbar_pr_v_pt)
    close(waitbar_pr_v_pt);
end

% Plot the results
fprintf('   - Plotting PRX vs. PTX\n');
plotPRXvsPTX(PTX_dBm_domain, PRX_total_dBm_domain, PRX_LOS_dBm_domain);



%% PATH LOSS MODEL
fprintf('\nFitting Path Loss Model\n');

% Define simulation parameters for path loss analysis
d_samp = 0.01;   % Sampling interval
d_local = 5.0;   % Must be greater than d_samp
d0 = 1;          % Reference distance for the model

% Define the distance range for the final path loss model
x_start_model = TX_pos(1) + 2 * d_local;
x_end_model   = L;

% Create a padded domain for simulation to ensure correct convolution at edges
x_start_padded = x_start_model - d_local;
x_end_padded   = x_end_model + d_local;

% We use a linear domain for the averaging. It cannot work with a log domain
distance_domain_padded = x_start_padded:d_samp:x_end_padded;

PRX = zeros(1, length(distance_domain_padded));

fprintf('   - Calculating instantaneous PRX for %d points\n', length(distance_domain_padded));

% Progress bar setup
queue_pl = parallel.pool.DataQueue;
waitbar_pl = waitbar(0, 'Path Loss - Calculating PRX');
set(waitbar_pl, 'UserData', 0);
afterEach(queue_pl, @(~) updateWaitbar(waitbar_pl, length(distance_domain_padded), 'Path Loss - Calculating PRX'));
tic;

% Calculate instantaneous power over the entire dense simulation domain
parfor i = 1:length(distance_domain_padded)
    d_i = distance_domain_padded(i);
    RX_pos_local = [d_i, 0];
    
    % Run ray tracing for this sample point
    [all_alphas_local, all_rays] = runRayTracing(walls, M, [0,0], RX_pos_local, params);
    
    isLOS = cellfun(@(r) strcmp(r.type, 'LOS'), all_rays);
 
    alpha_LOS = all_alphas_local(isLOS);

    if ~isempty(alpha_LOS)
        PRX_Friis(i) = PTX * sum(abs(alpha_LOS).^2);
    else
        PRX_Friis(i) = 0; 
    end
     

    if ~isempty(all_alphas_local)
        h_nb_local = sum(all_alphas_local);
        PRX(i) = PTX * abs(h_nb_local)^2;
    else
        PRX(i) = 0;
    end
    send(queue_pl, 1);
end

fprintf('      * Instantaneous PRX calculation took %.2f seconds.\n', toc);
if exist('waitbar_pl', 'var') && ishandle(waitbar_pl)
    close(waitbar_pl);
end

% Spatially average the power using 1D convolution
fprintf('   - Averaging instantaneous power using 1D convolution\n');
tic;
num_samples_local = round(d_local / d_samp);

if mod(num_samples_local, 2) == 0; num_samples_local = num_samples_local + 1; end % Ensure the kernel has a center
kernel = ones(1, num_samples_local) / num_samples_local; % Normalized 1D filter kernel

PRX_avg = conv(PRX, kernel, 'same');
fprintf('      * Convolution took %.2f seconds.\n', toc);

% Select only the valid, non-padded data for the model
valid_indices = find(distance_domain_padded >= x_start_model & distance_domain_padded <= x_end_model);
distance_domain_model = distance_domain_padded(valid_indices);

PRX_avg_model = PRX_avg(valid_indices);
PRX_Friis_model = PRX_Friis(valid_indices);

PRX_dBm = 10 * log10(PRX * 1000);
PRX_Friis_model_dBm = 10 * log10(PRX_Friis_model * 1000);
PRX_avg_model_dBm = 10 * log10(PRX_avg_model * 1000);

% Call the dedicated function with the CORRECT data
[n, L0_d0, sigma_L] = pathLoss(distance_domain_padded, PRX_dBm, distance_domain_model, PRX_Friis_model_dBm, PRX_avg_model_dBm, params, d0);



%% NARROWBAND HEATMAP
fprintf('\nGenerating 2D coverage heatmap and local averages\n');
TX_pos     = [10, 0];

x_start_model    = -20;
x_end_model      = 400;

y_start    = -w/2-30;
y_end      = w/2+30;

d_samp     = 0.05;
d_local    = 5.0;  % Cannot be smaller than d_samp

fprintf('   - Generating %.1fm x %.1fm grid with %.1fm sampling interval\n', ...
        (x_end_model - x_start_model), (y_end - y_start), d_samp);
        
num_x_points = round((x_end_model - x_start_model) / d_samp) + 1;
num_y_points = round((y_end - y_start) / d_samp) + 1;
RX_x_coordinates = linspace(x_start_model, x_end_model, num_x_points);
RX_y_coordinates = linspace(y_start, y_end, num_y_points);

% Run Simulation for each point on the grid
PRX_dBm = zeros(num_y_points, num_x_points);
fprintf('      * Generating Heatmap and Averaged Power Data \n');

% Progress bar setup
queue_hm = parallel.pool.DataQueue;
num_iterations = num_x_points;
waitbar_hm = waitbar(0, '2D Heatmap');
set(waitbar_hm, 'UserData', 0);
afterEach(queue_hm, @(~) updateWaitbar(waitbar_hm, num_iterations, '2D Heatmap'));

tic;
% First, calculate the instantaneous power P_RX for the entire grid
parfor i = 1:num_x_points
    % Create temporary column variables to avoid slicing errors in parfor
    temp_PRX_dBm_col = zeros(num_y_points, 1);
    RX_y_coordinates_temp = RX_y_coordinates;

    for j = 1:num_y_points
        RX_pos = [RX_x_coordinates(i), RX_y_coordinates_temp(j)];
        d_fixed = norm(RX_pos - TX_pos);
        
        % Avoid calculating at the transmitter's exact location
        if d_fixed < d_samp
            temp_PRX_dBm_col(j) = NaN;
            continue;
        end
        
        % Run ray tracing for this sample point
        [all_alphas, ~] = runRayTracing(walls, M, TX_pos, RX_pos, params);
        
        % calculate power if rays are found
        if ~isempty(all_alphas)
            h_nb = sum(all_alphas);
            PRX = PTX * abs(h_nb)^2;
            temp_PRX_dBm_col(j) = 10 * log10(PRX * 1000);
        else
            temp_PRX_dBm_col(j) = -Inf;
        end
    end
    
    % Assign the entire columns at once
    PRX_dBm(:, i) = temp_PRX_dBm_col;

    % Send a message to the queue to update the progress bar
    if ~isempty(queue_hm)
        send(queue_hm, 1);
    end
end

% Second, calculate the average power <PRX> by filtering the instantaneous results
PRX = 10.^((PRX_dBm - 30) / 10);
PRX(isnan(PRX)) = 0; % Handle NaN for filtering

num_samples_local = round(d_local / d_samp);
if mod(num_samples_local, 2) == 0; num_samples_local = num_samples_local + 1; end % Ensure the kernel has a center

kernel = ones(num_samples_local) / (num_samples_local^2); % Normalized filter kernel
PRX_avg_model = conv2(PRX, kernel, 'same');

PRX_avg_dBm = 10 * log10(PRX_avg_model * 1000);
PRX_avg_dBm(PRX_avg_model == 0) = -Inf;

fprintf('      * Data generation took %.2f seconds.\n', toc);

% Cleanup waitbar
if exist('waitbar_hm', 'var') && ishandle(waitbar_hm)
    close(waitbar_hm);
end

% Plotting Results
fprintf('   - Plotting heatmaps\n');
% Plot Instantaneous Power
figure('Name', 'Instantaneous Power Heatmap', 'NumberTitle', 'off');
plotHeatmap(gca, RX_x_coordinates, RX_y_coordinates, PRX_dBm, TX_pos, walls, sens_dBm);
title('Instantaneous Power $P_{RX}$', 'Interpreter', 'latex');

% Plot Averaged Power
% Trick to make the heatmap display well
walls(1).coordinates = walls(1).coordinates + d_local/2 * [[1  1]; [1  1]];
walls(2).coordinates = walls(2).coordinates + d_local/2 * [[1 -1]; [1 -1]];

figure('Name', 'Averaged Power Heatmap', 'NumberTitle', 'off');
plotHeatmap(gca, RX_x_coordinates, RX_y_coordinates, PRX_avg_dBm, TX_pos, walls, sens_dBm);
title('Averaged Power $\langle P_{RX} \rangle$', 'Interpreter', 'latex');



%% Function to update the waitbar
function updateWaitbar(h_bar, total_iterations, title_str)
    current_count = get(h_bar, 'UserData') + 1;
    set(h_bar, 'UserData', current_count);
    fraction_done = current_count / total_iterations;
    waitbar(fraction_done, h_bar, sprintf('%s  %.1f%%', title_str, fraction_done * 100));
end
