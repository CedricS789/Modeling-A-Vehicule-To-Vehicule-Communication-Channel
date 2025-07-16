%% V2V NARROWBAND CHANNEL SIMULATION - Parrallelism application for speed
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
params.Z_0 = 377;
params.R_a = 73.1;
params.PTX = 0.1; 
params.PRX_sens_dBm = -1000;
params.PTX_dBm = 10 * log10(params.PTX * 1000);
params.lambda = params.c / params.fc;

PTX = params.PTX;
sens_dBm = params.PRX_sens_dBm;
lambda = params.lambda;

M = 1;                     % Maximum number of reflections to consider
w = 20;
L = 10000e3;                % Length of wall in meters
eps_r = 4;                  % Relative permittivity building walls

walls(1).coordinates = [[0, w/2];  [L, w/2]];  walls(1).eps_r = eps_r;
walls(2).coordinates = [[0, -w/2]; [L, -w/2]];  walls(2).eps_r = eps_r;



%% NARROWBAND LOS ANALYSIS
d = 120;
TX_pos = [0, 0];
RX_pos = [d, 0];

% Run ray tracing for LOS path : M = 0
[~, LOS_rays_data] = runRayTracing(walls, 0, TX_pos, RX_pos, params);
fprintf('\nPerforming LOS analysis for d = %.1fm \n', d);

if ~isempty(LOS_rays_data)
    LOS_ray = LOS_rays_data{1};
    LOS_delay = LOS_ray.distance_total / params.c;
    h_nb_LOS = LOS_ray.alpha_n;
    
    % Calculate received power using the channel transfer function
    PRX_LOS = PTX * abs(h_nb_LOS)^2;
    PRX_LOS_dBm = 10 * log10(PRX_LOS * 1000);
    
    fprintf('   - For d = %.1fm:\n', d);
    fprintf('      * tau_LOS = %.3es\n', LOS_delay);
    fprintf('      * Narrowband Gain |h_NB| = %.3e, Angle = %.2f°\n', abs(h_nb_LOS), rad2deg(angle(h_nb_LOS)));
    fprintf('      * Received Power (LOS): %.2fdBm (should match Friis formula)\n', PRX_LOS_dBm);
else
    fprintf('   - No LOS path found for d = %.1fm.\n', d);
end



%% NARROWBAND PRX vs Distance and K-factor
fprintf('\nPerforming full multipath channel analysis\n');

% Analyze MPCs for one distance
fprintf('   - Analyzing MPCs for d = %.1fm\n', d);
[all_alphas, all_rays] = runRayTracing(walls, M, TX_pos, RX_pos, params);

% Plot the ray paths for the distance
fprintf('   - Plotting ray-tracing visualization\n');
plotRays(walls, TX_pos, RX_pos, all_rays, M);

% Display properties of each found ray
for i = 1:length(all_rays)
    ray = all_rays{i};
    fprintf('      * Ray %2d: Type = %-7s          d_%d = %7.2fm          |alpha_%d| = %.4e         arg(alpha_%2d) = %7.2f° \n', ... 
        i, ray.type, i, ray.distance_total, i, abs(ray.alpha_n), i, rad2deg(angle(ray.alpha_n)));
end

% Calculate total received power from all paths
h_nb_total = sum(all_alphas);
PRX_total = PTX * abs(h_nb_total)^2;
PRX_total_dBm = 10 * log10(PRX_total * 1000);
fprintf('\n   - Total Narrowband Gain at %.1fm: |h_NB| = %.3e\n', d, abs(h_nb_total));
fprintf('   - Total Received Power at %.1fm: PRX = %.2fdBm\n', d, PRX_total_dBm);

% Simulate over a range of distances
fprintf('   - Simulating over distances\n');
distances_domain = logspace(0, log10(L), 50000); % points from 1m to 10^
PRX_LOS_dBm_domain = zeros(1, length(distances_domain));
PRX_total_dBm_domain = zeros(1, length(distances_domain));
K_factor_dB_domain = zeros(1, length(distances_domain));

queue_1 = parallel.pool.DataQueue;
waitbar_1 = waitbar(0, 'P_{RX} vs. Distance');
set(waitbar_1, 'UserData', 0);
afterEach(queue_1, @(~) updateWaitbar(waitbar_1, length(distances_domain), 'P_{RX} vs. Distance'));

tic;
parfor i = 1:length(distances_domain)
    d_loop = distances_domain(i);
    current_tx_pos = [0, 0];
    current_RX_pos = [d_loop, 0];
    [all_alphas, all_rays] = runRayTracing(walls, M, current_tx_pos, current_RX_pos, params);
    if isempty(all_alphas)
        PRX_total_dBm_domain(i) = -Inf;
        PRX_LOS_dBm_domain(i) = -Inf;
        K_factor_dB_domain(i) = -Inf;
        send(queue_1, 1);
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
    send(queue_1, 1);
end
fprintf('      * PRX vs Distance simulation took %.2f seconds.\n', toc);

% Cleanup waitbar
if exist('waitbar_1', 'var') && ishandle(waitbar_1)
    close(waitbar_1);
end

fprintf('      * Data computation complete.\n');

% Plotting Results
fprintf('   - Plotting results\n');
plotPRXvsDistance(distances_domain, PRX_total_dBm_domain, PRX_LOS_dBm_domain);
plotKFactor(distances_domain, K_factor_dB_domain);



%% PATH LOSS MODEL
fprintf('\nFitting Path Loss Model\n');

% Define simulation parameters for path loss analysis
path_loss_distances = logspace(log10(50), log10(1000), 50); % 50 points from 50m to 1km
PRX_avg_W = zeros(1, length(path_loss_distances));
d_local = 5.0;      % 5m local area for averaging
d_samp = 0.5;       % Sampling interval for averaging
num_samples_local = round(d_local / d_samp);

fprintf('   - Calculating <PRX> for %d distance points...\n', length(path_loss_distances));
% Progress bar setup
queue_pl = parallel.pool.DataQueue;
waitbar_pl = waitbar(0, 'Path Loss Model');
set(waitbar_pl, 'UserData', 0);
afterEach(queue_pl, @(~) updateWaitbar(waitbar_pl, length(path_loss_distances), 'Path Loss Model'));
tic;

parfor i = 1:length(path_loss_distances)
    d_center = path_loss_distances(i);
    
    % Define the sample points for the local area along the x-axis
    local_x_coords = linspace(d_center - d_local/2, d_center + d_local/2, num_samples_local);
    
    PRX_local_samples_W = zeros(1, num_samples_local);
    
    for k = 1:num_samples_local
        RX_pos_local = [local_x_coords(k), 0]; % Assume movement along x-axis
        
        % Run ray tracing for this sample point
        [alphas_local, ~] = runRayTracing(walls, M, [0,0], RX_pos_local, params);
        
        if ~isempty(alphas_local)
            h_nb_local = sum(alphas_local);
            PRX_local_samples_W(k) = params.PTX * abs(h_nb_local)^2;
        else
            PRX_local_samples_W(k) = 0;
        end
    end
    
    % Calculate the average power <PRX> in Watts
    PRX_avg_W(i) = mean(PRX_local_samples_W);
    send(queue_pl, 1);
end
fprintf('      * <PRX> calculation took %.2f seconds.\n', toc);
if exist('waitbar_pl', 'var') && ishandle(waitbar_pl)
    close(waitbar_pl);
end

% --- Path Loss Model Fitting ---
fprintf('   - Fitting model to data...\n');

% 1. Convert <PRX> data to L0 data points
PRX_avg_dBm = 10 * log10(PRX_avg_W * 1000);
L0_data_dB = params.PTX_dBm + 2 * params.G_dBi - PRX_avg_dBm;

% 2. Perform linear regression on L0 vs log10(d)
p_coeffs = polyfit(log10(path_loss_distances), L0_data_dB, 1);
slope = p_coeffs(1);
intercept = p_coeffs(2);

% 3. Find the path loss exponent 'n'
n = slope / 10;
fprintf('      * Path Loss Exponent (n) = %.2f\n', n);

% 4. Find the reference path loss L0(d0) at d0 = 1m
d0 = 1;
L0_d0 = intercept + slope * log10(d0); % Since d0=1, log10(d0)=0, L0_d0 = intercept
fprintf('      * Reference Path Loss L0(d0=1m) = %.2f dB\n', L0_d0);

% 5. Plot the results
fprintf('   - Plotting path loss model fit\n');
figure('Name', 'Path Loss Model Fit', 'NumberTitle', 'off');
plot(log10(path_loss_distances), L0_data_dB, 'bo', 'DisplayName', 'Simulated L_0 Data');
hold on;
fitted_L0 = polyval(p_coeffs, log10(path_loss_distances));
plot(log10(path_loss_distances), fitted_L0, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Fitted Model (n=%.2f)', n));
grid on;
xlabel('log_{10}(Distance [m])');
ylabel('Basic Path Loss L_0 [dB]');
title('Path Loss Model Fit for V2V Urban Canyon');
legend('show');
hold off;



%% NARROWBAND PRX vs PTX
fprintf('\nAnalyzing Received Power vs. Transmitted Power\n');
d = 1000;
% Set a fixed distance between Tx and Rx
TX_pos = [0, 0];
RX_pos = [d, 0];
fprintf('   - Simulating for a fixed distance d = %.1fm\n', d);

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
queue_3 = parallel.pool.DataQueue;
waitbar_3 = waitbar(0, 'P_{RX} vs. P_{TX}');
set(waitbar_3, 'UserData', 0);
afterEach(queue_3, @(~) updateWaitbar(waitbar_3, length(PTX_dBm_domain), 'P_{RX} vs. P_{TX}'));

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
    send(queue_3, 1);
end
fprintf('      * PRX vs PTX simulation took %.2f seconds.\n', toc);

% Cleanup waitbar
if exist('waitbar_3', 'var') && ishandle(waitbar_3)
    close(waitbar_3);
end

% Plot the results
fprintf('   - Plotting PRX vs. PTX\n');
plotPRXvsPTX(PTX_dBm_domain, PRX_total_dBm_domain, PRX_LOS_dBm_domain);



%% NARROWBAND HEATMAP
fprintf('\nGenerating 2D coverage heatmap and local averages\n');
TX_pos     = [10, 0];
x_start    = -20;
x_end      = 500;
y_start    = -w/2-30;
y_end      = w/2+30;
d_samp     = 0.05;
d_local    = 5.0;  % Cannot be smaller than d_samp
fprintf('   - Generating %.1fm x %.1fm grid with %.1fm sampling interval\n', ...
        (x_end - x_start), (y_end - y_start), d_samp);
        
num_x_points = round((x_end - x_start) / d_samp) + 1;
num_y_points = round((y_end - y_start) / d_samp) + 1;
RX_x_coordinates = linspace(x_start, x_end, num_x_points);
RX_y_coordinates = linspace(y_start, y_end, num_y_points);

% Run Simulation for each point on the grid
PRX_dBm = zeros(num_y_points, num_x_points);
fprintf('      * Generating Heatmap and Averaged Power Data \n');

% Progress bar setup
queue_2 = parallel.pool.DataQueue;
num_iterations = num_x_points;
waitbar_2 = waitbar(0, '2D Heatmap');
set(waitbar_2, 'UserData', 0);
afterEach(queue_2, @(~) updateWaitbar(waitbar_2, num_iterations, '2D Heatmap'));

tic;
% First, calculate the instantaneous power P_RX for the entire grid
parfor i = 1:num_x_points
    % Create temporary column variables to avoid slicing errors in parfor
    temp_PRX_dBm_col = zeros(num_y_points, 1);
    RX_y_coordinates_temp = RX_y_coordinates;

    for j = 1:num_y_points
        RX_pos = [RX_x_coordinates(i), RX_y_coordinates_temp(j)];
        dist = norm(RX_pos - TX_pos);
        
        % Avoid calculating at the transmitter's exact location
        if dist < d_samp
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
    if ~isempty(queue_2)
        send(queue_2, 1);
    end
end

% Second, calculate the average power <PRX> by filtering the instantaneous results
PRX = 10.^((PRX_dBm - 30) / 10);
PRX(isnan(PRX)) = 0; % Handle NaN for filtering

num_samples_local = round(d_local / d_samp);
if mod(num_samples_local, 2) == 0; num_samples_local = num_samples_local + 1; end % Ensure the kernel has a center

kernel = ones(num_samples_local) / (num_samples_local^2); % Normalized filter kernel
PRX_avg = conv2(PRX, kernel, 'same');

PRX_avg_dBm = 10 * log10(PRX_avg * 1000);
PRX_avg_dBm(PRX_avg == 0) = -Inf;

fprintf('      * Data generation took %.2f seconds.\n', toc);

% Cleanup waitbar
if exist('waitbar_2', 'var') && ishandle(waitbar_2)
    close(waitbar_2);
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
