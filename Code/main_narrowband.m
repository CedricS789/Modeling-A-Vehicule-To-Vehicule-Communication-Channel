%% V2V NARROWBAND CHANNEL SIMULATION - Parallelism Implementation
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
params.PRX_sens_dBm = -70;
params.PTX_dBm = 10 * log10(params.PTX * 1000);
params.Gain = params.Z0 / (pi * params.Ra);
params.G_dBi = 10 * log10(params.Gain);
params.lambda = params.c / params.fc;

PTX = params.PTX;
PRX_sens_dBm = params.PRX_sens_dBm;
lambda = params.lambda;

M = 10;
w = 20;
L = 1000;
eps_r = 4;

TX_pos = [0, 0];
walls(1).coordinates = [[0, w/2];  [L, w/2]];  walls(1).eps_r = eps_r;
walls(2).coordinates = [[0, -w/2]; [L, -w/2]];  walls(2).eps_r = eps_r;

d_fixed = 100;
RX_pos = [d_fixed, 0];


%% NARROWBAND LOS ANALYSIS
[~, all_rays] = runRayTracing(walls, 0, TX_pos, RX_pos, params);
fprintf('\nPerforming LOS analysis for d = %.1fm \n', d_fixed);

if ~isempty(all_rays)
    LOS_ray = all_rays{1};
    LOS_delay = LOS_ray.tau_n;
    h_nb_LOS = LOS_ray.alpha_n;
    
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

fprintf('   - Analyzing MPCs for d = %.1f m, M = %.0f reflections\n', d_fixed, M);
[all_alphas, all_rays] = runRayTracing(walls, M, TX_pos, RX_pos, params);

fprintf('   - Plotting ray-tracing visualization\n');
plotRays(walls, TX_pos, RX_pos, all_rays, M);

for i = 1:length(all_rays)
    ray = all_rays{i};
    fprintf('      * Ray %2d: Type = %-7s          d_%.2d = %1.2f m          tau_%.2d = %5.2f ns          theta_%.2d = %5.2f°          gamma_tot_%.2d = %9.2e          |alpha_%.2d| = %8.4e         arg(alpha_%.2d) = %8.2f° \n', ... 
        i, ray.type, i, ray.distance_total, i, ray.tau_n*1e9, i, ray.theta_n, i, ray.gamma_tot_n, i, abs(ray.alpha_n), i, rad2deg(angle(ray.alpha_n)));
end

h_nb_total = sum(all_alphas);
PRX_total = PTX * abs(h_nb_total)^2;
PRX_total_dBm = 10 * log10(PRX_total * 1000);
fprintf('\n   - Total Narrowband Gain at %.1fm: |h_NB| = %.3e\n', d_fixed, abs(h_nb_total));
fprintf('   - Total Received Power at %.1fm: PRX = %.2fdBm\n', d_fixed, PRX_total_dBm);



%% Full Channel - NARROWBAND PRX vs Distance and K-factor
fprintf('   - Simulating over distances\n');
distance_domain = logspace(0, log10(L), 50000);
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

    h_nb = sum(all_alphas);
    PRX_total = PTX * abs(h_nb)^2;
    PRX_total_dBm_domain(i) = 10 * log10(PRX_total * 1000);
    
    isLOS = cellfun(@(r) strcmp(r.type, 'LOS'), all_rays);
 
    alpha_LOS = all_alphas(isLOS);

    if isempty(alpha_LOS)
        alpha_LOS = 0; 
    end
    
    P_LOS = PTX * sum(abs(alpha_LOS).^2); 
    P_NLOS = sum(PTX * abs(all_alphas(~isLOS)).^2);
    

    PRX_LOS_dBm_domain(i) = 10 * log10(P_LOS * 1000);
    if P_NLOS > 0
        K_factor = P_LOS / P_NLOS;
        K_factor_dB_domain(i) = 10 * log10(K_factor);
    else
        K_factor_dB_domain(i) = Inf;
    end
    send(queue_p_v_d, 1);
end
fprintf('      * PRX vs Distance simulation took %.2f seconds.\n', toc);

if exist('waitbar_p_v_d', 'var') && ishandle(waitbar_p_v_d)
    close(waitbar_p_v_d);
end

fprintf('      * Data computation complete.\n');

fprintf('   - Plotting results\n');
plotPRXvsDistance(distance_domain, PRX_total_dBm_domain, PRX_LOS_dBm_domain);
plotKFactor(distance_domain, K_factor_dB_domain);



%% NARROWBAND PRX vs PTX
fprintf('\nAnalyzing Received Power vs. Transmitted Power\n');
d_fixed = 1000;

TX_pos = [0, 0];
RX_pos = [d_fixed, 0];
fprintf('   - Simulating for a fixed distance d = %.1fm\n', d_fixed);

[all_alphas, all_rays] = runRayTracing(walls, M, TX_pos, RX_pos, params);

h_nb_total = sum(all_alphas);

isLOS = cellfun(@(r) strcmp(r.type, 'LOS'), all_rays);
alpha_LOS = all_alphas(isLOS);
if isempty(alpha_LOS)
    alpha_LOS = 0;
end

PTX_dBm_domain = -70:1:10; 
PRX_total_dBm_domain = zeros(1, length(PTX_dBm_domain));
PRX_LOS_dBm_domain = zeros(1, length(PTX_dBm_domain));

queue_pr_v_pt = parallel.pool.DataQueue;
waitbar_pr_v_pt = waitbar(0, 'P_{RX} vs. P_{TX}');
set(waitbar_pr_v_pt, 'UserData', 0);
afterEach(queue_pr_v_pt, @(~) updateWaitbar(waitbar_pr_v_pt, length(PTX_dBm_domain), 'P_{RX} vs. P_{TX}'));

tic;
parfor i = 1:length(PTX_dBm_domain)
    current_PTX_dbm = PTX_dBm_domain(i);
    current_PTX = 10^((current_PTX_dbm - 30) / 10);

    PRX_total = current_PTX * abs(h_nb_total)^2;
    PRX_total_dBm_domain(i) = 10 * log10(PRX_total * 1000);

    PRX_los = current_PTX * abs(alpha_LOS)^2;
    PRX_LOS_dBm_domain(i) = 10 * log10(PRX_los * 1000);
    send(queue_pr_v_pt, 1);
end
fprintf('      * PRX vs PTX simulation took %.2f seconds.\n', toc);

if exist('waitbar_pr_v_pt', 'var') && ishandle(waitbar_pr_v_pt)
    close(waitbar_pr_v_pt);
end

fprintf('   - Plotting PRX vs. PTX\n');
plotPRXvsPTX(PTX_dBm_domain, PRX_total_dBm_domain, PRX_LOS_dBm_domain);



%% PATH LOSS MODEL
fprintf('\nFitting Path Loss Model\n');

d_samp = 0.01;
d_local = 5.0;
d0 = 100;

x_start_model = TX_pos(1) + 2 * d_local;
x_end_model   = L;

x_start_padded = x_start_model - d_local;
x_end_padded   = x_end_model + d_local;

distance_domain_padded = x_start_padded:d_samp:x_end_padded;

PRX = zeros(1, length(distance_domain_padded));

fprintf('   - Calculating instantaneous PRX for %d points\n', length(distance_domain_padded));

queue_pl = parallel.pool.DataQueue;
waitbar_pl = waitbar(0, 'Path Loss - Calculating PRX');
set(waitbar_pl, 'UserData', 0);
afterEach(queue_pl, @(~) updateWaitbar(waitbar_pl, length(distance_domain_padded), 'Path Loss - Calculating PRX'));
tic;

parfor i = 1:length(distance_domain_padded)
    d_i = distance_domain_padded(i);
    RX_pos_local = [d_i, 0];
    
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

fprintf('   - Averaging instantaneous power using 1D convolution\n');
tic;
num_samples_local = round(d_local / d_samp);

if mod(num_samples_local, 2) == 0; num_samples_local = num_samples_local + 1; end
kernel = ones(1, num_samples_local) / num_samples_local;

PRX_avg = conv(PRX, kernel, 'same');
fprintf('      * Convolution took %.2f seconds.\n', toc);

valid_indices = find(distance_domain_padded >= x_start_model & distance_domain_padded <= x_end_model);
distance_domain_model = distance_domain_padded(valid_indices);

PRX_avg_model = PRX_avg(valid_indices);
PRX_Friis_model = PRX_Friis(valid_indices);

PRX_dBm = 10 * log10(PRX * 1000);
PRX_Friis_model_dBm = 10 * log10(PRX_Friis_model * 1000);
PRX_avg_model_dBm = 10 * log10(PRX_avg_model * 1000);

[n, L0_d0, sigma_L] = pathLoss(distance_domain_padded, PRX_dBm, distance_domain_model, PRX_Friis_model_dBm, PRX_avg_model_dBm, params, d0);



%% NARROWBAND HEATMAP
fprintf('\nGenerating 2D coverage heatmap and local averages\n');
TX_pos     = [10, 0];

x_start_model    = -20;
x_end_model      = 400;

y_start    = -w/2-30;
y_end      = w/2+30;

d_samp     = 0.05;
d_local    = 5.0;

fprintf('   - Generating %.1fm x %.1fm grid with %.1fm sampling interval\n', ...
        (x_end_model - x_start_model), (y_end - y_start), d_samp);
        
num_x_points = round((x_end_model - x_start_model) / d_samp) + 1;
num_y_points = round((y_end - y_start) / d_samp) + 1;
RX_x_coordinates = linspace(x_start_model, x_end_model, num_x_points);
RX_y_coordinates = linspace(y_start, y_end, num_y_points);

PRX_dBm = zeros(num_y_points, num_x_points);
fprintf('      * Generating Heatmap and Averaged Power Data \n');

queue_hm = parallel.pool.DataQueue;
num_iterations = num_x_points;
waitbar_hm = waitbar(0, '2D Heatmap');
set(waitbar_hm, 'UserData', 0);
afterEach(queue_hm, @(~) updateWaitbar(waitbar_hm, num_iterations, '2D Heatmap'));

tic;
parfor i = 1:num_x_points
    temp_PRX_dBm_col = zeros(num_y_points, 1);
    RX_y_coordinates_temp = RX_y_coordinates;

    for j = 1:num_y_points
        RX_pos = [RX_x_coordinates(i), RX_y_coordinates_temp(j)];
        d_fixed = norm(RX_pos - TX_pos);
        
        if d_fixed < d_samp
            temp_PRX_dBm_col(j) = NaN;
            continue;
        end
        
        [all_alphas, ~] = runRayTracing(walls, M, TX_pos, RX_pos, params);
        
        if ~isempty(all_alphas)
            h_nb = sum(all_alphas);
            PRX = PTX * abs(h_nb)^2;
            temp_PRX_dBm_col(j) = 10 * log10(PRX * 1000);
        else
            temp_PRX_dBm_col(j) = -Inf;
        end
    end
    
    PRX_dBm(:, i) = temp_PRX_dBm_col;

    if ~isempty(queue_hm)
        send(queue_hm, 1);
    end
end

PRX = 10.^((PRX_dBm - 30) / 10);
PRX(isnan(PRX)) = 0;

num_samples_local = round(d_local / d_samp);
if mod(num_samples_local, 2) == 0; num_samples_local = num_samples_local + 1; end

kernel = ones(num_samples_local) / (num_samples_local^2);
PRX_avg_model = conv2(PRX, kernel, 'same');

PRX_avg_dBm = 10 * log10(PRX_avg_model * 1000);
PRX_avg_dBm(PRX_avg_model == 0) = -Inf;

fprintf('      * Data generation took %.2f seconds.\n', toc);

if exist('waitbar_hm', 'var') && ishandle(waitbar_hm)
    close(waitbar_hm);
end

fprintf('   - Plotting heatmaps\n');
figure('Name', 'Instantaneous Power Heatmap', 'NumberTitle', 'off');
plotHeatmap(gca, RX_x_coordinates, RX_y_coordinates, PRX_dBm, TX_pos, walls, PRX_sens_dBm);
title('Instantaneous Power $P_{RX}$', 'Interpreter', 'latex');

walls(1).coordinates = walls(1).coordinates + d_local/2 * [[1  1]; [1  1]];
walls(2).coordinates = walls(2).coordinates + d_local/2 * [[1 -1]; [1 -1]];

figure('Name', 'Averaged Power Heatmap', 'NumberTitle', 'off');
plotHeatmap(gca, RX_x_coordinates, RX_y_coordinates, PRX_avg_dBm, TX_pos, walls, PRX_sens_dBm);
title('Averaged Power $\langle P_{RX} \rangle$', 'Interpreter', 'latex');



%% Function to update the waitbar
function updateWaitbar(h_bar, total_iterations, title_str)
    current_count = get(h_bar, 'UserData') + 1;
    set(h_bar, 'UserData', current_count);
    fraction_done = current_count / total_iterations;
    waitbar(fraction_done, h_bar, sprintf('%s  %.1f%%', title_str, fraction_done * 100));
end
