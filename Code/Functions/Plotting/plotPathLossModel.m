function [n_pl, sigma_L, PL_d0] = plotPathLossModel(distances, Prx_total_dBm, sim_params)
% PLOTPATHLOSSMODEL - Calculates and plots the large-scale path loss model.
%
% INPUTS:
%   distances     - Vector of distances for the x-axis.
%   Prx_total_dBm - Vector of total received power in dBm.
%   sim_params    - Struct with simulation parameters.
%
% OUTPUTS:
%   n_pl    - The calculated path loss exponent.
%   sigma_L - The calculated shadowing standard deviation in dB.
%   PL_d0   - The path loss at the reference distance d0=1m.

    % --- 1. Local Area Averaging ---
    window_size_m = 5;
    avg_spacing = mean(diff(distances));
    window_size_samples = round(window_size_m / avg_spacing);
    if window_size_samples < 2, window_size_samples = 2; end
    Prx_avg_dBm = movmean(Prx_total_dBm, window_size_samples, 'omitnan');

    % --- 2. Path Loss Calculation ---
    PL_instantaneous = sim_params.P_TX_dBm - Prx_total_dBm;
    PL_avg = sim_params.P_TX_dBm - Prx_avg_dBm;
    
    % --- 3. Fit Path Loss Model ---
    d0 = 1;
    log_dist = log10(distances / d0);
    
    [~, idx_d0] = min(abs(distances - d0));
    PL_d0 = PL_avg(idx_d0);
    
    valid_indices = ~isnan(PL_avg);
    p = polyfit(log_dist(valid_indices), PL_avg(valid_indices) - PL_d0, 1);
    n_pl = p(1) / 10;
    
    PL_model = PL_d0 + 10 * n_pl * log_dist;

    % --- 4. Calculate Variability (Shadowing) ---
    sigma_L = std(PL_avg - PL_model, 'omitnan');

    % --- 5. Plotting ---
    figure('Name', 'Path Loss Model', 'NumberTitle', 'off', 'Position', [200 200 800 600]);
    semilogx(distances, PL_instantaneous, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Instantaneous Path Loss');
    hold on;
    semilogx(distances, PL_avg, 'b-', 'LineWidth', 2, 'DisplayName', 'Local Average Path Loss');
    semilogx(distances, PL_model, 'r--', 'LineWidth', 2.5, 'DisplayName', sprintf('Fitted Model (n=%.2f)', n_pl));
    grid on;
    title('Large-Scale Path Loss Model');
    xlabel('Separation Distance d (m)');
    ylabel('Path Loss (dB)');
    legend('show', 'Location', 'southeast');
    xlim([1, 1000]);
    
    dim = [.2 .5 .3 .3];
    str = {sprintf('Path Loss Exponent n = %.2f', n_pl), sprintf('Shadowing Std Dev \\sigma_L = %.2f dB', sigma_L)};
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black');
    hold off;
end
