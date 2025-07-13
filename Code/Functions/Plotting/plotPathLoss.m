function [path_loss_exponent, shadowing_std_dev, pl_ref] = plotPathLoss(distances, total_rx_power_dBm, sim_params)
% plotPathLoss - Calculates, fits, and plots the large-scale path loss model.
%
% This function processes the received power data to extract key large-scale
% fading parameters: the path loss exponent (n) and the shadowing standard
% deviation (sigma).
%
% INPUTS:
%   distances          - Vector of distances for the x-axis.
%   total_rx_power_dBm - Vector of total received power in dBm.
%   sim_params         - Struct with simulation parameters.
%
% OUTPUTS:
%   path_loss_exponent - The calculated path loss exponent, n.
%   shadowing_std_dev  - The shadowing standard deviation, sigma_L (in dB).
%   pl_ref             - The path loss at the reference distance, d0=1m.

    % --- 1. Calculate Path Loss ---
    % Path Loss (dB) = Transmit Power (dBm) - Received Power (dBm)
    path_loss_instantaneous_dB = sim_params.transmit_power_dBm - total_rx_power_dBm;

    % --- 2. Spatially Average the Path Loss ---
    % This is done to remove the small-scale fading effects and reveal the
    % underlying large-scale trend (path loss + shadowing).
    window_size_m = 5; % Averaging window of 5 meters
    avg_spacing = mean(diff(distances));
    window_size_samples = max(2, round(window_size_m / avg_spacing));
    path_loss_averaged_dB = movmean(path_loss_instantaneous_dB, window_size_samples, 'omitnan');

    % --- 3. Fit the Path Loss Model ---
    % Model: PL(d)[dB] = PL(d0)[dB] + 10 * n * log10(d/d0)
    ref_distance_m = 1;
    log_distances = log10(distances / ref_distance_m);
    
    % Find the path loss at the reference distance d0 from the averaged data.
    [~, ref_idx] = min(abs(distances - ref_distance_m));
    pl_ref = path_loss_averaged_dB(ref_idx);
    
    % Use polyfit to find the slope of the line, which gives us 'n'.
    % We fit (PL_avg - PL_d0) vs log_dist. The slope is 10*n.
    valid_indices = ~isnan(path_loss_averaged_dB);
    p = polyfit(log_distances(valid_indices), path_loss_averaged_dB(valid_indices) - pl_ref, 1);
    path_loss_exponent = p(1) / 10;
    
    % Reconstruct the fitted model for plotting.
    path_loss_model_dB = pl_ref + 10 * path_loss_exponent * log_distances;

    % --- 4. Calculate Shadowing Standard Deviation ---
    % This is the standard deviation of the difference between the actual
    % (averaged) path loss and the value predicted by our model.
    shadowing_std_dev = std(path_loss_averaged_dB - path_loss_model_dB, 'omitnan');

    % --- 5. Plotting ---
    figure('Name', 'Large-Scale Path Loss Model', 'NumberTitle', 'off', 'Position', [200 200 800 600]);
    semilogx(distances, path_loss_instantaneous_dB, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Instantaneous Path Loss');
    hold on;
    semilogx(distances, path_loss_averaged_dB, 'b-', 'LineWidth', 2, 'DisplayName', 'Spatially Averaged Path Loss');
    semilogx(distances, path_loss_model_dB, 'r--', 'LineWidth', 2.5, 'DisplayName', sprintf('Fitted Model (n = %.2f)', path_loss_exponent));
    
    grid on;
    title('Large-Scale Path Loss Analysis');
    xlabel('Separation Distance, d (m)');
    ylabel('Path Loss (dB)');
    legend('show', 'Location', 'southeast');
    xlim([min(distances), max(distances)]);
    
    % Add an annotation box with the calculated parameters.
    annotation_str = {sprintf('Path Loss Exponent n = %.2f', path_loss_exponent), ...
                      sprintf('Shadowing Std. Dev. \\sigma_L = %.2f dB', shadowing_std_dev)};
    annotation('textbox', [0.2 0.6 0.3 0.2], 'String', annotation_str, ...
               'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black');
    hold off;
end
