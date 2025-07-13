function plotCellRangeAnalysis(distances, path_loss_exponent, pl_ref, shadowing_std_dev, sim_params)
% plotCellRangeAnalysis - Calculates and plots the cell range for different reliabilities.
%
% This function visualizes how the communication range is affected by the
% required service reliability, considering the fade margin needed to
% overcome shadowing.
%
% INPUTS:
%   distances          - Vector of distances for the x-axis.
%   path_loss_exponent - The path loss exponent, n.
%   pl_ref             - Path loss at the reference distance d0.
%   shadowing_std_dev  - The shadowing standard deviation in dB.
%   sim_params         - Struct with simulation parameters.

    % --- 1. Reconstruct the Mean Path Loss Model ---
    ref_distance_m = 1;
    path_loss_model_dB = pl_ref + 10 * path_loss_exponent * log10(distances / ref_distance_m);

    % --- 2. Define Reliability Levels and Calculate Fade Margins ---
    reliabilities = [0.50, 0.90, 0.99];
    colors = {'g', 'm', 'c'};
    
    figure('Name', 'Cell Range and Reliability Analysis', 'NumberTitle', 'off', 'Position', [350 350 800 600]);
    semilogx(distances, path_loss_model_dB, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Mean Path Loss Model');
    hold on;
    
    max_pl_values = zeros(size(reliabilities));
    for i = 1:length(reliabilities)
        R = reliabilities(i);
        
        % Calculate the required fade margin (M) for reliability R.
        % This is based on the properties of the Gaussian distribution (erfcinv).
        fade_margin_dB = -shadowing_std_dev * sqrt(2) * erfcinv(2 * R);
        
        % The maximum allowable path loss for this reliability.
        max_allowable_pl = sim_params.transmit_power_dBm - sim_params.receiver_sensitivity_dBm - fade_margin_dB;
        max_pl_values(i) = max_allowable_pl;
        
        % Invert the path loss formula to find the cell range.
        cell_range_m = ref_distance_m * 10^((max_allowable_pl - pl_ref) / (10 * path_loss_exponent));
        
        % Plot the max allowable path loss line.
        plot([1, 1e4], [max_allowable_pl, max_allowable_pl], '--', 'Color', colors{i}, 'LineWidth', 2, ...
             'DisplayName', sprintf('%.0f%% Reliability (Fade Margin = %.1f dB)', R*100, fade_margin_dB));
             
        % Plot a vertical line indicating the calculated cell range.
        plot([cell_range_m, cell_range_m], [0, max_allowable_pl], ':', 'Color', colors{i}, 'LineWidth', 1.5);
        text(cell_range_m, 25, sprintf(' d_{max} = %.0f m', cell_range_m), 'Rotation', 90, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    
    grid on;
    title('Communication Range vs. Service Reliability');
    xlabel('Separation Distance, d (m)');
    ylabel('Path Loss (dB)');
    legend('show', 'Location', 'southeast');
    xlim([1, 1000]);
    
    % Adjust y-axis to fit all plotted lines.
    ylim_max = max([max(path_loss_model_dB), max(max_pl_values)]) + 10;
    ylim_min = min([min(path_loss_model_dB), min(max_pl_values)]) - 10;
    ylim([ylim_min, ylim_max]);
    
    hold off;
end
