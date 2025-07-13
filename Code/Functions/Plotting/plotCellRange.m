function plotCellRange(distances, n_pl, PL_d0, sigma_L, sim_params)
% PLOTCELLRANGE - Calculates and plots the cell range for different reliabilities.
%
% INPUTS:
%   distances  - Vector of distances for the x-axis.
%   n_pl       - The path loss exponent.
%   PL_d0      - Path loss at the reference distance d0.
%   sigma_L    - The shadowing standard deviation in dB.
%   sim_params - Struct with simulation parameters.

    % --- 1. Reconstruct Path Loss Model ---
    d0 = 1;
    PL_model = PL_d0 + 10 * n_pl * log10(distances / d0);

    % --- 2. Calculate Max Allowable Path Loss for Reliabilities ---
    reliabilities = [0.50, 0.95, 0.99];
    colors = {'g', 'm', 'c'};
    
    figure('Name', 'Cell Range Analysis', 'NumberTitle', 'off', 'Position', [250 250 800 600]);
    semilogx(distances, PL_model, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Path Loss Model');
    hold on;
    
    all_L_max = zeros(size(reliabilities));
    for i = 1:length(reliabilities)
        R = reliabilities(i);
        M = -sigma_L * sqrt(2) * erfcinv(2 * R);
        
        L_max = sim_params.P_TX_dBm - sim_params.P_sens_dBm - M;
        all_L_max(i) = L_max;
        
        cell_range = d0 * 10^((L_max - PL_d0) / (10 * n_pl));
        
        plot([1, 1000], [L_max, L_max], '--', 'Color', colors{i}, 'LineWidth', 2, 'DisplayName', sprintf('%.0f%% Reliability (M=%.1f dB)', R*100, M));
        plot([cell_range, cell_range], [0, L_max], ':', 'Color', colors{i}, 'LineWidth', 1.5);
        text(cell_range, 20, sprintf('d = %.0f m', cell_range), 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    
    grid on;
    title('Cell Range for Different Service Reliabilities');
    xlabel('Separation Distance d (m)');
    ylabel('Path Loss (dB)');
    legend('show', 'Location', 'southeast');
    xlim([1, 1000]);
    
    ylim_max = max([max(PL_model), max(all_L_max)]) + 5;
    ylim_min = min([min(PL_model), min(all_L_max)]) - 5;
    ylim([ylim_min, ylim_max]);
    hold off;
end
