function [n, L0_d0_dB] = pathLoss(distance_domain_padded, PRX_dBm, distance_domain_model, PRX_avg_model_dBm, params, d0)
% This function performs three actions:
% 1. Generates a plot of instantaneous vs. spatially averaged power.
% 2. Fits a path loss model to the averaged data to find 'n' and 'L0'.
% 3. Generates a plot of the path loss model fit.
%
% INPUTS:
%   distance_domain_padded  - Vector of distances for the full simulation (m).
%   prx_instantaneous_dBm   - Vector of instantaneous received power (dBm).
%   distance_domain_model   - Vector of distances for the model fit (m).
%   prx_avg_model_dBm       - Vector of spatially averaged received power (dBm).
%   params                  - Struct with simulation parameters.
%   d0                      - Reference distance for the model (m).
%
% OUTPUTS:
%   n        - The calculated path loss exponent.
%   L0_d0_dB - The calculated reference path loss at d0 (dB).

    % Path Loss Model Fitting
    fprintf('   - Fitting model to data\n');
    
    PTX_dBm = params.PTX_dBm;
    G_dBi = params.G_dBi;
    
    % The prx_avg_model_dBm is already in dBm, so we can use it directly
    L0_d_dB_data = PTX_dBm + 2 * G_dBi - PRX_avg_model_dBm;
    L_d_dB_data = L0_d_dB_data - 2 * G_dBi;

    % Perform linear regression on L0 vs log10(d)
    p_coeffs = polyfit(log10(distance_domain_model), L0_d_dB_data, 1);
    slope = p_coeffs(1);
    intercept = p_coeffs(2);

    % Find the path loss exponent n
    n = slope / 10;
    fprintf('      * Path Loss Exponent (n) = %.2f\n', n);

    % Find the reference path loss L0(d0) at the specified d0
    L0_d0_dB = intercept + slope * log10(d0);
    fprintf('      * Reference Path Loss L0(d0=%.1fm) = %.2f dB\n', d0, L0_d0_dB);    

    % Plot Path Loss Model Fit
    fprintf('   - Plotting path loss model fit\n');
    figure('Name', 'Path Loss Model Fit', 'NumberTitle', 'off');
    
    % Use semilogx for a logarithmic distance axis
    semilogx(distance_domain_model, L_d_dB_data, 'LineWidth', 2, 'DisplayName', '$L_{data}(d)$');
    hold on;

    % Calculate the y-values for the fitted line
    L0_d_dB_fitted = polyval(p_coeffs, log10(distance_domain_model));
    L_d_dB_fitted = L0_d_dB_fitted - 2 * G_dBi;

    semilogx(distance_domain_model, L_d_dB_fitted, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('$L(d)$  (n=%.2f)', n));
    
    grid on;
    grid minor;
    title('Path Loss Model Fit for V2V Urban Canyon', 'FontSize', 20, 'Interpreter', 'latex');
    xlabel('$d$ (m)', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('Path Loss [dB]', 'FontSize', 18, 'Interpreter', 'latex');
    legend('show', 'Location', 'best', 'FontSize', 18, 'Interpreter', 'latex');
    axis tight;
    hold off;

    % Plot Instantaneous vs. Averaged Power
    PRX_avg_dBm_fit = PTX_dBm - L_d_dB_fitted;

    fprintf('   - Plotting instantaneous vs. averages power\n');
    figure('Name', 'Instantaneous vs. Averages Power', 'NumberTitle', 'off');
    
    % Plot instantaneous power
    semilogx(distance_domain_padded, PRX_dBm, 'DisplayName', 'Instantaneous Power $P_{RX}$');
    hold on;
    
    % Plot <PRX> power on top
    semilogx(distance_domain_model, PRX_avg_model_dBm, 'LineWidth', 2, 'DisplayName', 'Averaged Power $\langle P_{RX} \rangle$');
    
    % Plot <<PRX>>
    semilogx(distance_domain_model, PRX_avg_dBm_fit, 'r-', 'LineWidth', 2, 'DisplayName', '$\ll P_{RX} \gg$');

    grid on;
    grid minor;
    title('Instantaneous and Averaged Power vs. Distance', 'FontSize', 20, 'Interpreter', 'latex');
    xlabel('$d$ (m)', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('Received Power (dBm)', 'FontSize', 18, 'Interpreter', 'latex');
    legend('show', 'Location', 'northeast', 'FontSize', 18, 'Interpreter', 'latex');
    axis tight;
    hold off;
end