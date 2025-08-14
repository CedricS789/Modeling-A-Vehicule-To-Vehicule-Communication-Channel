function [n, L0_d0_dB, sigma_L] = pathLoss(distance_domain_padded, PRX_dBm, distance_domain_model, PRX_Friis_model_dBm, PRX_avg_model_dBm, params, d0)
% INPUTS:
%   distance_domain_padded  - Vector of distances for the full simulation (m).
%   PRX_dBm                 - Vector of instantaneous received power (dBm).
%   distance_domain_model   - Vector of distances for the model fit (m).
%   PRX_avg_model_dBm       - Vector of spatially averaged received power (dBm).
%   params                  - Struct with simulation parameters.
%   d0                      - Reference distance for the model (m).
%
% OUTPUTS:
%   n        - The calculated path loss exponent.
%   L0_d0_dB - The calculated reference path loss at d0 (dB).

    % Path Loss Model Fitting
    fprintf('   - Fitting model to data\n');
    
    PTX_dBm = params.PTX_dBm;
    G_dBi   = params.G_dBi;
    
    % The PRX_avg_model_dBm is already in dBm
    L_d_dB_data  = PTX_dBm - PRX_avg_model_dBm;
    L0_d_dB_data = L_d_dB_data + 2 * G_dBi;

    % Perform linear regression on L0 vs log10(d)
    p_coeffs  = polyfit(log10(distance_domain_model), L0_d_dB_data, 1);
    slope     = p_coeffs(1);
    intercept = p_coeffs(2);

    % Path loss exponent n
    n = slope / 10;
    fprintf('      * Path Loss Exponent n = %.2f\n', n);

    % Reference path loss L0(d0) at the specified d0
    L0_d0_dB = intercept + slope * log10(d0);
    fprintf('      * Reference Path Loss L0(d0=%.1fm) = %.2f dB\n', d0, L0_d0_dB);    

    % Shadowing std dev
    L0_d_dB_fitted = polyval(p_coeffs, log10(distance_domain_model));
    L_d_dB_fitted  = L0_d_dB_fitted - 2 * G_dBi;

    shadowing_error = L_d_dB_data - L_d_dB_fitted;
    sigma_L         = std(shadowing_error);
    fprintf('      * Shadowing Standard Deviation sigma_L = %.2f dB\n', sigma_L);

    % Plot Path Loss Model Fit
    fprintf('   - Plotting path loss model fit\n');
    figure('Name', 'Path Loss Model Fit', 'NumberTitle', 'off');
    
    % Use semilogx for a logarithmic distance axis
    semilogx(distance_domain_model, L_d_dB_data, 'b-', 'LineWidth', 1.2, ...
        'DisplayName', '$$L_{\mathrm{data}}(d)$$');
    hold on;

    % Fitted curve
    semilogx(distance_domain_model, L_d_dB_fitted, 'r-', 'LineWidth', 2, ...
        'DisplayName', sprintf('$$L_{\\mathrm{fitted}}(d),\\; n=%.2f$$', n));
    
    grid on; grid minor;
    title('Path Loss Model Fit for V2V Urban Canyon', 'FontSize', 20, 'Interpreter', 'latex');
    xlabel('$$d\\;[\\mathrm{m}]$$', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$$\\mathrm{Path\\ Loss}\\;[\\mathrm{dB}]$$', 'FontSize', 18, 'Interpreter', 'latex');
    legend('show', 'Location', 'best', 'FontSize', 18, 'Interpreter', 'latex');
    axis tight;
    hold off;

    % Plot Instantaneous vs. Averaged Power
    PRX_avg_dBm_fit = PTX_dBm - L_d_dB_fitted;

    fprintf('   - Plotting instantaneous vs. averages power\n');
    figure('Name', 'Instantaneous vs. Averages Power', 'NumberTitle', 'off');
    
    % Instantaneous power
    semilogx(distance_domain_padded, PRX_dBm, ...
        'DisplayName', '$$P_{\\mathrm{RX}}\\;\\text{(instantaneous)}$$');
    hold on;
    
    % <PRX> power
    semilogx(distance_domain_model, PRX_avg_model_dBm, 'b-', 'LineWidth', 1.2, ...
        'DisplayName', '$$\\langle P_{\\mathrm{RX}} \\rangle$$');
    
    % Friis model
    semilogx(distance_domain_model, PRX_Friis_model_dBm, '--', 'LineWidth', 1.2, ...
        'DisplayName', '$$P_{\\mathrm{RX,\\ Friis}}$$');

    % <<PRX>>
    semilogx(distance_domain_model, PRX_avg_dBm_fit, 'r-', 'LineWidth', 2, ...
        'DisplayName', '$$\\ll P_{\\mathrm{RX}} \\gg$$');
    
    grid on; grid minor;
    title('Instantaneous and Averaged Power vs. Distance', 'FontSize', 20, 'Interpreter', 'latex');
    xlabel('$$d\\;[\\mathrm{m}]$$', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$$\\mathrm{Received\\ Power}\\;[\\mathrm{dBm}]$$', 'FontSize', 18, 'Interpreter', 'latex');
    legend('show', 'Location', 'northeast', 'FontSize', 18, 'Interpreter', 'latex');
    axis tight;
    hold off;

    % Plot Shadowing Error
    fprintf('   - Plotting shadowing error\n');
    figure('Name', 'Shadowing Error', 'NumberTitle', 'off');
    
    semilogx(distance_domain_model, shadowing_error, 'k-', 'LineWidth', 1.2, ...
        'DisplayName', '$$L_{\\mathrm{data}}-L_{\\mathrm{fitted}}$$'); 
    hold on;

    % Zero line (dashed), guard if yline is unavailable
    if exist('yline','file')
        yline(0,'--','HandleVisibility','off');
    else
        xx0 = [min(distance_domain_model) max(distance_domain_model)];
        plot(xx0, [0 0], 'k--', 'HandleVisibility', 'off');
    end

    % Bold red line at +sigma_L with label using $$...$$
    if exist('yline','file')
        yline(sigma_L, 'r-', 'LineWidth', 3, 'DisplayName', '$$\sigma_L$$');
    else
        xxS = [min(distance_domain_model) max(distance_domain_model)];
        plot(xxS, [sigma_L sigma_L], 'r-', 'LineWidth', 3, 'DisplayName', '$$\sigma_L$$');
    end

    % Place the red label slightly above the line near the right side (log scale)
    x_min = min(distance_domain_model(:));
    x_max = max(distance_domain_model(:));
    x_text = 10^( log10(x_min) + 0.88*(log10(x_max) - log10(x_min)) );
    text(x_text, sigma_L, sprintf('$$\\sigma_L=%.2f\\,\\mathrm{dB}$$', sigma_L), ...
        'Interpreter','latex','Color','r','FontWeight','bold', ...
        'HorizontalAlignment','right','VerticalAlignment','bottom');

    grid on; grid minor;
    % Title WITHOUT the numeric value
    title('Shadowing Error', 'FontSize', 18, 'Interpreter', 'latex');
    xlabel('$$d\\;[\\mathrm{m}]$$', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('$$\\mathrm{Error}\\;[\\mathrm{dB}]$$', 'FontSize', 16, 'Interpreter', 'latex');
    legend('show', 'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
    hold off;
end