function [n, L0_d0_dB, sigma_L] = pathLoss(distance_domain_padded, PRX_dBm, distance_domain_model, PRX_Friis_model_dBm, PRX_avg_model_dBm, params, d0)

    fprintf('   - Fitting model to data\n');
    
    PTX_dBm = params.PTX_dBm;
    G_dBi   = params.G_dBi;
    
    L_d_dB_data  = PTX_dBm - PRX_avg_model_dBm;
    L0_d_dB_data = L_d_dB_data + 2 * G_dBi;
    
    log_dist = log10(distance_domain_model);
    log_dist_mean = mean(log_dist);
    p_coeffs  = polyfit(log_dist - log_dist_mean, L0_d_dB_data, 1);
    slope     = p_coeffs(1);
    intercept_centered = p_coeffs(2);
    
    n = slope / 10;
    fprintf('      * Path Loss Exponent n = %.2f\n', n);
    
    L0_d_dB_fitted = slope * (log_dist - log_dist_mean) + intercept_centered;
    
    L0_d0_dB = slope * (log10(d0) - log_dist_mean) + intercept_centered;
    fprintf('      * Reference Path Loss L0(d0=%.1fm) = %.2f dB\n', d0, L0_d0_dB);    
    
    shadowing_error = L_d_dB_data - L0_d_dB_fitted;
    sigma_L         = std(shadowing_error);
    fprintf('      * Shadowing Standard Deviation sigma_L = %.2f dB\n', sigma_L);
    
    fprintf('   - Plotting path loss model fit\n');
    figure('Name', 'Path Loss Model Fit', 'NumberTitle', 'off');
    
    L_d_dB_fitted_total = L0_d_dB_fitted - 2 * G_dBi;
    
    semilogx(distance_domain_model, L_d_dB_data, 'b-', 'LineWidth', 1.2, 'DisplayName', '$L_{\mathrm{data}}(d)$');
    hold on;
    semilogx(distance_domain_model, L_d_dB_fitted_total, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('$L_{\\mathrm{fitted}}(d), n=%.2f$', n));
    
    grid on; grid minor;
    title('Path Loss Model Fit for V2V Urban Canyon', 'FontSize', 20, 'Interpreter', 'latex');
    xlabel('$d\;[\mathrm{m}]$', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$\mathrm{Path\ Loss}\;[\mathrm{dB}]$', 'FontSize', 18, 'Interpreter', 'latex');
    legend('show', 'Location', 'best', 'FontSize', 18, 'Interpreter', 'latex');
    axis tight;
    hold off;
    
    fprintf('   - Plotting instantaneous vs. averages power\n');
    figure('Name', 'Instantaneous vs. Averages Power', 'NumberTitle', 'off');
    
    PRX_avg_dBm_fit = PTX_dBm - L_d_dB_fitted_total;
    
    semilogx(distance_domain_padded, PRX_dBm, 'DisplayName', '$P_{\mathrm{RX}}$');
    hold on;
    semilogx(distance_domain_model, PRX_avg_model_dBm, 'b-', 'LineWidth', 1.2, 'DisplayName', '$\langle P_{\mathrm{RX}} \rangle$');
    semilogx(distance_domain_model, PRX_Friis_model_dBm, '--', 'LineWidth', 1.2, 'DisplayName', '$P_{\mathrm{RX,\ Friis}}$');
    semilogx(distance_domain_model, PRX_avg_dBm_fit, 'r-', 'LineWidth', 2, 'DisplayName', '$\ll P_{\mathrm{RX}} \gg$');
    
    grid on; grid minor;
    title('Power vs. Distance', 'FontSize', 20, 'Interpreter', 'latex');
    xlabel('$d\;[\mathrm{m}]$', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$\mathrm{Received\ Power}\;[\mathrm{dBm}]$', 'FontSize', 18, 'Interpreter', 'latex');
    legend('show', 'Location', 'northeast', 'FontSize', 18, 'Interpreter', 'latex');
    axis tight;
    hold off;
    
    fprintf('   - Plotting shadowing error\n');
    figure('Name', 'Shadowing Error', 'NumberTitle', 'off');
    
    semilogx(distance_domain_model, shadowing_error, 'Color', 'black', 'LineWidth', 1.2, 'DisplayName', '$L_{0,\mathrm{data}}-L_{0,\mathrm{fitted}}$'); 
    hold on;

    if exist('yline','file')
        yline(0,'--','HandleVisibility','off');
    else
        xx0 = [min(distance_domain_model) max(distance_domain_model)];
        plot(xx0, [0 0], 'k--', 'HandleVisibility', 'off');
    end

    if exist('yline','file')
        yline(sigma_L, 'r-', 'LineWidth', 3, 'DisplayName', '$\sigma_L$');
    else
        xxS = [min(distance_domain_model) max(distance_domain_model)];
        plot(xxS, [sigma_L sigma_L], 'r-', 'LineWidth', 3, 'DisplayName', '$\sigma_L$');
    end
    
    x_min = min(distance_domain_model(:));
    x_max = max(distance_domain_model(:));
    x_text = 10^( log10(x_min) + 0.1*(log10(x_max) - log10(x_min)) ); 
    text(x_text, sigma_L, sprintf('$\\sigma_L=%.2f \\mathrm{dB}$', sigma_L), 'Interpreter','latex','Color','r','FontWeight','bold', 'HorizontalAlignment','left','VerticalAlignment','bottom', 'FontSize', 16);
    
    grid on; grid minor;
    title('Shadowing Error', 'FontSize', 18, 'Interpreter', 'latex');
    xlabel('$d\;[\mathrm{m}]$', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('$\mathrm{Error}\;[\mathrm{dB}]$', 'FontSize', 16, 'Interpreter', 'latex');
    legend('show', 'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex');
    hold off;
    
    fprintf('   - Calculating fade margin and cell range\n');
    
    reliabilities = [0.50, 0.95, 0.99];
    PRX_sens_dBm = params.PRX_sens_dBm;
    
    L_max = PTX_dBm - PRX_sens_dBm;
    fprintf('      * Max allowed total path loss L_max = %.2f dB\n', L_max);
    
    fprintf('      --------------------------------------------------\n');
    fprintf('      | Reliability | Fade Margin (M) | Cell Range (R) |\n');
    fprintf('      --------------------------------------------------\n');
    
    for i = 1:length(reliabilities)
        reliability = reliabilities(i);
        
        if reliability == 0.5
            M = 0;
        else
            outage_prob = 1 - reliability;
            M = sigma_L * sqrt(2) * erfcinv(2 * outage_prob);
        end
        
        L0_target_at_edge = (L_max - M) + 2 * G_dBi;
        
        log10_R = (L0_target_at_edge - intercept_centered) / slope + log_dist_mean;
        R = 10^(log10_R);
        
        fprintf('      |    %2.0f %%     |     %.2f dB     |    %.2f m    |\n', reliability*100, M, R);
    end
    fprintf('      --------------------------------------------------\n');
end