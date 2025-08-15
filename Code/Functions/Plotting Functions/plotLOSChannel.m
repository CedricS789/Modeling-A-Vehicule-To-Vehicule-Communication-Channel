function plotLOSChannel(ray, params)
    d_tau = params.resolution;
    L = params.Ltaps;
    B = params.BRF;

    tau_LOS   = ray.tau_n;
    alpha_LOS = ray.alpha_n;

    figure('Name','Impulse Response h(\tau) - LOS','NumberTitle','off');
    dt_s   = max(1e-9, tau_LOS/200);
    t_s    = linspace(tau_LOS - 20*dt_s, tau_LOS + 20*dt_s, 401);
    t_ns   = t_s * 1e9;
    tau_ns = tau_LOS * 1e9;
    h_abs = zeros(size(t_ns));
    [~, k] = min(abs(t_ns - tau_ns));
    h_abs(k) = abs(alpha_LOS);
    stem(t_ns, h_abs, 'filled', 'LineWidth', 2, 'DisplayName', '$h_{\mathrm{LOS}}(\tau)$');
    hold on; 
    xline(tau_ns, 'k--', 'HandleVisibility','off'); 
    grid on;
    title('Impulse Response for LOS Channel', 'Interpreter','latex', 'FontSize',20);
    xlabel('$\tau$ (ns)', 'Interpreter','latex', 'FontSize',18);
    ylabel('$|h(\tau)|$', 'Interpreter','latex', 'FontSize',18);
    lg = legend('show','Location','best','FontSize',16); 
    set(lg,'Interpreter','latex'); 
    axis tight; 
    hold off;

    f = linspace(-B/2, B/2, 2001);
    Hf = alpha_LOS .* exp(-1j*2*pi*f*tau_LOS);
    figure('Name','Frequency Response H(f) - LOS','NumberTitle','off');
    plot(f/1e6, abs(Hf), 'LineWidth', 2, 'DisplayName', '$|H_{\mathrm{LOS}}(f)|$');
    grid on;
    title('Frequency Response Magnitude for LOS Channel', 'Interpreter','latex', 'FontSize',16);
    xlabel('$f$ (MHz)', 'Interpreter','latex', 'FontSize',14);
    ylabel('$|H(f)|$', 'Interpreter','latex', 'FontSize',14);
    lg2 = legend('show','Location','best','FontSize',14); 
    set(lg2,'Interpreter','latex'); 
    axis tight;

    l = 0:L;
    tau_l = l * d_tau;
    h_tdl = alpha_LOS .* sinc(B * (tau_LOS - tau_l));
    figure('Name','Tapped Delay Line Model - LOS Channel','NumberTitle','off');
    stem(tau_l*1e9, abs(h_tdl), 'filled', 'LineWidth', 2, 'DisplayName', '$h_{\mathrm{TDL}}(\tau)$');
    grid on;
    title('Tapped Delay Line Model - LOS Channel', 'Interpreter','latex', 'FontSize',16);
    xlabel('$\tau$ (ns)', 'Interpreter','latex', 'FontSize',14);
    ylabel('$|h_{\mathrm{TDL}}|$', 'Interpreter','latex', 'FontSize',14);
    lg3 = legend('show','Location','best','FontSize',14); 
    set(lg3,'Interpreter','latex'); 
    axis tight;
end