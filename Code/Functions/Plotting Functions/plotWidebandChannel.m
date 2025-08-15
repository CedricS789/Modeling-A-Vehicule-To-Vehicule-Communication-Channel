function plotWidebandChannel(all_rays, params)
    B   = params.BRF;
    dT  = params.resolution;
    L   = params.Ltaps;

    N = length(all_rays);
    all_tau_n   = zeros(1,N);
    all_alpha_n = zeros(1,N);
    for n = 1:N
        all_tau_n(n)   = all_rays{n}.tau_n;
        all_alpha_n(n) = all_rays{n}.alpha_n;
    end

    figure('Name','Physical Impulse Response h(\tau) - Full Channel','NumberTitle','off');
    stem(all_tau_n*1e9, abs(all_alpha_n), 'filled', 'LineWidth', 2, ...
         'Marker','o','MarkerSize',5, 'DisplayName', '$h(\tau)$');
    grid on;
    title('Physical Impulse Response - Full Channel', 'Interpreter','latex', 'FontSize',16);
    xlabel('$\tau$ (ns)', 'Interpreter','latex', 'FontSize',14);
    ylabel('$|h(\tau)|$', 'Interpreter','latex', 'FontSize',14);
    lg1 = legend('show','Location','northeast','FontSize',14);
    set(lg1,'Interpreter','latex');
    axis tight;

    N_fft = 4096;
    f = linspace(-B/2, B/2, N_fft);
    Hf = zeros(1, N_fft);
    for n = 1:N
        Hf = Hf + all_alpha_n(n) * exp(-1j * 2 * pi * f * all_tau_n(n));
    end

    figure('Name','Frequency Response H(f) - Full Channel','NumberTitle','off');
    plot(f/1e6, 20*log10(abs(Hf)), 'LineWidth', 2, 'DisplayName', '$|H(f)|$');
    grid on;
    title('Frequency Response Magnitude - Full Channel', 'Interpreter','latex', 'FontSize',16);
    xlabel('$f$ (MHz)', 'Interpreter','latex', 'FontSize',14);
    ylabel('$|H(f)|$ (dB)', 'Interpreter','latex', 'FontSize',14);
    lg2 = legend('show','Location','northeast','FontSize',14);
    set(lg2,'Interpreter','latex');
    axis tight;

    l = 0:L;
    tau_dom = l * dT;
    h_TDL = zeros(1, L+1);
    for li = 1:(L+1)
        for n = 1:N
            h_TDL(li) = h_TDL(li) + all_alpha_n(n) * sinc(B * (all_tau_n(n) - tau_dom(li)));
        end
    end

    figure('Name','Tapped Delay Line Model - Full Channel','NumberTitle','off');
    stem(tau_dom*1e9, abs(h_TDL), 'filled', 'LineWidth', 2, ...
         'Marker','o','MarkerSize',5, 'DisplayName', '$h_{\mathrm{TDL}}(\tau)$');
    grid on;
    title('Tapped Delay Line Model - Full Channel', 'Interpreter','latex', 'FontSize',16);
    xlabel('$\tau$ (ns)', 'Interpreter','latex', 'FontSize',14);
    ylabel('$|h_{\mathrm{TDL}}|$', 'Interpreter','latex', 'FontSize',14);
    lg3 = legend('show','Location','northeast','FontSize',14);
    set(lg3,'Interpreter','latex');
    axis tight;
end