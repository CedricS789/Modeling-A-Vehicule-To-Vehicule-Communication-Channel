function plotWidebandChannel(rays, params)
    % Parameters
    B   = params.BRF;         % Hz
    dT  = params.resolution;  % Tap spacing (s)
    L   = params.Ltaps;       % Number of taps

    % Extract path delays and gains
    N = length(rays);
    all_tau_n   = zeros(1,N);
    all_alpha_n = zeros(1,N);
    for n = 1:N
        all_tau_n(n)   = rays(n).tau_n;
        all_alpha_n(n) = rays(n).alpha_n;
    end

    % Tap delays
    l = 0:L;
    tau_dom = l * dT;

    % h_TDL
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
    lg = legend('show','Location','best','FontSize',14); 
    set(lg,'Interpreter','latex'); 
    axis tight;

    % h(tau)
    t_start = 0;
    t_end = max(tau_dom) + 5*dT;
    Nt   = 10001;
    t = linspace(t_start, t_end, Nt);

    h_tau = zeros(1, Nt);
    for li = 1:(L+1)
        h_tau = h_tau + h_TDL(li) * sinc(B * (t - tau_dom(li)));
    end

    figure('Name','Physical Impulse Response h(\tau) - Full Channel','NumberTitle','off');
    plot(t*1e9, abs(h_tau), 'LineWidth', 2, 'DisplayName', '$|h(\tau)|$');
    grid on;
    title('Physical Impulse Response - Full Channel', 'Interpreter','latex', 'FontSize',16);
    xlabel('$\tau$ (ns)', 'Interpreter','latex', 'FontSize',14);
    ylabel('$|h(\tau)|$', 'Interpreter','latex', 'FontSize',14);
    lg2 = legend('show','Location','best','FontSize',14); 
    set(lg2,'Interpreter','latex'); 
    axis tight;

    % H(f)
    Nf = 10001;
    f = linspace(-B/2, B/2, Nf);
    Hf = zeros(1, Nf);
    for k = 1:Nf
        for li = 1:(L+1)
            Hf(k) = Hf(k) + h_TDL(li) * exp(-1j*2*pi * f(k) * tau_dom(li));
        end
    end

    figure('Name','Frequency Response H(f) - Full Channel','NumberTitle','off');
    plot(f/1e6, abs(Hf), 'LineWidth', 2, 'DisplayName', '$|H(f)|$');
    grid on;
    title('Frequency Response Magnitude - Full Channel', 'Interpreter','latex', 'FontSize',16);
    xlabel('$f$ (MHz)', 'Interpreter','latex', 'FontSize',14);
    ylabel('$|H(f)|$', 'Interpreter','latex', 'FontSize',14);
    lg3 = legend('show','Location','best','FontSize',14); 
    set(lg3,'Interpreter','latex'); 
    axis tight;
end
