function plotKFactor(distances, K_factor_dB)
    figure('Name', 'Rician K-Factor vs. Distance', 'NumberTitle', 'off');
    
    semilogx(distances, K_factor_dB, 'k-', 'LineWidth', 2);
    
    grid on;
    grid minor;
    
    title('$K = f(d)$', 'FontSize', 20, 'Interpreter', 'latex');
    xlabel('$d$ (m)', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$K$-factor (dB)', 'FontSize', 18, 'Interpreter', 'latex');
    
    axis tight;
end