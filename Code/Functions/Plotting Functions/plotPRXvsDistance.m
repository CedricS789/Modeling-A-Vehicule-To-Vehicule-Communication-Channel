function plotPRXvsDistance(distances, total_RX_power_dBm, los_power_dBm)
    figure('Name', 'Received Power vs. Distance', 'NumberTitle', 'off');
    
    semilogx(distances, total_RX_power_dBm, 'b-', 'LineWidth', 2, 'DisplayName', 'Total Received Power');
    hold on;
    
    semilogx(distances, los_power_dBm, 'r--', 'LineWidth', 2, 'DisplayName', 'Friis Formula');
    
    grid on;
    title('$P_{RX} = f(d)$', 'FontSize', 20, 'Interpreter', 'latex');
    xlabel('$d$ (m)', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$P_{RX}$ (dBm)', 'FontSize', 18, 'Interpreter', 'latex');
    legend('show', 'Location', 'northeast', 'FontSize', 18, 'Interpreter', 'latex');
    axis tight;
    hold off;
end