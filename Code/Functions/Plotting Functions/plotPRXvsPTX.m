function plotPRXvsPTX(ptx_dBm_domain, prx_total_dBm_domain, prx_los_dbm)
    figure('Name', 'Received Power vs. Transmitted Power', 'NumberTitle', 'off');
    
    plot(ptx_dBm_domain, prx_total_dBm_domain, 'b-', 'LineWidth', 2, 'DisplayName', 'Total Received Power');
    hold on;
    
    plot(ptx_dBm_domain, prx_los_dbm, 'r--', 'LineWidth', 2, 'DisplayName', 'Friis Formula');
    
    grid on;
    title('$P_{RX} = f(P_{TX})$', 'FontSize', 20, 'Interpreter', 'latex');
    xlabel('$P_{TX}$ (dBm)', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$P_{RX}$ (dBm)', 'FontSize', 18, 'Interpreter', 'latex');
    legend('show', 'Location', 'best', 'FontSize', 18, 'Interpreter', 'latex');
    axis tight;
    hold off;
    
end