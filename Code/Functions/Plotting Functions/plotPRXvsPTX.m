function plotPRXvsPTX(ptx_dBm_domain, prx_total_dBm_domain, prx_los_dbm)
% Plots received power vs. transmitted power.
% This function plots the total received power from all multipath components
% and compares it against the Friis transmission formula.
%
% INPUTS:
%   ptx_dbm         - Vector of transmitted power values for the x-axis.
%   total_prx_dbm   - Vector of total received power in dBm (multipath).
%   los_prx_dbm     - Vector of LOS-only received power in dBm (Friis).

    figure('Name', 'Received Power vs. Transmitted Power', 'NumberTitle', 'off');
    
    % Plot the total received power 
    plot(ptx_dBm_domain, prx_total_dBm_domain, 'b-', 'LineWidth', 2, 'DisplayName', 'Total Received Power');
    hold on;
    
    % Plot Friis
    plot(ptx_dBm_domain, prx_los_dbm, 'r--', 'LineWidth', 2, 'DisplayName', 'Friis Formula');
    
    grid on;
    title('$P_{RX} = f(P_{TX})$', 'FontSize', 20, 'Interpreter', 'latex');
    xlabel('$P_{TX}$ (dBm)', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel('$P_{RX}$ (dBm)', 'FontSize', 18, 'Interpreter', 'latex');
    legend('show', 'Location', 'best', 'FontSize', 18, 'Interpreter', 'latex');
    axis tight;
    hold off;
    
end
