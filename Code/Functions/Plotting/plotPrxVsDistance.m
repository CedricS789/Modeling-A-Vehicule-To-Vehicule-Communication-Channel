function plotPrxVsDistance(distances, Prx_total_dBm, Prx_friis_dBm)
% PLOTPRXvsDISTANCE - Plots received power vs. distance and compares to Friis.
%
% INPUTS:
%   distances      - Vector of distances for the x-axis.
%   Prx_total_dBm  - Vector of total received power in dBm.
%   Prx_friis_dBm  - Vector of received power from Friis formula in dBm.

    figure('Name', 'Received Power vs. Distance', 'NumberTitle', 'off', 'Position', [100 100 800 600]);
    semilogx(distances, Prx_total_dBm, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Total Received Power (Full Channel)');
    hold on;
    semilogx(distances, Prx_friis_dBm, 'r--', 'LineWidth', 2, 'DisplayName', 'LOS Power (Friis Formula)');
    grid on;
    title('Received Power vs. Distance');
    xlabel('Separation Distance d (m)');
    ylabel('Received Power (dBm)');
    legend('show', 'Location', 'northeast');
    ylim([-100, -20]);
    xlim([1, 1000]);
    hold off;
end
