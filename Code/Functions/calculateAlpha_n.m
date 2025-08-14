function alpha_n = calculateAlpha_n(ray_data, params)
% INPUTS:
%   ray_data    - A struct containing the ray's properties: dist, gamma_tot_n.
%   params      - A struct with simulation parameters: fc, c, Z0, Ra.
%
% OUTPUTS:
%   alpha_n     - The complex channel gain coefficient for the ray.

    d_n = ray_data.distance_total;
    gamma_tot_n = ray_data.gamma_tot_n;
    fc = params.fc;
    c = params.c;
    Z0 = params.Z0;
    Ra = params.Ra;
    lambda = params.lambda;
    
    tau_n = d_n / c;
    phase_shift = exp(-1j * 2 * pi * fc * tau_n);
    
    amplitude = (lambda * Z0) / (4 * pi^2 * Ra * d_n);
    
    alpha_n = 1j * amplitude * phase_shift * gamma_tot_n;
end
