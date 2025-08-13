function alpha_n = calculateAlpha_n(ray_data, params)
% INPUTS:
%   ray_data    - A struct containing the ray's properties: dist, gamma_tot_n.
%   params      - A struct with simulation parameters: fc, c, Z_0, R_a.
%
% OUTPUTS:
%   alpha_n     - The complex channel gain coefficient for the ray.

    d_n = ray_data.distance_total;
    gamma_tot_n = ray_data.gamma_tot_n;
    fc = params.fc;
    c = params.c;
    Z_0 = params.Z_0;
    R_a = params.R_a;
    lambda = params.lambda;
    
    tau_n = d_n / c;
    phase_shift = exp(-1j * 2 * pi * fc * tau_n);
    
    amplitude = (lambda * Z_0) / (4 * pi^2 * R_a * d_n);
    
    alpha_n = 1j * amplitude * phase_shift * gamma_tot_n;
end
