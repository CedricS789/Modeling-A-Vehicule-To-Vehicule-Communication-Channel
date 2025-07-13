function alpha_n = calculateAlpha_n(ray_data, params)
% calculateAlpha_n - Calculates the complex channel gain for a single ray.
%
% This function implements the formula for the complex channel gain (alpha_n),
% which combines the effects of path loss, phase shift, and reflections.
%
% INPUTS:
%   ray_data   - A struct containing the ray's properties (dist, gamma_prod).
%   params - A struct with simulation parameters (fc, c, Z_0, R_a).
%
% OUTPUTS:
%   alpha_n    - The complex channel gain coefficient for the ray.

    % Unpack necessary parameters
    d_n = ray_data.total_distance;
    gamma_prod = ray_data.gamma_prod;
    fc = params.fc;
    c = params.c;
    Z_0 = params.Z_0;
    R_a = params.R_a;
    lambda = params.lambda;
    
    % Time delay (phase shift) due to path length.
    tau_n = d_n / c;
    phase_shift = exp(-1j * 2 * pi * fc * tau_n);
    
    % Amplitude reduction based on a form of the Friis formula.
    % This term accounts for free-space path loss and antenna characteristics.
    amplitude = (lambda * Z_0) / (4 * pi^2 * R_a * d_n);
    
    % Combine all effects: amplitude, phase, and reflections.
    alpha_n = 1j * amplitude * phase_shift * gamma_prod;
end
