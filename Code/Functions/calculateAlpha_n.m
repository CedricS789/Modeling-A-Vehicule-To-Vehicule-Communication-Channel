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
    path_distance = ray_data.dist;
    reflection_product = ray_data.gamma_prod;
    fc = params.fc;
    c = params.c;
    Z_0 = params.Z_0;
    R_a = params.R_a;
    
    % --- Calculation Steps ---
    
    % 1. Wavelength of the carrier signal.
    wavelength = c / fc;
    
    % 2. Time delay (phase shift) due to path length.
    time_delay = path_distance / c;
    phase_shift_term = exp(-1j * 2 * pi * fc * time_delay);
    
    % 3. Amplitude reduction based on a form of the Friis formula.
    % This term accounts for free-space path loss and antenna characteristics.
    amplitude_term = (wavelength * Z_0) / (4 * pi^2 * R_a * path_distance);
    
    % 4. Combine all effects: amplitude, phase, and reflections.
    % The leading '1j' term is related to the radiation characteristics of a dipole.
    alpha_n = 1j * amplitude_term * phase_shift_term * reflection_product;
end
