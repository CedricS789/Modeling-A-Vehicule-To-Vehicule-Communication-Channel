function alpha_n = calculateAlpha_n(ray_data, simulation_parameters)
% calculateAlpha_n - Calculates the complex channel gain for a single ray.
%
% This function implements the formula for the complex channel gain, 'alpha_n',
% which represents the change in amplitude and phase that a signal experiences
% as it travels along a specific ray.
%
% INPUTS:
%   ray_data             - A struct containing the ray's properties.
%   simulation_parameters - A struct with simulation parameters.
%
% OUTPUTS:
%   alpha_n                 - The complex channel gain coefficient for the ray.

    d_n = ray_data.dist;
    gamma_prod = ray_data.gamma_prod;
    fc = simulation_parameters.fc;
    c = simulation_parameters.c;
    Z_0 = simulation_parameters.Z_0;
    R_a = simulation_parameters.R_a;
    
    wavelength = c / fc;
    time_delay = d_n / c;
    
    phase_shift_term = exp(-1j * 2 * pi * fc * time_delay);
    amplitude_term = (wavelength * Z_0) / (4 * pi^2 * R_a * d_n);
    
    alpha_n = 1j * amplitude_term * phase_shift_term * gamma_prod;
end
