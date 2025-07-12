%% MAIN SIMULATION SCRIPT - V2V STREET CANYON (Simplified)
% =========================================================================
% This script sets up and runs the ray-tracing simulation for a simplified
% urban street canyon scenario, represented by two long, continuous walls.
%
% 1. It defines the physical parameters for V2V communication.
% 2. It generates a simple environment with two parallel walls.
% 3. It places the TX and RX vehicles in the middle of the street.
% 4. It calls the ray-tracing engine to compute the paths.
% 5. It calls the plotting function to visualize the results.
% 6. It displays a summary of the key physical parameters for each found ray.
% =========================================================================
clear; clc; close all;
addpath('Functions')

% =================================================================
% 1: DEFINE SIMULATION PARAMETERS
% =================================================================

% --- Physical Channel Parameters ---
sim_params.fc = 5.9e9;          % Carrier frequency in Hz
sim_params.c = 3e8;             % Speed of light in m/s
sim_params.Z_0 = 377;           % Impedance of free space in Ohms
sim_params.R_a = 73.1;          % Radiation resistance for a half-wave dipole (Ohms)
sim_params.B_RF = 100e6;        % RF bandwidth in Hz
sim_params.P_TX = 0.1;          % Transmit power in Watts
sim_params.G_TX = 1.7;          % Transmit antenna gain
sim_params.G_RX = 1.7;          % Receive antenna gain
sim_params.l_area = 5;          % length of the local areas in meters
sim_params.lamda = sim_params.c/sim_params.fc;              % Wavelength in meters
sim_params.Z_0 = 377;                                       % Impedance of free space in Ohms
sim_params.he = sim_params.lamda/pi;                        % Transversal effective height in the horizontal plane (tetha=pi/2) in meters
sim_params.beta = 2*pi*sim_params.fc/sim_params.c;          % Wave number in rad/m

% --- Ray Tracing Configuration ---
k_max = 3; % The maximum order of reflections to compute.

% =================================================================
% 2: GENERATE THE ENVIRONMENT LAYOUT
% =================================================================

% --- Street and Building Geometry ---
w = 20;     % Total width of the street (2w = 20m)
L = 500;  % Total length of the simulated street area
eps_r = 4;      % Relative permittivity of the buildings

% --- Vehicle Positions ---
dist = 150; % Separation distance 'd' between TX and RX in meters
tx_pos = [-dist/2, 0]; % TX is at -d/2, in the middle of the lane
rx_pos = [dist/2, 0];   % RX is at +d/2, in the middle of the lane

% --- Create two simple walls for the street canyon ---
walls = [];
% Top wall
walls(1).coords = [-L, w/2; L, w/2];
walls(1).eps_r = eps_r;
% Bottom wall
walls(2).coords = [-L, -w/2; L, -w/2];
walls(2).eps_r = eps_r;


% =================================================================
% 3: EXECUTE THE RAY-TRACING CALCULATION
% =================================================================
fprintf('Starting ray-tracing calculation for up to %d reflections...\n', k_max);
[alphas, rays] = runRayTracing(walls, k_max, tx_pos, rx_pos, sim);
fprintf('Calculation complete. Found %d valid propagation rays.\n', length(alphas));

% =================================================================
% 4: VISUALIZE THE RESULTS
% =================================================================
fprintf('Generating plot of the environment and found rays...\n');
plotRays(walls, tx_pos, rx_pos, rays, k_max);
% Adjust plot limits to focus on the area around the vehicles
xlim([tx_pos(1)-50, rx_pos(1)+50]);
ylim([-w, w]);
fprintf('Simulation finished.\n');

% =================================================================
% 5: DISPLAY RAY DATA SUMMARY
% =================================================================
fprintf('\n--- Ray Data Summary ---\n');
for i = 1:length(rays)
    r = rays{i};
    g_mag = abs(r.gain);
    g_phase = rad2deg(angle(r.gain));
    
    fprintf('Ray %2d: Type = %-7s | Distance = %7.2f m | Gain Mag = %.4e | Gain Phase = %7.2f deg\n', ...
            i, r.type, r.dist, g_mag, g_phase);
end
fprintf('------------------------\n');
