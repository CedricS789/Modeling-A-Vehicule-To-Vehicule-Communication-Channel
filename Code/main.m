% In your main.m script

%% Step 1 - System Parameters
clear; clc; close all;
sim_params.fc = 5.9e9;          % Carrier frequency in Hz
sim_params.c = 3e8;             % Speed of light in m/s
sim_params.Z_0 = 377;           % Impedance of free space in Ohms
sim_params.R_a = 73.1;          % Radiation resistance for a half-wave dipole (Ohms)
sim_params.B_RF = 100e6;       % RF bandwidth in Hz
sim_params.P_TX = 0.1;         % Transmit power in Watts
sim_params.G_TX = 1.7;         % Transmit antenna gain
sim_params.G_RX = 1.7;         % Receive antenna gain
sim_params.eps_r = 1;          % Relative permittivity of the  buildgins
sim_params.l_area = 5;         % length of the local areas in meters
sim_params.lamda = sim_params.c/sim_params.fc;       % Wavelength in meters
sim_params.Z_0 = 377;          % Impedance of free space in Ohms
sim_params.he = sim_params.lamda/pi;      % Transversal effective height in the horizontal plane (tetha=pi/2) in meters
sim_params.beta = 2*pi*sim_params.fc/sim_params.c;   % Wave number in rad/m

%% Step 3 - Full Channel, Narrowband Analysis
% Define the environment
w = 10; % Half the street width (total width 2w = 20m)
d = 100;
building_eps_r = 4; % Relative permittivity of buildings

% Define walls as a struct array
walls(1).coords = [0 w; d w];   walls(1).eps_r = building_eps_r; % Top wall
walls(2).coords = [0 -w; d -w]; walls(2).eps_r = building_eps_r; % Bottom wall

% Define TX and RX positions
tx_pos = [d/10, 0];
rx_pos = [d - d/10, 0];

% Set max number of reflections
k_max = 3;

% --- Run the Ray Tracing Simulation ---
[alphas, paths_data] = ray_tracing_v2(walls, k_max, tx_pos, rx_pos, sim_params);

% --- Display Results ---
fprintf('Found %d paths.\n', length(alphas));
disp('Complex Gain Coefficients (alphas):');
disp(alphas.'); % Display as a column vector

% The narrowband channel transfer function is the sum of all alphas
h_NB = sum(alphas);
fprintf('\nNarrowband Channel Transfer Function (h_NB): %.4e + %.4ei\n', real(h_NB), imag(h_NB));
fprintf('Magnitude |h_NB|^2: %.4e\n', abs(h_NB)^2);