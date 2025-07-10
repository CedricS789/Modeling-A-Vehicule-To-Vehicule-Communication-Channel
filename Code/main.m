%% Step 1 - Theoretical preliminaries
clear; clc; close all;
fc = 5.9e9;         % Carrier frequency in Hz
c = 3e8;            % Speed of light in m/s
B_RF = 100e6;       % RF bandwidth in Hz
P_TX = 0.1;         % Transmit power in Watts
G_TX = 1.7;         % Transmit antenna gain in dBi
G_RX = 1.7;         % Receive antenna gain in dBi
eps_r = 1;          % Relative permittivity of the  buildgins
l_area = 5;         % length of the local areas in meters
lamda = c/fc;       % Wavelength in meters
Z_0 = 377;          % Impedance of free space in Ohms
he = lamda/pi;      % Transversal effective height in the horizontal plane (tetha=pi/2) in meters
beta = 2*pi*fc/c;   % Wave number in rad/mµ

%% Step 2 - LOS channel
Fs = fc * 3;                % Sampling frequency in Hz
d_0 = 1;                    % Distance LOS in meters  
tau_0 = d_0/c;              % Delay of the LOS path in seconds
tau = 0:1/Fs:2*tau_0;       % Time vector in seconds
h_tau = zeros(size(tau));   % Initialize impulse response vector
[~, i] = min(abs(tau - tau_0));
h_tau(i) = lamda/(4*pi*d_0) * sqrt(G_TX*G_RX) * exp(-1j*beta*d_0);
stem(tau, abs(h_tau), 'LineWidth', 1.5);