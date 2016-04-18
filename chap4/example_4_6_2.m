% Establish a DMT symbol period:
T = 250e-6;

% Establish the number of used subchannels:
N_bar = 256;  % Used subchannels
N = 2 * N_bar % FFT size and number of real dimensions

% Establish the subchannel bandwidth
delta_f = 4.3125e3;

% This determines the necessary minimum sampling frequency:
fs = N * delta_f;
Ts = 1/fs;

% Which then implies the cyclic prefix length available:
% Since T/(N + nu) = Ts
nu = (T - N*Ts) / Ts

% Transmission bandwidth
B = fs/2