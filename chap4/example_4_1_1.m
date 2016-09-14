%% Example 4.1.1 - Introduction to multi-channel modulation
clearvars, clc
addpath(genpath('../lib'))

% Parameters
N       = 8;    % FFT size and the number of real dimensions
Ex      = N;    % Total energy budget for each OFDM symbol

% Derivations
N_bar   = N/2 + 1;        % Available subchannels (with DC and Nyquist)

% Noise energy per dimension (N0/2):
sigma_n = 0.181;

% Vector of used real dimensions
nDim = [ 1 2 2 2 0 ];

% Energy per subchannel
% At this point, water-filling has not been presented, so the energy and
% bit load are simply "guessed". An equal amount of energy is allocated to
% each available (used) dimension.
Ex_bar = Ex / sum(nDim); % Energy per dimension
En     = Ex_bar * nDim;  % Energy per subchannel
En_bar = En ./ nDim;     % Energy per dimension per subchannel

% Arbitrary bit load aiming at 1 bit/dimension:
bn     = [ 1.6 3 2.4 1 0 ];
bn_bar = bn ./ nDim;

% Double check the total number of bits, raw and per dimension:
b     = bn_bar(1) + 2*sum(bn_bar(2:4));
b_bar = b / N;

fprintf('Arbitrary bit load at: %g bit/dimension\n\n', b_bar);

%% Channel

h = [0.9 1];

% Channel frequency response
H = fft(h, N);

%% SNR and Q func arg

% Signal to noise ratio
SNR_n = Ex_bar * abs(H).^2 / sigma_n;

% Argument of Q function

argQ = sqrt( (3 ./ ((2.^(2*bn_bar)) - 1)) .* SNR_n(1:N_bar));

display('Arg Q-func (dB)');
20*log10(argQ)

% Compare these with the SNRmfb performance (for binary modulation it is 10
% dB)
