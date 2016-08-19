function [ rx_data, Z ] = dmtFreqPrecReceiver(y, demodulator, modChoice,...
    scale, FEQ, b_bar, dmin, usedTones)
% DMT Receiver for the Frequency-domain Precoded Transmissions
% ------------------------------------------------------------------------
%
%   Applies the mod-M operator to the FEQ output and then applies
%	regular decision.
%
%   Inputs:
% y             - Non-extended time-domain vectors (before FFT and FEQ)
%				  Dimensions: N x nSymbols
% demodulator   - Cell array with unique demodulators
% modChoice     - Vector with the modulator choices for each subchannel
% scale         - Scaling factors for each constellation
% FEQ           - FEQ Taps (Hermitian vector)
% b_bar         - Vector of bits/dimension for each subchannel
% dmin          - Vector of minimum distances for each subchannel
% usedTones     - Vector indicating one side of the DFT tones that are used
%
%   Outputs:
% rx_data       - Decision output data

%% Initialization
% Infer number of symbols and FFT size
nSymbols   = size(y, 2);
N          = size(y, 1);

% Preallocate
rx_data = zeros(N/2 + 1, nSymbols);

% Initialize full vector of bits per dimension:
b_bar_padded            = zeros(N/2 +1, 1);
% Set the actual number of bits/dim corresponding to the used tones:
b_bar_padded(usedTones) = b_bar;

% Initialize full vector of minimum distances:
D            = zeros(N/2 +1, 1);
D(usedTones) = dmin;

%% Receiver

% FFT
Y = (1/sqrt(N)) * fft(y, N);

% FEQ - One-tap Frequency Equalizer
Z = diag(FEQ) * Y(usedTones, :);

% Modulo operation
Z = moduloOperator(Z, b_bar, dmin);

% Constellation decoding (decision)
for iModem = 1:length(demodulator)
    iTones = (modChoice == iModem);
    % Demodulate
    rx_data(iTones, :) = ...
        demodulator{iModem}.demodulate(...
        diag((1./scale(iTones))) * Z(iTones, :));
end

end

