function [ rx_data, Z ] = dmtFreqPrecReceiver(y, demodulator, modChoice,...
    scale, FEQ, modOrder, dmin, usedTones, usedTonesHerm)
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
% modOrder      - Vector of modulation orders for each subchannel
% dmin          - Vector of minimum distances for each subchannel
% usedTones     - Vector indicating one side of the DFT tones that are used
% usedTonesHerm - Hermitian vector indicating all used DFT tones
%
%   Outputs:
% rx_data       - Decision output data

%% Initialization
% Infer number of symbols and FFT size
nSymbols   = size(y, 2);
N          = size(y, 1);

% Preallocate
rx_data = zeros(N/2 + 1, nSymbols);

% Initialize full vector of modulation orders:
M            = ones(N/2 +1, 1); % Note an unitary order corresponds to b=0
% Set the actual order corresponding to the used tones:
M(usedTones) = modOrder;

% Initialize full vector of minimum distances:
D            = zeros(N/2 +1, 1);
D(usedTones) = dmin;

%% Receiver

% FFT
Y = (1/sqrt(N)) * fft(y, N);

% FEQ - One-tap Frequency Equalizer
Z = diag(FEQ) * Y(usedTones, :);

% Modulo operation
for iSym=1:nSymbols
    Z(:,iSym) = moduloOperator(Z(1:N/2+1, iSym), ...
        M, D, 'Hermitian');
end

% Constellation decoding (decision)
for iModem = 1:length(demodulator)
    iTones = (modChoice == iModem);
    % Demodulate
    rx_data(iTones, :) = ...
        demodulator{iModem}.demodulate(...
        diag((1./scale(iTones))) * Z(iTones, :));
end

end

