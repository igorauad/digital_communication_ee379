function [ rx_sym, Z ] = dmtFreqPrecReceiver(y, demodulator, modChoice, scale, FEQ, M, D)
% DMT Receiver for the Frequency-domain Precoded Transmissions
% ------------------------------------------------------------------------
% dmtFreqPrecReceiver(y, demodulator, modChoice, scale, FEQ)
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
%
%   Outputs:
% rx_sym        - Decision output data

% Infer number of symbols and FFT size
nSymbols   = size(y, 2);
N          = size(y, 1);

% Preallocate
rx_sym = zeros(N/2 + 1, nSymbols);

% FFT
Y = (1/sqrt(N)) * fft(y, N);

% FEQ - One-tap Frequency Equalizer
Z = diag(FEQ) * Y;

% Modulo operation
for iSym=1:nSymbols
    Z(:,iSym) = moduloOperator(Z(1:N/2+1, iSym), ...
        M, D, 'Hermitian');
end

% Constellation decoding (decision)
for iModem = 1:length(demodulator)
    iTones = (modChoice == iModem);
    % Demodulate
    rx_sym(iTones, :) = ...
        demodulator{iModem}.demodulate(...
        diag((1./scale(iTones))) * Z(iTones, :));
end

end

