function [ rx_data, Z ] = dmtTdDfeReceiver(y, modulator, demodulator, ...
    modChoice, scale, TD, FEQ, usedTones)
% DMT Time-domain Decision-feedback Receiver
% ------------------------------------------------------------------------
% dmtTdDfeReceiver(y, modulator, demodulator, modChoice, scale, TD, FEQ)
%
%   Re-generates ISI through decision-feedback and outputs the decisions.
%
%   Inputs:
% y             - Non-extended time-domain vectors (before FFT and FEQ)
%                 Dimensions: N x nSymbols
% modulator     - Cell array with unique modulators
% demodulator   - Cell array with unique demodulators
% modChoice     - Vector with the modulator choices for each subchannel
% scale         - Scaling factors for each constellation
% TD            - Structure with the Time-domain Precoder/Equalizer design
% FEQ           - FEQ Taps (Hermitian vector)
% usedTones     - Vector indicating the DFT tones that are used
%
%   Outputs:
% rx_data       - Decision output data

% Infer number of symbols and FFT size
nSymbols   = size(y, 2);
N          = size(y, 1);

% Preallocate
rx_data = zeros(N/2 + 1, nSymbols);

for iSymbol = 1:nSymbols
    % ISI equalization, except for the first symbol:
    if (iSymbol > 1)
        y(TD.isi.significantRows, iSymbol) = ...
            y(TD.isi.significantRows, iSymbol) ...
            - TD.isi.Wisi * x_tx(TD.isi.significantCols);
    end

    % FFT
    Y = fft(y(:, iSymbol)) / sqrt(N);

    % FEQ + Detection:
    Z = diag(FEQ) * Y(usedTones, :);

    % Erase "sliced" symbol from previous iteration
    X_tx_p = zeros(N/2 + 1, 1);

    for iModem = 1:length(demodulator)
        iTones = (modChoice == iModem);
        % Demodulate
        rx_data(iTones, iSymbol) = ...
            demodulator{iModem}.demodulate(...
            (1./scale(iTones)) .* Z(iTones));
        % And modulate again to reconstruct the Tx symbol ("sliced"
        % symbol)
        X_tx_p(usedTones(iTones)) = scale(iTones) .* ...
            modulator{iModem}.modulate(...
            rx_data(iTones, iSymbol));
    end

    % Re-generate the original freq-domain Tx symbol:
    X_tx = [X_tx_p; conj(flipud(X_tx_p(2:N/2)))];

    % x_tx is the estimated time domain original Tx symbol,
    % which when multiplied by Wisi yields the ISI due to the
    % precoded symbol. NOTE: Wisi includes the precoding matrix
    x_tx = real(ifft(X_tx) * sqrt(N)); % Orthonormal FFT
end

end

