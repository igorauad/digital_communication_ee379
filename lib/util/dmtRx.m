function [rx_data] = dmtRx(y, dmt)
% Random DMT Data Generation
% [rx_data] function = dmtRx(dmt)
%
% Inputs:
%  y       -> Received (frame aligned) sequence
%
% Outputs:
%  rx_data -> Bits decoded per subchannel

% Equalizer Types
EQ_FREQ_PREC = 2;
EQ_TIME_PREC = 3;

Nfft     = dmt.Nfft;
N_subch  = dmt.N_subch;
nSymbols = dmt.nSymbols;
nu       = dmt.nu;
tau      = dmt.tau;

% Preallocate
rx_data    = zeros(N_subch, nSymbols);

%% Serial to Parallel

y_sliced = reshape(y, Nfft + nu, nSymbols);

%% Extension removal

y_no_ext = y_sliced(nu + 1:end, :);

%% Equalization
switch (dmt.equalizer)
    case EQ_TIME_PREC % Time-domain ISI DFE
        [ rx_data, Z ] = dmtTdDfeReceiver(y_no_ext, ...
            dmt.modulator, ...
            dmt.demodulator, ...
            dmt.modem_n, ...
            dmt.scale_n, ...
            dmt.Precoder, ...
            dmt.FEQ_n, ...
            dmt.iTonesTwoSided);

    case EQ_FREQ_PREC % DMT with additional modulo operation
        [ rx_data, Z ] = dmtFreqPrecReceiver(y_no_ext,...
            dmt.demodulator, ...
            dmt.modem_n, ...
            dmt.scale_n, ...
            dmt.FEQ_n(1:N_subch), ...
            dmt.b_bar_n, ...
            dmt.dmin_n, ...
            dmt.iTones);

    otherwise
        % FFT
        Y = (1/sqrt(Nfft)) * fft(y_no_ext, Nfft);

        % FEQ - One-tap Frequency Equalizer
        Z = diag(dmt.FEQ_n) * Y(dmt.iTonesTwoSided, :);

        %% Constellation decoding (decision)

        % Iterate over the distinct demodulators
        for iModem = 1:length(dmt.demodulator)
            iSubChs = (dmt.modem_n == iModem); % Loaded subchannels
            % Demodulate
            rx_data(iSubChs, :) = ...
                dmt.demodulator{iModem}.demodulate(...
                diag(1./dmt.scale_n(iSubChs)) * Z(iSubChs, :));
        end
end
end