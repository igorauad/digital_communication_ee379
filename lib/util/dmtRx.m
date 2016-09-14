function [rx_data] = dmtRx(y, dmt)
% Random DMT Data Generation
% [rx_data] function = dmtRx(dmt)
%
% Inputs:
%  y       -> Received (non-aligned) sequence
%
% Outputs:
%  rx_data -> Bits decoded per subchannel

% Equalizer Types
EQ_TEQ       = 1;

% Parameters
Nfft     = dmt.Nfft;
N_subch  = dmt.N_subch;
nSymbols = dmt.nSymbols;
nu       = dmt.nu;
n0       = dmt.n0;

% Preallocate
rx_data    = zeros(N_subch, nSymbols);

%% Time-domain Equalization
switch (dmt.equalizer)
    case EQ_TEQ
        z = conv(dmt.w, y);
    otherwise
        z = y;
end

%% Frame Synchronization
% Note: synchronization introduces a phase shift that should be taken
% into account in the FEQ.

nRxSamples = (Nfft+nu)*nSymbols;
y_sync     = z((n0 + 1):(n0 + nRxSamples));

%% Serial to Parallel

y_sliced = reshape(y_sync, Nfft + nu, nSymbols);

%% Extension removal

y_no_ext = y_sliced(nu + 1:end, :);

%% FFT
Y = (1/sqrt(Nfft)) * fft(y_no_ext, Nfft);

%% FEQ - One-tap Frequency Equalizer
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