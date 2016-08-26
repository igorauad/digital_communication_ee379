function [u, x] = dmtTx(tx_data, dmt)
% Random DMT Symbol generation
% [u] = dmtTx(tx_data, dmt)
%
% Inputs:
%  tx_data -> Bits allocated per subchannel
%  dmt     -> Struct with DMT parameters
%
% Output:
%  u       -> Extended/windowed modulated sequence
%  x       -> IDFT output (non-extended modulated sequence)

% Equalizer Types
EQ_FREQ_PREC = 2;
EQ_TIME_PREC = 3;

Nfft     = dmt.Nfft;
nSymbols = dmt.nSymbols;
nu       = dmt.nu;
tau      = dmt.tau;

% Erase any previous entries within the symbols
X  = zeros(Nfft, nSymbols);

% Iterate over the distinct modulators
for iModem = 1:length(dmt.modulator)
    M = dmt.modulator{iModem}.M;       % Modulation order
    iSubChs = (dmt.modem_n == iModem); % Loaded subchannels
    % Constellation Encoding
    X(dmt.iTones(iSubChs), :) = diag(dmt.scale_n(iSubChs)) * ...
        dmt.modulator{iModem}.modulate(tx_data(iSubChs, :));
end

% Hermitian symmetry
X(Nfft/2 + 2:Nfft, :) = flipud(conj( X(2:Nfft/2, :)));

%% Per-tone Precoder
if (dmt.equalizer == EQ_FREQ_PREC)
    X = precodeFreqDomain(X, dmt.Precoder, dmt.b_bar_n, ...
        dmt.dmin_n, dmt.iTones);
end

x = sqrt(Nfft) * ifft(X, Nfft);

if (dmt.equalizer == EQ_TIME_PREC)
    x = precodeTimeDomain( x, dmt.Precoder );
end

%% Cyclic extension -> Windowing + overlap -> Parallel to serial
if (dmt.windowing)
    x_ext = [x(Nfft-nu+1:Nfft, :); x; x(1:tau,:)];
    x_ce = windowAndOverlap(x_ext, dmt.window, ...
        Nfft, nu, tau);
    u = x_ce(:);
else
    x_ext = [x(Nfft-nu+1:Nfft, :); x];
    u = x_ext(:);
end

end