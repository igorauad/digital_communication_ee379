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

Nfft     = dmt.Nfft;
nSymbols = dmt.nSymbols;
nu       = dmt.nu;

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

% IFFT
x = sqrt(Nfft) * ifft(X, Nfft);

%% Cyclic extension -> Parallel to serial

x_ext = [x(Nfft-nu+1:Nfft, :); x];
u = x_ext(:);


end