function [bn, En, SNR_n, Rb, n_loaded] = dmtLoading(dmt, gn)
% Computes the bit-loading and adapt it if necessary

Ex_budget          = dmt.Ex_budget;
N                  = dmt.N;
Nfft               = dmt.Nfft;
nDim               = dmt.nDim;
N_subch            = dmt.N_subch;
equalizer          = dmt.equalizer;
Tsym               = dmt.Tsym;
gap_db             = dmt.gap_db;
max_load           = dmt.max_load;

gap = 10^(gap_db/10); % Gap in linear scale

%% Discrete-loading: Levin Campello Rate Adaptive

% Rate-adaptive Levin-Campello loading:
[En, bn] = DMTLCra(...
    gn(1:N_subch),...
    Ex_budget,...
    N, gap_db, ...
    max_load, ...
    dmt.dim_per_subchannel);

% Residual unallocated energy
fprintf('Unallocated energy:      \t %g\n', Ex_budget - sum(En));

% Energy per real dimension
En_bar = En ./ dmt.dim_per_subchannel;

%% Bit loading computations

% Save a vector with the index of the subchannels that are loaded
n_loaded = dmt.iTones(bn ~= 0);

% Total bits per dimension:
b_bar_discrete = 1/nDim*(sum(bn));

% SNRdmt from the number of bits per dimension
SNRdmt    = gap*(2^(2*b_bar_discrete)-1);
SNRdmt_db = 10*log10(SNRdmt);

% SNR on each tone, per real dimension:
SNR_n     = En_bar .* gn(1:N_subch);

% Bit rate
Rb = sum(bn) / Tsym;

fprintf('b_bar:                    \t %g bits/dimension', b_bar_discrete)
fprintf('\nBit rate:               \t %g mbps\n', Rb/1e6);
fprintf('Multi-channel SNR (SNRdmt): \t %g dB\n', ...
    SNRdmt_db);
end