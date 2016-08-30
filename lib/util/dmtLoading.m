function [bn, En, SNR_n, Rb, n_loaded] = dmtLoading(dmt, gn)
% Computes the bit-loading and adapt it if necessary

% Equalizer Types
EQ_FREQ_PREC = 2;
EQ_TIME_PREC = 3;

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

%% Loading adaptation to avoid power penalties in full ICPD mitigation
% In case the full ICPD equalizers are used, either the frequency-domain or
% the time-domain precoder, power increase can occur. This power penaly
% shall be pre-compensated in the energy budget that is passed to the bit
% loader. The main difficulty in this process, however, is that the power
% increase itself depends on the energy load. Our approach is to first
% compute an initial energy load (the one from previous section), assuming
% initially no power increase due to precoding. Then, based on the computed
% energy load, the total energy after precoding is computed and compared to
% the original budget. The reciprocal of the factor by which the energy
% increases due to precoding is used to reduce the budget. In the end, the
% bit/energy load is re-computed.

if (equalizer == EQ_TIME_PREC || equalizer == EQ_FREQ_PREC)

    fprintf('\n------ Energy-load adaptation for the ICPD Precoder -----\n\n');

    % Full Hermitian En_bar vector for the Levin-Campello energy load
    En_bar_herm = zeros(Nfft, 1);
    En_bar_herm(dmt.iTonesTwoSided(1:N_subch)) = En_bar;
    En_bar_herm(dmt.iTonesTwoSided(N_subch+1:end)) = ...
        fliplr(En_bar);

    % Average transmit energy per symbol after precoding
    if (equalizer == EQ_TIME_PREC)
        Ex_precoded = real(trace(dmt.Precoder.ici.W * ...
            diag(En_bar_herm) * dmt.Precoder.ici.W'));
    else
        Ex_precoded = real(trace(dmt.Precoder.W * ...
            diag(En_bar_herm) * dmt.Precoder.W'));
    end

    % By how much the precoded energy exceeds the budget:
    Ex_budget_excess = Ex_precoded / Ex_budget;

    % Reduce the budget by the following amount
    budget_red_factor = 1/Ex_budget_excess;

    % Rate-adaptive Levin-Campello loading:
    [En, bn] = DMTLCra(...
        gn(1:N_subch),...
        budget_red_factor * Ex_budget,...
        N, gap_db, ...
        max_load, ...
        dmt.dim_per_subchannel);

    % Residual unallocated energy
    fprintf('Unallocated energy:      \t %g\n', Ex_budget - sum(En));

    % Energy per real dimension
    En_bar = En ./ dmt.dim_per_subchannel;

    fprintf('Energy budget was reduced by %.2f %%\n', ...
        100*(1 - budget_red_factor));
end

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