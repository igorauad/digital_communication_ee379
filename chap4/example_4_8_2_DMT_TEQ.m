%% Example 4.8.2 - Full-band TEQ
clearvars, clc
addpath(genpath('../lib'))

% Parameters
N       = 128;     % FFT size and the number of used real dimensions
nu      = 3;       % Prefix
nDim    = N + nu;  % Total number of real dimensions per DMT symbol
gap_db  = 8.8;     % SNR gap to capacity (dB)
L       = 1;       % Oversampling factor
Ex      = nDim;    % Total energy budget for each DMT symbol
% Note: since repetition of samples in the prefix increases energy,
% differently to VC, here Ex is not directly passed as the water-fill
% energy. Instead, Ex is scaled by a factor of (N/nDim), so that the
% effective transmit energy in the end is equal to Ex.

nSymbols = 100;     % Number of transmit symbols

% Derivations
gap     = 10^(gap_db/10); % Gap in linear scale

% Noise energy per dimension (N0/2):
sigma_n_sq = 0.1;

% Energy per real dimension:
Ex_bar = Ex / nDim;
% Note: only N real dimensions are effectively are used for data, namely
% the DC and Nyquist real dimensions + (N/2)-1 complex dimensions, but all
% N + nu real dimensions are loaded with energy. In another words, since
% repetition of samples in the prefix increases energy, differently to VC,
% Ex_bar here considers the redundant real dimensions, so that the
% effective transmit energy in the end is equal to Ex.


%% Pulse Response

p = [-.729 .81 -.9 2 .9 .81 .729];
% Notes: According to the model developed in Chap4 of EE379, p(t), the
% pulse response, corresponds to the combination between the post-DAC
% filter, the channel impulse response and the pre-ADC filtering.

% Pulse frequency response
H = fft(p, N);

% Channel-description matrix
P = convmtx(p, N);
% P is an N x (N + nu) real matrix

% Pulse response length
Lh = length(p);

% Circulant matrix
P_circ = P(:, 1:N);
P_circ(:,1:nu) = P_circ(:,1:nu) + P(:, N+1:N+nu);

% Matched-filter Bound
SNRmfb = (Ex_bar * norm(p).^2) / sigma_n_sq;
fprintf('SNRmfb:                  \t %g dB\n\n', 10*log10(SNRmfb))

%% Performance with DFE scheme in a single-carrier system
Nf = 14;
Nb = 6;

fprintf('---------------------------------------\n');
fprintf('Single-carrier DFE performance:\n\n');

% The FIR Equalizer can be obtained using the DFE program:
[SNR_mmse_unbiased_db,w_t,opt_delay] = dfsecolorsnr(...
    L,...
    p,...
    Nf,...
    Nb,...
    -1,... % choose the best delay delta
    Ex_bar,...
    (L * sigma_n_sq)*[1; zeros(Nf*L-1,1)]);
% Save the optimum delay

% Expected Performance
SNR_fir_mmse_le_unbiased = 10^(SNR_mmse_unbiased_db/10);
SNR_fir_mmse_le_biased   = SNR_fir_mmse_le_unbiased + 1;
gamma_mmse_le            = 10*log10(SNRmfb / SNR_fir_mmse_le_unbiased);

fprintf('Biased MMSE SNR:\t %g dB\n',...
    10*log10(SNR_fir_mmse_le_biased));
fprintf('Unbiased MMSE SNR:\t %g dB\n',...
    10*log10(SNR_fir_mmse_le_unbiased));
fprintf('MMSE gap to SNRmfb:\t %g dB\n', gamma_mmse_le);
fprintf('Complexity:        \t %d MACs per sample\n', ...
    Nf + Nb);
fprintf('\n');

%% Water filling

fprintf('---------------------------------------\n');
fprintf('Non-equalized DMT:\n\n');

% SNR for unitary energy transmit symbol:
gn = (abs(H).^2) / sigma_n_sq;

[bn_bar, ~, usedTones] = waterFilling(gn, Ex_bar*N, N, gap);
% Recall the budget "Ex" must be reduced by a factor N/(N + nu), to account
% for the fact that repetition in the cyclic prefix will increase the
% energy (in VC there is no repetition, but simply zer-padding). So here
% the SNR is expected to be inferior than the one in VC, which is the price
% paid for complexity reduction.
dim_per_subchannel = [1 2*ones(1, N/2-1) 1 2*ones(1, N/2-1)];

% Number of bits per dimension
b_bar = (1/nDim)*(sum(bn_bar));
fprintf('b_bar:                  \t %g bits/dimension\n', b_bar)
% For gap=0 and N->+infty, this should be the channel capacity in bits per
% real dimension.

% Corresponding SNR:
spectral_efficiency = 2*b_bar; % b/2D (bits per 2 dimensions)
% Minimum SNR to satisfy the Shannon limit:
SNR_min = 2^spectral_efficiency - 1;
% The normalized signal-to-noise ratio (gap to capacity) is defined as:
%   gap = SNR / SNR_min
% Thus, given the gap that was chosen for the water-filling computation,
% the SNR can be computed as (de-normalization):
SNRdmt = 10*log10(gap * SNR_min);

% Number of used tones, according to the water-filling:
N_star = length(usedTones);

fprintf('Multi-channel SNR (SNRdmt):\t %g dB\n', SNRdmt)
fprintf('\nNote:\n');
fprintf('Pe is not guaranteed due to ISI!\n');
fprintf('\n');
%% Now evaluate the TEQ
nTaps      = 11; % Number of taps
filtertype = 1;  % FIR
delta      = 10;

fprintf('---------------------------------------\n');
fprintf('DMT preceded by a TEQ\n\n');

[w,b,SNRteq,bias] = ...
    teq(p, nTaps, nu, delta, sigma_n_sq, Ex_bar, filtertype);

fprintf('New SNRmfb (TEQ):         \t %g dB\n\n', 10*log10(SNRteq))

% Now compute the water-filling solution for the target pulse response
%
% Notes:
%   # 1) The water-filling solution assumes no ISI/ICI. This is safe
%   provided that the TEQ constrains the pulse response energy to a portion
%   that can be covered by the guard band.
%   # 2) Instead of computing the water-fill solution with the N0/2 noise
%   variance per dimension, the error power should be used.
%
unbiased_error_energy_per_dim = ( (norm(p)^2) * Ex_bar) / SNRteq;
% New effective channel:
H_teq = fft(b, N);
% New Unitary-energy SNR:
gn_teq = (abs(H_teq).^2) / unbiased_error_energy_per_dim;
% Water-filling:
[bn_bar_teq, ~, usedTones] = waterFilling(gn_teq, Ex_bar*N, N, gap);

% Number of bits per dimension
b_bar_teq = (1/nDim)*(sum(bn_bar_teq));
fprintf('b_bar (TEQ):               \t %g bits/dimension\n', b_bar_teq)

spectral_efficiency_teq = 2*b_bar_teq; % b/2D (bits per 2 dimensions)
SNR_min_teq = 2^spectral_efficiency_teq - 1;
SNRteq = 10*log10(gap * SNR_min_teq);
fprintf('Multi-channel SNR (SNRdmt_teq):\t %g dB\n', SNRteq)
fprintf('Complexity:               \t %g MACs per sample\n', nTaps + ...
(N/nDim) * log2(N));

fprintf('\n');
%% Finally, compare the performance of a DMT with N = 1024 and no TEQ
% Note the way to increase performance by increasing N is to keep the
% sampling frequency constant, so that the tone spacing is reduced. By
% doing so, the transmit energy, being Px*Tsym, has to increase, because
% Tsym = (N + nu) * Ts and N increases. Furthermore, this guarantees that
% no other cyclic prefix is required if it is already chosen sufficiently,
% because the channel is sampled by the same frequency. In our case, we
% will slightly increase the CP to make it sufficient.

fprintf('---------------------------------------\n');
fprintf('DMT without TEQ, but with N=1024\n\n');
% Previous power, assuming Ts = 1 for simplicity
Px = Ex/nDim;
% New dimensions
N    = 1024;
nu   = 6;       % Sufficient cyclic prefix
nDim = N + nu;
% Now, let's maintain the same power and ajust the energy for the new Tsym
Ex     = Px * nDim;
Ex_bar = Ex / nDim;

% Pulse frequency response
H = fft(p, N);

% SNR for unitary energy transmit symbol:
gn = (abs(H).^2) / sigma_n_sq;

% Water-filling
[bn_bar, ~, usedTones] = waterFilling(gn, Ex_bar*N, N, gap);

% Bits per dimension
b_bar = (1/nDim)*(sum(bn_bar));
fprintf('b_bar (TEQ):               \t %g bits/dimension\n', b_bar)

% Spectral efficiency:
spectral_efficiency = 2 * b_bar; % b/2D (bits per 2 dimensions)

% SNR
SNR_min = 2^spectral_efficiency - 1;
SNRdmt = 10*log10(gap * SNR_min);

fprintf('Multi-channel SNR (SNRdmt):\t %g dB\n', SNRdmt)
fprintf('Complexity:               \t %g MACs per sample\n', ...
    (N/nDim) * log2(N));
fprintf('\nNote:\n');
fprintf('nu = %d now. There is no ISI and, thus, Pe is guaranteed!\n', ...
    nu);