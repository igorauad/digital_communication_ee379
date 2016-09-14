%% Example 4.6.1 - DMT
% By introducing the cyclic prefix, a DMT transmitter is similar to a VC
% transmitter, but with the SVD replaced by the eigendecomposition and with
% the decomposition unitary matrices being related to the DFT/IDFT matrix,
% which avoids the need for knowing the channel impulse response in order
% to design transmit/receive basis functions.

clearvars, clc
addpath(genpath('../lib'))

% Parameters
N       = 8;       % FFT size and the number of subchannels
nu      = 1;       % Cyclic Prefix Length
nDim    = N + nu;  % Number of real dimensions
gap_db  = 0;       % SNR gap to capacity (dB)
Ex      = nDim;    % Total energy budget for each DMT symbol
% Note the goal is to have Ex_bar = 1 for comparison with other invocations
% of the example.

% Derivations
gap     = 10^(gap_db/10); % Gap in linear scale

% Noise energy per dimension (N0/2):
sigma_n = 0.181;

% Energy per real dimension:
Ex_bar = Ex / nDim
% Note: only N real dimensions are effectively are used for data, but all N
% + nu real dimensions are loaded with energy. In another words, since
% repetition of samples in the prefix increases energy, differently to VC,
% Ex_bar here considers the redundant real dimensions, so that the
% effective transmit energy in the end is equal to Ex.

%% Channel

h = [0.9 1];

% Channel frequency response
H = fft(h, N);

% Channel-description matrix
P = convmtx(h, N);
% P is an N x (N + nu) real matrix

% Channel impulse response length:
Lh = length(h);

% Circulant matrix
P_circ = P(:, 1:N);
P_circ(:,1:nu) = P_circ(:,1:nu) + P(:, N+1:N+nu);

SNRmfb = (Ex_bar * norm(h).^2)/sigma_n

%% Eigendecomposition, DFT and Channel Response comparison

[M,S] = eig(P_circ);

% By comparing M to the DFT matrix or S to the channel response, it is
% possible to check that they contain the same elements. The difference is
% that the eigendecomposition sorts differently.

%% Water filling

% SNR for unitary energy transmit symbol:
gn = (abs(H).^2) / sigma_n;

[bn_bar, En_bar, usedTones] = waterFilling(gn, Ex_bar*N, N, gap)
dim_per_subchannel = [1 2*ones(1, N/2-1) 1 2*ones(1, N/2-1)]

% Number of bits per dimension
b_bar = (1/nDim)*(sum(bn_bar));
% Corresponding SNR:
SNRdmt = 10*log10(gap*(2^(2*b_bar)-1));
% Note, there is no need to account for the number of dimensions in each
% subchannel. The Hermitian symmetry in the bn_bar vector already takes
% into account a complex subchannel twice.

% Similarly to the case of VC, the energy per dimension passed to the
% water-filler is the one computed dividing the total budget by the number
% of used dimensions (N), which excludes the "wasted" dimensions. Likewise,
% the number of dimensions passed is N, not N + nu, because the extra "nu"
% real dimensions cannot be used to carry information. The difference here
% with respect to VC is that the budget "Ex" must be reduced by a factor
% N/(N + nu), to account for the fact that repetition in the cyclic prefix
% will increase the energy (in VC there is no repetition, but simply
% zer-padding). So here the SNR is expected to be inferior than the one in
% VC, which is the price paid for complexity reduction.

% Number of used tones, according to the water-filling:
N_star = length(usedTones);

%% SNR

% SNR per subchannel
SNRn = En_bar(usedTones) .* gn(usedTones)

% Multi-channel SNR:
SNRm_u = ((prod(1 + SNRn/gap)^(1/(nDim))) - 1) * gap

fprintf('Multi-channel SNR: \t %g dB\n\n', 10*log10(SNRm_u))
% Note the SNR is lower than the one obtained for VC

%% Capacity

% Capacity per subchannel, per dimension
cn = 0.5 * log2(1 + SNRn);

% Note: in this case each subchannel has a single dimension

% Multi-channel capacity, per dimension:
c = sum(cn) / nDim

%% Discrete-loading
% Levin Campello Rate Adaptive

[En, bn] = DMTLCra(gn(1:N/2 + 1), Ex_bar*N, N, gap)

% calculate b_bar (bits per dimension)
b_bar=1/nDim*(sum(bn));
% SNRdmt from the number of bits per dimension
SNRdmt=10*log10(gap*(2^(2*b_bar)-1));

fprintf('\nDiscrete Loading:\n\n');
fprintf('Multi-channel SNR: \t %g dB\n\n', 10*log10(SNRm_u))
% Note the SNR is lower than the one obtained for VC
