%% Example 4.5.1 - Vector Coding
clearvars, clc
addpath(genpath('../lib'))

% Parameters
N       = 8;       % FFT size and the number of subchannels
nu      = 1;       % Prefix
nDim    = N + nu;  % Number of real dimensions
gap_db  = 0;       % SNR gap to capacity (dB)
Ex      = nDim;    % Total energy budget for each VC symbol

% Derivations
gap     = 10^(gap_db/10); % Gap in linear scale

% Noise energy per dimension (N0/2):
sigma_n = 0.181;

% Energy per used dimension (only N are used, although there are N + nu)
Ex_bar = Ex / N


%% Channel

h = [0.9 1];

% Channel frequency response
H = fft(h, N);

% Channel-description matrix
P = convmtx(h, N);


%% Vector Coding
% SVD of the channel matrix leads to the matrices that contain the transmit
% basis vectors and receive "matched filters" for each subchannel.

[F,S,M] = svd(P);

% In the end, the response at each subchannel is given by the singular
% values in "S"

% Note the problem is that the transmiter and receiver need to know the
% channel response.

%% Water filling

gn = (diag(S).^2) / sigma_n;

[bn_bar, En_bar, usedTones] = waterFilling(gn, Ex, N, gap);
% Note the energy per dimension passed to the water-filler is the one
% computed by dividing the total budget by the number of used dimensions
% (N), which excludes the "wasted" dimensions. Likewise, the number of
% dimensions passed is N, not N + nu.

% Number of used tones, according to the water-filling:
N_star = length(usedTones);

%% Transmit energy analys

% In VC, the energy is allocated to N symbols composing vector X of
% dimensions N x 1. This vector then is zero-padded and left-multiplied by
% M from the SVD yielding an (N + nu) x 1 vector with the transmit samples.
% Since M is unitary, the original allocated energy is preserved. In
% contrast, in DMT the energy is allocated to N + nu samples, since nu
% consists necessarily of sample repetitions, which increase the transmit
% energy (think of a non unitary matrix transformation).

% Double check that M is unitary
norm(M)

% Double check the transmit energy
x = M * [ sqrt(En_bar).' ; zeros(nu, 1)];
sum(abs(x).^2)
% Should be equal to
Ex

%% SNR

% SNR per subchannel
SNRn = En_bar(usedTones) .* gn(usedTones).'

% Multi-channel SNR:
SNRm_u = ((prod(1 + SNRn/gap)^(1/(N + nu))) - 1) * gap

fprintf('Multi-channel SNR: \t %g dB \n\n', 10*log10(SNRm_u))

%% Capacity

% Capacity per subchannel, per dimension
cn = 0.5 * log2(1 + SNRn)

% Note: in this case each subchannel has a single dimension (but then all N
% subchannels are considered)

% Multi-channel capacity, per dimension:
c = sum(cn) / nDim

% Note #1: Although the energy distribution per subchannel was computed
% dividing the total energy budget over the number of available subchannels
% (excluding the prefix), the capacity per dimension takes the prefix into
% account in order to allow a fair comparison with other systems. Or, it
% can be interpreted as follows: the extra dimensions from the guard band
% do not require extra energy (or energy allocation) in VC, because they
% are zeroed, but they do occupy some sample periods in transmission, which
% reduces the capacity.

% Note #2: Here the results are given per dimension. The sampling frequency
% and the corresponding transmission bandwith define the number of degrees
% of freedom (dimensions) per second. Since each sample out of the VC
% transmitter implies a single dimension, a sampling frequency fs directly
% maps to dimensions/second. Thus, for example, if fs = 1~MHz, the system
% capacity  ould approach:
%
%       (1.55 bits/dimension)*(1e6 dimensions / second) = 1.55 Mbps.
%
% However, note 1.55 bits/dimension is only valid in the limit for
%   N -> +infty
