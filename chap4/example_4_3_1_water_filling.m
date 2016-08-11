%% Example 4.3.1 - Water filling
clearvars, clc
addpath(genpath('../lib'))

% Parameters
N       = 8;    % FFT size and the number of dimensions
gap_db  = 0;    % SNR gap to capacity (dB)
Ex      = N;    % Total energy budget for each OFDM symbol


% Derivations
N_bar   = N/2 + 1;        % Used subchannels (including DC and Nyquist)
Ex_bar  = Ex / N;         % Energy per dimension
gap     = 10^(gap_db/10); % Gap in linear scale

% Noise energy per dimension (N0/2):
sigma_n = 0.181;

%% Channel

h = [0.9 1];

H = fft(h, N);

figure
plot(10*log10(abs(H).^2))
xlabel('Subcarrier (n)')
title('Channel Response')
grid on

%% SNRn when the transmitter applies unit energy to the subchannel

gn = abs(H).^2 / sigma_n;

%% Compute bit-loading using waterfilling
% This distributes the energy budget Ex over the subchannels with the
% objective of maximizing the data-rate. As a result, a non-integer number
% of bits is obtained as well as the corresponding energy allocated to each
% dimension of the subchannels.

[bn_bar, En_bar] = waterFilling(gn, Ex, N, gap)

% Sanity check. Is the sum of the energy in each dimension less than the
% total?
fprintf('Total allocated energy: \t %g \n', sum(En_bar));

% Ideally they should be equal, because this is the constraint in the
% optimization.

% Total number of bits per dimension
b_bar = sum(bn_bar) / N;

figure
stem(bn_bar(1:N_bar), 'linewidth', 1.2)
xlabel('Subchannel n')
ylabel('$b_n$ (dB)')
grid on

%% SNR per subchannel

SNRn = En_bar .* gn;

figure
stem(10*log10(SNRn))
xlabel('Subchannel n')
ylabel('$SNR_n$ (dB)')
grid on

%% Resulting multi-channel SNR
% "a single SNR measure that characterizes the set of subchannels by an
% equivalent single AWGN that achieves the same data rate"

SNRm_u = ((prod(1 + SNRn/gap)^(1/N)) - 1) * gap

% Alternative way of computing the number of bits per dimension:
b_bar_2 = 0.5 * log2(1 + SNRm_u/gap)

%% Margin

% First, note the difference between the gap and the margin definitions:
%
%   The gap is obtained by noting that any spectral efficiency should be
%   lower than the channel capacity, i.e. rho < C = log2(1+SNR). The
%   inequality can be stated as:
%
%       SNR > 2^rho - 1
%
%   Which leads to the conclusion that "2^rho - 1" is a lower bound on the
%   SNR.

% Channel Capacity
C  = log2(1 + SNRm_u); % Bits per 2 dimensions (b/2D)
% The acual spectral efficiency must be lower than C:
spec_effiency =  2*b_bar; % Bits per 2 dimensions (b/2D)
% SNR lower bound:
SNR_minimum = 2^spec_effiency - 1;
% Gap to capacity:
gap_computed = SNRm_u / SNR_minimum;
% See Chap. 4 of Forney's course on Principles of Digital Communications II
% at MIT OCW

% However, sometimes the spectral effiency is arbitrarily reduced such that
% more margin is left for the errors in the communication. By recalling
% that spectral efficiency is R / W, one concludes normally this implies a
% lower data rate.
%   The constrained spectral effiency is given by:
%       2 * b_bar_constrained = log2(1 + SNR/(gap * gamma))
%   Hence, gamma (the margin) can be obtained by:
%       gamma = (SNR/gap) / (2^(2*b_bar_constrained) - 1)
%   Finally, note (SNR/gap) is the SNR lower bound of Shannon, which is
%   equivalent to 2^(2*b_bar_unconstrained) - 1.

% Definition in (4.7):
gamma_m = (SNRm_u/gap) / (2^(2*b_bar) - 1);

% In dB:
gamma_m_db = 10*log10(gamma_m)

% Note: margin is expected to be very low, since we are maximizing
% data-rate at its expense.