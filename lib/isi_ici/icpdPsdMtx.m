function [S_icpd, S_post, S_pre] = icpdPsdMtx(Hisi, HpreIsi, Ex_arg, Ndft)
% Insuficient Cyclic Prefix Distortion PSD Based on ISI/ICI Matrices
% -------------------------------------------------------------------------
%   [ S_icpd, S_post, S_pre ] = icpdPsd(Hisi, HpreIsi, Ex_arg, N)
%
%   Inputs
% Hisi     -> Post-cursor ISI Matrix
% HpreIsi  -> Pre-cursor ISI Matrix
% Ex_arg   -> Average energy per dimension or input autocorrelation matrix
% Ndft     -> DFT Size
%
%   Notes:
% 1) The ICI matrices could also be passed in place of the ISI matrices.
%    Result is equivalent.
%
% 2) The DFT Size N can be higher than the dimensions of Hisi and HpreIsi
% if a higher DFT resolution is sought. In this case, these matrices must
% be zero-padded.
%
% 3) When the argument "Ex_arg" is a scalar, it is interpreted as the
% average energy per dimension. In this case, it is assumed that the
% transmit signal is uncorrelated, such that its autocorrelation matrix is
% Ex_bar * eye(N). Otherwise (if "Ex_arg" is a matrix), it is interpreted
% as the autocorrelation matrix.
%
% 4) Since the ISI PSD is identical to the ICI PSD, the ISI matrices (pre
% and post cursor) are sufficient for the computation.
%
% 5) PSD scaling factors given zero-padding by L
%
%   A factor of (L = Ndft/N) is used below to scale the PSDs. We can
% understand it by analyzing three possible spectral densities: Watts per
% Hz, watts per normalized frequency and u.e. (unit of energy) per
% dimension. In the first, by increasing N (DFT length) by L, we diminsh
% the resolution in Hertz (=fs/N) by a factor of (1/L) and increase the
% number of density values by a factor of L, so the total energy remains
% roughly the same if the spectral density values are simply interpolated
% between the original DFT points. The same is valid for the case of Watts
% per normalized frequency, the resolution (=1/N) reduces by a factor of L
% and the number of measurements increase by L, so the spectral density
% values remain roughly the same near the original values. In the latter
% case of u.e. per dimension, zero padding by L does not increase the
% number of dimensions of the transmission system, but only the number of
% points computed in the DFT. Each DFT point, therefore, represents 1/L
% dimensions, and since the number of measurements increase by L,
% necessarily the spectral density values must stay the same. In all cases,
% we see that the spectral density must stay the same.
%
% Before proceeding with an example, it is useful to note that the PSD
% computed from mean(|DFT|^2) is related to the density of energy per units
% of normalized frequency (resolution of 1/Ndft). The normalized DFT
% squared is related to the density of energy per tone of the DFT. From
% either of them, in case we want the energy per each of the N degrees of
% freedom, since each DFT tone represents (1/L) dimensions, respectively,
% we need to scale the former computation by (L/Ndft) and the latter by L,
% so that obtain the same total energy is obtained.
%
% Suppose a white noise at 1 u.e./dim is analyzed. Suppose then 8 degrees
% of freedom in the multicarrier system (N = 8 subchannels), so the total
% energy would be 8. The normalized frequency resolution is 1/N = 1/8, so
% mean(|DFT|^2) at each subchannel must be 8. In contrast, the mean of the
% normalized DFT squared should be 1 at each subchannel, such that total
% energy is 8. Now, if we zero pad the sequence to length 16, mean(|DFT|^2)
% continues to be 8, but now we have 16 tones in the DFT and a resolution
% of 1/16 in the normalized frequency, so the total energy continues to be
% 8. Using the normalized DFT, the normalized DFT squared leads to 0.5 in
% each tone. If we want the spectral density in terms of energy per
% dimension, the normalized DFT must be scaled by L = 2, such that with a
% resolution of 1/L dimensions we obtain the same total energy of 8.
%% Initialization

N = size(Hisi,1); % Number of dimensions in the transmission system

if (Ndft < N)
    error('DFT length should be higher than the number of subchannels');
end

if (numel(Ex_arg) > 1)
    if (size(Ex_arg, 1) ~= Ndft)
        error('Input autocorrelation should have dimensions of the DFT');
    end
    Rxx = Ex_arg;
else
    Rxx = Ex_arg * eye(Ndft);
end

% Unitary DFT Matrix
Q = (1/sqrt(Ndft)) * fft(eye(Ndft));

%% Zero-pad the ISI Matrices if necessary

if (Ndft > N)
    nPaddedZeros = Ndft - N;
    Hisi    = padarray(Hisi,nPaddedZeros*[1 1], 'post');
    HpreIsi = padarray(HpreIsi,nPaddedZeros*[1 1], 'post');
end

%% Post-cursor PSD

L = (Ndft/N); % PSD scaling factor

% By explicit writing the formal equation presented in the paper, the
% post-cursor ICPD PSD would be written as:

%   S_post = real(2 * L * diag(Q * Hisi * Rxx * Hisi' * Q'));

% However, it is faster to implement it relying on the FFT algorithm:
tmp = Hisi * Rxx * Hisi';
tmp = fft(tmp, Ndft);
tmp = fft(conj(tmp), Ndft, 2); % Same as " tmp = fft(tmp', Ndft)' "
S_post = real(2 * L * diag((1/Ndft) * tmp));


%% Pre-cursor PSD

% Again, the formal equation presented in the paper is:

% S_pre  = real(2 * L * diag(Q * HpreIsi * Rxx * HpreIsi' * Q'));

% But it is faster to implement it using the FFT algorithm:
tmp = HpreIsi * Rxx * HpreIsi';
tmp = fft(tmp, Ndft);
tmp = fft(conj(tmp), Ndft, 2); % Same as " tmp = fft(tmp', Ndft)' "
S_pre = real(2 * L * diag((1/Ndft) * tmp));


%% Total ICPD PSD:

S_icpd = S_post + S_pre;
end

