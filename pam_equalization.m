% PAM Equalization
%
%   Author: Igor Freire
%
%   Based on the material for Stanford EE 379A - Digital Communication:
%   Signal Processing, by John M. Cioffi
%
% The models used in this script are mostly based on Figure 3.7 and the
% conventions of Table 3.1 are used throughout.
%
% Tips to understand conversion factors of T (Table 3.1):
%   - If a function has unitary energy in continuous time, then its
%   discrete-time equivalent has norm^2 of 1/Ts (inverse of the sampling
%   period.
%

clearvars, close all;
clc

% Parameters
N           =   1;      % Dimensions per symbol
nBits       =   4e5;    % Number of transmit symbols
debug       =   0;      % Enable plots and extra information
rollOff     =   0.1;    % Roll-off factor
L           =   4;      % Oversampling (support only for integer values)
W           =   5e3;    % Nominal bandwith (Hz)
N_T         =   10;     % Raised-cosine group-delay in symbols
Px          =   1e-3;   % Transmit Power (W)
N0_over_2   =   1e-13;  % Noise PSD (W/Hz/dim) and variance per dimension
M           =   16;
ideal_chan  =   0;
en_noise    =   1;
equalizer   =   2;      % 0) no equalizer; 1) FIR ZF-LE; 2) FIR MMSE-LE
% FIR ZF-LE
% FIR MMSE-LE
Nf = 10;

% Derived computations:
b        = log2(M);         % Bits per symbol
Rsym     = 2 * W;           % Symbol rate (real degrees of freedom per sec)
Tsym     = 1 / Rsym;        % Symbol Period
Fs       = Rsym * L;        % Symbol Rate
Ts       = 1 / Fs;          % Sampling period
nSymbols = ceil(nBits / b); % Number of Tx symbols


fprintf('Baud rate:\t%g symbols/s\n', Rsym);
fprintf('Data rate:\t%g kb/s\n', Rsym * b / 1e3);
fprintf('Fs:       \t%g Hz\n', Fs);

if (~en_noise)
    fprintf('\nBackground Noise Disabled\n');
end
%% Constellation and Energy/Power computations

Ex     = Px * Tsym; % Average energy of a constellation
Ex_bar = Ex / N;    % Energy per dimension

% Scale factor for the PAM constellation to present average energy of "Ex":
Scale = modnorm(pammod(0:(M-1), M), 'avpow', Px);
% Why Px, instead of Ex? Because the energy within symbol-spaced samples of
% the transmit signal is actually Tsym * E{ |x_k|^2 }, due to the fact that
% the sincs in the sampling theorem sinc-interpolation formula are not
% orthonormal, but only orthogonal. The former implies the samples must be
% scaled by sqrt(Tsym), so a factor of Tsym appears multiplying the norm.
% Since
%
%   Tsym * E{ |x_k|^2 } = Ex                                        (1)
%
% it follows that:
%
%   E{ |x_k|^2 } = Ex/Tsym = Px                                     (2)
%
% Q.E.D
%


% Noise energy per dimensions
%
%   It depends on the Rx filter, which is different for each receiver (a
%   receiver in this script is characterized by a given equalizer).
%
% - Matched Filter
%   For the matched filter receiver, the noise energy per dimension remains
% N0/2 regardless of oversampling.
%
% - MMSE Receiver:
%
%   For the MMSE receiver, in the presence of oversampling, the Rx filter
% is assumed to be a "brick-wall" filter of bandwidth "L" times larger than
% the nominal bandwidth, but with the same conventional magnitude sqrt(T)
% that preserves the spectrum within the central period.
%   Thus, the analog filter energy becomes L, rather than unitary, so that
% the noise energy per dimension becomes (N0/2) * L. In contrast, the Tx
% signal energy is assumed to be contained within -1/T to 1/T and, thus,
% does not change. As a result, the SNRmfb is reduced by a factor of L.
switch (equalizer)
    case 2
        noise_en_per_dim = L * N0_over_2;
    otherwise
        noise_en_per_dim = N0_over_2;
end

%% Generate random symbols

tx_decSymbols = randi(M, nSymbols, 1) - 1; % Decimals from 0 to M-1
tx_binSymbols = de2bi(tx_decSymbols);      % Corresponding binaries
tx_bitStream  = tx_binSymbols(:);          % Bitstream

%% Modulation

unscaled_signals = real(pammod(tx_decSymbols, M));
tx_signals = Scale * unscaled_signals;

%% Tx pulse shaping filter

if (L > 1)
    % Apply a square-root raised cosine pulse shaping filter:
    htx    = rcosine(1, L, 'sqrt', rollOff, N_T);
    % Energy of the continuous-time transmit basis function:
    E_htx  = Ts * sum(abs(htx).^2);
    % Normalize for unitary energy (in continuous-time):
    htx    = htx * (1/sqrt(E_htx));
    % Ts * sum(abs(htx).^2) now is unitary
else
    % Without oversampling, pulse shaping (other than the T-spaced sinc)
    % can not be applied.
    htx = 1/sqrt(Tsym);
    % Note: this is the same as sampling (1/sqrt(T))*sinct(t/T) at t=kT.
    % All samples, except the one for k=0, are zero.
    % Note "Ts * sum(abs(htx).^2)" is unitary
end

% Filter response
if (debug)
    figure
    freqz(htx)
    title('Tx Filter')
end

%% Baseband channel response

if (ideal_chan)
    h = 1;
else
    h = [0.9 1];
end
% Note: the energy is not very important here, because there is not too
% much to do about the channel attenuation/gain anyway.

%% Pulse response

p = Ts * conv(h, htx);
% Note: p (therefore h and htx) are sampled with Ts = Tsym/L (oversampled)

% Pulse norm:
norm_p_sq = Ts * norm(p)^2;
norm_p = sqrt(norm_p_sq);
% Note: SNR_{MFB} has to be found using p(t) before anti-aliasing filter.

% Unitary-energy (in continuous-time) pulse response:
phi_p = p / norm_p;

fprintf('\n--------- MFB ---------\n');
SNRmfb = Ex_bar * norm_p_sq / noise_en_per_dim;
fprintf('\nSNRmfb:   \t %g dB\n', 10*log10(SNRmfb))
% NNUB based on the SNRmfb
Pe = 2 * (1 - 1/M) * qfunc(sqrt(3*SNRmfb / (M^2 - 1)));
fprintf('Pe (NNUB):\t %g\n', Pe);


%% Waveform generation - upsample and filter

signals_up          = zeros(1,nSymbols*L);
signals_up(1:L:end) = tx_signals;

% Shaped waveform:
tx_waveform = Ts * conv(htx, signals_up(:));

if (debug)
   % To understand the following, consult page 26, chap 9 of Gallager's
   % book on Digital Comm I.
   fprintf('\n--- Energy/Power Measurements ---\n');
   % Due to the invariance of the inner product, the average transmit
   % energy (given the basis are orthonormal) should be close to Ex:
   tx_avg_energy = Tsym * mean(abs(tx_signals).^2);
   fprintf('Measured average Tx energy:\t %g\n', tx_avg_energy);
   fprintf('Spec average Tx energy (Ex):\t %g\n', Ex);
   % Upsampled sequence average energy
   tx_avg_energy_sampled = Ts * mean(abs(signals_up).^2);
   fprintf('Average sample energy (Es):\t %g\n', tx_avg_energy_sampled);
   fprintf('Observe Es = Ex/L\n');
   fprintf('--\n');
   % The transmit power is equivalent to the mean in the transmit signal
   % sequence.
   Ex_over_Tsym = tx_avg_energy / Tsym;
   Es_over_Ts = tx_avg_energy_sampled / Ts;
   fprintf('Ex/Tsym:\t %g\n', Ex_over_Tsym);
   fprintf('Es/Ts:  \t %g\n', Es_over_Ts);
   fprintf('Spec Px:\t %g\n', Px);
end

%% Transmission through channel

% Receive signal past channel, but pre noise:
rx_pre_noise = Ts * conv(h, tx_waveform);

% AWGN:
noise = sqrt(N0_over_2) * randn(size(rx_pre_noise));
% There is a related note pointed at the document below about the
% factor of Ts that scales the noise variance.
%
% http://web.stanford.edu/group/cioffi/ee379a/extra/Hoi_problem3.34.zip
%
% Another related note can be found in Robert Gallager's material for the
% Principles of Digital Communications I course, Chaper 9, footnote 25.
% Rephrased slightly (for compatibility with the nomenclature in this
% script), it can be understood as follows:
%   The sinc functions in the orthogonal expansion (sampling theorem) have
%   energy Ts, so the variance of each real and imaginary coefficient in
%   the noise expansion must be scaled down by (Ts) from the noise energy
%   N0/2 per degree of freedom to compensate the sinc scaling. In the end,
%   iid variables with N0/2 variance are obtained. TO-DO: how to adapt the
%   above?
if (debug)
    fprintf('\n--- Noise Power Measurements ---\n');
    % Compare with Mathworks results
    AWGN = comm.AWGNChannel;
    AWGN.NoiseMethod = 'Signal to noise ratio (Es/No)';
    AWGN.EsNo = 10*log10((Ex/(2*N0_over_2)));
    AWGN.SignalPower = Px;
    AWGN.SamplesPerSymbol = L;
    noise_mtwks = AWGN.step(zeros(size(rx_pre_noise))); % AWGN channel
    fprintf('Nominal N0/2:\t %g\n', N0_over_2);
    fprintf('Measured noise variance per real dim:\t %g\n', ...
        Ts * var(noise));
    fprintf('Mathworks noise variance per real dim:\t %g\n', ...
        Ts * var(noise_mtwks));
end

% Rx waveform
if (en_noise)
    rx_waveform = rx_pre_noise + noise;
else
    rx_waveform = rx_pre_noise;
end

if (debug)
   pwelch(tx_waveform,[],[],[],Fs,'twosided');
   figure
   pwelch(noise,[],[],[],Fs,'twosided');
end

%% Equalization

switch (equalizer)
    case 1
        error('FIR ZF is not implemented yet');

    case 2
        fprintf('\n------- FIR MMSE -------\n');
        % First, MMSE is fractionally spaced and incorporates both
        % matched filtering and equalization.
        % Secondly, it is preceded by an anti-alias filter whose gain is
        % sqrt(T) from -l/T to l/T (a support of 1/Ts). The filter in
        % time-domain is given by:
        %
        %   (sqrt(Tsym)/Ts) * sinc(t/Ts) = (L/sqrt(Tsym)) * sinc(t/Ts)
        %
        % When sampled at t = kTs, since sinc(t/Ts) is 1 for k=0 and 0
        % elsewhere, the discrete sequence is [ L/sqrt(Tsym) ]
        hrx = L/sqrt(Tsym);
        % Note: Ts * norm(hrx).^2 = L (filter energy is L).

        % Combined pulse response + anti-aliasing filter:
        p_tilde = Ts * conv(p, hrx);

        % Anti-aliasing receive filtering:
        rx_waveform = Ts * conv(rx_waveform, hrx);

        % First split p into blocks (columns) of L samples each. This is
        % the effect of a clockwise commutator in the L-branch
        % decomposition.

        % Number of blocks:
        n_blocks = 1 + ceil((length(p_tilde)-1)/L);
        % Last block may require zero-padding
        if (mod(length(p_tilde)-1, L) ~= 0)
            n_missing_samples = L - mod(length(p_tilde)-1, L);
        else
            n_missing_samples = 0;
        end

        % Zero-pad and split:
        p_0 = [p_tilde(1); zeros(L-1,1)];
        p_split = [p_0 flipud(reshape([p_tilde(2:end) ...
            zeros(1,n_missing_samples)], L, n_blocks-1))];
        % Note: flipud is used to distribute the p_tilde samples as a CW
        % commutator over the L branches

        nu = size(p_split,2) - 1;   % Pulse response dispersion
        delta = round((Nf + nu)/2); % Equalized system delay

        nRows = Nf * L;
        nCol = Nf + nu;
        % Preallocate
        P = zeros(nRows, nCol);

        nZerosRight = nCol - (nu + 1);
        nZerosLeft = 0;

        nBlocks = Nf;

        for iBlock = 1:nBlocks
            P((iBlock-1)*L + 1 : (iBlock*L),:) = ...
                [zeros(L, nZerosLeft) p_split zeros(L, nZerosRight)];
            nZerosLeft = nZerosLeft + 1;
            nZerosRight = nZerosRight - 1;
        end

        R_Yx = Ex_bar * P * [zeros(delta,1); 1; zeros(Nf + nu - delta - 1, 1)];
        R_YY = (Ex_bar * (P * P')) + (noise_en_per_dim * eye(Nf*L));
        % MMSE-LE FIR Equalizer:
        w = (R_YY\R_Yx)';
        % Alternatively, the FIR Equalizer can be obtained using the
        % DFE program:
        if (debug)
            fprintf('\n-- EE379A DFE Implementation --\n');
            [SNR_mmse,w_t,opt_delay]=dfsecolorsnr(L,p_tilde,Nf,0,delta,...
                Ex,noise_en_per_dim*[1; zeros(Nf*L-1,1)]);
            figure
            plot(w, '--', 'linewidth', 1.2)
            hold on
            plot(w_t, 'r')
            title('FIR MMSE-LE');
            legend('Our', 'EE379A Program');
            fprintf('Unbiased MMSE SNR:\t %g dB\n',SNR_mmse);
        end

        % Same as (inv(R_YY)*R_Yx)';
        z = conv(w, rx_waveform);
        % Note: the Ts factor is not necessary here, since there is not
        % turning back to analog domain again past this point

        % Skip MMSE filter delay and Acquire a window with nSymbols * L
        % samples. Again, recall nu and Nf are given in terms of T-spaced
        % symbols, not samples, so multiplication by L is required.
        z = z( (delta*L + 1) : (delta + nSymbols)*L );
        % Down-sample
        z_k = z(1:L:nSymbols*L).';
        % Normalize:
        z_k = z_k / (sqrt(L) * norm(w));
        % Note: w filter is fractionally spaced, so its "downsampled"
        % version will present a norm scaled by sqrt(L).

        % Performance
        error_power = Ex_bar - w * R_Yx;
        SNR_fir_mmse_le_biased = Ex_bar / error_power;
        fprintf('Biased MMSE SNR:\t %g dB\n',...
            10*log10(SNR_fir_mmse_le_biased));
        SNR_fir_mmse_le_unbiased = SNR_fir_mmse_le_biased - 1;
        fprintf('Unbiased MMSE SNR:\t %g dB\n',...
            10*log10(SNR_fir_mmse_le_unbiased));
        % NNUB based on the SNR_fir_mmse_le_unbiased
        Pe = 2 * (1 - 1/M) * qfunc(sqrt(3*SNR_fir_mmse_le_unbiased / (M^2 - 1)));
        fprintf('Pe (NNUB):      \t %g\n', Pe);
        gamma_mmse_le = 10*log10(SNRmfb / SNR_fir_mmse_le_unbiased);
        fprintf('MMSE gap to SNRmfb:\t %g dB\n', gamma_mmse_le);

    otherwise
        % Matched filter receiver
        hrx  = conj(fliplr(htx));
        [~, n0] = max(conv(p, hrx));
        % Matched filtering
        y = Ts * conv(rx_waveform, hrx);
        % Symbol timing synchronization
        % Acquire a window with nSymbols * L samples
        y_s = y(n0:n0 + nSymbols*L - 1);
        % Followed by downsampling (when L > 1)
        % T-spaced received symbols:
        y_k = y_s(1:L:end);
        % Equalized are equal to rx symbols (i.e. no equalization)
        z_k = y_k;
end

% Receiver Gain control (Assumed to be known)
z_k_unscaled = z_k / (Scale * norm_p * Ts);

% if (debug)
    figure
    plot(Tsym*(0:nSymbols-1), unscaled_signals, 'o')
    xlabel('Tempo (s)')
    ylabel('Amplitude')
    grid on
    hold on
    plot(Tsym*(0:nSymbols-1), z_k_unscaled, 'ro')
    legend('Original','Equalized')
% end

%% Decision
rx_decSymbols = pamdemod(z_k_unscaled, M);
% Filter NaN
rx_decSymbols(isnan(rx_decSymbols)) = 0;
rx_binSymbols = de2bi(rx_decSymbols);
rx_bitstream  = rx_binSymbols(:);

%% Symbol error
fprintf('\n----- Performance -----\n');

[nSymErrors, SER, symErrArray] = symerr(tx_decSymbols, rx_decSymbols(:));

fprintf('\nSymbol errors:\t %g\n',nSymErrors);
fprintf('SER:     \t %g\n', SER);

%% Bit error

[nBitErrors, BER, errArray] = biterr(tx_bitStream, rx_bitstream);

fprintf('\nBit errors:\t %g\n',nBitErrors);
fprintf('BER:      \t %g\n', BER);