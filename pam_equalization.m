clearvars, close all;
clc

% Parameters
N           =   1;      % Dimensions per symbol
nBits       =   4e6;    % Number of transmit symbols
debug       =   0;      % Enable plots and extra information
rollOff     =   0.2;    % Roll-off factor
L           =   4;      % Oversampling (support only for integer values)
Fs          =   1e3;
N_T         =   10;     % Raised-cosine group-delay in symbols
Px          =   1e-3;   % Transmit Power (W)
N0_over_2   =   1e-9;   % Noise PSD (W/Hz) and variance per dimension
M           =   4;
ideal_chan  =   0;
en_noise    =   1;
equalizer   =   2;      % 0) no equalizer; 1) FIR ZF-LE; 2) FIR MMSE-LE
% FIR ZF-LE
% FIR MMSE-LE
Nf = 10;

% Derived computations:
b        = log2(M);         % Bits per symbol
Ts       = 1 / Fs;          % Sampling period
Rsym     = Fs / L;          % Symbol Rate
Tsym     = 1 / Rsym;        % Symbol Period
nSymbols = ceil(nBits / b); % Number of Tx symbols
Wb       = Rsym/2;          % Nominal bandwith (half the symbol rate)
Wb_norm  = Wb / (Fs/2);     % Normalized nominal bandwidth

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
Scale  = sqrt(Ex) * modnorm(pammod(0:(M-1), M),'avpow',1);

% Noise

noise_en_per_dim = L * N0_over_2;
% Remember that in the presence of oversampling, the receiver analog filter 
% is assumed to be a "brick-wall" filter of bandwidth "l" times larger than 
% the nominal bandwidth, but with the same conventional magnitude sqrt(T).
% Thus, the analog filter energy becomes l, rather than unitary. Then, the
% noise energy per dimension becomes (N0/2) * l.

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
    E_htx  = sum(abs(htx).^2);      % pulse energy
    htx    = htx * (1/sqrt(E_htx)); % normalize for unitary energy
else
    % Without oversampling, pulse shaping (other than the T-spaced sinc)
    % can not be applied.
    htx = 1;
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

%% Pulse response

p = conv(h, htx);
% Note: p (therefore h and htx) are sampled with L*Rsym (oversampled)

% Combined response with the gain of the anti-alias receive filter
p_tilde = sqrt(Tsym) * p;
% It is the convolution between the pulse response and the anti-alias 
% filter.

% Pulse norm:
norm_p_sq = norm(p_tilde)^2;
norm_p = norm(p_tilde);
% Note: SNR_{MFB} has to be found using p(t) before anti-aliasing filter,
% or using p_tilde with a slightly different formula (the one above).

% Unitary-energy pulse response
phi_p = p / norm_p;

fprintf('\n--------- MFB ---------\n');
SNRmfb = Ex_bar * norm_p_sq / noise_en_per_dim;
fprintf('\nSNRmfb:   \t %g dB\n', 10*log10(SNRmfb))
% NNUB based on the SNRmfb
Pe = 2 * (1 - 1/M) * qfunc(sqrt(3*SNRmfb / (M^2 - 1)));
fprintf('Pe (NNUB):\t %g\n', Pe);

%% Matched Rx filter

hrx  = conj(fliplr(phi_p));     % matched filter

[~, n0] = max(conv(p, hrx));

%% Combined Response
% The response that satisfies Nyquist Criterion
%
% $$q(t) = \text{sinc}\left(\frac{t}{T}\right)$$

q = conv(phi_p, hrx);

if (debug)
    N_q = length(q);
    t = [-(N_q/2 - 0.5):1:(N_q/2)]*Ts;
    figure
    plot(Tsym * q)
    hold on
    plot(sinc(t/Tsym))
    title('Nyquist Pulse vs. Pulse Response autocorrelation');
    legend('Autocorr','Ideal Nyquist')
end

%% Waveform generation - upsample and filter

signals_up          = zeros(1,nSymbols*L);
signals_up(1:L:end) = tx_signals;

% Shaped waveform:
tx_waveform = conv(p, signals_up(:));

%% Transmission through channel

% Pre-noise rx waveform
rx_no_noise = tx_waveform;

% AWGN:
noise = sqrt(noise_en_per_dim/Tsym) * randn(size(rx_no_noise));
% Scale the noise_en_per_dim because it is being generated before the
% receive "analog" filter, whose energy is (1/Tsym)

% Rx waveform
if (en_noise)
    rx_waveform = rx_no_noise + noise;
else
    rx_waveform = rx_no_noise;
end

rx_waveform = Tsym * rx_waveform;
% Scale by Tsym as in Table 3.1
%
% Ts could be used as well, but then the anti-imaging filter gain of the
% interporlator should be taken into account (gain of L).

%% Equalization

switch (equalizer)
    case 1
        error('FIR ZF is not implemented yet');

    case 2
        fprintf('\n--------- MMSE --------\n');
        % Anti-alias receive filter gain of sqrt(Tsym) is incorporated
        rx_waveform = (1/sqrt(Tsym)) * rx_waveform;
        % This is the same as the convolution with the sinc filter whose
        % gain is sqrt(Tsym)

        % IMPORTANT: MMSE is fractionally spaced and incorporates both
        % matched filtering and equalization

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
        [SNR_mmse,w_t,opt_delay]=dfsecolorsnr(L,p_tilde,Nf,0,delta,...
            Ex,noise_en_per_dim*[1; zeros(Nf*L-1,1)])
        
        % Same as (inv(R_YY)*R_Yx)';
        z = conv(w, rx_waveform);

        % Skip MMSE filter delay and Acquire a window with nSymbols * L
        % samples. Again, recall nu and Nf are given in terms of T-spaced
        % symbols, not samples, so multiplication by L is required.

        z = z( (delta*L + 1) : (delta + nSymbols)*L );
        % Down-sample
        z_k = z(1:L:nSymbols*L).';
        % Remove
        z_k = z_k * (1/norm(w));

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
        % Matched filtering
        y = conv(rx_waveform, hrx);  % matched filtering
        % Symbol timing synchronization
        % Acquire a window with nSymbols * L samples
        y_s = y(n0:n0 + nSymbols*L - 1);
        % Followed by downsampling
        % T-spaced received symbols:
        y_k = y_s(1:L:end);
        % Equalized are equal to rx symbols (i.e. no equalization)
        z_k = y_k;
end

% Receiver Gain control (Assumed to be known)
z_k_unscaled = z_k / (Scale * norm_p);

if (debug)
    figure
    stem(Tsym*(0:nSymbols-1), unscaled_signals)
    xlabel('Tempo (s)')
    ylabel('Amplitude')
    grid on
    hold on
    stem(Tsym*(0:nSymbols-1), z_k_unscaled, 'r')
    legend('Original','Equalized')
end
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