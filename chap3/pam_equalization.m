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
% Tips to properly use the conversion factors of T from Table 3.1:
%   - If a function has unitary energy in continuous time, then its
%   discrete-time equivalent has norm^2 of 1/Ts (inverse of the sampling
%   period).
%   - Be very attentive to orthogonal vs. orthonormal expansions. The
%   factor of T above only applies in the specific case of the
%   sinc-interpolation formula, which is an orthogonal expansion. For
%   orthonormal expansions, it should not be used.
%

clearvars, close all;
addpath(genpath('../lib'))

% Parameters
N           =   1;      % Dimensions per symbol
nBits       =   1e3;    % Number of transmit symbols
debug       =   0;      % Enable plots and extra information
rollOff     =   0.1;    % Roll-off factor
L           =   2;      % Oversampliqng (support only for integer values)
W           =   1e3;    % Nominal bandwith (Hz)
N_T         =   10;     % Raised-cosine group-delay in symbols
Px          =   2e3;    % Transmit Power (W)
N0_over_2   =   0.181;  % Noise PSD (W/Hz/dim) and variance per dimension
M           =   4;
ideal_chan  =   0;
en_noise    =   1;
equalizer   =   2;      % 0) no equalizer; 1) FIR MMSE-DFE; 2) FIR MMSE-LE
% MMSE Parameters (if used)
Nf = 3;
% MMSE-DFE:
Nb = 2;

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
Scale = modnorm(pammod(0:(M-1), M), 'avpow', Ex);
% - Why Ex, when the modnorm criteria is "avpow"?
%   First of all, most DSP authors consider the average signal power as
% mean(abs(x).^2). Naturally, this is the convention adopted here within
% "modnorm" function. Note, however, that digital communications authors
% seldom forget about continuous-time, so they reconcile the two "worlds".
% For PAM, each constellation symbol is a coefficient that multiplies a
% unitary-energy transmit basis function, so by taking the integral in CT,
% it is easy to show that each Tsym-spaced pulses will have energy equal
% to |x_k|^2, where x_k is a symbol taken from the constellation alphabet
% (thought as an impulse in CT). Therefore, clearly "mean(abs(x).^2)" in
% this case is the average energy, not the average power. Finally, since
% mean(abs(x).^2) is the criteria that is adopted in modnorm, we pass Ex as
% argument, because:
%
%   E{ |x_k|^2 } = Ex
%

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
    % Ts * sum(abs(htx).^2) now is 1
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
    h = (1/Ts)*[0.9 1];
end
% Note: the energy is not very important here, because there is not too
% much to do about the channel attenuation/gain anyway.

%% Pulse response

p = Ts * conv(h, htx);
% Note: p (therefore h and htx) are sampled with Ts = Tsym/L (oversampled)

% Pulse norm:
norm_p_sq = Ts * norm(p)^2;
norm_p = sqrt(norm_p_sq);

% Unitary-energy (in continuous-time) pulse response:
phi_p = p / norm_p;

% Combined transmit -> matched filter
q = Ts * conv(phi_p, conj(fliplr(phi_p)));
[q_max,i_q0] = max(q);
if(q_max-1 > 1e-8)
    warning('q(t) peak is not unitary.');
end

fprintf('\n------------------ MFB ------------------\n\n');
% Note: SNR_{MFB} has to be found using p(t) before the anti-aliasing
% filter (see the answer to Exercise 3.34 in HW7).
SNRmfb = Ex_bar * norm_p_sq / N0_over_2;
fprintf('SNRmfb:   \t %g dB\n', 10*log10(SNRmfb))
% Average number of nearest neighbors:
Ne = 2 * (1 - 1/M);
% NNUB based on the SNRmfb
Pe = Ne * qfunc(sqrt(3*SNRmfb / (M^2 - 1)));
fprintf('Pe (NNUB):\t %g\n', Pe);

% Also consider distortion
if (~ideal_chan)
    fprintf('\n--------- ISI Characterization ----------\n\n');
    % Maximum value for |x_k|
    x_abs_max = max(abs(pammod(0:(M-1), M)));
    % Mean-square distortion - Eq. (3.33)
    D_ms = Ex * norm_p_sq * (sum(abs(q).^2) - 1);
    % -1 in the summation removes the magnitude of q_0
    % From (1.216):
    d_min = sqrt((12 * Ex) / (M^2 - 1));
    % Then, from (3.34):
    Pe = Ne * qfunc((norm_p * d_min) / (2 * sqrt(N0_over_2 + D_ms)));
    % Prints
    fprintf('Mean-Square Distortion:\t %g db/Hz\n', 10*log10(D_ms));
    fprintf('Pe (NNUB):             \t %g\n', Pe);
end

%% Equalizer Design

% Define equalizers and receive filters
switch (equalizer)
    case {1,2}
        % First of all, MMSE is fractionally spaced and incorporates both
        % matched filtering and equalization.
        % Secondly, it is preceded by an anti-alias filter whose gain is
        % sqrt(Tsym) over the entire band from -L*pi/Tsym to L*pi/Tsym,
        % namely over -pi/Ts to pi/Ts. Since a regular unitar-energy filter
        % over this band would present magnitude sqrt(Ts) in the freq.
        % domain, in order for the filter to have magnitude sqrt(Tsym), it
        % has to be scaled by a factor of "sqrt(Tsym/Ts)", which is the
        % same as "sqrt(L)".
        %   Therefore, the filter in time-domain is given by:
        %
        %   sqrt(L) * [ (1/sqrt(Ts)) * sinc(t/Ts) ]
        %
        % When sampled at t = kTs, since sinc(t/Ts) is 1 for k=0 and 0
        % elsewhere, the discrete sequence is [ L/sqrt(Tsym) ]
        hrx = sqrt(L/Ts);
        % Note: Ts * norm(hrx).^2 = L (filter energy is L).

        % Combined pulse response + anti-aliasing filter:
        p_tilde = Ts * conv(p, hrx);
        % Note: the Ts * convolution here is the same effect as:
        %   Ts * sqrt(L/Ts)
        %   = sqrt(Ts) * sqrt(L)
        %   = sqrt(L * Ts)
        %   = sqrt(Tsym)
        %
        %   So, the same effect could be obtained if the pulse response was
        %   multiplied by sqrt(Tsym) directly, as in the solution for
        %   Exercise 3.18

        %
        % FIR design
        %
        nu = ceil((length(p_tilde)-1)/L);  % Pulse response dispersion

        if (equalizer == 2)
            fprintf('\n-------------- FIR MMSE-LE --------------\n\n');
            Nb = 0;
        else
            fprintf('\n-------------- FIR MMSE-DFE -------------\n\n');
        end

        % The FIR Equalizer can be obtained using the DFE program:
        [SNR_mmse_unbiased_db,w_t,opt_delay]=dfsecolorsnr(...
            L,...
            p_tilde,...
            Nf,...
            Nb,...
            -1,... % choose the best delay delta
            Ex,...
            (L * N0_over_2)*[1; zeros(Nf*L-1,1)]);
        % Save the optimum delay
        delta = opt_delay;
        % Notes:
        %   #1: The last argument is the noise autocorrelation vector
        %       (one-sided).
        %   #2: For the MMSE receiver, in the presence of oversampling, the
        %       Rx filter is assumed to be preceded by a "brick-wall"
        %       anti-alising filter with cuttof at fs/2, or, equivalently,
        %       at 2*pi*fs/2 = pi/Ts, but with the same conventional
        %       magnitude sqrt(Tsym) that preserves the spectrum within the
        %       central period. Thus, the analog filter energy becomes L,
        %       rather than unitary. Consequently, the noise energy per
        %       dimension becomes (N0/2) * L. In contrast, the Tx signal
        %       energy is assumed to be contained within -1/Tsym to 1/Tsym
        %       and, thus, does not change.
        w = w_t(1:(Nf*L));
        b = -w_t((Nf*L) + 1:end);
        % The dfecolor function can adjust the length of b from the
        % specified value. Hence, we must update Nb to the actual length
        % that was designed:
        Nb = length(b);

        % Expected Performance
        SNR_fir_mmse_le_unbiased = 10^(SNR_mmse_unbiased_db/10);
        SNR_fir_mmse_le_biased   = SNR_fir_mmse_le_unbiased + 1;
        fprintf('Biased MMSE SNR:\t %g dB\n',...
            10*log10(SNR_fir_mmse_le_biased));
        fprintf('Unbiased MMSE SNR:\t %g dB\n',...
            10*log10(SNR_fir_mmse_le_unbiased));
        % NNUB based on the SNR_fir_mmse_le_unbiased
        Pe = 2 * (1 - 1/M) * ...
            qfunc(sqrt(3*SNR_fir_mmse_le_unbiased / (M^2 - 1)));
        fprintf('Pe (NNUB):      \t %g\n', Pe);
        gamma_mmse_le = 10*log10(SNRmfb / SNR_fir_mmse_le_unbiased);
        fprintf('MMSE gap to SNRmfb:\t %g dB\n', gamma_mmse_le);

        % Factor to remove bias:
        unbiasing_factor = SNR_fir_mmse_le_biased / SNR_fir_mmse_le_unbiased;

    otherwise
        % Matched filter receiver
        hrx  = conj(fliplr(htx));
        % Notes "Ts * conv(htx, hrx)" should be delta_k (for t = kTsym)

        % For the matched filter receiver, the factor that removes the bias
        % prior to the decision device is the reciprocal of ||p||. See,
        % e.g., the solution for exercise 3.4, which considers the
        % conventional matched filter receiver.
        unbiasing_factor = (1/norm_p);
end

%% Waveform generation - upsample and filter

signals_up          = zeros(1,nSymbols*L);
signals_up(1:L:end) = tx_signals;
% When upsampling a sequence and low-pass filtering to remove images, the
% ideal LPF should have normalized bandwidth pi/L and be of unitary energy.
% Hence, its gain has to be sqrt(L).

% Shaped waveform:
tx_waveform = conv(htx, signals_up(:));
% IMPORTANT: this is the only convolution that is not actually in the
% sense of Table 3.1. Convolution here is just an artifact for producing
% the orthogonal expansion that produces PAM modulated signal. Hence, the
% factor of Ts is not required.

if (debug)
    % To understand the following, consult page 26, chap 9 of Gallager's
    % book on Digital Comm I.
    fprintf('\n------- Energy/Power Measurements -------\n\n');
    % Due to the invariance of the inner product, the average transmit
    % energy (given the basis are orthonormal) should be close to Ex:
    tx_avg_energy = mean(abs(tx_signals).^2);
    % Note the above does not have the Ts factor, because "tx_signals"
    % multiply orthonormal basis already.
    fprintf('Measured average Tx energy:\t %g\n', tx_avg_energy);
    fprintf('Spec average Tx energy (Ex):\t %g\n', Ex);
    % Upsampled sequence average energy
    tx_total_energy = Ts * norm(tx_waveform).^2;
    tx_avg_energy_sampled = tx_total_energy / length(tx_waveform);
    % In contrast to "tx_avg_energy", the above comes from samples of the
    % sinc-interpolation formula, which is not an orthonormal expansion, but
    % orthogonal. Thus, the Ts is required.
    fprintf('Average sample energy (Es):\t %g\n', tx_avg_energy_sampled);
    fprintf('Observe that Es = Ex/L\n');
    fprintf('    (Ex/L):                \t %g\n', Ex/L);
    fprintf('\n');
    % The transmit power is equivalent to the mean in the transmit signal
    % sequence.
    Ex_over_Tsym = tx_avg_energy / Tsym;
    Es_over_Ts = tx_avg_energy_sampled / Ts;
    fprintf('Ex/Tsym:\t %g\n', Ex_over_Tsym);
    fprintf('Es/Ts:  \t %g\n', Es_over_Ts);
    fprintf('Spec Px:\t %g\n', Px);
    fprintf('\nTotal Energy\n');
    fprintf('%.E symbols require:\t %g u.e\n', nSymbols, Ex * nSymbols);
    fprintf('Measured energy:\t %g u.e\n', tx_total_energy);
end

%% Transmission through channel

% Receive signal past channel, but pre noise:
rx_pre_noise = Ts * conv(h, tx_waveform);

% AWGN:
noise = sqrt(N0_over_2/Ts) * randn(size(rx_pre_noise));
% An explanation can be found in Robert Gallager's material for the
% Principles of Digital Communications I course, Chaper 9, footnote 25.
% Rephrased slightly (for compatibility with the nomenclature in this
% script), it can be understood as follows:
%   The sinc functions in the orthogonal expansion (sampling theorem) have
%   energy Ts, so the variance of each real and imaginary coefficient in
%   the noise expansion must be scaled down by (Ts) from the noise energy
%   N0/2 per degree of freedom to compensate the sinc scaling. In the end,
%   iid variables with N0/2 variance are obtained.

if (debug)
    fprintf('\n-------- Noise Power Measurements -------\n\n');
    fprintf('Nominal N0/2:\t %g\n', N0_over_2);
    fprintf('Measured noise variance per real dim:\t %g\n', ...
        Ts * var(noise));
end

%% Anti-aliasing LPF at Rx

% Rx waveform - After anti-aliasing receive filtering:
if (en_noise)
    rx_waveform = Ts * conv(rx_pre_noise + noise, hrx);
else
    rx_waveform = Ts * conv(rx_pre_noise, hrx);
end

if (debug)
    pwelch(tx_waveform,[],[],[],Fs,'twosided');
    plot_obj = findobj(gcf);
    alllines=findall(plot_obj,'Type','line');
    set(alllines,'Color','red');
    hold on
    pwelch(noise,[],[],[],Fs,'twosided');
    hold on
    pwelch(rx_pre_noise,[],[],[],Fs,'twosided');
    plot_obj = findobj(gcf);
    plot_obj(end-2).Color = [1 0.5 0];
    legend('Tx-Pre-Chan', 'Noise', 'Tx-Post-Chan')
end

%% Equalize the received samples

% Preallocate decoded symbols
rx_decSymbols = zeros(nSymbols, 1);

switch (equalizer)
    case 1 % MMSE-DFE
        % Feed-forward section
        z = conv(w, rx_waveform);
        % Note: the Ts factor is not necessary here, since it is not
        % considered in the derivation. Furthermore, convolution here is
        % used solely to implement the inner product of (3.344).

        % Skip MMSE filter delay and Acquire a window with nSymbols * L
        % samples. Again, recall nu and Nf are given in terms of T-spaced
        % symbols, not samples, so multiplication by L is required.
        z = z( (delta*L + 1) : (delta + nSymbols)*L );
        % Down-sample to obtain a symbol-spaced sequence
        z_k = z(1:L:nSymbols*L).';

        % Since the DFE involves feedback, decision has to be made
        % iteratively
        z_dec = zeros(nSymbols, 1);
        b_mem = zeros(length(b) - 1);
        % Feedback loop
        for k = 1 : nSymbols

            if (k > 1)
                % Filter past "sliced" symbols (after decision)"
                [b_out, b_mem] = filter(b, 1, z_dec(k-1), b_mem);

                % Symbol after feedback filter
                z_prime_k = z_k(k) - b_out(1);
            else
                % For the first-ever received symbol, there is no feedback
                z_prime_k = z_k(k);
            end

            % Remove bias before decision and unscale back to original
            % constellation:
            rx_decSymbols(k) = ...
                pamdemod(z_prime_k * unbiasing_factor / Scale, M);
            % Re-modulate to obtain the corresponding symbols for the
            % decision:
            z_dec(k) = Scale * pammod(rx_decSymbols(k), M);
        end

    case 2 % MMSE-LE
        % Feed-forward section
        z = conv(w, rx_waveform);
        % Note: the Ts factor is not necessary here, since it is not
        % considered in the derivation. Furthermore, convolution here is
        % used solely to implement the inner product of (3.292).

        % Skip MMSE filter delay and Acquire a window with nSymbols * L
        % samples. Again, recall nu and Nf are given in terms of T-spaced
        % symbols, not samples, so multiplication by L is required.
        z = z( (delta*L + 1) : (delta + nSymbols)*L );
        % Down-sample
        z_k = z(1:L:nSymbols*L).';

    otherwise
        % Cursor for Symbol timing synchronization
        [~, n0] = max(conv(p, hrx));
        % Acquire a window with nSymbols * L samples
        y_s = rx_waveform(n0:n0 + nSymbols*L - 1);
        % Followed by downsampling (when L > 1)
        % T-spaced received symbols:
        y_k = y_s(1:L:end);
        % There is no equalizer, so:
        z_k = y_k;
end

% Remove bias:
z_k = z_k * unbiasing_factor;

% Scale back to "standard" pam constellation with d=2 for comparison:
z_k_unscaled = z_k / Scale;

if (debug)
    figure
    stem(Tsym*(0:nSymbols-1), unscaled_signals, 'o')
    xlabel('Time (s)', 'FontSize', 12)
    ylabel('Amplitude', 'FontSize', 12)
    grid on
    hold on
    stem(Tsym*(0:nSymbols-1), z_k_unscaled, 'ro', 'MarkerSize', 10)
    legend('Tx','Rx')
    set(gca, 'YTick', [-(M - 1):2:(M)])
end

%% Decision

% Except for the DFE, perform decision in the entire vector at once
if (equalizer ~= 1)
    rx_decSymbols = pamdemod(z_k_unscaled, M);
end
% Filter NaN
rx_decSymbols(isnan(rx_decSymbols)) = 0;
rx_binSymbols = de2bi(rx_decSymbols);
rx_bitstream  = rx_binSymbols(:);

%% Symbol error
fprintf('\n-------------- Performance --------------\n');

[nSymErrors, SER, symErrArray] = symerr(tx_decSymbols, rx_decSymbols(:));

fprintf('\nSymbol errors:\t %g\n',nSymErrors);
fprintf('SER:     \t %g\n', SER);

%% Bit error

[nBitErrors, BER, errArray] = biterr(tx_bitStream, rx_bitstream);

fprintf('\nBit errors:\t %g\n',nBitErrors);
fprintf('BER:      \t %g\n', BER);