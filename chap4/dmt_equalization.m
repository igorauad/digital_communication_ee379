%% DMT Equalization Analysis through Monte-Carlo Simulation
clearvars, clc

% Parameters
debug       = 1;       % Enable debug information
L           = 1;       % Oversampling (support only for integer values)
W           = 1e3;     % Nominal bandwith (Hz)
Px          = 1e-3;    % Transmit Power (W)
N0_over_2   = 1e-13;   % Noise PSD (W/Hz/dim) and variance per dimension
N           = 1024;    % FFT size and the number of used real dimensions
nu          = 24;      % Prefix
nDim        = N + nu;  % Total number of real dimensions per DMT symbol
gap_db      = 8.8;     % SNR gap to capacity (dB)
delta_f     = 1e3;     % Subchannel bandwidth
L           = 1;       % Oversampling Ratio
nSymbols    = 100;     % Number of transmit symbols
loading     = 1;       % 0 - Water-fill; 1 - Discrete (LC Rate Adaptive)

% Derived computations:
Fs        = N * delta_f;
Ts        = 1 / Fs;
gap       = 10^(gap_db/10); % Gap in linear scale
Tsym      = (N + nu) * Ts;  % Symbol Period
Rsym      = 1 / Tsym;       % DMT Symbol rate (real dimensions per sec)
Ex        = Px * Tsym;      % Average DMT symbol energy
% Recall due to repetition of samples in the prefix, the energy budget used
% by the water-fill contraint must be Ex*(N/nDim), so that the effective
% transmit energy in the end is equal to Ex.
Ex_bar    = Ex / N; % Energy per used dimension
% Used dimensions: DC + Nyquist + (N/2)-1 complex dimensions

%% Pulse Response

p = [-.729 .81 -.9 2 .9 .81 .729];
% Note: According to the model developed in Chap4 of EE379, p(t), the pulse
% response, corresponds to the combination between the post-DAC filter, the
% channel impulse response and the pre-ADC filtering.

% Pulse frequency response
H = fft(p, N);

% Pulse response length
Lh = length(p);

% Matched-filter Bound
SNRmfb = (Ex_bar * norm(p).^2) / N0_over_2;
fprintf('SNRmfb:                  \t %g dB\n\n', 10*log10(SNRmfb))

%% Water filling

% SNR for unitary energy transmit symbol:
gn = (abs(H).^2) / N0_over_2;

[bn_bar, En_bar, usedTones] = waterFilling(gn, Ex_bar*(N/nDim), N, gap);
dim_per_subchannel = [1 2*ones(1, N/2-1) 1 2*ones(1, N/2-1)];
unusedTones = setdiff(1:N, usedTones);

% Number of bits per dimension
b_bar = (1/nDim)*(sum(bn_bar));
fprintf('\nb_bar:                  \t %g bits/dimension\n\n', b_bar)
% For gap=0 and N->+infty, this should be the channel capacity per real
% dimension.

% Corresponding SNR:
SNRdmt = 10*log10(gap*(2^(2*b_bar)-1));

% Number of used tones, according to the water-filling:
N_star = length(usedTones);

fprintf('Multi-channel SNR (SNRdmt):\t %g dB\n\n', SNRdmt)

%% Discrete-loading: Levin Campello Rate Adaptive

[En_discrete, bn_discrete] = DMTLCra(...
    gn(1:N/2 + 1),...
    Ex_bar*(N/nDim),...
    N, gap);

% Energy per real dimension
En_bar_discrete = [En_discrete, flipud(conj(En_discrete(2:N/2)))] ...
    ./ dim_per_subchannel;

% Total bits per dimension:
b_bar_discrete = 1/nDim*(sum(bn_discrete));
% SNRdmt from the number of bits per dimension
SNRdmt_discrete = 10*log10(gap*(2^(2*b_bar_discrete)-1));

fprintf('Multi-channel SNR with discrete Loading: \t %g dB\n\n', ...
    SNRdmt_discrete);

% Compare water-filling and discrete-loading
if (debug && loading)
    stem(bn_bar(1:N/2+1) .* dim_per_subchannel(1:N/2+1))
    hold on
    stem(bn_discrete, 'g')
    legend('Water-filling', 'Discrete Loading')
    xlabel('Subchannel');
    ylabel('Bits');
    set(gca,'XLim',[0 N/2+1]);
end
%% Monte-carlo

% Constants
maxNumErrs = 100;
maxNumBits = 1e12;

% Preallocate
X          = zeros(N, nSymbols);
tx_symbols = zeros(N/2 + 1, nSymbols);
Scale_n    = zeros(N/2 + 1, nSymbols); % Constellation scaling factors

numErrs = 0; numBits = 0; results=zeros(3,1);

% Sys Objects
BitError = comm.ErrorRate;

% Constantes

% Normalized FFT Matrix
Q = (1/sqrt(N))*fft(eye(N));

%% Discrete-loading of bits and energy

% Modulation order on each subchannel
if (loading == 0)
    % Non-optimal discrete loading
    modOrder = 2.^floor(bn_bar(1:N/2+1) .* dim_per_subchannel(1:N/2+1));
else
    modOrder = 2.^bn_discrete;
end

% Scale factors for each subchannel n to achieve the En average energy
for k = 1:(N/2 + 1)
    if (En_bar(k) > 0 & modOrder(k) > 1)
        if (dim_per_subchannel(k) == 2)
            if (loading == 0)
                Scale_n(k) = modnorm(qammod(0:(modOrder(k)-1),...
                    modOrder(k)),...
                    'avpow', En_bar(k));
            else
                Scale_n(k) = modnorm(qammod(0:(modOrder(k)-1),...
                    modOrder(k)),...
                    'avpow', 0.5 * En_discrete(k));
            end
        else
            if (loading == 0)
                Scale_n(k) = modnorm(pammod(0:(modOrder(k)-1),...
                    modOrder(k)),...
                    'avpow', En_bar(k));
            else
                Scale_n(k) = modnorm(pammod(0:(modOrder(k)-1),...
                    modOrder(k)),...
                    'avpow', En_discrete(k));
            end
        end
    end
end

%% Iterative Transmissions

while ((numErrs < maxNumErrs) && (numBits < maxNumBits))
    % Random Symbol generation
    for k = 1:(N/2 + 1)
        tx_symbols(k, :) = randi(modOrder(k), 1, nSymbols) - 1;
    end

    %% Constellation Encoding
    for k = 1:(N/2 + 1)
        if (dim_per_subchannel(k) == 2)
            X(k, :) = Scale_n(k) * qammod(tx_symbols(k, :), modOrder(k));
        else
            X(k, :) = Scale_n(k) * pammod(tx_symbols(k, :), modOrder(k));
        end
    end

    % Hermitian symmetry
    X(N/2+2:end, :) = flipud( conj( X(2:N/2, :) ) );

    %% Modulation

    x = Q' * X;

    %% Cyclic extension

    x_ext = [x(N-nu+1:N, :); x];

    %% Parallel to serial

    u = x_ext(:);

    if (debug)
        % Note: "u" should become samples leaving the DAC. In that case,
        % they would repreent coefficients of the sampling theorem's sinc
        % interpolation formula, which is an orthogonal (non-normal)
        % expansion. However, note x comes from an orthonormal expansion,
        % which is the normalized IDFT. Hence, the energy in x at this
        % point is still given simply by:
        tx_total_energy = norm(u).^2;
        % A Ts factor should multiply the norm if u was a vector of samples
        % out of the DAC, but in this case there would be a scaling factor
        % introduced by the DAC anti-imaging LPF. Both would cancel each
        % other.
        fprintf('Transmit energy per symbol:\t %g\n', ...
            tx_total_energy / nSymbols);
        fprintf('Nominal transmit energy per symbol:\t %g\n',Ex);
    end

    %% Channel

    y = conv(u, p);
    % Note:
    %   In contrast to the derivation of Chapter 3, here the scaling of Ts
    %   is not used. The reason is that "u" comes from an orthonormal
    %   expansion (the normalized IFFT) and p also satisfies the inner
    %   product invariance, due to the more involved explanation in the
    %   sequel.
    %
    %   The text mentions that "the notational use of P for the channel
    %   matrix suggests that any anti-alias analog filters at transmitter
    %   and receiver have been convolved with the channel impulse response
    %   h(t) and included in the discrete-time response of the matrix
    %   channel". The model is given in (4.185).
    %
    %   So let us further interpret that: once we have a measurement of the
    %   channel impulse response, we have samples and, therefore,
    %   coefficients of an orthogonal expansion. Thus, convolution or inner
    %   product (for energy computation) can not be applied directly,
    %   only with a Ts factor in front. Nonetheless, in practice the
    %   samples would be obtained by first passing the signal through an
    %   anti-alias LPF, ideally with unitary-energy. Such filter
    %   effectively scales each "pre-ADC" sample by sqrt(Ts). To understand
    %   that, note a unitary energy anti-alias filter has continuous-time
    %   response:
    %
    %       1/sqrt(Ts) * sinc(t/Ts)
    %
    %   Then, by sampling at t = kTs, the only non-zero value of the
    %   sampled sequence is "1/sqrt(Ts)". Finally, by convolving the
    %   "pre-ADC" samples (from an orthogonal expansion) with the receive
    %   filter samples, one obtains:
    %
    %       h = Ts * conv(rx_filter, h)
    %         = Ts * (1/sqrt(Ts)) * h
    %         = sqrt(Ts) * h_pre_adc
    %
    %   Essentially, when we sample the channel, we get h, instead of
    %   h_pre_adc. The channel impulse response energy can be computed by
    %   either "Ts * norm(h_pre_adc)^2" or "norm(h)^2", because both are
    %   equivalent. Similarly, the sampled response "h" can be used
    %   directly in the convolution, without a Ts factor in front.

    % Add noise
    y = y + (sqrt(N0_over_2) * randn(length(y),1));
    % Important considerations:
    %
    % First, recall the noise continuous-time PSD coincides with the noise
    % energy per dimension. Second, remember that the sinc functions in the
    % orthogonal expansion of the sampling theorem have energy 1/2W, so the
    % variance of each real and imaginary coefficient in the noise
    % expansion must be scaled up by 2W from the noise energy N0/2 per
    % degree of freedom. Since AWGN is flat, multiplication of N0/2
    % (two-sided PSD) by 2W yields the total transmit power. Hence, if "y"
    % consisted of samples, the target variance for the "randn" sequence
    % would be the noise power N0_over_2 * 2W. However, the catch here is
    % that y does not represent the samples

    %% Synchronization

    nRxSamples = (N+nu)*nSymbols;
    y_sync     = y(1:nRxSamples);

    %% Slicing

    y_sliced = reshape(y_sync, N + nu, nSymbols);

    %% Extension removal

    y_no_ext = y_sliced(nu + 1:end, :);

    %% Demodulation

    Y = Q * y_no_ext;

    %% FEQ - One-tap Frequency Equalizer

    H_freq = fft(p, N);
    FEQ    = 1 ./ H_freq;

    Z = diag(FEQ) * Y;

    %% Constellation decoding (decision)
    rx_symbols = zeros(N/2+1, nSymbols);

    for k = 1:(N/2 + 1)
        if (dim_per_subchannel(k) == 2)
            rx_symbols(k, :) = qamdemod((1/Scale_n(k)) * Z(k, :), ...
                modOrder(k));
        else
            rx_symbols(k, :) = pamdemod((1/Scale_n(k)) * Z(k, :), ...
                modOrder(k));
        end
    end

    results = BitError.step(tx_symbols(:), rx_symbols(:)); % Update BER
    numErrs = results(2);
    numBits = results(3);

    if (debug)
        fprintf('SER:\t%g\t', results(1));
        fprintf('nErrors:\t%g\t', results(2));
        fprintf('nSymbols:\t%g\n', results(3));
    end
end

fprintf('SER:          \t %g\n', results(1));
fprintf('Total errors: \t %g\n', results(2));
fprintf('Total symbols:\t %g\n', results(3));
