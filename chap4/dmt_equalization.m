%% DMT Equalization Analysis through Monte-Carlo Simulation
clearvars, clc
addpath('../lib/')

%% Debug levels
debug               = 1;  % Enable debug information
debug_constellation = 0;  % Debug a certain subchannel constellation
debug_tone          = 16; % Tone whose constellation is debugged
debug_Pe            = 1;  % Debug error probabilities
debug_loading       = 0;  % Debug bit loading
debug_tx_energy     = 0;  % Debug transmit energy
debug_teq           = 0;  % Debug TEQ design

%% Parameters
alpha      = 1;         % Increase FFT size by this factor preserving Fs
% Note: this is useful to evaluate the DMT performance as N -> infty
L          = 1;         % Oversampling (support only for integer values)
Px         = 1e-3;      % Transmit Power (W)
N0_over_2  = 1e-10;     % Noise PSD (W/Hz/dim) and variance per dimension
N          = 128;       % FFT size and the number of used real dimensions
nu         = 8;         % Cyclic Prefix Length
nDim       = N + nu;    % Total number of real dimensions per DMT symbol
gap_db     = 8.8;       % SNR gap to capacity (dB)
delta_f    = 51.75e3;   % Subchannel bandwidth
nSymbols   = 1e3;       % Number of DMT symbols per transmission iteration
max_load   = inf;        % Maximum allowed bit load for each subchannel
equalizer  = 0;         % 0 - None; 1) MMSE-TEQ
% MMSE-TEQ Parameters
teqType    = 0;         % 0 - MMSE; 1 - SSNR; 2 - GeoSNR
% Monte-Carlo Parameters
maxNumErrs   = 100;
maxNumDmtSym = 1e12;

%% Derived computations:

% Number of used real dimensions and tone spacing adjusted by alpha
N         = N*alpha;
delta_f   = delta_f/alpha;
% Fs does not change with alpha.

% Sampling frequency and FFT size, based on the oversampling ratio:
Nfft      = L * N;
Fs        = Nfft * delta_f;
% delta_f does not change with L.

% Note the difference between alpha and L. Variable "alpha" is used to
% increase the density of the channel partitioning, namely reduce the tone
% spacing solely by increasing the number of used dimensions, while
% preserving the sampling frequency. In contrast, "L", the oversampling
% ratio, is used to increase the sampling frequency without altering the
% number of used dimensions, but only the FFT size. The latter naturally
% has to change according to "L" because the tone spacing must be
% preserved.

Ts        = 1 / Fs;
gap       = 10^(gap_db/10); % Gap in linear scale
Tofdm     = 1 / delta_f;    % OFDM symbol duration (without CP)
Tsym      = Tofdm + nu*Ts;  % Cyclic-prefixed multicarrier symbol period
Rsym      = 1 / Tsym;       % DMT Symbol rate (real dimensions per sec)
Ex        = Px * Tsym;      % Average DMT symbol energy
% Recall due to repetition of samples in the prefix, the energy budget used
% by the water-fill contraint must be Ex*(N/(N + nu)), so that the
% effective transmit energy in the end is equal to Ex.
Ex_bar    = Ex / nDim;      % Energy per real dimension

%% Constants

% TEQ criterion
TEQ_MMSE    = 0;
TEQ_SSNR    = 1;
TEQ_GEO_SNR = 2;

% Normalized FFT Matrix
Q = (1/sqrt(Nfft))*fft(eye(Nfft));

% Dimensions per subchannel:
dim_per_subchannel = [1 2*ones(1, N/2-1) 1 2*ones(1, N/2-1)];

% Vector of used subchannels among the NFFT indices
used_fft_indices = [1:(N/2), (Nfft - N/2 + 1):Nfft].';
%% Pulse Response

p = [-.729 .81 -.9 2 .9 .81 .729];
% Note: According to the model developed in Chap4 of EE379, p(t), the pulse
% response, corresponds to the combination between the post-DAC filter, the
% channel impulse response and the pre-ADC filtering.

% Pulse frequency response
H = fft(p, Nfft);
% Store only the response at the used indices of the FFT
H = H(used_fft_indices);

% Pulse response length
Lh = length(p);

% Matched-filter Bound
SNRmfb = (Ex_bar * norm(p).^2) / N0_over_2;
fprintf('SNRmfb:    \t %g dB\n\n', 10*log10(SNRmfb))

%% Equalizers
% The TEQ is designed before loading, by assuming flat input spectrum.
% However, since bit loading alters the energy allocation among
% subchannels, the TEQ is redesigned after bit loading.

switch (equalizer)
    case 1

fprintf('\n-------------------- MMSE-TEQ Design ------------------- \n\n');

        if (nu >= (Lh - 1))
            error('MMSE-TEQ is unecessary. CP is already sufficient!');
        end

        % Search optimum length and delay for the TEQ design
        [nTaps, delta] = optimizeTeq(teqType, p, nu, L, N0_over_2, ...
            Ex_bar, Nfft, debug_teq);
        fprintf('Optimal Equalizer Length:\t %d\n', nTaps);
        fprintf('Optimal Delay:           \t %d\n', delta);

        % Design final TEQ
        switch (teqType)
            case TEQ_MMSE
                [w, SNRteq] = ...
                    mmse_teq(p, L, delta, floor(nTaps/L), nu, Ex_bar, ...
                    N0_over_2, debug_teq);
                fprintf('New SNRmfb (TEQ):\t %g dB\n', 10*log10(SNRteq))
            case TEQ_SSNR
                w = ssnr_teq(p, L, delta, floor(nTaps/L), nu, debug_teq);
        end


        if(~isreal(w))
            warning('MMSE-TEQ designed with complex taps');
        end

        % New effective channel:
        p_eff = conv(p,w);

        % Shortening SNR:
        ssnr = teqSSNR( w, p, delta, nu );
        fprintf('SSNR:\t %g dB\n', 10*log10(ssnr));

        % Notes:
        %   # 1) The water-filling solution assumes no ISI/ICI. This is
        %   safe provided that the TEQ constrains the pulse response energy
        %   to a portion that can be covered by the guard band. However,
        %   with windowing+overlap this fails.
        %   # 2) Note the feed-foward TEQ at the receiver shapes the
        %   spectrum of the noise, so the noise PSD becomes |H_w|^2 * N0/2,
        %   where H_w is given below:
        H_w = fft(w, Nfft);
        %   # 3) Meanwhile, the transmit signal is subject to the
        %   compounded response of the channel + equalizer, namely the
        %   effective response:
        H_eff = fft(p_eff, Nfft);
        %   The gain to noise ratio at the receiver becomes:
        gn_teq = (abs(H_eff).^2)./(N0_over_2 * abs(H_w).^2);
        % Store only the used tones
        gn_teq = gn_teq(used_fft_indices);

    otherwise
        % When MMSE-TEQ is not used, at least a bias should be generically
        % defined:
        bias = 1;
end

%% 1-tap Frequency Equalizer

% Frequency domain response
switch (equalizer)
    case 1
        % Use the effective pulse response in the FEQ
        H_freq = fft(p_eff, Nfft);
    otherwise
        H_freq = fft(p, Nfft);
end

% Store only the response at the used indices of the FFT
H_freq = H_freq(used_fft_indices);

% Cursor
switch (equalizer)
    case 1
        % MMSE-TEQ Chosen Delay
        n0 = delta;
        % The cursor considers the MMSE-TEQ delay.
    otherwise
        [~, iMax] = find(abs(p)>1e-5, 1, 'first');
        n0 = iMax - 1;

end

% Corresponding phase shift due to cursor
phaseShift = exp(1j*2*pi*(n0/Nfft)*(used_fft_indices.' - 1));

% Include the unbiasing factor into the FEQ. Note bias=1 when MMSE-TEQ
% is not used.
FEQ    = (1/bias) * (1 ./ (H_freq .* phaseShift));

%% SNR for unitary energy transmit symbol

switch (equalizer)
    case 1
        % For the MMSE-TEQ, the effective pulse responsive becomes the
        % result of the convolution between the actual pulse response and
        % the feed-forward equalizer.
        gn = gn_teq;
        Ex_red_factor = 1;
    otherwise
        gn = (abs(H).^2) / N0_over_2;
        Ex_red_factor = 1;
end

%% Water filling

fprintf('\n--------------------- Water Filling -------------------- \n\n');

% Water-filling:
[bn_bar, En_bar, usedTones] = waterFilling(gn, Ex_bar, N, gap);

% Bits per subchannel
bn = bn_bar(1:N/2+1) .* dim_per_subchannel(1:N/2+1);

% Number of bits per dimension
b_bar = (1/nDim)*(sum(bn_bar));
fprintf('\nb_bar:                  \t %g bits/dimension\n', b_bar)
% For gap=0 and N->+infty, this should be the channel capacity per real
% dimension.

% Corresponding multi-channel SNR:
SNRdmt = 10*log10(gap*(2^(2*b_bar)-1));
% SNR at each tone, per dimension:
SNR_n = En_bar .* gn;
% Normalized SNR on each tone, per dimension (should approach the target
% gap):
SNR_n_norm = SNR_n ./ (2.^(2*bn_bar) - 1);

fprintf('Multi-channel SNR (SNRdmt):\t %g dB\n', SNRdmt)

if (equalizer == 1)
    fprintf('Note: shortened response was used for water-filling.\n');
end

%% Discrete-loading: Levin Campello Rate Adaptive

fprintf('\n------------------ Discrete Loading -------------------- \n\n');

% Rate-adaptive Levin-Campello loading:
noDcNyquist = 1; % Avoid loading DC and Nyquist subchannels
[En_discrete, bn_discrete] = DMTLCra(...
    gn(1:N/2 + 1),...
    Ex_bar,...
    N, gap_db, ...
    max_load,...
    noDcNyquist);

% Save a vector with the index of the subchannels that are loaded
n_loaded = find(bn_discrete ~= 0);

% Energy per real dimension
En_bar_lc = [En_discrete, fliplr(conj(En_discrete(2:N/2)))] ...
    ./ dim_per_subchannel;
% Bits per subchannel per dimension
bn_bar_lc = [bn_discrete(1:N/2+1), fliplr(bn_discrete(2:N/2))] ./ ...
    dim_per_subchannel;

% Total bits per dimension:
b_bar_discrete = 1/nDim*(sum(bn_discrete));

% SNRdmt from the number of bits per dimension
SNRdmt_discrete    = gap*(2^(2*b_bar_discrete)-1);
SNRdmt_discrete_db = 10*log10(SNRdmt_discrete);
% SNR on each tone, per dimension:
SNR_n_lc           = En_bar_lc .* gn; % SNR per real dimension
% Normalized SNR on each tone, per dimension (should approach the gap)
SNR_n_norm_lc      = SNR_n_lc ./ (2.^(2*bn_bar_lc) - 1);

% Bit rate
Rb = sum(bn_discrete) / Tsym;


fprintf('b_bar:                    \t %g bits/dimension', b_bar_discrete)
fprintf('\nBit rate:               \t %g mbps\n', Rb/1e6);
fprintf('Multi-channel SNR (SNRdmt): \t %g dB\n', ...
    SNRdmt_discrete_db);

% Compare water-filling and discrete-loading
if (debug && debug_loading)
    figure
    plot(bn(1:N/2+1), ...
        'linewidth', 1.1)
    hold on
    plot(bn_discrete, 'g')
    legend('Water-filling', 'Discrete Loading')
    xlabel('Subchannel');
    ylabel('Bits');
    set(gca,'XLim',[1 N/2+1]);
end

%% Channel Capacity
% Channel capacity is computed considering the SNR that results from LC
% discrete loading.

fprintf('\n------------------ Channel Capacity -------------------- \n\n');

% Capacity per real dimension
cn = 0.5 * log2(1 + SNR_n_lc);
% Multi-channel capacity, per dimension:
c = sum(cn) / nDim;
% Note #1: for the capacity computation, all real dimensions are
% considered, including the overhead. See the example of (4.208)
% Note #2: the actual capacity is only obtained for N -> infty, so the
% above is only an approximation.

fprintf('capacity:               \t %g bits/dimension', c)
fprintf('\nBit rate:               \t %g mbps\n', c * Rsym * (N+nu) /1e6);

%% Analysis of the Error Probability per dimension
% Comparison between the water-filling and the discrete loading in terms of
% the nearest-neighbors union bound (also known as Union Bound Estimate).
% Variable "Ne" is used below to denote the average number of Nearest
% Neighbors per real dimension.

fprintf('\n----------------- Error Probabilities ------------------ \n\n');

% Preallocate
Pe_bar_n    = zeros(N/2 + 1, 1);
Pe_bar_n_lc = zeros(N/2 + 1, 1);

for k = n_loaded
    if (dim_per_subchannel(k) == 2)
        % QAM Nearest-neighbors Union Bound assuming QAM-SQ constellations
        % for both even and odd loads. "Hybrid QAM" constellations are used
        % as "Square" for odd loads.

        % Water-filling (with fractional load):
        if ((mod(bn(k),2) ~= 0))
            % For odd b, Hybrid QAM is used:
            M = 2^bn(k);
            Ne = 2*(1 - sqrt(2/M) + 1/(2*M));
            Pe_bar_n(k) = Ne * qfunc(sqrt( 3*SNR_n_norm(k)) );
            % Note: the argument here should be actually:
            %   sqrt((6 * SNR_n(k)) / (2*M -1)),
            % which is approximately equivalent to sqrt( 3*SNR_n_norm(k))
            % for sufficiently large M. The approximation won't be very
            % tight for low M. Since for WF any fractional M is allowed,
            % using the actual value will lead to higher than expected
            % error probabilities, so we simply use the approximation. For
            % LC the Pe computation will be tighter.
        else
            % QAM-SQ
            Ne = 2 * (1 - 1/(2^bn_bar(k)));
            Pe_bar_n(k) = Ne * qfunc(sqrt( 3*SNR_n_norm(k)) );
        end

        % Levin-Campello (LC):
        if ((mod(bn_discrete(k),2) ~= 0))
            % For odd b, Hybrid QAM is used
            M = 2^bn_discrete(k);
            Ne = 2*(1 - sqrt(2/M) + 1/(2*M));
            Pe_bar_n_lc(k) = Ne * qfunc(sqrt((6 * SNR_n_lc(k))/(2*M -1)));
            % Note the aforementioned q-function argument for Hybrid QAM
        else
            % QAM-SQ
            Ne = 2 * (1 - 1/(2^bn_bar_lc(k)));
            Pe_bar_n_lc(k) =  Ne * qfunc(sqrt( 3*SNR_n_norm_lc(k) ));
        end

    else
        % PAM Nearest-neighbors Union Bound

        % Water-filling (with fractional load):
        Ne = 2 * (1 - 1/(2^bn_bar(k)));
        Pe_bar_n(k) = Ne * qfunc(sqrt( 3*SNR_n_norm(k)) );

        % Levin-Campello (LC):
        Ne = 2 * (1 - 1/(2^bn_bar_lc(k)));
        Pe_bar_n_lc(k) = Ne * qfunc(sqrt( 3*SNR_n_norm_lc(k) ));
    end
end

if (debug && debug_Pe)
    figure
    plot(Pe_bar_n, 'linewidth', 1.1)
    hold on
    plot(Pe_bar_n_lc, 'r')
    legend('Water-filling', 'Levin-Campello')
    title('Pe per dimension on each subchannel')
    set(gca,'XLim',[1 N/2+1]);
end

fprintf('Approximate NNUB Pe per dimension:\n');
fprintf('Fractional-load (WF):\t %g\n', mean(Pe_bar_n(n_loaded), ...
    'omitnan'));
fprintf('Discrete-load (LC)  :\t %g\n', mean(Pe_bar_n_lc(n_loaded), ...
    'omitnan'));

%% Modulators

% Modulation order on each subchannel
modOrder = 2.^bn_discrete;

oneDimModOrders = [modOrder(1), modOrder(N/2 + 1)];
twoDimModOrders = modOrder(2:N/2);
oneDim_const_orders = unique(oneDimModOrders(oneDimModOrders~=1));
twoDim_const_orders = unique(twoDimModOrders(twoDimModOrders~=1));

%Preallocate modems
modulator = cell(length(twoDim_const_orders), 1);
demodulator = cell(length(twoDim_const_orders), 1);

% Configure 2-dimensional modems for each distinct bit loading:
for i = 1:length(twoDim_const_orders)
    M = twoDim_const_orders(i);

    if (mod(log2(M),2) ~= 0)
        modulator{i} = modem.genqammod('Constellation', ...
            qamHybridConstellation(M));
        demodulator{i} = modem.genqamdemod('Constellation', ...
            qamHybridConstellation(M));
    else
        modulator{i} = modem.qammod('M', M, 'SymbolOrder', 'Gray');
        demodulator{i} = modem.qamdemod('M', M, 'SymbolOrder', 'Gray');
    end
end

for l = 1:length(oneDim_const_orders)
    i = i + 1;
    M = oneDim_const_orders(l);
    modulator{i} = modem.pammod('M', M, 'SymbolOrder', 'Gray');
    demodulator{i} = modem.pamdemod('M', M, 'SymbolOrder', 'Gray');
end

%% Look-up table for each subchannel indicating the corresponding modem

modem_n = zeros(N/2 + 1, 1);

for k = n_loaded % Iterate over loaded subchannels
    if (dim_per_subchannel(k) == 2)
        iModem = find (twoDim_const_orders == modOrder(k));
        if (iModem)
            modem_n(k) = iModem;
        end
    else
        iModem = find (oneDim_const_orders == modOrder(k));
        if (iModem)
            modem_n(k) = length(twoDim_const_orders) + iModem;
        end

    end
end

%% Energy loading (constellation scaling factors) and minimum distances

[Scale_n, dmin_n] = dmtSubchanScaling(modulator, modem_n, En_discrete);

%% Monte-carlo

fprintf('\n---------------------- Monte Carlo --------------------- \n\n');

% Preallocate
X          = zeros(Nfft, nSymbols);
tx_data    = zeros(N/2 + 1, nSymbols);
rx_data    = zeros(N/2 + 1, nSymbols);
sym_err_n  = zeros(N/2 + 1, 1);

numErrs = 0; numDmtSym = 0;

% Sys Objects
BitError = comm.ErrorRate;

%% Iterative Transmissions

iTransmission = 0;

while ((numErrs < maxNumErrs) && (numDmtSym < maxNumDmtSym))
    iTransmission = iTransmission + 1;

    % Random Symbol generation
    for k = n_loaded
        tx_data(k, :) = randi(modOrder(k), 1, nSymbols) - 1;
    end

    %% Constellation Encoding
    for k = n_loaded
        if (modem_n(k) > 0)
            X(k, :) = Scale_n(k) * ...
                modulator{modem_n(k)}.modulate(tx_data(k, :));
        end
    end

    % Hermitian symmetry
    X(Nfft - N/2 + 2:Nfft, :) = flipud( conj( X(2:N/2, :) ) );

    %% Modulation

    x = sqrt(Nfft) * ifft(X, Nfft);

    %% Cyclic extension

    x_ext = [x(Nfft-nu+1:Nfft, :); x];

    %% Parallel to serial

    u = x_ext(:);

    if (debug && debug_tx_energy)
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
        fprintf('Tx Energy p/ Sym:\t%g\t', ...
            tx_total_energy / nSymbols);
        fprintf('Nominal:\t%g\t',Ex);
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

    %% Time-domain Equalization
    switch (equalizer)
        case 1
            z = conv(w, y);
        otherwise
            z = y;
    end

    %% Synchronization
    % Note: synchronization introduces a phase shift that should be taken
    % into account in the FEQ.

    nRxSamples = (Nfft+nu)*nSymbols;
    y_sync     = z((n0 + 1):(n0 + nRxSamples));

    %% Slicing

    y_sliced = reshape(y_sync, Nfft + nu, nSymbols);

    %% Extension removal

    y_no_ext = y_sliced(nu + 1:end, :);

    %% Regular Demodulation (without decision feedback)

    % FFT
    Y = (1/sqrt(Nfft)) * fft(y_no_ext, Nfft);

    % FEQ - One-tap Frequency Equalizer
    Z = diag(FEQ) * Y(used_fft_indices, :);

    %% Constellation decoding (decision)

    for k = n_loaded
        if (modem_n(k) > 0)
            rx_data(k, :) = demodulator{modem_n(k)}.demodulate(...
                (1/Scale_n(k)) * Z(k, :));
        end
    end

    % Symbol error count
    sym_err_n = sym_err_n + symerr(tx_data, rx_data, 'row-wise');
    % Symbol error rate per subchannel
    ser_n     = sym_err_n / (iTransmission * nSymbols);
    % Per-dimensional symbol error rate per subchannel
    ser_n_bar = ser_n ./ dim_per_subchannel(1:N/2+1).';

    % Preliminary results
    numErrs   = sum(sym_err_n);
    numDmtSym = iTransmission * nSymbols;

    fprintf('Pe_bar:\t%g\t', mean(ser_n_bar(n_loaded)));
    % Note: consider only the loaded subchannels in the above
    fprintf('nErrors:\t%g\t', numErrs);
    fprintf('nDMTSymbols:\t%g\n', numDmtSym);


end

%% Constellation plot for debugging
if (debug && debug_constellation && modem_n(debug_tone) > 0)
    k = debug_tone;
    figure
    if (dim_per_subchannel(k) == 2)
        plot(Z(k, :), 'o')
        hold on
        plot(Scale_n(k) * ...
            modulator{modem_n(k)}.modulate(0:modOrder(k) - 1),...
            'ro', 'MarkerSize', 8, 'linewidth', 2)
    else
        plot(Z(k, :) + j*eps, 'o')
        hold on
        plot(Scale_n(k) * ...
            modulator{modem_n(k)}.modulate(0:modOrder(k) - 1) ...
            + j*eps, 'ro', 'MarkerSize', 8, 'linewidth', 2)
    end
    legend('Rx', 'Tx')
    title(sprintf('Tone: %d', debug_tone));
end

%% Results
fprintf('\n----------------------- Results ------------------------ \n\n');
fprintf('Pe_bar:       \t %g\n', mean(ser_n_bar));

if (debug && debug_Pe)
    figure
    stem(ser_n_bar)
    hold on
    stem(Pe_bar_n, 'g')
    hold on
    stem(Pe_bar_n_lc, 'r')
    title('Results: Pe per dimension')
    xlabel('Subchannel (n)')
    ylabel('$\bar{Pe}(n)$')
    legend('Measured','WF','LC')
    set(gca,'XLim',[1 N/2+1]);
end
