%% DMT Equalization Analysis through Monte-Carlo Simulation
clearvars, clc
addpath(genpath('../lib'))

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
tau        = 8;         % Cyclic Suffix
windowing  = 0;         % Activate Lcs windowing + Overlap
gap_db     = 8.8;       % SNR gap to capacity (dB)
delta_f    = 51.75e3;   % Subchannel bandwidth
nSymbols   = 1e3;       % Number of DMT symbols per transmission iteration
max_load   = inf;        % Maximum allowed bit load for each subchannel
equalizer  = 0;         % 0 - None; 1) TEQ; 2) Cheong; 3) Time Domain
dcNyquist  = 0;         % Flag to enable loading DC and Nyquist subchan
% MMSE-TEQ Parameters
teqType    = 0;         % 0 - MMSE; 1 - SSNR; 2 - GeoSNR
% Monte-Carlo Parameters
maxNumErrs   = 100;
maxNumDmtSym = 1e12;
% Channel
channelChoice = 0;

%% Derived computations:

% Number of used real dimensions and tone spacing adjusted by alpha
N         = N*alpha;
delta_f   = delta_f/alpha;
% Fs does not change with alpha.

% Sampling frequency and FFT size, based on the oversampling ratio:
Nfft      = L * N;
Fs        = Nfft * delta_f;
% delta_f does not change with L.

% Total number of real dimensions per DMT symbol
nDim       = Nfft + nu;

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
% Due to repetition of samples in the prefix, the energy budget is not
% entirely used for data transmission. Hence, the water-fill contraint must
% be designed as the energy budget discounted by the energy allocated to
% the cylic prefix. By setting the contraint to Ex*(Nfft/(Nfft + nu)), the
% effective transmit energy in the end is equal to Ex, as desired.
Ex_budget = Ex*(Nfft/(Nfft + nu)); % Energy budget passed to the WaterFill
Ex_bar    = Ex / nDim;      % Energy per real dimension

%% DMT Struct
% Create an object with all the parameters

dmtObj = [];

% Fixed Parameters
dmtObj.nSymbols          = nSymbols;
dmtObj.Nfft              = Nfft;
dmtObj.N                 = N;
dmtObj.L                 = L;
dmtObj.nu                = nu;
dmtObj.nDim              = nDim;
dmtObj.N0_over_2         = N0_over_2;
dmtObj.Ex_bar            = Ex_bar;
dmtObj.Ex_budget         = Ex_budget;
dmtObj.Tsym              = Tsym;
dmtObj.gap_db            = gap_db;
dmtObj.max_load          = max_load;

%% Constants
POST_PRE_ICPD_FLAG = 0;

% TEQ criterion
TEQ_MMSE    = 0;
TEQ_SSNR    = 1;

% Equalizer Types
EQ_NONE      = 0;
EQ_TEQ       = 1;
EQ_FREQ_PREC = 2;
EQ_TIME_PREC = 3;

% Normalized FFT Matrix
Q = (1/sqrt(Nfft))*fft(eye(Nfft));

%% Look-up tables for the indexes and number of dimensions of subchannels

[subCh_tone_index, subCh_tone_index_herm] = ...
    dmtSubchIndexLookup(N, Nfft, L, dcNyquist);

[dim_per_subchannel, dim_per_dft_tone] = ...
    dmtSubchDimLookup(Nfft, subCh_tone_index);

% Number of available subchannels
N_subch  = length(subCh_tone_index);

% Copy to DMT Object
dmtObj.iTonesTwoSided     = subCh_tone_index_herm;
dmtObj.iTones             = subCh_tone_index;
dmtObj.dim_per_subchannel = dim_per_subchannel;
dmtObj.N_subch            = N_subch;

%% Pulse Response

switch (channelChoice)
    case 0
        p = [-.729 .81 -.9 2 .9 .81 .729];
        % Note: According to the model developed in Chap4 of EE379, p(t),
        % the pulse response, corresponds to the combination between the
        % post-DAC filter, the channel impulse response and the pre-ADC
        % filtering.
    case 1
        % D2-H2
        N_old = N;
        load('/Users/igorfreire/Documents/Lasse/gfast_simulator/Channel_Model/data/all_models_106mhz.mat');
        N = N_old;
        iModel = 2;
        H_herm = [H(:,iModel,iModel); conj(flipud(H(2:end-1,iModel,iModel)))];
        if (any(imag(ifft(H_herm, N)) > 1e-8))
            warning('H is not Hermitian symmetric');
        end
        clear H
        h = real(ifft(H_herm));

        p = truncateCir(h).';
    case 2
        % D2-H1
        load('/Users/igorfreire/Documents/Lasse/gfast_simulator/Cables/D2-H1.mat');
        h = real(ifft(H));

        if (any(imag(ifft(H, N)) > 0))
            warning('H is not Hermitian symmetric');
        end

        p = truncateCir(h).';
    case 3
        p = [1e-5 1e-5 .91 -.3 .2 .09 .081 .0729];
end

if (length(p) > Nfft)
    warning('Pulse response longer than Nfft');
end

% Pulse response length
Lh = length(p);

% Matched-filter Bound
SNRmfb = (Ex_bar * norm(p).^2) / N0_over_2;
fprintf('SNRmfb:    \t %g dB\n\n', 10*log10(SNRmfb))

%% Windowing

if (windowing)
    dmtWindow = designDmtWindow(Nfft, nu, tau);
else
    % When windowing is not used, the suffix must be 0
    tau = 0;
end

% Windowing
dmtObj.windowing = windowing;
dmtObj.tau       = tau;
if (windowing)
    dmtObj.window = dmtWindow;
end

%% Cursor
% Except for when the TEQ is used, the cursor corresponds to the index of
% the peak in the CIR.

if (equalizer ~= EQ_TEQ)
    [~, iMax] = max(abs(p));
    n0 = iMax - 1;
end

%% Equalizers
% The TEQ is designed before loading, by assuming flat input spectrum.
% However, since bit loading alters the energy allocation among
% subchannels, the TEQ is redesigned after bit loading.

switch (equalizer)
    case EQ_TEQ
fprintf('\n-------------------- MMSE-TEQ Design ------------------- \n\n');

        if (nu >= (Lh - 1))
            error('MMSE-TEQ is unecessary. CP is already sufficient!');
        end

        % Search optimum length and delay for the TEQ design
        [nTaps, delta] = optimizeTeq(teqType, p, nu, L, N0_over_2, ...
            Ex_bar, Nfft, debug_teq);
        fprintf('Chosen Equalizer Length:\t %d\n', nTaps);
        fprintf('Chosen Delay:           \t %d\n', delta);

        % Design final TEQ
        Nf = floor(nTaps/L);
        switch (teqType)
            case TEQ_MMSE
                [w, SNRteq] = ...
                    mmse_teq(p, L, delta, Nf, nu, Ex_bar, N0_over_2, ...
                    debug_teq);
                fprintf('New SNRmfb (TEQ):\t %g dB\n', 10*log10(SNRteq))
            case TEQ_SSNR
                w = ssnr_teq(p, L, delta, Nf, nu, debug_teq);
        end

        if(~isreal(w))
            warning('MMSE-TEQ designed with complex taps');
        end

        % Shortening SNR:
        ssnr_w = ssnr( w, p, delta, nu );
        fprintf('SSNR:\t %g dB\n', 10*log10(ssnr_w));

        % Cursor based on the TEQ delay.
        n0 = delta;

    case EQ_FREQ_PREC
fprintf('\n------------------- Freq DMT Precoder ------------------ \n');
        FreqPrecoder = dmtFreqPrecoder(p, Nfft, nu, tau, n0, windowing);
    case EQ_TIME_PREC
fprintf('\n------------------- Time DMT Precoder ------------------ \n\n');
        TimePrecoder = dmtTimePrecoder(p, n0, nu, tau, Nfft, windowing);
end

% Copy cursor to DMT Object
dmtObj.n0 = n0;

% Equalization/Precoding
dmtObj.equalizer         = equalizer;
switch equalizer
    case EQ_TEQ
        dmtObj.w        = w;
        dmtObj.teqTaps  = nTaps;
        dmtObj.teqDelta = delta;
        dmtObj.teqNf    = Nf;
        dmtObj.teqType  = teqType;
    case EQ_FREQ_PREC
        dmtObj.Precoder = FreqPrecoder;
    case EQ_TIME_PREC
        dmtObj.Precoder = TimePrecoder;
end

%% Effective pulse response

switch (equalizer)
    case EQ_TEQ
        % New effective channel:
        p_eff = conv(p,w);
    otherwise
        p_eff = p;
end

%% Frequency Equalizer

FEQn = dmtFEQ(p_eff, dmtObj);

%% Gain-to-noise Ratio

[ gn ] = dmtGainToNoise(p_eff, dmtObj);

%% Water filling

fprintf('\n--------------------- Water Filling -------------------- \n\n');

% Water-filling:
[bn_bar_wf, En_bar_wf] = waterFilling(gn, Ex_budget, N, gap);

% Residual unallocated energy
fprintf('Unallocated energy:      \t  %g\n', Ex_budget - sum(En_bar_wf));

% Bits per subchannel
bn_wf = bn_bar_wf(1:N_subch) .* dim_per_subchannel;
% Number of bits per dimension
b_bar_wf = (1/nDim)*(sum(bn_bar_wf));
fprintf('b_bar:                  \t %g bits/dimension\n', b_bar_wf)
% For gap=0 and N->+infty, this should be the channel capacity per real
% dimension.

% Corresponding multi-channel SNR:
SNRdmt_wf = 10*log10(gap*(2^(2*b_bar_wf)-1));
% SNR at each tone, per dimension:
SNR_n_wf = En_bar_wf .* gn;
% Normalized SNR on each tone, per dimension (should approach the target
% gap):
SNR_n_norm_wf = SNR_n_wf ./ (2.^(2*bn_bar_wf) - 1);

fprintf('Multi-channel SNR (SNRdmt):\t %g dB\n', SNRdmt_wf)

if (equalizer == EQ_TEQ)
    fprintf('Note: shortened response was used for water-filling.\n');
end

%% Loading

fprintf('\n------------------ Discrete Loading -------------------- \n\n');

[bn, En, SNR_n, n_loaded] = dmtLoading(dmtObj, gn);

% Number of subchannels that are loaded
N_loaded = length(n_loaded);
% Dimensions in each loaded subchannel
dim_per_loaded_subchannel = dim_per_dft_tone(n_loaded);

% Bits per subchannel per dimension
bn_bar = bn ./ dim_per_subchannel;

%% Compare water-filling and discrete-loading
if (debug && debug_loading)
    figure
    plot(subCh_tone_index, bn_wf, ...
        'linewidth', 1.1)
    hold on
    plot(subCh_tone_index, bn, 'g')
    legend('Water-filling', 'Discrete Loading')
    xlabel('Subchannel');
    ylabel('Bits');
    grid on
    title('Bit loading')
end

%% Channel Capacity
% Channel capacity is computed considering the SNR that results from LC
% discrete loading.

fprintf('\n------------------ Channel Capacity -------------------- \n\n');

% Capacity per real dimension
cn_bar = 0.5 * log2(1 + SNR_n);
% Capacity per subchannel
cn = cn_bar .* dim_per_subchannel;
% Multi-channel capacity, per dimension:
c = sum(cn) / nDim;
% Note #1: for the capacity computation, all real dimensions are
% considered, including the overhead. See the example of (4.208)
% Note #2: the actual capacity is only obtained for N -> infty, so the
% above is only an approximation.

fprintf('capacity:               \t %g bits/dimension', c)
fprintf('\nBit rate:               \t %g mbps\n', c * Rsym * nDim /1e6);

%% Error Probability per dimension
% Nearest-neighbors union bound probability of error.

fprintf('\n----------------- Error Probabilities ------------------ \n\n');

% Levin-Campello:
Pe_bar_n = dmtPe(bn, SNR_n, dim_per_subchannel);

% NNUB Pe per dimension:
Pe_bar = mean(Pe_bar_n, 'omitnan');

fprintf('Approximate NNUB Pe per dimension:\n');
fprintf('For Discrete-load (LC)  :\t %g\n', Pe_bar);

%% Modulators

% Modulation order on each subchannel
modOrder = 2.^bn;

[modulator, demodulator] = dmtGenerateModems(modOrder, dim_per_subchannel);

%% Look-up table for each subchannel indicating the corresponding modem

modem_n = dmtModemLookUpTable(modOrder, dim_per_subchannel);

%% Energy loading (constellation scaling factors) and minimum distances

[scale_n, dmin_n] = dmtSubchanScaling(modulator, modem_n, ...
                    En, dim_per_subchannel);

%% Prepare DMT Struct for Monte-Carlo

% Mod/Demod Objects
dmtObj.modulator         = modulator;
dmtObj.demodulator       = demodulator;

% Bit-loading dependent parameters for each subchannel
dmtObj.b_bar_n           = bn_bar;          % Bit load per dimension
dmtObj.modem_n           = modem_n;         % Constellation scaling
dmtObj.scale_n           = scale_n;         % Subchannel scaling factor
dmtObj.dmin_n            = dmin_n;          % Minimum distance
dmtObj.FEQ_n             = FEQn;            % FEQ

%% Traing Loading based on modulated sequence

fprintf('\n------------------ Trained Loading -------------------- \n\n');

% Random DMT Data divided per subchannel
[tx_data] = dmtRndData(dmtObj);

% DMT Modulation
[u, x] = dmtTx(tx_data, dmtObj);

% Transmit Autocorrelation

% Input Autocorrelation based on actual transmit data
[rxx, ~] = xcorr(x(:), Nfft-1, 'unbiased');

[ dmtObj, bn, En, SNR_n, n_loaded, p_eff ] = dmtTrainining(p, ...
    dmtObj, rxx );

% Number of subchannels that are loaded
N_loaded = length(n_loaded);
% Dimensions in each loaded subchannel
dim_per_loaded_subchannel = dim_per_dft_tone(n_loaded);

% Update the probability of error
Pe_bar_n = dmtPe(bn, SNR_n, dim_per_subchannel);
Pe_bar = mean(Pe_bar_n, 'omitnan');

% Print the results of the new bit-load
fprintf('Pe_bar (LC)  :\t %g\n', Pe_bar);

%% Monte-carlo

fprintf('\n---------------------- Monte Carlo --------------------- \n\n');

% Preallocate
sym_err_n  = zeros(N_loaded, 1);
numErrs = 0; numDmtSym = 0;

% Sys Objects
BitError = comm.ErrorRate;

%% Iterative Transmissions

iTransmission = 0;
iReTraining   = 0;

while ((numErrs < maxNumErrs) && (numDmtSym < maxNumDmtSym))
    iTransmission = iTransmission + 1;

    %% Random DMT Data divided per subchannel
    [tx_data] = dmtRndData(dmtObj);

    %% DMT Modulation
    [u, x] = dmtTx(tx_data, dmtObj);

    %% Transmit Autocorrelation

    % Input Autocorrelation based on actual transmit data
    [rxx_current, ~] = xcorr(x(:), Nfft-1, 'unbiased');

    % Average the autocorrelation
    rxx = (rxx_current + rxx)/2;

    %% Debug Tx Energy

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
        % Design energy is the energy budget plus the excess energy in the
        % prefix
        fprintf('Design value:\t%g\t', Ex_budget*(Nfft + nu)/Nfft);
    end

    %% Channel
    [y] = dmtChannel(u, p, dmtObj);

    %% DMT Receiver
    [rx_data] = dmtRx(y, dmtObj);

    %% Error results

    % Symbol error count
    sym_err_n = sym_err_n + symerr(tx_data(bn > 0, :), ...
                                  rx_data(bn > 0, :), 'row-wise');
    % Symbol error rate per subchannel
    ser_n     = sym_err_n / (iTransmission * nSymbols);
    % Per-dimensional symbol error rate per subchannel
    ser_n_bar = ser_n ./ dim_per_loaded_subchannel.';

    % Preliminary results
    numErrs   = sum(sym_err_n);
    numDmtSym = iTransmission * nSymbols;

    fprintf('Pe_bar:\t%g\t', mean(ser_n_bar));
    % Note: consider only the loaded subchannels in the above
    fprintf('nErrors:\t%g\t', numErrs);
    fprintf('nDMTSymbols:\t%g\n', numDmtSym);

    %% Re-training of the bit-loading (and possibly the equalizer)
    % The initial bit-loading can often be innacurate, mostly due to the
    % ICPD that is initially computed assuming the input is perfectly
    % uncorrelated (equivalently, the energy load is flat). During
    % show-time, we can compute the actual correlation of the transmit
    % signal and use it to compute a more accurate ICPD PSD. Using this
    % PSD, in turn, we update the bit-load and restart the transmission. In
    % case the MMSE-LE TEQ is adopted, it can also be updated using the
    % actual estimated input autocorrelation Rxx.

    % If the error is too high, bit-loading shall be re-trained
    if ((mean(ser_n_bar) > 3 * Pe_bar) && ...
          (iTransmission * nSymbols) > 1e4)
        % Number of transmitted symbols is checked to avoid triggering
        % re-training when the SER has been measured for a short period

        % Keep track of how many times re-training was activated
        iReTraining = iReTraining + 1;

        fprintf('\n## Re-training the ICPD PSD and the bit-loading...\n');

        [ dmtObj, bn, En, SNR_n, n_loaded, p_eff ] = dmtTrainining(p, ...
            dmtObj, rxx );

        % Number of subchannels that are loaded
        N_loaded = length(n_loaded);
        % Dimensions in each loaded subchannel
        dim_per_loaded_subchannel = dim_per_dft_tone(n_loaded);

        % Update the probability of error
        Pe_bar_n = dmtPe(bn, SNR_n, dim_per_subchannel);
        Pe_bar = mean(Pe_bar_n, 'omitnan');

        % Print the results of the new bit-load
        fprintf('Pe_bar (LC)  :\t %g\n', Pe_bar);
        fprintf('## Restarting transmission...\n\n');

        % Finally, reset the SER computation:
        sym_err_n  = zeros(N_loaded, 1);
        numErrs = 0; numDmtSym = 0; iTransmission = 0;
    end

    %% Constellation plot for debugging
    if (debug && debug_constellation && modem_n(debug_tone) > 0 ...
        && iTransmission == 1)
        k = debug_tone;

        viewConstellation(Z, scale_n(k) * ...
                modulator{modem_n(k)}.modulate(0:modOrder(k) - 1), k);
    end

end

%% Results
fprintf('\n----------------------- Results ------------------------ \n\n');
Pe_bar_measured =  mean(ser_n_bar);
fprintf('Pe_bar:       \t %g\n', Pe_bar_measured);

if (debug && debug_Pe)
    figure
    semilogy(n_loaded, ser_n_bar, 's')
    hold on
    semilogy(subCh_tone_index, Pe_bar_n, 'r*')
    title('Measured SER per dimension vs. nominal $\bar{Pe}$', ...
        'Interpreter', 'latex')
    xlabel('Subchannel (n)')
    ylabel('$\bar{Pe}(n)$', 'Interpreter', 'latex')
    legend('Measured','Theoretical Approx')
    grid on
end
