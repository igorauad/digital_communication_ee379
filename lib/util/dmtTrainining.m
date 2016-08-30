function [ Dmt, bn, SNR_n, Rb, n_loaded ] = dmtTrainining(p, Dmt, rxx )
% Training of DMT FEQ and Bit Loading
%
% Inputs
%  p     -> Original channel pulse response
%  dmt   -> DMT Object
%  rxx   -> Tx Unbiased Autocorrelation Vector

%% Constants
% TEQ criterion
TEQ_MMSE    = 0;

% Equalizer Types
EQ_TEQ       = 1;

%% Parameters
Nfft               = Dmt.Nfft;          % FFT size
equalizer          = Dmt.equalizer;     % Equalizer Type
teqType            = Dmt.teqType;       % Type of TEQ
L                  = Dmt.L;             % Oversampling ratio
nu                 = Dmt.nu;            % Cyclic prefix length
N0_over_2          = Dmt.N0_over_2;     % Avg Noise Energy per dim

if (equalizer  == EQ_TEQ)
    nTaps = Dmt.teqTaps;
    delta = Dmt.teqDelta;
    Nf    = Dmt.teqNf;
    w     = Dmt.w;
end

%% Effective channel pulse response

switch (equalizer)
    case EQ_TEQ
        % New effective channel:
        p_eff = conv(p, w);
    otherwise
        p_eff = p;
end

% Pulse response length
Lh = length(p);

%% For an MMSE_TEQ, jointly design the TEQ
if (equalizer == EQ_TEQ && teqType == TEQ_MMSE)

    % Autocorrelation matrix with the appropriate dimensions
    Rxx = toeplitz(rxx(Nfft:Nfft + floor(nTaps/L)*L + (Lh-1) - 1));

    % Re-design the TEQ
    [w, SNRteq] = ...
        mmse_teq(p, L, delta, Nf, nu, Rxx, N0_over_2);
    fprintf('New SNRmfb (TEQ):\t %g dB\n', 10*log10(SNRteq))

    % Shortening SNR:
    ssnr_w = ssnr( w, p, delta, nu );
    fprintf('SSNR:\t %g dB\n', 10*log10(ssnr_w));
    if(~isreal(w))
        warning('MMSE-TEQ designed with complex taps');
    end

    % Update the effective channel:
    p_eff = conv(p,w);

    % Update the FEQ
    FEQn = dmtFEQ(p_eff, Dmt);

    % Update equalizers in the DMT Object
    Dmt.w      = w;    % TEQ
    Dmt.FEQ_n  = FEQn; % FEQ
end

%% Nfft x Nfft Autocorrelation Matrix
Rxx = toeplitz(rxx(Nfft:end));

%% Update the gain-to-noise ratio:
gn = dmtGainToNoise(p_eff, Dmt, Rxx);

%% Re-compute the bit loading
[bn, En, SNR_n, Rb, n_loaded] = dmtLoading(Dmt, gn);

% Bits per subchannel per dimension
bn_bar = bn ./ Dmt.dim_per_subchannel;

%% Update the vector of modulation orders
modOrder = 2.^bn;
% Update modem objects
[modulator, demodulator] = dmtGenerateModems(modOrder, ...
    Dmt.dim_per_subchannel);
% Re-generate modem look-up table
modem_n = dmtModemLookUpTable(modOrder, Dmt.dim_per_subchannel);
% Re-generate the subchannel scaling factors
[scale_n, dmin_n] = dmtSubchanScaling(modulator, modem_n, ...
    En, Dmt.dim_per_subchannel);

%% Update DMT Object
Dmt.modulator    = modulator;
Dmt.demodulator  = demodulator;
Dmt.b_bar_n      = bn_bar;        % Bit load per dimension
Dmt.modem_n      = modem_n;       % Constellation scaling
Dmt.scale_n      = scale_n;       % Subchannel scaling factor
Dmt.dmin_n       = dmin_n;        % Minimum distance

end

