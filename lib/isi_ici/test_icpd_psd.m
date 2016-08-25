clearvars, clc, close all
addpath(genpath('../'))

%% Parameters
M         = 16; % Constellation order
Nfft      = 128;
Nfft_psd  = 512;
nu        = 16;
tau       = 8;
n0        = 12;
nSymbols  = 1e3;
windowing = 1;
flat      = 0;
% CIR:
h = [randn(1,10) 0.01 0.1 1 exp(-0.15*(1:50))].';

%% Derived parameters
symbolSize  = Nfft + nu;
Lh          = length(h);
L           = Lh - 1;

% Start/end index boundaries after symbol timing recovery at the receiver:
rxStart = (n0 + 1);
rxEnd   = nSymbols*symbolSize + rxStart - 1;

% Number of samples affected by post-cursor ICPD:
if (windowing)
    delta = L - (nu - tau) - n0
else
    delta = L - nu - n0
    % Force "tau" (CS length) to 0 when windowing is not used
    tau = 0;
end

% Window applied to CP and CS:
w = designDmtWindow(Nfft, nu, tau);

if (delta <= 0)
    warning('Parameters yield absent post-cursor ICPD already');
end

%% Random energy loading for each subchannel

% Average energy per dimension
if (flat)
    En = [0 1*ones(1, Nfft/2 - 1) 0]; % Subchannel energy
else
    % Note a large portion of the subchannels is purposely not loaded
    % (energy set to zero). This disturbs the correlation matrix with
    % respect to the assumption that it is equivalent to Ex_bar *
    % eye(Nfft). In this case, the ICPD PSD computed using the measured Rxx
    % (the one using the icpdPsdMtx function) tends to be more accurate.
    En = [0 10*randi(100, 1, Nfft/8) ...
        zeros(1, 3*Nfft/8 - 1) 0]; % Subchannel energy
end
scale_qam = zeros(Nfft, 1);
for i = 1:length(En)
    if (En(i) > 0)
        scale_qam(i) = modnorm(qammod(0:M-1, M), 'avpow',...
            En(i)/2); % half for each dimension
    end
end
scale_qam(Nfft/2 + 2:Nfft, :) = flipud( scale_qam(2:Nfft/2) );

%% Energy computations

% Total Energy
% The total energy must be the total loaded energy plus the amount of
% extra energy due to the CP.
Ex = sum(En)*(1 + nu/Nfft);

% Energy Per dimension
Ex_bar = Ex/(Nfft + nu);

%% ISI/ICI Matrices + Circulant Matrix

[ Hisi, Hici, Hcirc, HpreIsi, HpreIci ] = ...
    dmtIsiIciMatrices(h, n0, nu, tau, Nfft, windowing);

%% Transmission with random energy load

% Random data
tx_data = randi(M, Nfft/2 - 1, nSymbols + 1) - 1;

% Preallocate DMT Symbols
X       = zeros(Nfft, nSymbols + 1);

% Hermitian-symmetric random symbols:
X(2:Nfft/2, :)        = qammod(tx_data, M);
X(Nfft/2 + 2:Nfft, :) = flipud( conj( X(2:Nfft/2, :) ) );

X = diag(scale_qam) * X;

% IFFT
x = sqrt(Nfft) * ifft(X, Nfft);

if (windowing)
    x_ce = [x(end -nu + 1:end,:); x; x(1:tau,:)];

    x_ce = windowAndOverlap(x_ce, w, Nfft, nu, tau);
else
    x_ce = [x(end -nu + 1:end,:); x];
end

x_stream = x_ce(:);

%% Received sequence

y_stream = conv(x_stream, h);

y_stream = y_stream(rxStart:rxEnd);

y_ce = reshape(y_stream, symbolSize, nSymbols);

y = y_ce((nu + 1):(nu + Nfft), :);

%% Test Formulas

% Throw away last symbol in x and y:
x_withLastSymbol = x;
x = x(:, 1:end-1);

y_ideal = Hcirc * x;

%% Post-cursor

y_postIci = Hici * x;
y_postIsi = Hisi * circshift(x(:,1:end-1), -tau);
y_postIsi = [zeros(Nfft,1), y_postIsi];

y_post_icpd = y_postIci + y_postIsi;

%% Pre-cursor

% Pre-cursor ICI
y_preIci = HpreIci * x;

% Pre-cursor ISI
%Multiply by "future" symbols:
y_preIsi = HpreIsi * circshift(x_withLastSymbol(:, 2:end), nu);

y_pre_icpd = y_preIci + y_preIsi;

%% Total ICPD
y_icpd = y - y_ideal;

%% Ex_bar approximation

display('Measured Ex_bar')
mean(abs(x_ce(:)).^2)
display('Nominal Ex_bar')
Ex_bar

%% Measured PSDs
% Frequency is not explicitly given, so it is assumed to be 1. PSD for
% two-sided representation is divided by 2*pi, because of the axis.

% Post-cursor ICI
[Pici_ici, W]       = pwelch(y_postIci(:), ...
    Nfft_psd, Nfft_psd/8, Nfft_psd, 'twosided');
% Post-cursor ISI
[Pisi_isi, ~]       = pwelch(y_postIsi(:), Nfft_psd, Nfft_psd/8, W);

% Post-cursor ICPD
[Ppost_post, ~]     = pwelch(y_post_icpd(:), Nfft_psd, Nfft_psd/8, W);

% Pre-cursor ICI
[Ppreici_preici, ~] = pwelch(y_preIci(:), Nfft_psd, Nfft_psd/8, W);

% Pre-cursor ISI
[Ppreisi_preisi, ~] = pwelch(y_preIsi(:), Nfft_psd, Nfft_psd/8, W);

% Pre-cursor ICPD
[Ppre_pre, ~]       = pwelch(y_pre_icpd(:), Nfft_psd, Nfft_psd/8, W);

% Total ICPD
[Picpd_icpd, ~]     = pwelch(y_icpd(:), Nfft_psd, Nfft_psd/8, W);

% Transmit signal
[Pxx, ~]            = pwelch(x(:), Nfft_psd, Nfft_psd/8, W);

%% Computed PSDs

[ S_icpd, S_post, S_pre ] = ...
    icpdPsd(h, Nfft, Nfft_psd, nu, tau, n0, Ex_bar, w);

%% Compare ICPD PSD computed using matrix formulation
% IMPORTANT: the third argument of "icpdPsdMtx" can be a scalar (Ex_bar) or
% a matrix (Rxx). When Ex_bar is passed, the function assumes the input is
% uncorrelated and has Rxx = Ex_bar * eye(Nfft). Our goal here is to
% contrast the ICPD PSD computed based on the assumption of
% uncorrelatedness and the PSD computed with the correlation matrix Rxx
% that is computed based on the transmit data.

% Measured Autocorrelation
[r, l] = xcorr(x(:), Nfft_psd-1, 'unbiased');

% Measured Autocorrelation matrix
Rxx = toeplitz(r(Nfft_psd:end));
[ S_icpd2, S_post2, S_pre2 ] = icpdPsdMtx(Hisi, HpreIsi, Rxx, Nfft_psd);

figure
plot(10*log10(abs(S_icpd)))
hold on
plot(10*log10(abs(S_icpd2)), '--r')
legend('Using Summation', 'Using Matrices')
title('PSD Computation Using Matrices')

%% Transmit PSDs

L         = Nfft_psd / Nfft; % Oversampling
En_herm   = [En fliplr(conj(En(2:end-1)))];
En_interp = filter(ones(L,1), 1, upsample(En_herm, L));
En_bar    = En_interp/2; % Per dimension

figure
plot(W, 10*log10(En_bar), 'g', 'linewidth', 1.1)
hold on
plot(W, 10*log10(Pxx*2*pi), 'b--')
legend('Computed', 'Measured')
title('Transmit Signal PSD')

%% Plot Post-cursor PSDs

figure
plot(W, 10*log10(S_post/2), 'g', 'linewidth', 1.1)
hold on
plot(W, 10*log10(S_post2/2), 'y', 'linewidth', 1.1)
hold on
plot(W, 10*log10(Pici_ici*2*pi), 'b--')
hold on
plot(W, 10*log10(Pisi_isi*2*pi), 'r--')
legend('Computed ISI/ICI', 'Computed ISI/ICI (w/ Rxx)',...
    'Measured ICI', 'Measured ISI')
title('Post-cursor ICI/ISI')

figure
plot(W, 10*log10(S_post), 'g', 'linewidth', 1.1)
hold on
plot(W, 10*log10(S_post2), 'y', 'linewidth', 1.1)
hold on
plot(W, 10*log10(Ppost_post*2*pi), 'b--')
legend('Computed', 'Computed (w/ Rxx)', 'Measured')
title('Total Post-cursor ICPD')

%% Plot Pre-cursor PSDs

figure
plot(W, 10*log10(S_pre/2), 'g', 'linewidth', 1.1)
hold on
plot(W, 10*log10(S_pre2/2), 'y', 'linewidth', 1.1)
hold on
plot(W, 10*log10(Ppreici_preici*2*pi), 'b--')
hold on
plot(W, 10*log10(Ppreisi_preisi*2*pi), 'r--')
legend('Computed ISI/ICI', 'Computed ISI/ICI (w/ Rxx)',...
    'Measured ICI', 'Measured ISI')
title('Pre-cursor ICI/ISI')

figure
plot(W, 10*log10(S_pre), 'g', 'linewidth', 1.1)
hold on
plot(W, 10*log10(S_pre2), 'y', 'linewidth', 1.1)
hold on
plot(W, 10*log10(Ppre_pre*2*pi), 'b--')
legend('Computed', 'Computed (w/ Rxx)', 'Measured')
title('Total Pre-cursor ICPD')

%% Plot Total ICPD PSD

% Old computation
POST_PRE_ICPD_FLAG = 1;
icpd_psd = interferencePsd(h(:), Nfft_psd, nu, tau, n0, ...
    Ex, POST_PRE_ICPD_FLAG, windowing);

figure
plot(W, 10*log10(icpd_psd), 'ks')
hold on
plot(W, 10*log10(S_icpd), 'g', 'linewidth', 1.1)
hold on
plot(W, 10*log10(S_icpd2), 'y', 'linewidth', 1.1)
hold on
plot(W, 10*log10(Picpd_icpd*2*pi), 'b--')
legend('Old Computation', 'New Computation', ...
    'New Computation (w/ Rxx)', 'Measured')
title('Total ICPD')