clearvars, clc
addpath(genpath('../'))

%% Parameters
Nfft      = 4;
nu        = 2;
tau       = 2;
n0        = 1;
nSymbols  = 2;
windowing = 0;
% CIR:
h = [0 0.1 1 0.9 0.3 0.2 0.1].';

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
dmtWindow = designDmtWindow(Nfft, nu, tau);

if (delta <= 0)
    warning('Parameters yield absent post-cursor ICPD already');
end
%% Circulant channel matrix

h_padded = [h;  zeros(Nfft - Lh, 1)];
Hcirc = toeplitz(h_padded, ...
    flipud(circshift(h_padded,-1)));

%% Transmission

x = rand(Nfft, nSymbols + 1);

if (windowing)
    x_ce = [x(end -nu + 1:end,:); x; x(1:tau,:)];

    x_ce = windowAndOverlap(x_ce, dmtWindow, Nfft, nu, tau);
else
    x_ce = [x(end -nu + 1:end,:); x];
end

x_stream = x_ce(:);

y_stream = conv(x_stream, h);

y_stream = y_stream(rxStart:rxEnd);

y_ce = reshape(y_stream, symbolSize, nSymbols);

y = y_ce((nu + 1):(nu + Nfft), :);

%% Test Formulas

% Throw away last symbol in x and y:
x_withLastSymbol = x;
x = x(:, 1:end-1);

y_ideal = Hcirc * circshift(x, -n0);

%% Post-cursor

% Preallocate ISI/ICI Matrices
Hisi = zeros(Nfft,Nfft);
Hici = zeros(Nfft,Nfft);

% Window
W_t = diag(dmtWindow(end-delta+1:end));

% ISI Core Matrix
Ht = toeplitz([h(end) zeros(1, L-(nu-tau)-n0-1)], flipud(h((nu - tau + n0 + 2):end)));
Ht_windowed = Ht * W_t;
% IMPORTANT: the above works both when windowing is used and when it is not
% used. The only condition is that tau must be set to 0 when not using
% windowing, first because it is really 0, second because W_t becomes an
% identity.

% Post-cursor ISI
Hisi(1:(delta), (Nfft - L + nu + 1):(Nfft + tau - n0) ) = Ht_windowed;
if (n0 < tau)
    % When n0 < tau, the matrix has to be circularly shifted
    shiftBy = tau - n0;
    Hisi = Hisi(:,shiftBy+1:end);
    Hisi = circshift(Hisi, [0, shiftBy]);
end

% Post-cursor ICI
Hici(1:(delta), (Nfft - L + 1):(Nfft -nu + tau - n0)) = -Ht_windowed;

y_postIci = Hici * circshift(x, -n0);
y_postIsi = Hisi * circshift(x(:,1:end-1), -n0);
y_postIsi = [zeros(Nfft,1), y_postIsi];

disp('Test removing post-cursor ICPD:')
y - (y_ideal + y_postIci + y_postIsi)

%% Pre-cursor

% Common matrix
H_tilde_t = toeplitz(h(1:n0), [h(1) zeros(1, n0 - 1)]);

W_tilde_t = diag(dmtWindow(1:n0));

H_tilde_t_windowed = H_tilde_t * W_tilde_t;

% Pre-cursor ICI

%ICI Matrix
HpreIci = zeros(Nfft, Nfft);
HpreIci(end-n0+1:end,end-n0+1:end) = - H_tilde_t_windowed;

y_preIci = HpreIci * circshift(x, -n0);

% Pre-cursor ISI
HpreIsi = zeros(Nfft, Nfft);
HpreIsi(end-n0+1:end,(Nfft - nu - n0 + 1):(Nfft - nu)) = H_tilde_t_windowed;

% Multiply by "future" symbols:
y_preIsi = HpreIsi * circshift(x_withLastSymbol(:, 2:end), -n0);

disp('Test removing pre-cursor ICPD:')
y - (y_ideal + y_preIci + y_preIsi)

%% Complete ICPD Removal
disp('Test removing post and pre-cursor ICPD:')
y - (y_ideal + y_postIci + y_postIsi + y_preIci + y_preIsi)
