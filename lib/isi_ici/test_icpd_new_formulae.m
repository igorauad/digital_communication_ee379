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
h = [0 0.1 1 0.9 0.3 0.2 0.1 0.1];
h = h(:);

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

h_padded         = [h;  zeros(Nfft - Lh, 1)]; % Zero-padded CIR to length N
h_padded_shifted = circshift(h_padded, -n0);  % Circular-shifted by -n0

% Circulant matrix:
Hcirc = toeplitz(h_padded_shifted, ...
    flipud(circshift(h_padded_shifted,-1)));

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

y_ideal = Hcirc * x;

%% Post-cursor

% Preallocate ISI/ICI Matrices
Hisi = zeros(Nfft,Nfft);
Hici = zeros(Nfft,Nfft);

% Core Convolution (Toeplitz) Matrix
Ht = toeplitz([h(end) zeros(1, L-(nu-tau)-n0-1)],...
               flipud(h((nu - tau + n0 + 2):end)));

% Window Diagonal Matrix
W_w = diag(dmtWindow(end-delta+1:end));

% Windowed Core Matrix
Ht_windowed = Ht * W_w;
% IMPORTANT: the above works both when windowing is used and when it is not
% used. The only condition is that tau must be set to 0 when not using
% windowing, first because it is indeed 0 in thia case, second because W_t
% becomes an identity.

% Post-cursor ISI Matrix
Hisi(1:delta, (Nfft - delta + 1):Nfft ) = Ht_windowed;
% Row indexes: reason is that the first delta samples are affected by ISI
%
% Column indexes: reason is that the multiplied x column-vectors will be
% rotated by tau upwards, which will leave the last samples of the vector
% being the samples from the overlapped suffix, which are exactly the
% "first" to cause ISI. So the convolution matrix must multiply the last
% delta elements of the vectors.

% Post-cursor ICI Matrix
Hici(1:delta, (Nfft - L + n0 + 1):(Nfft -nu + tau)) = -Ht_windowed;
% Row indexes: reason is that the first delta samples are affected by ICI
%
% Column indexes: reason is more complex, related to the "missing samples".
% The first column of the conv matrix has x[n0] in its first row and (n0 +
% nu) samples available in the rows beneath. However, tau of these samples
% will be naturally affected by ICI due to windowing. Therefore, to
% complete the L + 1 rows with non-zero values, the missing samples will be
% from index (n0 - L) mod-N, equivalent to x[N + (n0 -L)] (because n0 < L),
% up to x[N -(nu - tau) - 1], where N -(nu - tau) - 1 < N (namely the index
% will be within [0,N-1]) because (nu > tau).

% Post-cursor ISI
y_postIsi = Hisi * circshift(x(:,1:end-1), -tau);
y_postIsi = [zeros(Nfft,1), y_postIsi]; % First symbol does not suffer ISI

% Post-cursor ICI
y_postIci = Hici * x;

disp('Test removing post-cursor ICPD:')
y - (y_ideal + y_postIci + y_postIsi)

%% Pre-cursor

% Preallocate Pre-cursor ISI/ICI Matrices
HpreIci = zeros(Nfft, Nfft);
HpreIsi = zeros(Nfft, Nfft);

if (n0 > 0)
    % Core Convolution (Toeplitz) Matrix
    H_tilde_t = toeplitz(h(1:n0), [h(1) zeros(1, n0 - 1)]);

    % Window Matrix
    W_tilde_w = diag(dmtWindow(1:n0));

    % Windowed Core Matrix
    H_tilde_t_windowed = H_tilde_t * W_tilde_w;

    % ISI Matrix
    HpreIsi(Nfft-n0+1:end,1:n0) = H_tilde_t_windowed;
    % Row indexes: last n0 samples are affected by pre-cursor ISI
    %
    % Column indexes: the samples x[n] causing pre-cursor ISI always start
    % at the CP of the succeding DMT symbol. Since we circularly shift the
    % symbol vectors by nu downwards, the first sample of the vectors will
    % become x[N-nu], exactly the first of the sequence that causes
    % pre-cursor ISI. Therefore, we need to place the convolution matrix at
    % the first n0 columns, which will span samples x[N-nu] to
    % x[((n0 -nu -1))mod-N].

    %ICI Matrix
    HpreIci = -HpreIsi;
    % Row indexes: last n0 samples are affected by pre-cursor ICI
    %
    % Column indexes: the samples x[n] that will be missing in the
    % convolution summation at the last n0 outputs y[n] are always the
    % first n0 samples x[0] to x[n0 - 1].

end

% Pre-cursor ICI
y_preIci = HpreIci * x;

% Pre-cursor ISI
%Multiply by "future" symbols:
y_preIsi = HpreIsi * circshift(x_withLastSymbol(:, 2:end), nu);

disp('Test removing pre-cursor ICPD:')
y - (y_ideal + y_preIci + y_preIsi)

%% Complete ICPD Removal
disp('Test removing post and pre-cursor ICPD:')
y - (y_ideal + y_postIci + y_postIsi + y_preIci + y_preIsi)

%% Test DMT ISI/ICI Matrix function

fprintf('\n---------------------------------------------------------\n');
fprintf('Testing Matrices designed by \"dmtIsiIciMatrices():\"\n\n');

[ Hisi2, Hici2, Hcirc2, HpreIsi2, HpreIci2 ] = ...
    dmtIsiIciMatrices(h, n0, nu, tau, Nfft, windowing);

if (any(Hisi(:) ~= Hisi2(:)))
    fprintf('Hisi\t Error (different)\n');
else
    fprintf('Hisi\t Passed (identical)\n');
end

if (any(Hici(:) ~= Hici2(:)))
    fprintf('Hici\t Error (different)\n');
else
    fprintf('Hici\t Passed (identical)\n');
end

if (any(Hcirc(:) ~= Hcirc2(:)))
    fprintf('Hcirc\t Error (different)\n');
else
    fprintf('Hcirc\t Passed (identical)\n');
end

if (any(HpreIsi(:) ~= HpreIsi2(:)))
    fprintf('HpreIsi\t Error (different)\n');
else
    fprintf('HpreIsi\t Passed (identical)\n');
end

if (any(HpreIci(:) ~= HpreIci2(:)))
    fprintf('HpreIci\t Error (different)\n');
else
    fprintf('HpreIci\t Passed (identical)\n');
end
