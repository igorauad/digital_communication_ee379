clear all
clc

%% Parameters
Nfft = 4;
Lcp = 2;
Lcs = 2;
n0 = 1;
nSymbols = 2;

h = [0 1 0.9 0.3].';

%% Derived parameters
symbolSize = Nfft + Lcp;
Lh = length(h);
L = Lh - 1;
windowStart = (n0 + 1);
windowEnd   = nSymbols*symbolSize + windowStart - 1;

% Window:
dmtWindow = designDmtWindow(Nfft, Lcp, Lcs);

nAffected = L - (Lcp - Lcs) - n0; % Number of affected samples

%% Transmission

x = rand(Nfft, nSymbols + 1);

x_ce = [x(end -Lcp + 1:end,:); x; x(1:Lcs,:)];

x_ce = windowAndOverlap(x_ce, dmtWindow, Nfft, Lcp, Lcs);

x_stream = x_ce(:);

y_stream = conv(x_stream, h);

y_stream = y_stream(windowStart:windowEnd);

y_ce = reshape(y_stream, symbolSize, nSymbols);

y = y_ce((Lcp + 1):(Lcp + Nfft), :);

h_padded = [h;  zeros(Nfft - Lh, 1)];
Hcirc = toeplitz(h_padded, ...
    flipud(circshift(h_padded,-1)));

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
W_t = diag(dmtWindow(end-nAffected+1:end));

Ht = toeplitz([h(end) zeros(1, L-(Lcp-Lcs)-n0-1)], flipud(h((Lcp - Lcs + n0 + 2):end)));

Ht_windowed = Ht * W_t;

% Post-cursor ISI
if (n0 >= Lcs)
    Hisi(1:(nAffected), (Nfft - L + Lcp + 1):(Nfft + Lcs - n0) ) = Ht_windowed;
else
    shiftBy = Lcs - n0;
    Hisi(1:(nAffected), (Nfft - L + Lcp + 1):(Nfft + Lcs - n0) ) = Ht_windowed;
    Hisi = Hisi(:,shiftBy+1:end);
    Hisi = circshift(Hisi, [0, shiftBy]);
end

% Post-cursor ICI
Hici(1:(nAffected), (Nfft - L + 1):(Nfft -Lcp + Lcs - n0)) = -Ht_windowed;

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
HpreIsi(end-n0+1:end,(Nfft - Lcp - n0 + 1):(Nfft - Lcp)) = H_tilde_t_windowed;

% Multiply by "future" symbols:
y_preIsi = HpreIsi * circshift(x_withLastSymbol(:, 2:end), -n0);

disp('Test removing pre-cursor ICPD:')
y - (y_ideal + y_preIci + y_preIsi)

disp('Test removing post and pre-cursor ICPD:')
y - (y_ideal + y_postIci + y_postIsi + y_preIci + y_preIsi)
