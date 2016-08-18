function [ Hisi, Hici, Hcirc, HpreIsi, HpreIci ] = ...
    dmtIsiIciMatrices(p, n0, nu, tau, N, windowing)
% Compute the time-domain ISI anc ICI matrices
% ---------------------------------------------
% [ Hisi, Hici, Hcirc ] = dmtIsiIciMatrices(h, n0, nu, tau, N, ...
%                                           preCursor, windowing)
%
%   Computes the ISI and ICI matrices such that the post-cursor ISI and ICI
%   are given by:
%       y_isi     = Hisi * circshift(x, -tau),
%       y_postIci = Hici * x,
%   and the pre-cursor ISI/ICI are given by:
%       y_tilde_ici = HpreIci * x;
%       y_tilde_isi = HpreIsi * circshift(x, nu);
%
%   IMPORTANT: note that the circular shifts applied to x for post-cursor
%   and pre-cursor ISI are assumed to be different.
%
%   Input Parameters
% p                 Channel Pulse Response
% n0                CIR delay (index of its peak)
% nu                Cyclic Prefix Length
% tau               Cyclic Suffix length
% N                 FFT size
% preCursor         When pre-cursor ICPD should be accounted
% windowing         When spectrum shaping windowing should be accounted
%
% Outpu:
%   Hisi   ->  Time-domain ISI matrix
%   Hici   ->  Time-domain ICI matrix
%   Hcirc  ->  Circulant matrix for the Ideal channel

%% Verification of Outputs

% Default
computeHcirc   = (nargout > 2);
computeHpreIsi = (nargout > 3);
computeHpreIci = (nargout > 4);

%% Initialization

% Force column vector:
p = p(:);

% Force tau to 0 when windowing is not used and all derivations will
% continue to hold.
if (~windowing)
    tau = 0;
end

% Channel length and dispersion, respectively:
Lh = length(p);
L  = Lh-1;

% Number of samples affected by post-cursor ICPD
delta = L - (nu - tau) - n0;

%% Circulant Matrix

% CIR zero padded to a vector of length N
h_padded = [p;  zeros(N - Lh, 1)];

% Note h has to be a column vector
Hcirc = toeplitz(h_padded, flipud(circshift(h_padded, -1)));

%% Windowing Sequence

% Window
dmtWindow = designDmtWindow(N, nu, tau);

%% Post-cursor ICPD

% Preallocate
Hici = zeros(N,N); % Post-cursor ICI Matrix
Hisi = zeros(N,N); % Post-cursor ISI Matrix

% Check whether post-cursor ICPD effectively occurs
if (delta <= 0)
    warning('Post-cursor ICPD does not occur');
else

    % Core Convolution (Toeplitz) Matrix
    Ht = toeplitz([p(end) zeros(1, delta-1)], ...
        flipud(p((nu - tau + n0 + 2):end)));

    % Window Diagonal Matrix
    W_w = diag(dmtWindow(end-delta+1:end));

    % Windowed Core Matrix
    Ht_windowed = Ht * W_w;

    % Post-cursor ISI Matrix
    Hisi(1:delta, (N - delta + 1):N ) = Ht_windowed;
    % It is assumed that the column-vector DMT symbols will be circularly
    % shifted by -tau (upwards) before multiplying Hisi.

    % Post-cursor ICI Matrix
    Hici(1:delta, (N - L + n0 + 1):(N -nu + tau)) = -Ht_windowed;
    % It is assumed that the DMT symbols that multiply Hici are not
    % shifted.

end
%% Pre-cursor ICI

% Check Pre-Cursor ICPD Energy

% Check the energy prior to the cursor
preCursorEnergy = sum(abs(p(1:(n0-1))).^2);
if (preCursorEnergy >= 0.1 * norm(p)^2)
    warning('Pre-cursor ICPD is not irrelevant\n');
    assert(delta > 0, ...
        'Post-cursor ICPD does not occur');
end

% Preallocate
HpreIci = zeros(N, N);
HpreIsi = zeros(N, N);

% Check whether pre-cursor ICPD effectively occurs
if (n0 > 0)
    % Core Convolution (Toeplitz) Matrix
    H_tilde_t = toeplitz(p(1:n0), [p(1) zeros(1, n0 - 1)]);

    % Window Matrix
    W_tilde_w = diag(dmtWindow(1:n0));

    % Windowed Core Matrix
    H_tilde_t_windowed = H_tilde_t * W_tilde_w;

    if (computeHpreIsi) % Pre-cursor ICI is also mitigated
        % ISI Matrix
        HpreIsi(N-n0+1:end,1:n0) = H_tilde_t_windowed;
        % It is assumed that the column-vector DMT symbols will be circularly
        % shifted by nu (downwards) before multiplying HpreIsi.
    end

    if (computeHpreIci)
        % ICI Matrix
        HpreIci = -HpreIsi;
        % It is assumed that the DMT symbols that multiply HpreIci are not
        % shifted.
    end
end

end

