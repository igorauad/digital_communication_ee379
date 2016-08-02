function [ Hisi, Hici, Hcirc ] = dmtIsiIciMatrices(p, n0, nu, tau, N, preCursor, windowing)
% Compute the time-domain ISI anc ICI matrices
% ---------------------------------------------
% [ Hisi, Hici, Hcirc ] = dmtIsiIciMatrices(h, n0, nu, tau, N, ...
%                                           preCursor, windowing)
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

%% Post-cursor ICPD

% Preallocate
Hici = zeros(N,N); % Post-cursor ICI Matrix
Hisi = zeros(N,N); % Post-cursor ISI Matrix

% Check whether post-cursor ICPD effectively occurs
assert(delta > 0, 'Post-cursor ICPD does not occur');

% Window
dmtWindow = designDmtWindow(N, nu, tau);

% Post-cursor ISI core matrix
Ht = toeplitz([p(end) zeros(1, L-(nu-tau)-n0-1)], ...
    flipud(p((nu - tau + n0 + 2):end)));

% Window
W_w = diag(dmtWindow(end-delta+1:end));

% Post-cursor ISI Matrix
if (n0 >= tau)
    Hisi( 1:(delta), (N - L + nu + 1):(N + tau - n0)) = Ht*W_w;
else
    % When n0 < tau, the matrix has to be circularly shifted
    shiftBy = tau - n0;
    Hisi(1:(delta), (N - L + nu + 1):(N + tau - n0) ) = Ht*W_w;
    Hisi = Hisi(:,shiftBy+1:end);
    Hisi = circshift(Hisi, [0, shiftBy]);
end

% Post-cursor ICI Matrix
Hici(1:(delta),(N - L + 1):(N -nu + tau - n0)) = -Ht*W_w;

%% Pre-cursor ICI

% Check Pre-Cursor ICPD Energy

% Check the energy prior to the cursor
preCursorEnergy = sum(abs(p(1:(n0-1))).^2);
if (preCursorEnergy >= 0.1 * norm(p)^2)
    warning('Pre-cursor ICPD is not irrelevant\n');
    assert(delta > 0, ...
        'Post-cursor ICPD does not occur');
end

if (preCursor) % Pre-cursor ICI is also mitigated
    % Preallocate
    HpreIci = zeros(N, N);

    H_tilde_t = toeplitz(p(1:n0), [p(1) zeros(1, n0 - 1)]);

    W_tilde_w = diag(dmtWindow(1:n0));

    H_tilde_t_windowed = H_tilde_t * W_tilde_w;

    HpreIci(end-n0+1:end,end-n0+1:end) = - H_tilde_t_windowed;

    % Add pre-cursor ICI matrix to the post-cursor ICI matrix:
    Hici = Hici + HpreIci;
end

end

