function [ S_icpd, S_post, S_pre ] = ...
    icpdPsd(h, Nfft, Nfft_psd, nu, tau, n0, Ex_bar, w )
% Insuficient Cyclic Prefix Distortion PSD
% ---------------------------------------------------
%   [ S_icpd, S_post, S_pre ] = icpdPsd(h, Nfft, Nfft_psd, nu, tau,
%                                       n0, Ex_bar, w )
%
%
%   Input Arguments:
%       w   ->  Window applied to transmit DMT symbols. It can be either a
%       vector or a scalar. In case it is a vector, w is the actual window
%       vector. In case w is a scalar, it indicates whether windowing
%       should be applied or not, so that the window is designed internally
%       when applicable.
%
%   Notes: 1) "Nfft" corresponds to the FFT size used in the DMT system,
%   whereas Nfft_psd is the one used for the spectral density estimation.

%% Initialization

% Channel Memory
L = length(h) - 1;

% Number of samples affected by post-cursor ICPD:
delta = L - (nu - tau) - n0;

% Define a unitary window when not defined
if (nargin < 8)
    w = ones(Nfft + nu + tau, 1);
elseif (numel(w) == 1)
    w = designDmtWindow(Nfft, nu, tau);
end

%% Post-cursor PSD

% Preallocate
S_isi = zeros(Nfft_psd, 1); % Post-cursor ISI PSD

% The noise can be interepreted as the sum of all CIR "tails", such as
% [h(end)], [h(end-1:end)], [h(end-2:end)] and so forth, up to
% [h(end-delta+1:end)];
for iRow = 1:delta
    % CIR Tail
    h_tail = h(end -delta + iRow:end);
    % Accumulate the corresponding power spectral lines:
    S_isi  = S_isi + ...
        abs(w(Nfft + nu + tau + 1 - iRow)).^2 * ...
        abs(fft(h_tail(:), Nfft_psd, 1)).^2;
end
% Finally, scale:
S_isi = (Ex_bar/Nfft) * S_isi;

% Since the ICI has the same PSD, the total post-cursor ICPD becomes:
S_post = 2*S_isi;

%% Pre-cursor PSD

% Preallocate
S_pre_isi = zeros(Nfft_psd, 1);

% The noise can be interepreted as the sum of all CIR "heads", such as
% [h(1)], [h(1:2)], [h(1:3)] and so forth, up to [h(1:n0)];
for iRow = 1:n0
    h_head = h(1:iRow);
    S_pre_isi = S_pre_isi + ...
        abs(w(n0 + 1 - iRow)).^2 * ...
        abs(fft(h_head(:), Nfft_psd, 1)).^2;
end

% Finally, scale:
S_pre_isi = (Ex_bar/Nfft) * S_pre_isi;

% Since the ICI has the same PSD, the total pre-cursor ICPD becomes:
S_pre = 2*S_pre_isi;

%% Total ICPD PSD:

S_icpd = S_post + S_pre;
end

