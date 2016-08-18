function [ X ] = precodeFreqDomain( X, Precoder, modOrder, dmin, usedTones )
%   Precode frequency domain symbols
% --------------------------------------------------
% [ x ] = precodeTimeDomain( x_nce, params )
%       Applies frequency-domain precoding on DMT symbols in order to
%       provide ICI mitigation.
%
% Input:
%   X         ->    Frequency-domain DMT symbols
%   Precoder  ->    Struct with the time-domain precoding matrix
%   modOrder  ->    Vector of modulation orders for each subchannel
%   dmin      ->    Vector of minimum distances for each subchannel
%   usedTones ->    Vector indicating the DFT tones that are used
%
% Output:
%   x         ->    Precoded time-domain DMT symbols

%% Initialization

% Infer number of symbols and FFT size
nSymbols   = size(X, 2);
N          = size(X, 1);

% Initialize full vector of modulation orders:
M            = ones(N/2 +1, 1); % Note an unitary order corresponds to b=0
% Set the actual order corresponding to the used tones:
M(usedTones) = modOrder;

% Initialize full vector of minimum distances:
D            = zeros(N/2 +1, 1);
D(usedTones) = dmin;

%% Precode

% First symbol (no feedback term)
X(:,1) = Precoder.W * moduloOperator(X(1:N/2+1,1), M, ...
    D, 'Hermitian');

% Remaining symbols
for sym_index = 2:nSymbols
    % Current symbol + feedback term due to previous precoded symbol
    X_prev_prec = X(:,sym_index) + ...
        Precoder.B * sqrt(N) * ifft(X(:,sym_index-1),N);
    % Precode:
    X(:,sym_index) = Precoder.W * moduloOperator(X_prev_prec(1:N/2+1), ...
        M, D, 'Hermitian');
end

end

