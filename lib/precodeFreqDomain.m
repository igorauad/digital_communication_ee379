function [ X ] = precodeFreqDomain( X, Precoder, M, D )
%   Precode frequency domain symbols
% --------------------------------------------------
% [ x ] = precodeTimeDomain( x_nce, params )
%       Applies frequency-domain precoding on DMT symbols in order to
%       provide ICI mitigation.
%
% Input:
%   X         ->    Frequency-domain DMT symbols
%   Precoder  ->    Struct with the time-domain precoding matrix
%   M         ->    Vector of constellation orders for each subchannel
%   D         ->    Vector of minimum distances for each subchannel
%
% Output:
%   x         ->    Precoded time-domain DMT symbols

% Infer number of symbols and FFT size
nSymbols   = size(X, 2);
N          = size(X, 1);

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

