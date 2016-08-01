function [ x_out ] = windowAndOverlap( x_ce, dmtWindow, N, nu, tau )
% Spectrum shaping through windowing and overlap between CP and CS
% -------------------------------------------------------------------------
%   Window DMT symbols with the given window and apply overlap between
%   cyclic suffix samples and cyclic prefix samples. The result is the
%   matrix containing the cyclic prefixed DMT symbols whose cyclic prefix
%   contain the overlap from the previous cyclic suffix. In another words,
%   x_out has dimensios (N + nu) x nSymbols, whereas the input x_ce has
%   dimensions (N + nu + tau) x nSymbols.
%
%       Input:
%   x_ce       --> Cyclic extended DMT symbols (N + nu + tau samples each)
%   dmtWindow  --> Window to be applied on each DMT symbol
%   N          --> FFT size
%   nu         --> Cyclic Prefix Length
%   tau        --> Cyclic Suffix Length
%
%       Output:
%   x_out      --> Windowed and overlapped time-domain DMT symbols

% Apply window:
x_windowed = bsxfun(@times, x_ce, dmtWindow);

% Overlap
x_windowed(1:tau,:) = x_windowed(1:tau,:) + ...
    [zeros(tau,1), x_windowed((end - tau + 1):end, 1:(end-1))];

% Cut cyclic suffix copy (the other copy was already
% summed in the "Overlap" step):
x_out = x_windowed(1:(N+nu),:);


end

