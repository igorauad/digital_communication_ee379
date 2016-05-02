function [ window ] = designDmtWindow( N, nu, tau )
% Window for the cyclic extended DMT symbols
% --------------------------------------------
% [ window ] = designDmtWindow.m( N, nu, tau )
%
%   Designs a window with size N + nu + tau (cyclic extended DMT symbol
% length). The window portion that multiplies the N information samples is
% flat and equal to 1, while, the tapered part of the window, which
% multiplies the Cyclic Prefix and the Cyclic Suffix, uses a raised-cosine
% function.
%
%   Input Parameters
% N                     FFT size (non-extended DMT symbol length)
% nu                    Cyclic prefix length
% tau                   Cyclic suffix length
%
%   Output:
% window                Window sequence
%

windowIndexes = (-tau:1:(tau-1)).';
raisedCosineWindow = (1/2)*(1 + cos(pi*(windowIndexes/tau)));
lowerWindow = raisedCosineWindow(1:tau); % Length of lowerWindow is Lcs
upperWindow = raisedCosineWindow(tau+1:end);

% Put all parts together:
window = [lowerWindow; ones(N + nu - tau,1); upperWindow];

end

