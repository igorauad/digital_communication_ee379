function [ x ] = precodeTimeDomain( x, Precoder )
%   Precode time domain symbols
% --------------------------------------------------
% [ x ] = precodeTimeDomain( x_nce, params )
%       Applies time-domain precoding on DMT symbols in order to provide
%   ICI mitigation.
%
% Input:
%   x         ->    Non-cyclic-extended DMT symbols
%   Precoder  ->    Struct with the time-domain precoding matrix
%
% Output:
%   x         ->    Precoded time-domain DMT symbols

% Precoding matrix for the current line:
Wt    = Precoder.ici.Wt;
iRows = Precoder.ici.significantRows;
iCols = Precoder.ici.significantCols;

% Precode the symbols:
x(iRows, :) = x(iRows, :) + Wt*x(iCols, :);

end

