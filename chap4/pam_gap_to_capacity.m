% A script to understand the uncoded PAM GAP to capacity
clearvars, clc

% First question is, what is the target Pe?
Pe = 1e-6;

N     = 1;           % dimensions
b     = [1 2 3 4 5]; % bits
b_bar = b / N;       % bits per dimension
M     = 2.^b         % modulation order
rho   = 2*b_bar;     % Spectral efficiency

% Now discover the SNR required to achieve the target Pe for each spectral
% effiency. Recall that, for PAM, the NNUB is exact and given by:
%
%   Pe_bar = 2*(1 - 1/M) qfunc( sqrt( 3 * SNRnormalized ) ),
%
% where SNRnormalized (the gap to capacity) is given by "SNR/(2^rho - 1)".
% Let us first discover the required argument for the Q function:
arg_qfunc = qfuncinv(Pe ./ ( 2*( 1 - (1./M) )) );
%
% Now, since the argument is equivalent to sqrt( 3 * SNRnormalized ), it is
% possible to obtain the gap for each case:
SNRnormalized = (arg_qfunc.^2) / 3;

% Finally, express it in dB:
fprintf('Gap to Capacity:\n');
10*log10(SNRnormalized)

% Note the gap is not constant due to the "2*(1 - 1/M)" term in Pe