function [] = viewConstellation(Z, Xconst, k)
% Shows Rx symbols in the original constellation
%
% Inputs
% Z -> Received Freq. Domain Symbols
% X -> Original Freq. Domain Symbols
% k -> Tone to be inspected

N = size(Z, 1);

    figure
    if (k == 1 || k == N/2)
    	% One-dimensional subchannel
    	plot(Z(k, :) + j*eps, 'o')
        hold on
        plot(Xconst + j*eps, 'ro', 'MarkerSize', 8, 'linewidth', 2)
    else
    	% Two-dimensional subchannel
        plot(Z(k, :), 'o')
        hold on
        plot(Xconst, 'ro', 'MarkerSize', 8, 'linewidth', 2)
    end
    legend('Rx', 'Tx')
    title(sprintf('Tone: %d', k));

end