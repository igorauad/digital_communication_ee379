function [ psd ] = checkInterferencePsd(Y, H, X, N, fs, PLOT_PSD)
%checkInterferencePsd Estimates the interference PSD
%   This script estimates the interference (ISI/ICI) power spectral density
%   based on the information of the received symbols and the originally
%   transmitted symbols. Given the channel response, it is possible to
%   calculate the ideally received symbols based on the originally transmit
%   symbols. The difference between the actual received symbols and the
%   ideally received symbols is used to estimate the PSD.
%
%   Note 1: The PSD is estimated using Welch`s periodogram method.
%
%   Note 2: If the received symbols Y is disturbed by background noise, the
% estimated PSD will not be solely the interference PSD (ISI+ICI), but the
% combination of interference and noise.
%
%   Input parameters:
% Y                         Received frequency domain symbols
% H                         Channel frequency response
% X                         Transmit symbols
% fs                        Sampling frequency
% PLOT_PSD                  Flag to activate plotting
%
%   Output parameters
% psd                       Estimated PSD

    % Ideally received symbols:
    Y_ideal = bsxfun(@times,H,X);

    % Compute interference (or the combination of interference + background
    % noise):
    Y_interference = Y - Y_ideal;

    % Time domain interference:
    y_interference = ifft(Y_interference,N);

    % Estimate PSD using Welch's periodogram method:
    [psd, freq] = pwelch(y_interference, [], [], N/2, fs, 'twosided');

    % Plot
    if (PLOT_PSD)
        figure
        hold on

        plot(freq/1e6, 10*log10(psd))
        xlabel('Frequency (MHz)')
        f_max = max(freq);
        set(gca,'XTick',[0:(f_max/8):f_max]*1e-6)
        set(gca,'Xlim',[0 ((f_max/2)*1e-6)])
        ylabel('PSD (dB/Hz)')
        grid on

        h_fig_ck = gcf;
    end

end

