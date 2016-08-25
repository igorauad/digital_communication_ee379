function [Precoder] = dmtFreqPrecoder(h, N, nu, tau, n0, windowing)
%assembleCheongPrecoder Assembles the Cheong [1] Precoding Matrices for
%Insufficient Cyclic Prefix Distortion mitigation
%
%   Input Parameters
% h                 Channel Impulse Response (CIR)
% N                 FFT size
% nu                Cyclic Prefix Length
% tau               Cyclic Suffix length
% n0                CIR delay (index of its peak)
% windowing         Indicates whether windowing+overlap is used
%
%   Output parameters
%
% W                 Preocoding matrix on direct path (after modulo
%                   operation), which is used for ICI mitigation
% B                 Precoding matrix on feedback path, which is used for
%                   ISI mitigation
%
% [1] Kok-Wui Cheong and J. M. Cioffi, "Precoder for DMT with insufficient
% [cyclic prefix," Communications, 1998. ICC 98. Conference Record. 1998
% IEEE [International Conference on, Atlanta, GA, 1998, pp. 339-343 vol.1.

    % Ensure the CIR is a column vector
    h = h(:);

    fprintf('\n\nCheong Precoder\n')

    Q = (1/sqrt(N))*fft(eye(N));

    [ Hisi, Hici, Hcirc ] = dmtIsiIciMatrices(h,...
            n0, nu, tau, N, windowing);

    % Assume that the time-domain DMT vectors multiplying the ISI matrix
    % will be circularly shifted by -tau in the time-domain.

    % Compute Ideal Freq. Domain Channel Matrix (CIR FFT):
    %
    % The diagonal Matrix could be computed by (Q * Hcirc * Q'), but
    % instead it is optimized to:
    tmp  = fft(Hcirc, N);
    H = (1/N) * fft(conj(tmp), N, 2);
    % Note: CIR is assumed to already contain a phase shift

    % Feed-forward Precoder (Q * W_common * Q')
    W_common = (Hcirc+Hici)\Hcirc;

    tmp  = fft(W_common, N);
    W = (1/N) * fft(tmp', N)';

    % Feedback equalizer

    % Explicit equation:
    % Q = (1/sqrt(N)) * fft(eye(N));
    % B = - (Hcirc * Q')\circshift(Hisi, [0 tau]);
    % which differs from Cheong's formulation in terms of the circular
    % shift due to tau

    % Equivalent implementation using the FFT algorithm:
    tmp  = (1/sqrt(N)) * fft(Hcirc', N)';   % (Hcirc * Q')
    B    = - tmp \ circshift(Hisi, [0 tau]);

    % Save for external usage
    Precoder.W = W;
    Precoder.B = B;
end

