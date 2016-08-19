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

    % Phase shift (in freq. domain) due to circular shift in time domain:
    phaseShift = sparse(diag(exp(1j*2*pi*(tau/N)*(0:N-1))));
    % Assume that the time-domain DMT vectors multiplying the ISI matrix
    % will be circularly shifted by -tau in the time-domain.

    % Ideal Freq. Domain Channel Matrix (CIR FFT with a phase shift):
    H = Q * Hcirc * Q'; % Diagonal Matrix with the Channel Gains
    % Feed-forward Precoder

    % Precoder
    W_common = (Hcirc+Hici)\Hcirc;
    W = Q * W_common * Q';

    % Feedback equalizer
    B = - H\(Q*Hisi*Q'*phaseShift*Q);

    clear Hcirc Hici Hisi Ht Phase_shift

    Precoder.W = W;
    Precoder.B = B;
    % Row norms:
    Precoder.wk = sum(abs(W).^2,2).';

end

