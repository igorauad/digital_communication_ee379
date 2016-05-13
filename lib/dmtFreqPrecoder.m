function [Precoder] = dmtFreqPrecoder(h, N, Lcp, Lcs, n0, WINDOWING)
%assembleCheongPrecoder Assembles the Cheong [1] Precoding Matrices for
%Insufficient Cyclic Prefix Distortion mitigation
%
%   Input Parameters
% h                 Channel Impulse Response (CIR)
% N                 FFT size
% Lcp               Cyclic Prefix Length
% Lcs               Cyclic Suffix length
% n0                CIR delay (index of its peak)
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

    Lh = length(h);
    h = h(:);

    fprintf('\n\nCheong Precoder\n')

    Q = (1/sqrt(N))*fft(eye(N));

    [ Hisi, Hici, Hcirc ] = dmtIsiIciMatrices(h,...
            n0, Lcp, Lcs, N, 1, WINDOWING);

    % Phase shift (in freq. domain) due to circular shift in time domain:
    phaseShift = sparse(diag(exp(1j*2*pi*(n0/N)*(0:N-1))));

    % Ideal Freq. Domain Channel Matrix (CIR FFT with a phase shift):
    H = Q*Hcirc*Q'*phaseShift; % Diagonal Matrix with the Channel Gains
    Hdiag = diag(H);
    % Feed-forward Precoder

    % Precoder
    W_common = (Hcirc+Hici)\Hcirc;
    W_common = circshift(W_common,[n0 n0]);
    W = Q * W_common * Q';

    % Feedback equalizer
    B = - H\(Q*Hisi*Q'*phaseShift*Q);
    
    clear Hcirc Hici Hisi Ht Phase_shift

    Precoder.W = W;
    Precoder.B = B;
    % Row norms:
    Precoder.wk = sum(abs(W).^2,2).';

end

