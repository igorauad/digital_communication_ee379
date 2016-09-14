function [w, SSNR] = ssnr_teq(h, l, delta, Nf, nu, debug)
% SSNR TEQ design based on Melsa's paper in [1]
%   Designs an optimal equalizer such that the energy of the shortened
%   impulse response (SIR) outside the cyclic prefix (window) is minimum
%   with respect to the energy of the SIR within the prefix. Noise is
%   neglected in the optimization problem.
%
%   Inputs:
% h     -> Channel impulse response
% l     -> Oversampling ratio
% delta -> Delay
% Nf    -> Number of output "symbols"
% nu    -> Cyclic prefix length, equivalent to the TIR memory
% debug -> Enable debug plots
%
%   Output:
% w     -> Equalizer
% SSNR  -> Shortening SNR
%
% [1] P. J. W. Melsa, R. C. Younce and C. E. Rohrs, "Impulse response
%     shortening for discrete multitone transceivers," in IEEE Transactions
%     on Communications, vol. 44, no. 12, pp. 1662-1672, Dec 1996.

if (nargin > 5)
    if (debug)
        show_plots = 1;
    else
        show_plots = 0;
    end
else
    show_plots = 0;
end

%% Constants

t = Nf * l; % Length of the equalizer (in Melsa's notation)
M = length(h); % Impulse response length (in Melsa's notation)

%% Error Checking

if (t > nu)
    warning('Matrix B will not be full rank');
    warning('Cholesky decomposition will fail');
end

%% Design

H = convmtx(h.', t);
% Dimensions (M + t - 1) x t

% Matrix that contains the taps of the SIR value within the window of
% interest (containing nu + 1) consecutive samples:
win_indexes = delta + 1 : delta + nu + 1;
H_win = H(win_indexes, :);
% Has dimensions (nu + 1) x t
% Assuming t < nu, H_win is a skinny matrix

H_wall = [ H(1:delta, :); H(delta + nu + 2:end, :)];
% Has dimensions (M + t -nu -2) x t
% Assuming M > nu + 2, H_wall is also a skinny matrix.

% Symmetric and positive semidefinite matrices of (5) and (6) in [1]:
A = H_wall' * H_wall;
B = H_win' * H_win;
% Both are t x t matrices.
% Matrix B must be invertible, but this is not necessarily true.
%
% Note:
% x' * A' * A * x = (A * x)' * (A * x)
% Since A*x is a (t x 1) column vector, the above is the norm of the
% vector, which is necessarily greater than or equal to 0.

% Cholesky decomposition:
[Q, Lambda] = eig(B);
sqrt_B = Q*sqrt(Lambda);
% Note 1: B must be full rank and positive semidefinite matrix for the
% Cholesky decomposition to exit.
%
% Note 2: The sqrt(B) could also be computed by the function
% chol(B,'lower'), however this function returns an error when B is not
% full rank, so the aproach used above was preferred.

C = inv(sqrt_B) * A * inv(sqrt_B');

[Eigvec_c, Eigval_c] = eig(C);
% Find the minimum eigenvalue:
[~, i_opt] = min(diag(Eigval_c));
% Find the eigenvector corresponding to the minimum eigen value.
l_min = Eigvec_c(:,i_opt); % Note this is b_opt' (transpose)

% Optimum equalizer:
w = inv(sqrt_B') * l_min;
w = w.'; % Output as row vector

% Shortening SNR:
SSNR = ssnr( w, h, delta, nu );

%% Debug Plots

if (show_plots)
    sir = conv(h, w); % Shortened impulse response

    figure
    plot(sir)
    xlabel('Index')
    title('Shortened impulse response');
    ylabel('Amplitude')
    % Set ticks marking the beginning and end of the TIR "window"
    set(gca, 'XTick', sort([delta + 1, delta + nu + 1]));
    grid on
end

end

