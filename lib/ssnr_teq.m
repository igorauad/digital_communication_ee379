function [w, SNR, bias] = ssnr_teq(h, l, delta, Nf, nu, debug)
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

H_wall = [ H(1:delta, :); H(delta + nu + 2:end, :)];
% Has dimensions (M + t -nu -2) x t

% Symmetric and positive semidefinite matrices of (5) and (6) in [1]:
A = H_wall' * H_wall;
B = H_win' * H_win;

% [Q, Lambda] = eig(B);
% sqrt_B = Q*sqrt(Lambda);

% Cholesky decomposition:
sqrt_B = chol(B,'lower'); % sqrt_B * sqrt_B' = B
% Note: B must be full rank and positive semidefinite matrix for the
% Cholesky decomposition to exit.

C = inv(sqrt_B) * A * inv(sqrt_B');

[Eigvec_c, Eigval_c] = eig(C);
% Find the minimum eigenvalue:
[~, i_opt] = min(diag(Eigval_c));
% Find the eigenvector corresponding to the minimum eigen value.
l_min = Eigvec_c(:,i_opt); % Note this is b_opt' (transpose)

% Optimum equalizer:
w = inv(sqrt_B') * l_min;
w = w.'; % Output as row vector

%% Debug Plots

if (show_plots)
    sir = conv(h, w); % Shortened impulse response

    figure
    plot(sir)
    xlabel('Index')
    title('Shortened impulse response');
    ylabel('Amplitude')
end

end

