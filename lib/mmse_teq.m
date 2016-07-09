function [w, b_opt] = mmse_teq(h, l, delta, Nf, Nb, Sx, sigma, debug)
% MMSE TEQ design based on Al-Dhahir's paper in [1]
%
%   Inputs:
% h     -> Channel impulse response
% l     -> Oversampling ratio
% delta -> Delay
% Nf    -> Number of output "symbols"
% Nb    -> Target impulse response (TIR) memory (length will be Nb + 1)
% Sx    -> Tx energy per dimension (equivalent to the PSD)
%
%   Optional Inputs
% Rxx -> autocorrelation of the transmit samples
%
%
% Notes:
% 1) Nf is stated as the number of "symbols", but in the case of OFDM, it
% can be interpreted as samples of the IFFT taken at every interval of l
% samples. In this case, it is less confusing to think only in terms of the
% number of samples l*Nf
%
% [1] Efficiently Computed Reduced-Parameter Input-Aided MMSE Equalizers
% for ML Detection: A Unified Approach

if (nargin > 7)
    if (debug)
        show_plots = 1;
    else
        show_plots = 0;
    end
else
    show_plots = 0;
end

%% Parameters
complexity       = 0; % 0 - Full ; 1 - Reduced
constraint       = 0; % 0 - UEC  ; 1 - UTC
noise_assumption = 0; % 0 - White; 1 - Colored, WSS;
input_assumption = 0; % 0 - White symbols; 1 - Colored, WSS;


%% Constants
nu       = length(h) - 1;
w_length = Nf * l; % Length of the equalizer

%% Design

% Channel matrix (possibly fractionally spaced):
H = toeplitz([h(1), zeros(1,w_length-1)]',[h,zeros(1,w_length-1)]);
% Dimensions: (Nf*l) x (Nf*l + nu - 1)


% Input autocorrelation
%
Rxx = Sx * eye( Nf*l + nu);

% Noise autocorrelation
switch (noise_assumption)
    case 0
        % White-noise auto-correlation matrix:
        Rnn = sigma * eye(w_length, w_length);
        % White noise implies the autocorrelation is non-zero only at lag
        % 0. Thus, the autocorrelation matrix is diagonal.
    case 1
        % If noise is a stationary process, then necessarily Rnn is a
        % Toeplitz square matrix, because each matrix diagonal corresponds
        % to a fixed interval between the instants where correlation is
        % measured. In this case, "sigma" must be a vector of length Nf*l
        % containing the autocorrelation values.
        Rnn = toeplitz(sigma);
end

% Input-output cross-correlation:
Rxy = Rxx * H';

% Output autocorrelation:
Ryy = H * Rxx * H' + Rnn;

% Matrix defined in Eq. (10):
if (l == 1)
    % For white input autocorrelation (only applicable for l = 1), i.e.,
    % Rxx = Sx * I, calculation of R_x_pipe_y can be simplified as shown in
    % Section III of [1]. Note it has a single inverse computation, rather
    % than three.
    R = (H * H') + (1/Sx)*Rnn;
    R_x_pipe_y = Sx*(eye(Nf*l + nu) - H' * inv(R) * H);
else
    R_x_pipe_y = inv(inv(Rxx) +  (H' * inv(Rnn) * H));
end

% Permutation/zero-padding matrix of Eq. (14) in [1]:
s = Nf*l + nu - delta - Nb - 1;
P = [zeros(Nb + 1, delta), eye(Nb + 1), zeros(Nb + 1, s)];
R_delta = P * R_x_pipe_y * P';

if (~constraint) % UEC
    % Solution is shown to be the eigenvector of R_delta that corresponds
    % to the minimum eigenvalue. Thus, find first eigenvectors and
    % eigenvalues:
    [B, Lambda] = eig(R_delta);
    % Find the minimum eigenvalue:
    [lambda_min, i_opt] = min(diag(Lambda));
    % Find the optimal target impulse response (TIR):
    b_opt = B(:,i_opt)';
    % Padded TIR:
    b_tilde_opt = [zeros(1, delta), b_opt, zeros(1, s)];
    % Then, obtain the equalizer from Eq. (21):
    w = b_tilde_opt * Rxy * inv(Ryy);
    % Find the MMSE:
    mmse = lambda_min;
end

%% Debug Plots

if (show_plots)
    sir = conv(h, w); % Shortened impulse response

    figure
    plot(sir)
    hold on
    plot(b_tilde_opt, '--r')
    xlabel('Index')
    legend('SIR', 'TIR')
    title('Shortened vs. target impulse response');
    ylabel('Amplitude')
end

end

