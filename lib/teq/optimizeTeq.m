function [ nTaps_o, delta_o ] = optimizeTeq(type, p, nu, l, sigma, Ex_bar, N, debug)
%optimizeTeq Searches TEQ length and delay parameters
%   TEQ performance is very sensitive to the chosen delay and TEQ length.
%   An analytical formulation for these parameters is not available, so the
%   approach is to perform a search for the optimal parameters. As
%   summarized in the dissertation [1], there are three main TEQ designs:
%   based on MMSE, the Shortening SNR (SSNR) and the Geometric SNR. The
%   three modes of operation are configurable in this script.
%
%   General (but loose) guidelines for the choices of delay and length are
%   presented in [3]. Nonetheless, specific conditions are applicable to
%   each type of TEQ design.
%
% Inputs:
%   type    =>  TEQ design criterion
%                   0) Minimization of MSE;
%                   1) Maximimization of Shortening SNR (SSNR);
%                   2) Maximimization of Geometric SNR (GeoSNR);
%   p       =>  Channel pulse response
%   l       ->  Oversampling factor
%   nu      =>  Memory of the target impulse response (b)
%   sigma   =>  Noise energy per dimension
%   Ex_bar  =>  Signal energy per dimension
%
%
%   References:
% [1] http://users.ece.utexas.edu/~bevans/projects/adsl/index.html
% [2] N. Al-Dhahir and J. Cioffi, "Optimum finite-length equalization for
%     multicarrier transceivers," Global Telecommunications Conference,
%     1994. GLOBECOM '94. Communications: The Global Bridge., IEEE,
%     San Francisco, CA, 1994, pp. 1884-1888 vol.3.
% [3] R. K. Martin and C. R. Johnson, "Adaptive equalization: transitioning
%     from single-carrier to multicarrier systems," in IEEE Signal
%     Processing Magazine, vol. 22, no. 6, pp. 108-122, Nov. 2005.
% [4] P. J. W. Melsa, R. C. Younce and C. E. Rohrs, "Impulse response
%     shortening for discrete multitone transceivers," in IEEE Transactions
%     on Communications, vol. 44, no. 12, pp. 1662-1672, Dec 1996.
% [5] Efficiently Computed Reduced-Parameter Input-Aided MMSE Equalizers
%     for ML Detection: A Unified Approach

if (nargin > 7)
    if (debug)
        show_plots = 1;
    else
        show_plots = 0;
    end
else
    show_plots = 0;
end

%% Channel Length (considering energy)

chLenThresh = 5e-4;
chLen = find(abs(p) > chLenThresh, 1, 'last') - ...
    find(abs(p) > chLenThresh, 1, 'first');

%% Equalizer Length Range
% Equalizer length is advised in [3] to be 3-5 times the difference of
% channel and prefix lengths. In contrast, [4] states that typical DMT
% transceivers use an equalizer (there called SIRF) length which is shorter
% than the length of the CP. For the MMSE-TEQ, we opt to adapt the
% guidelines of [3] to only 1 to 3 times the difference of channel and
% prefix lengths. The reason is that 3-5 times the difference yields
% lengeths that are too long, so that the simulation becomes very slow. In
% the specific case of MaxSSNR TEQ, the equalizer length must strictly be
% less than or equal to the cyclic prefix length in order for one matrix of
% the derivation to be positive semidefinite.

switch (type)
    case {0,2} % MSE and GeoSNR
        maxTaps  = 4*(chLen - nu);
        minTaps  = 1*(chLen - nu);
        % Force at least a range of nu
        if (maxTaps < minTaps + nu)
            maxTaps = minTaps + nu;
        end
    case 1 % SSNR
        maxTaps  = nu;
        minTaps  = l;
    otherwise
        error('Non-available TEQ type');
end

% Avoid null number of symbol-spaced taps (floor of a div by l is used):
if (minTaps < l)
    minTaps = l;
end

% Vector of TEQ lengths:
nTap_vec = minTaps:maxTaps;

% Downsample to avoid a very long search
searchLength = 20;
step_size = round(length(nTap_vec)/searchLength);
if (step_size > 1)
    nTap_vec = downsample(nTap_vec, step_size);
else
    % Step size can be rounded to 0. In this case, it is forced to 1.
    step_size = 1;
end

fprintf('Optimizing TEQ length between %d and %d\n', minTaps, maxTaps);
fprintf('TEQ Length search step:\t %d\n', step_size);

%% Delay Range
% Desired delay is advised in [3] to be the location of nu-length window of
% channel with largest energy + half the length of the shortener.

win_energy = zeros(length(p) - nu, 1);
for i = 1:length(p) - nu
    win_energy(i) = norm(p(1+i : i+nu))^2;
end
% Maximum windowed energy
[~, iMax]  = max(win_energy);
% Location of the window center
iWinCenter = iMax + nu/2;

% Minimum delay in the range
delta_min  = iWinCenter - ceil(maxTaps/2);
if (delta_min < 0)
    delta_min = 0;
end

% Maximum delay in the range
delta_max  = iWinCenter + ceil(maxTaps/2);

% Vector of delays:
delta_vec  = delta_min:delta_max;

% Downsample to avoid a very long search
step_size = round(length(delta_vec)/searchLength);
if (step_size > 1)
    delta_vec = downsample(delta_vec, step_size);
else
    % Step size can be rounded to 0. In this case, it is forced to 1.
    step_size = 1;
end

fprintf('Optimizing delay between %d and %d\n', delta_min, delta_max);
fprintf('Delay search step:\t %d\n', step_size);

%% Initialize

% Constants
nTap_length   = length(nTap_vec);
nDelta_length = length(delta_vec);

% Preallocate
objective = NaN * ones(nTap_length, nDelta_length);

%% Search
for i = 1:nTap_length
    % Number of taps
    nTaps = nTap_vec(i);
    % Symbol-spaced length of the equalizer (actual length is Nf*l)
    Nf = floor(nTaps/l);

    % Search delay
    for k = 1:length(delta_vec)
        % Delay
        delta = delta_vec(k);

        switch (type)
            case 0 % MSE
                type_label = 'MSE';
                objective_label = '10\\log(SNRmfb)';
                % Consider the contraint presented in [5] for the delay
                % delta with respect to the equalizer length.
                if (delta <= Nf + (length(p) - 1) - nu - 1)
                    % Design the TEQ and evaluate the resulting SNRmfb:
                    [w, SNRteq] = ...
                        mmse_teq(p, l, delta, Nf, nu, Ex_bar, ...
                        sigma);
                    % Maximizes the TEQ MFB that considers the minimzed MSE
                    % in the denominator:
                    objective(i,k) = SNRteq;
                end
            case 1 % SSNR
                type_label = 'SSNR';
                objective_label = '10\\log(SSNR)';
                % Design the TEQ and evaluate the resulting SSNR:
                [w, SSNR] = ssnr_teq(p, l, delta, Nf, nu);
                % Maximizes the Shortening SNR
                objective(i,k) = SSNR;
            case 2 % Geo SNR
                % Maximizes the Geometric (multich-channel) SNR
                type_label = 'GeoSNR';
                % New effective channel (see Eq. 7 in [2]):
                B_i = fft(b, N);
                W_i = fft(w, N);
                % New Unitary-energy SNR:
                gn_teq = (abs(B_i).^2)./(sigma*(abs(W_i).^2));

                % Cost function
                L_b = sum(log(gn_teq))/N;

                objective(i,k) = L_b;
            otherwise
                error('Non-available TEQ type');
        end
    end
end

%% Delay and equalizer length choices

% Find the optimal combination of number of taps and delay:
maximum = max(objective(~isnan(objective)));

% Choose the equalizer that yields an objective value within 5% of the
% maximum, but within the minimum length and delay.
[i, k] = find(objective > 0.95*maximum);
nTaps_o = min(nTap_vec(i));
delta_o = min(delta_vec(k));

% TODO: consider the optimal tradeoff curve of the multiobjective
% optimization and avoid obtaining a significantly longer length at the
% expense of minor gain in the objective.

%% Results

if (show_plots)
    % 3-D plot of the cost function
    figure
    surf(delta_vec, nTap_vec, 10*log10(objective))
    xlabel('Delay (Delta)', 'Interpreter', 'latex')
    ylabel('TEQ Length (L)', 'Interpreter', 'latex')
    zlabel(objective_label, 'Interpreter', 'latex')
    title('TEQ Delay and Length Search', 'Interpreter', 'latex')
    title(sprintf('TEQ optimized in terms of %s', type_label));
end

end

