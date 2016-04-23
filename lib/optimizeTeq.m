function [ nTaps_o, delta_o ] = optimizeTeq(p, nu, sigma, Ex_bar, Lmax)
%optimizeTeq Searches the optimal TEQ length and delay
%   TEQ performance is very sensitive to the chosen delay and TEQ length.
%   An analytical formulation for these parameters is not available, so the
%   approach is to perform a search for the optimal parameters.
%
% Inputs:
%   p      =>    channel pulse response
%   nu     =>    the dispersion of the desirable pulse response (b)
%   sigma  =>    noise power in linear scale
%   Ex_bar =>    Tx Energy per dimension
%   Lmax   =>    Maximum number of taps

% Initialize
L_vec     = 1:Lmax; % Vector of TEQ lengths
delta_vec = 1:Lmax; % Vector of delays

% Constants
filtertype    = 1; % FIR Filter Type
nTap_length   = length(L_vec);
nDelta_length = length(delta_vec);

% Preallocate
SNRteq = NaN * ones(nTap_length, nDelta_length);

%% Search

for i = 1:nTap_length
    % Number of taps
    L = L_vec(i);
    % The delay has to be less than or equal to the number of taps
    iDeltaLength = sum(delta_vec <= L);

    % Search delay
    for k = 1:iDeltaLength
        % Delay
        delta = delta_vec(k);

        % Design the TEQ and evaluate the resulting SNRmfb:
        [~, ~, SNRteq(i,k), ~] = ...
                teq(p, L, nu, delta, sigma, Ex_bar, filtertype);

        % Prevent negative linear SNR
        if (SNRteq(i,k) < 0)
            SNRteq(i,k) = NaN;
        end
    end
end

%% Results

% 3-D plot of the SNRs
figure
surf(delta_vec, L_vec, 10*log10(SNRteq))
xlabel('Delta')
ylabel('L')
zlabel('SNR_{teq} (dB)')
title('TEQ Delay and Length Search')

% Find the optimal combination of number of taps and delay:
maxSnr = max(SNRteq(~isnan(SNRteq)));
[i, k] = find(SNRteq == maxSnr);
nTaps_o = L_vec(i);
delta_o = delta_vec(k);

end

