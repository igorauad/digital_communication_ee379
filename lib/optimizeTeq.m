function [ nTaps_o, delta_o ] = optimizeTeq(criterion, p, nu, sigma, Ex_bar, gap, N, Lmax, delta_min, delta_max)
%optimizeTeq Searches the optimal TEQ length and delay
%   TEQ performance is very sensitive to the chosen delay and TEQ length.
%   An analytical formulation for these parameters is not available, so the
%   approach is to perform a search for the optimal parameters. According
%   to the dissertation in [1], there are three main search criteria: the
%   MSE, the Shortening SNR (SSNR) and the Geometric SNR. The three modes
%   of operation are configurable in this script
%
% Inputs:
%   criterion   =>    Optimization criterion
%           0) Optimize using MSE;
%           1) Optimize using SSNR;
%           2) Optimize using Geo-SNR;
%   p           =>    channel pulse response
%   nu          =>    the dispersion of the desirable pulse response (b)
%   sigma       =>    noise power in linear scale
%   Ex_bar      =>    Tx Energy per dimension
%   Lmax        =>    Maximum number of taps
%
%
%   References:
% [1] http://users.ece.utexas.edu/~bevans/projects/adsl/index.html
% [2] N. Al-Dhahir and J. Cioffi, "Optimum finite-length equalization for
% multicarrier transceivers," Global Telecommunications Conference, 1994.
% GLOBECOM '94. Communications: The Global Bridge., IEEE, San Francisco,
% CA, 1994, pp. 1884-1888 vol.3.
%
if (delta_min > delta_max)
    error('Search range is null');
end
% Initialize
L_vec     = 1:Lmax; % Vector of TEQ lengths
delta_vec = delta_min:delta_max; % Vector of delays

% Constants
filtertype    = 1; % FIR Filter Type
nTap_length   = length(L_vec);
nDelta_length = length(delta_vec);

% Preallocate
objective = NaN * ones(nTap_length, nDelta_length);

gap_db = 10*log10(gap);

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
        [w, b, SNRteq, ~] = ...
            teq(p, L, nu, delta, sigma, Ex_bar, filtertype);

        switch (criterion)
            case 0 % MSE
                criterion_label = 'MSE';
                % Maximizes the TEQ MFB that considers the minimzed MSE in
                % the denominator:
                objective(i,k) = SNRteq;
            case 1 % SSNR
                criterion_label = 'SSNR';
                % Maximizes the Shortening SNR
                objective(i,k) = teqSSNR( w, p, delta, nu );
            case 2 % Geo SNR
                % Maximizes the Geometric (multich-channel) SNR
                criterion_label = 'GeoSNR';
                % New effective channel (see Eq. 7 in [2]):
                B_i = fft(b, N);
                W_i = fft(w, N);
                % New Unitary-energy SNR:
                gn_teq = (abs(B_i).^2)./(sigma*(abs(W_i).^2));

                % Cost function
                L_b = sum(log(gn_teq))/N;

                objective(i,k) = L_b;
            otherwise
                error('Non-available TEQ optimization criterion');
        end
    end
end

%% Results

% 3-D plot of the SNRs
figure
surf(delta_vec, L_vec, objective)
xlabel('Delay (Delta)')
ylabel('TEQ Length (L)')
zlabel('Objective')
title('TEQ Delay and Length Search')
title(sprintf('Optimized in terms of %s', criterion_label));

% Find the optimal combination of number of taps and delay:
maxSnr = max(objective(~isnan(objective)));
[i, k] = find(objective == maxSnr);
nTaps_o = min(L_vec(i));
delta_o = min(delta_vec(k));

end

