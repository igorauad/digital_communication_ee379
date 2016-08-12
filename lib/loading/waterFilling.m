function [ bn_bar, En_bar, iUsed ] = waterFilling(gn, Ex, N, gap)
%WATERFILLING Water-filling bit-loading solution
%   Finds the rational bit-load that satisfies the water-filling solution
%   subject to the constraint that Ex = Ex_bar*N, namely the transmit
%   energy is equal to the Energy per dimension * number of dimensions.
%
%   Inputs
% gn      -> SNRn when the transmitter applies unit energy per dimension
%            of the subchannel.
% Ex      -> Total energy budget
% N       -> Number of dimensions
% gap     -> Linear SNR gap to capacity
%
%   Outputs
% bn_bar  -> Optimal bit loading vector per dimension
% En_bar  -> Energy per dimension
% iUsed   -> Used dimensions
%
%
% Note: two dimensions can pertain to the same subchannel. For example, for
% DMT, the vector always has 2 real subchannels at DC and Nyquist and N/2 -
% 1 complex subchannels with two real dimensions.

%% Initialization

% Preallocate output vectors:
En_bar  = zeros(1, length(gn));
bn_bar  = zeros(1, length(gn));

% Initialize the vector of indexes corresponding to the used dimensions:
iUsed = find(gn ~= 0);

%% Run search for the water-fill level:

i = 0;
positiveEnergies = 0;

while (~positiveEnergies)
    gn_vector = gn(iUsed);

    N_i = N - i;

    % Use (4.44) to compute the constant that satisfy the water-filling
    % equation (4.24):
    K = (1 / N_i) * (Ex +  gap * sum( 1 ./ gn_vector ) );

    % Then, the corresponding energies for each dimension (Eq. 4.45):
    En_bar_i = K - gap*(1./gn_vector);

    % Find the one with the minimum energy
    [minEn, iMinEn] = min(En_bar_i);

    % Is the minimum energy negative?
    if (minEn < 0)
        % Remove the dimension if minimum energy from the vector of usable
        % dimensions:
        iUsed(iMinEn) = [];
        % Attempt one more time
        i = i + 1;
    else
        % All energies are positive -> constant K was found!
        positiveEnergies = 1;

        fprintf('Water-fill level:           \t %g \n', K);
    end

end

%% Optimal loading:

% Energies per per dimension (in this vector, a 0 means the dimension
% should not be used):
En_bar(iUsed) = En_bar_i;

% Aslanis formula (4.46) at the used dimensions:
b_bar_used = .5 * log2(K * gn(iUsed) / gap);

% Fill the complete output vector (where 0 means the dimension is not
% used):
bn_bar(iUsed) = b_bar_used;


end

