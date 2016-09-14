function [En,bn] = DMTLCra(gn, Ex, N, gap_db, max_load, n_dim)
%
% EE379C 2008 Spring
%
% Levin Campello?s Method with DMT
%
% Inputs
% gn            Unitary-energy SNRs (gain-to-noise ratios)
% Ex            Energy budget
% N             Number of used real dimensions
% gap_db is     Gap to capacity in dB
% max_load      Maximum allowed bit load for a given subchannel (optional)
% n_dim         Number of real dimensions for each subchannel (optional)
%
% Outputs
% En is the vector energy distribution (PAM or QAM) per subchannel
% bn is the vector bit distribution (PAM or QAM) per subchannel
%
% If vector "n_dim" is not passed as argument, then gn is inferred to come
% from the positive half of an Hermitian symmetric vector. In this case,
% the first and last bins are PAM (one-dimensional), while the rest are
% QAM.

%% Parameter Definitions

% Define the maximum  bit load
if (nargin < 5)
    max_load = inf;
end

% Define the number of real dimensions corresponding to each index of gn
if (nargin < 6)
    % If the number of real dimensions used in each index is not given,
    % define it according to an Hermitian symmetric vector:
    n_dim = [1; 2*ones(N/2 -1, 1); 1];
else
    if (length(n_dim) ~= length(gn))
        error('Vector with number of dimensions has innapropriate length');
    end
end

% Define sets containing the indexes of unused, one-dimensional and
% two-dimensional subchannels
unused_indexes  = find(n_dim == 0);
one_dim_indexes = find(n_dim == 1);
two_dim_indexes = find(n_dim == 2);

% Sanity check
if (sum(n_dim) > N)
    error('Number of used real dimensions exceeds the limit');
end

% Gap to capacity in linear scale
gap = 10^(gap_db/10);

%% Levin Campello Initialization

% Preallocate/initialize
En             = zeros(1, length(gn)); % Energy per subchannel
bn             = zeros(1, length(gn)); % Bits per subchannel
decision_table = zeros(1, length(gn)); % Decision table
E_so_far       = 0;                    % Energy used so far

% Decision Table Initialization

% Two-dimensional subchannels
decision_table(two_dim_indexes) = 2*gap./gn(two_dim_indexes);

% One-dimensional subchannels
decision_table(one_dim_indexes) = 3*gap./gn(one_dim_indexes);

% Unused subchannels
decision_table(unused_indexes)  = inf;

%% Levin Campello Loading
while(1)
    [y,index] = min(decision_table);
    if (isinf(y))
        warning('Consider pouring the remaining energy over subchannels');
        error('All suchannels are maximally loaded');
    end
    E_so_far = E_so_far+y;
    if E_so_far > Ex
        break;
    else
        En(index) = En(index) + y;
        bn(index) = bn(index) + 1;
        % Prevent constellations above a certain maximum load
        if( (bn(index) == max_load) )
            decision_table(index) = inf;
        end
        % Update the decision table, considering whether the index
        % corresponds to a 1-dimensional or 2-dimensional subchannel
        if (n_dim(index) == 1)
            decision_table(index) = 4 * decision_table(index);
        else
            % From 4.121, the incremental energy is 3dB for each bit.
            % However, note this only works when the gap for the two loads
            % (before and after the bit increment) are the same. Thus, care
            % must be taken for QAM, for which different gaps might be
            % applicable if QAM-Cross is used for odd number of bits and
            % QAM-SQ is used for even number of bits.
            decision_table(index) = 2*decision_table(index);
        end
    end
end

end