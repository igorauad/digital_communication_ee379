function [ Scale_n, dmin_n ] = dmtSubchanScaling(modulator, modChoice, ...
    En, n_dim )
% Compute the scaling factors for the constellations in each subchannel
% ------------------------------------------------------------------------
%   Compute the factor that scales the average symbol energy in each
%   constellation to the target average energy Ex. Also returns the
%   corresponding constellation minimum distance after scaling.
%
%   Inputs:
%       modulator   -   Cell array with the unique used modulators
%       modChoice   -   Vector with the modulator choices for each
%                       subchannel
%       En          -   Average symbol energy for each subchannel
%       n_dim       -   Number of dimensions for each subchannel
%
%
% Note 2-dimensional subchannels whose bit loading is 1 (i.e. M=2) use a
% QAM constellation equivalent to a 2-PAM rotated by 90 degrees, which
% implies the two dimensions are really used regardless.
%
% For the computation of d_min for each subchannel, note it is equivalent
% to 2*Scale_n(k) for QAM-SQ and sqrt(2)*2*Scale_n(k) for QAM-Hybrid.
%
% Finally, note scaling factors and corresponding minimum distances are
% computed only for the positive half of the spectrum and notice the fact
% that the Hermitian symmetry will double the energy in a certain
% subchannel  is already pre-compensated by allocating only half of the
% subchannel energy budget through the scaling factor.

if (length(n_dim) ~= length(En))
    error('Vector with number of dimensions must be of same length as En');
end

% Preallocate
Scale_n = zeros(length(En), 1);
dmin_n  = zeros(length(En), 1);

% Iterate over subchannels
for k = 1:length(En)
    % Check if the subchannel is loaded
    k_loaded = modChoice(k) > 0;
    %% One-dimensional Subchannels (if used)
    if (n_dim(k) == 1 && k_loaded)
        % Modulation order for the k-th subchannel
        M_k = modulator{modChoice(k)}.M;

        % The last argument should be the Energy per real dimension
        Scale_n(k) = modnorm(...
            modulator{modChoice(k)}.constellation,...
            'avpow', En(k));
        dmin_n(k) = 2*Scale_n(k);
    end
    %% Two-dimensional Subchannels
    if (n_dim(k) == 2 && k_loaded)
        % Modulation order and bit load for the k-th subchannel
        M_k = modulator{modChoice(k)}.M;
        b_k = log2(M_k);

        % The last argument should be the Energy per 2 dimensions.
        % However, since Hermitian symmetry will double the energy, it
        % is chosen as the energy per real dimension.
        Scale_n(k) = modnorm(...
            modulator{modChoice(k)}.constellation,...
            'avpow', 0.5 * En(k));

        % Compute the corresponding minimum distance
        if (mod(b_k, 2) ~= 0 && b_k ~= 1)
            % Odd bit load > 1: Hybrid QAM
            dmin_n(k) = sqrt(2)*2*Scale_n(k);
        elseif (mod(b_k, 2) ~= 0 && b_k == 1)
            % Odd bit load == 1: "rotated" PAM's
            dmin_n(k) = 2*Scale_n(k);
        else
            % SQ-QAM (even bit load)
            dmin_n(k) = 2*Scale_n(k);
        end
    end
end
end

