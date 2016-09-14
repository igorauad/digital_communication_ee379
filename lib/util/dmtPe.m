function [ Pe_bar_n ] = dmtPe( bn, SNR_n, dim )
% Computes the Probability of Error per Dimension for the DMT system
% -------------------------------------------------------------------------
%   Probability of error (Pe) per dimension (Pe_bar) is computed for each
%   subchannel of the DMT system. The probabilities consist of the
%   nearest-neighbor union bound (NNUB) probability of error (also known as
%   union bound estimate). For each subchannel, depending on the bit load,
%   different constellations can be used (Hybrid QAM, QAM-Square or PAM).
%   The following routine infers which constellation is used and computes
%   the NNUB accordingly.
%
%   Inputs:
%  bn         -> Bit load vector for each subchannel
%  SNR_n      -> Per-dimensional SNR for each subchannel
%  dim        -> Vector with the number of dimensions associated to each
%                subchannel

% Preallocate
Pe_bar_n    = zeros(length(dim), 1);

for k = 1:length(dim)
    % "iSubCh" iterates over the loaded subchannels, but does not give the
    % exact DFT tone index.

    M_n       = 2^bn(k);      % Modulation Order
    bn_bar    = bn(k)/dim(k); % Bits per dimension
    SNRnorm_n = SNR_n ./ (2.^(2*bn_bar) - 1);

    if (dim(k) == 2)
        % QAM Nearest-neighbors Union Bound assuming QAM-SQ constellations
        % for both even and odd loads. "Hybrid QAM" constellations are used
        % as "Square" for odd loads.

        % For odd b, Hybrid QAM is used
        % For even b, conventional square QAM (QAM-SQ) is used.
        if ((mod(bn(k),2) ~= 0))
            % Average Number of nearest neighbors
            Ne          = 2*(1 - sqrt(2/M_n) + 1/(2*M_n));
            % Pe per dimension:
            Pe_bar_n(k) = Ne * qfunc(sqrt((6 * SNR_n(k))/(2*M_n -1)));
            % Note the aforementioned q-function argument for Hybrid QAM
        else
            % QAM-SQ
            % Average Number of nearest neighbors
            Ne          = 2 * (1 - 1/(2^bn_bar));
            % Pe per dimension:
            Pe_bar_n(k) =  Ne * qfunc(sqrt( 3*SNRnorm_n(k) ));
        end

    else
        % PAM Nearest-neighbors Union Bound

        % Average Number of nearest neighbors
        Ne = 2 * (1 - 1/(2^bn_bar));
        % Pe per dimension:
        Pe_bar_n(k) = Ne * qfunc(sqrt( 3*SNRnorm_n(k) ));
    end
end

end

