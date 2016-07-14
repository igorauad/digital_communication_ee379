function qamConst_o = qamHybridConstellation(M)
% Generates a hybrid QAM square constellation with M symbols.
% If log2(M) is even, returns a regular SQ-QAM constellation.
%
%   A Hybrid QAM constellation for b bits per symbol is obtained by
%   designing a QAM-SQ constellation of b+1 bits per symbol and then
%   selecting every other symbol in the constellation.
%
%   While regular SQ-QAM constellation has average symbol energy:
%       Ex = d^2 * ((2^b - 1)/6)
%   The Hybrid QAM has average energy:
%       Ex = d^2 * ((2^(b+1) - 1)/6),
%   where d in this case is the original distance for the SQ-QAM (d=2 on
%   MATLAB) constellation of b+1 bits per symbol (not the minimum distance
%   in the Hybrid QAM constellation, which is sqrt(2)*d.
%   The average energy for the Hybrid QAM can be rewritten as:
%       Ex = d^2 * ((2*M - 1)/6),
%   and since d = dmin/sqrt(2), it is equivalently:
%       Ex = dmin^2 * ((2*M - 1)/12).
%   Therefore, the minimum distance is given by:
%       dmin = sqrt( (12 * Ex)/(2*M - 1) )
%            = sqrt( (24 * Ex_bar)/(2*M - 1) )
%   Then, the Q function argument "dmin/2*sigma" is equivalent to:
%       (1/2) * sqrt( (24 * Ex_bar)/(2*M - 1) ) / sigma
%       = sqrt( (6 * Ex_bar)/ ((2*M - 1)*sigma^2) )
%       = sqrt( (6 * SNR) / (2*M - 1))
%   Finally, considering the normalized SNR to be:
%       SNRnorm = SNR / (M - 1)
%    the Q-function argument is only aproximately sqrt(3*SNRnorm) for large
%    M, when the "-1" term disapears. This reveals that the gap for hybrid
%    and square QAM is approximately equal for large M. However, note that
%    the average number of nearest neighbors is different, so slightly
%    different Pe's are obtained, as explained in the sequel.
%
%       The probability of error is also different. Consider three
%   possibilities in terms of nearest neighbors: the corner points, which
%   have only a single nearest neighbor, the edge points, which have 2
%   nearest neighbors and the inner points, which have 4 nearest neighbors.
%   For every Hybrid QAM, there are:
%       - 2 corner points
%       - 4*(sqrt(M/2) - 1) edge points
%       - M - 2*sqrt(2*M) + 2 inner points
%
%   Hence, it can be shown that the average number of nearest neighbors is:
%
%       Ne = 4*(1 - sqrt(2/M) + 1/(2*M)),
%
%   which asymptotically (for large M) reaches the same Ne=4 of a regular
%   square QAM constellation. For M = 128, for example, Ne = 3.5156, and
%   for M = 32, Ne = 3.0625. For comparison, recall that for QAM-SQ:
%
%       Ne = 4*(1 - 1/sqrt(M)),
%
%   and for QAM-Cross:
%
%       Ne = 4*(1 - 1/sqrt(2*M))
%
%   By plotting the three Ne's, for Hybrid, Cross and SQ QAM, one can see
%   that the Ne for Hybrid QAM is always slightly lower than the Ne for
%   square QAM and that the Ne for Cross is the highest (always slightly
%   higher than the one for SQ).
%
%

b = log2(M);

if (mod(b,2) ~= 0)
    % For odd b, form a b+1 SQ constellation, then select only every other
    b_used = b + 1;
else
    b_used = b;
end

% In-phase/Quadrature bits
b_i = ceil(b_used/2); % number of bits for in-phase component
b_q = b_used-b_i;     % remaining bits are for quadrature component
M_i = 2^b_i;
M_q = 2^b_q;

% QAM is formed from a Cartesian product of PAM constellations:
pamConst_i = -(M_i-1):2:M_i-1; % In-phase PAM constellation
pamConst_q = -(M_q-1):2:M_q-1; % Quadrature PAM constellation

% Preallocate
qamConst = zeros(M_q, M_i);
qamConst_o = zeros(M_q, M_i/2);

% Full QAM constellation
for i=1:M_q
    for k=1:M_i
        qamConst(i,k) = pamConst_i(k) + 1j*pamConst_q(i);
    end
end

%% Filter full QAM constellation for odd b and output

% For odd b, select every other point:
if (mod(b,2) ~= 0)
    fixedIndexes = 2:2:M_i;
    for i = 1:M_i
        indexes = fixedIndexes - mod(i, 2);
        qamConst_o(i,:) = qamConst(i, indexes);
    end
    % Scale for a minimum distance of 2:
    qamConst_o = (1/sqrt(2)) * qamConst_o;
    qamConst_o = transpose(qamConst_o(:));
else
    qamConst_o = transpose(qamConst(:)); %make it a row vector
end