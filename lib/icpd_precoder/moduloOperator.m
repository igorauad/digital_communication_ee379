function [ Gamma_M ] = moduloOperator( X, M, D, flag )
%MODULO OPERATOR Non linear function defined on an M-ary (PAM or QAM) input
%   constellation with uniform spacing d.
%
% Inputs:
%   X - Vector with QAM or PAM symbols
%   M - Vector with the number of points in the constellation of each subchannel
%   D - Vector with the minimum distances for each subchannel
%   flag - Optional modes of operation
%        * "Hermitian" - A vector with N/2 + 1 should be passed to be
%                        interpreted as half of an Hermitian vector. The output
%                        is the Hermitian vector.
%
bypass = 0;
nSymbols = size(X, 2);

if (~bypass)
    M = M(:);
    D = D(:);

    if (nargin > 3 && strcmp(flag, 'Hermitian'))
        % When the "Hermitian" flag is asserted, assume all entries, except the
        % ones for DC and Nyquist are complex. For those, the operator should
        % be modulo sqrt(M) for each dimension.

        N = 2*(size(X,1) - 1); % FFT Size

        % Order per dimension
        M_bar = [M(1); sqrt(M(2:N/2)); M(N/2+1)];
        % Minimum distances per dimension are the same
    else
        N = size(X,1); % Number of rows in X
        % Assume all subsymbols are complex
        M_bar = sqrt(M);
    end

    Gamma_M_real = real(X) - diag(M_bar .* D) * ...
        floor(inv(diag(M_bar.*D)) * ...
        (real(X) + repmat(M_bar.*D/2, [1 nSymbols]) ));
    Gamma_M_imag = imag(X) - diag(M_bar .* D) * ...
        floor(inv(diag(M_bar.*D)) * ...
        (imag(X) + repmat(M_bar.*D/2, [1 nSymbols]) ));

    % For a single symbol (nSymbols = 1), the above is equivalent to:
    %     Gamma_M_real = real(X) - M_bar.*D.*floor((real(X)+ M_bar.*D/2)./(M_bar.*D));
    %     Gamma_M_imag = imag(X) - M_bar.*D.*floor((imag(X)+ M_bar.*D/2)./(M_bar.*D));
    Gamma_M = Gamma_M_real + 1j*Gamma_M_imag;

    % Important: The division by (M_bar.*D) can yield NaN for the
    % tones in which M or D is 0. Correct this by:
    Gamma_M(M_bar==0 | D==0) = 0;

    if (nargin > 3 && strcmp(flag, 'Hermitian'))
        Gamma_M = [Gamma_M; flipud(conj(Gamma_M(2:end-1,:)))];
    end
else
    Gamma_M = X;
    Gamma_M = [Gamma_M; flipud(conj(Gamma_M(2:end-1,:)))];
end



end

