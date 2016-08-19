function [ Gamma_M ] = moduloOperator( X, b_bar, D, flag )
% MODULO OPERATOR Non-linear function defined on an M-ary input
% constellation with uniform spacing d.
%
% Inputs:
%   X      -  Vector with QAM or PAM symbols
%   b_bar  -  Vector with the number of bits per dimension for each
%             subchannel
%   D      -  Vector with the minimum distances for each subchannel
%   flag   -  Optional modes of operation
%             * "Hermitian" - A vector with N/2 + 1 should be passed to be
%                             interpreted as half of an Hermitian vector.
%                             The output is the Hermitian vector.

%% Initialization

% Debug purpose:
bypass   = 0;           % Bypass the modulo operator

% Ensure that the arguments are used as column vectors:
b_bar    = b_bar(:);
D        = D(:);

% Constants
nSymbols = size(X, 2);
M_bar    = 2.^b_bar;    % Constellation order per dimension

% Used tones:
usedTones = M_bar~=0 & D~=0;

% Ignore the constellation orders and minimum distances of the tones that
% are not used
M_bar = M_bar(usedTones);
D     = D(usedTones);
% Important: In the ensuing formulation, there is an inversion
% inv(diag(M_bar.*D). If a vector with M_bar=0 or D=0 is adopted, then the
% result will yield a matrix NaNs. By extracting only the used tones, this
% problem is avoided.

% Preallocate the modulo-operated symbols
Gamma_M   = zeros(size(X));

%% Main
if (~bypass)

    % Apply modulo-M in the real part of the symbols
    Gamma_M_real = real(X(usedTones, :)) - diag(M_bar .* D) * ...
        floor(inv(diag(M_bar.*D)) * ...
        (real(X(usedTones, :)) + repmat(M_bar.*D/2, [1 nSymbols]) ));

    % Apply modulo-M in the imaginary part of the symbols
    Gamma_M_imag = imag(X(usedTones, :)) - diag(M_bar .* D) * ...
        floor(inv(diag(M_bar.*D)) * ...
        (imag(X(usedTones, :)) + repmat(M_bar.*D/2, [1 nSymbols]) ));

    % For a single symbol (nSymbols = 1), the above is equivalent to:
    %     Gamma_M_real = real(X) - ...
    %         M_bar.*D.*floor((real(X)+ M_bar.*D/2)./(M_bar.*D));
    %     Gamma_M_imag = imag(X) - ...
    %         M_bar.*D.*floor((imag(X)+ M_bar.*D/2)./(M_bar.*D));

    % Form the complex modulo-operated symbols:
    Gamma_M(usedTones, :) = Gamma_M_real + 1j*Gamma_M_imag;

    % In case the Hermitian flag is assert, output the derived Hermitian
    % vector:
    if (nargin > 3 && strcmp(flag, 'Hermitian'))
        Gamma_M = [Gamma_M; flipud(conj(Gamma_M(2:end-1,:)))];
    end
else
    % Bypass the symbols to the output
    Gamma_M = X;

    if (nargin > 3 && strcmp(flag, 'Hermitian'))
        Gamma_M = [Gamma_M; flipud(conj(Gamma_M(2:end-1,:)))];
    end
end



end

