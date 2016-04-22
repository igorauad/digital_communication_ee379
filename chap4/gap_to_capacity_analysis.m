% A script to understand the uncoded PAM GAP to capacity
clearvars, clc

N     = 1;           % dimensions
b     = 1:10;        % bits
b_bar = b / N        % bits per dimension
M     = 2.^b;        % modulation order
rho   = 2*b_bar;     % Spectral efficiency

% Define the target modulation:
encoder_type = 2; % 0 - PAM; 1 - QAM-SQ; 2 - QAM-Cross

% Define one of the two searches:
%   'Gap' - Find the Gap for a given Pe per dimension.
%   'Pe'  - Find the Pe per dimension for a given Gap.
search  = 'Pe';
% If 'Gap', first define the target Pe:
Pe_bar  = 1e-6;
% Else if 'Pe', first define the Gap in dB:
gap_db  = 8.8;

%% Solution

if (strcmp(search, 'Gap'))

    % Discover the SNR required to achieve the target Pe for each spectral
    % effiency. Recall that, for PAM, the NNUB is exact and given by:
    %
    %   Pe_bar = 2*(1 - 1/(2.^b_bar)) qfunc( sqrt( 3 * SNRnorm ) ),
    %
    % where the normalized SNR (SNRnorm), a.k.a. the gap to capacity, is
    % given by:
    %
    %   gap = SNRnorm = SNR/(2^rho - 1)
    %
    % For QAM with SQ constellation, the error probability is the same, but
    % the NNUB is not exact. Moreover, for QAM with cross constellation,
    % the NNUB Pe per dimension is different and given by:
    %
    %   Pe_bar = 2*(1 - 1/(sqrt(2)*M))*qfunc(sqrt((3/(31/32)) * SNRnorm)).
    %
    % Two steps are required:
    %
    % - First discover the required argument for the Q function.
    %
    % - Then, obtain SNRnorm by equating the Qfunc arg to its value.


    switch (encoder_type)
        case {0,1}
            % PAM and SQ-QAM
            arg_qfunc = qfuncinv(Pe_bar ./ ( 2*( 1 - (1./(2.^b_bar)))) );
            % The argument is equivalent to sqrt( 3 * SNRnorm )
            SNRnorm = (arg_qfunc.^2) / 3;
        case 2
            % Cross-QAM
            arg_qfunc = qfuncinv(Pe_bar ./ ...
                (2*( 1 - (1./(2.^(b_bar + 0.5))))));
            % The argument is equivalent to sqrt((3/(31/32)) * SNRnorm )
            SNRnorm = (arg_qfunc.^2) / (3/(31/32));
    end

    % Finally, express it in dB:
    fprintf('Gap to Capacity:\n');
    10*log10(SNRnorm)

    % Note the gap is not constant due to the "2*(1 - 1/M)" term in Pe

elseif (strcmp(search, 'Pe'))

    % Convert to linear scale:
    gap = 10^(gap_db/10);
    % Find the achievable Pe_bar

    switch (encoder_type)
        case {0,1}
            Pe_bar = 2 * (1 - (1./(2.^b_bar))) * qfunc( sqrt( 3 * gap ) )
        case 2
            Pe_bar = 2*(1 - (1./(sqrt(2)*M))) * ...
                qfunc(sqrt((3/(31/32)) * gap))
    end

else

    error('Search not supported');

end
