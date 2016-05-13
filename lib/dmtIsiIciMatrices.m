function [ Hisi, Hici, Hcirc ] = dmtIsiIciMatrices(p, n0, nu, tau, N, preCursor, windowing)
% Compute the time-domain ISI anc ICI matrices
% ---------------------------------------------
% [ Hisi, Hici, Hcirc ] = dmtIsiIciMatrices(h, n0, nu, tau, N, ...
%                                           preCursor, windowing)
%
%   Input Parameters
% p                 Channel Pulse Response
% n0                CIR delay (index of its peak)
% nu                Cyclic Prefix Length
% tau               Cyclic Suffix length
% N                 FFT size
% preCursor         When pre-cursor ICPD should be accounted
% windowing         When spectrum shaping windowing should be accounted
%
% Outpu:
%   Hisi   ->  Time-domain ISI matrix
%   Hici   ->  Time-domain ICI matrix
%   Hcirc  ->  Circulant matrix for the Ideal channel

%% Initialization

% Force column vector:
p = p(:);

% Channel length and dispersion, respectively:
Lh = length(p);
L  = Lh-1;

h_padded = [p;  zeros(N - Lh, 1)];

%% Computation

% Preallocate
Hici = zeros(N,N);
Hisi = zeros(N,N);

if (windowing)

    % Number of affected samples
    nAffected = L - (nu - tau) - n0;

    % Window
    dmtWindow = designDmtWindow(N, nu, tau);

    % Assert conditions
    %
    % The only situation in which the function should return with an error
    % is when the ICPD types that are to be mitigated do not exist. Namely,
    % when pre-cursor and post-cursor do not occurr and preCursor is
    % asserted, or when post-cursor does not occurr.
    %
    % For type 1, there are just 2 possibilities:
    % 1) Post-cursor occur
    % 2) Post-cursor does not occur
    %   Action: return Error
    %
    % For type 2, there are 4 possibilities:
    % 1) Pre-cursor does not occur, but post-cursor does
    %   Action: force change to type 1
    % 2) Pre-cursor occurs, but post-cursor does not
    %   Action: just warn the user that post-cursor was turned off
    % 3) Neither of the two occurr
    %   Action: return Error
    % 4) Both occurr
    if (preCursor)
        % Check the energy prior to the cursor
        preCursorEnergy = sum(abs(p(1:(n0-1))).^2);
        if (preCursorEnergy <= 1e-5)
            warning('Pre-cursor ICPD is irrelevant\n');
            warning('Mitigation changed to post-cursor only\n');
            % Filter #1 - change (downgrade) type:
            preCursor = 0;
            % Filter #3 - continue only if post-cursor ICPD in fact exists:
            assert(nAffected > 0, ...
                'Post-cursor ICPD does not occur');
        end
    else
        % Filter #2 for type 1
        assert(nAffected > 0, ...
            'Post-cursor ICPD does not occur');
    end

    % Filter #2 for type 2:
    % If the number of samples affected by Post-cursor ICPD is positive,
    % add the matrices, otherwise, turn post-cursor mitigation off.
    if (nAffected > 0)
        Ht = toeplitz([p(end) zeros(1, L-(nu-tau)-n0-1)], ...
            flipud(p((nu - tau + n0 + 2):end)));

        % Window
        W_w = diag(dmtWindow(end-nAffected+1:end));

        % Post-cursor ISI
        if (n0 >= tau)
            Hisi( 1:(nAffected), (N - L + nu + 1):(N + tau - n0)) = Ht*W_w;
        else
            shiftBy = tau - n0;
            Hisi(1:(nAffected), (N - L + nu + 1):(N + tau - n0) ) = Ht*W_w;
            Hisi = Hisi(:,shiftBy+1:end);
            Hisi = circshift(Hisi, [0, shiftBy]);
        end

        % Post-cursor ICI
        Hici(1:(nAffected),(N - L + 1): (N -nu + tau - n0)) = -Ht*W_w;
    else
        % The only way to reach the following warning is when type==2.
        warning('Post-cursor ICPD mitigation turned off for line %d\n', ...
            iLine);
        fprintf('Mitigating only pre-cursor ICPD\n');
    end

    % Pre-cursor ICI

    if (preCursor) % Pre-cursor ICI is also mitigated
        % Preallocate
        HpreIci = zeros(N, N);

        H_tilde_t = toeplitz(p(1:n0), [p(1) zeros(1, n0 - 1)]);

        W_tilde_w = diag(dmtWindow(1:n0));

        H_tilde_t_windowed = H_tilde_t * W_tilde_w;

        HpreIci(end-n0+1:end,end-n0+1:end) = - H_tilde_t_windowed;

        % Add pre-cursor ICI matrix to the post-cursor ICI matrix:
        Hici = Hici + HpreIci;
    end

else

    % Samples affected by ISI/ICI:
    nAffected = Lh - nu -n0 -1;
    assert(nAffected > 0, 'Channel length shorter than nu')

    %ISI Matrix
    Ht = toeplitz([p(end) zeros(1,Lh-nu-n0-2)], flipud(p((nu+n0+2):end)));

    % Check assumption of n0 > tau
    if (n0 >= tau)
        Hisi( 1:(nAffected),(N-(nAffected)-(n0-tau)+1):(N-(n0-tau)) ) = Ht;
    else
        shiftBy = tau - n0;
        Hisi( 1:(nAffected),(N-(nAffected)-(n0-tau)+1):(N-(n0-tau)) ) = Ht;
        Hisi = Hisi(:,shiftBy+1:end);
        Hisi = circshift(Hisi, [0, shiftBy]);
    end

    %ICI Matrix
    Hici(1:(nAffected),(N-(nAffected)+1-(nu+n0)):N-(nu+n0)) = -Ht;
    Hici = sparse(Hici);

end

%Circulant Matrix (note h has to be a column vector)
Hcirc = toeplitz(h_padded, flipud(circshift(h_padded, -1)));

end

