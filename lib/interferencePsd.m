function [ psd ] = interferencePsd(h, N, Lcp, Lcs, n0, Px, type, WINDOWING, fs)
% Calculates the ICP Distortion(ISI/ICI) PSD
% ---------------------------------------------------
% [ psd ] = calculateInterferencePsd(h, N, Lcp, Lcs, n0, Px, PLOT, txNo,...
% WINDOWING, fs, type)
%
% This script implements the ISI PSD Equation described on the internal
% report "Intersymbol and Intercarrier Interference on DMT Systems".
%
%   Input Parameters
% h                     - Channel Impulse Response
% N                     - FFT size
% Lcp                   - Cyclic prefix length
% Lcs                   - Cyclic suffix length
% n0                    - Rx Cursor (CIR peak index)
% Px                    - Tx signal average power
% WINDOWING             - Flag to indicate windowing is being used
% fs                    - Sampling Frequency
% type                  - Scenario for the PSD:
%                           - (0) Post-cursor + Pre-cursor ISI/ICI
%                           - (1) Post-cursor ISI/ICI
%                           - (2) Pre-cursor ISI/ICI
%                           - (3) Pre-cursor ISI
%
%   Output
% psd                   - The interference PSD

if (nargin < 8)
    WINDOWING = 0;
end

if (nargout > 0 || nargin <= 6)
    plotOverAnalogFreq = 0;
else
    plotOverAnalogFreq = 1;
end

% Ensure column-vector
h = h(:);

%% Definitions

Lh = length(h);
L = Lh - 1;
symbolSize = N + Lcp + Lcs;
nTxPower = size(Px,2); % Number of simulated Tx Power

% Window:
% The tapered part of the window uses a raised-cosine function.
if(WINDOWING)
    window = designDmtWindow( N, Lcp, Lcs );
end

%% Allocation
psd = zeros(N,nTxPower);
isi_summation = zeros(N,1); % Summation from psd formula
pre_isi_summation = zeros(N,1); % Summation from psd formula


%% Post-cursor ISI/ICI

if Lcp < Lh-1 % If post-cursor ISI and ICI occurs:

    if(WINDOWING)
        for phi = (n0 + (Lcp - Lcs) + 1) : L
            mu = phi - (n0 + (Lcp - Lcs) + 1);
            w_index =  symbolSize - mu;
            Hphi = fft( h( (phi+1) : end ), N, 1); %FFT of the channel tail
            % phi + 1 - to convert the MATLAB indexing
            isi_summation = isi_summation + ...
                (abs(window(w_index)).^2)*(abs(Hphi).^2);
            % The previous line partially implements equation (10) from
            % the reference paper.
        end
    else
        for phi = (n0 + Lcp + 1) : (Lh-1)
            Hphi = fft(h( (phi+1) : end ), N, 1); %FFT of the channel tail
            % phi + 1 - to convert the MATLAB indexing
            isi_summation = isi_summation + (abs(Hphi).^2);
            % The previous line partially implements equation (10) from
            % the reference paper.
        end
    end

end

%% Pre-cursor ISI/ICI

if (WINDOWING)
    % For the pre-cursor ISI/ICI calculation, when windowing+overlap is
    % used, the number of affected samples is the same as if Lcs = 0;
    Lcs = 0;
end

% If pre-cursor ISI and ICI occurs:
if n0 > Lcs
    for phi = (-Lcp) : (n0 - Lcs - Lcp -1) % should be another iterator
        mu = phi + Lcp;
        m = mod(phi,N);
        upperLimDft = n0 - Lcs -1 -mu + 1; % The last +1 is on MATLAB only
        if(WINDOWING)
            pre_isi_summation = pre_isi_summation + ...
             (abs(window(mu + 1))^2)*(abs(fft(h(1:upperLimDft), N, 1)).^2);
        else
            pre_isi_summation = pre_isi_summation + ...
                    abs(fft(h(1:upperLimDft), N, 1)).^2;
        end
    end
end

%% Additional processing

% Calculate the interference psd for each Tx power passed as input to
% the function
for ii = 1:nTxPower;
    % For the post-cursor ISI+ICI psd, multiply the CIR Tail FFTs
    % by the Tx signal variance
    SpostIsi = (Px(ii)) * isi_summation;

    % For the pre-cursor ISI+ICI psd, multiply the CIR Head FFTs
    % by the Tx signal variance as well
    SpreIsi = (Px(ii)) * pre_isi_summation;

    % The psd for the ICI is equal to the PSD for the ISI, both for
    % pre-cursor and post-cursor intereferences. Thus, the PSD for
    % ISI+ICI is twice the PSD for ISI alone:
    SpostIcpd = 2 * SpostIsi;
    SpreIcpd = 2 * SpreIsi;

    % Finally, divide by N to comply with the periodogram formula for
    % the classical FFT/IFFT transform convention:
    switch type
        case 0
            psd(:,ii) = (SpostIcpd + SpreIcpd);
        case 1
            psd(:,ii) = (SpostIcpd);
        case 2
            psd(:,ii) = (SpreIcpd);
        case 3
            psd(:,ii) = (SpreIsi);
    end

     psd(:,ii) =  psd(:,ii)/N;

end


%% Plot

% Plot intereference PSD for each Tx power on the same figure.
% Note the figure handler (h_fig) is used to plot on the same figure
if(nargout == 0)
    figure
    if (plotOverAnalogFreq)
        plot(10*log10( ( (psd/fs) )*1e3 ), 'r') % in dBm
    else
        plot(10*log10( ( (psd) )*1e3 ), 'r') % in dBm
    end
    hold on
    title('ICP Distortion PSD')
    xlabel('Subcarrier')
    ylabel('PSD (dBm/Hz)')
    set(gca, 'XTick', [0:N/8:(N-N/8) (N-1)])
    set(gca, 'XTickLabel', [0:N/8:(N-N/8) (N-1)])
    set(gca, 'XLim', [0 (N/2 - 1)])
    grid on
end

end

