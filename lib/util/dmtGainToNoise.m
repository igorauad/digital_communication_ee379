function [ gn ] = dmtGainToNoise(p, dmt, Rxx)
% Computes the Gain-to-Noise Ratio for each subchannel
%
% Inputs
%  p     -> Effective channel pulse response
%  dmt   -> DMT Object
%  Rxx   -> (Optional) Tx Autocorrelation Matrix

% Equalizer Types
EQ_NONE      = 0;
EQ_TEQ       = 1;
EQ_FREQ_PREC = 2;
EQ_TIME_PREC = 3;

% Parameters
Nfft      = dmt.Nfft;
nu        = dmt.nu;
n0        = dmt.n0;
tau       = dmt.tau;
Ex_bar    = dmt.Ex_bar;
N0_over_2 = dmt.N0_over_2;
windowing = dmt.windowing;
eqType    = dmt.equalizer;

if (nargin < 3)
    % When the function is invoked without the Rxx parameter,
    % it is implied that the transmit sequence is uncorrelated.
    Rxx = Ex_bar * eye(Nfft);
    fprintf('ICPD PSD computed assuming flat energy load\n');
end

%% Frequency domain response (use the effective pulse response in the FEQ)
H = fft(p, Nfft);
% Store only the response at the used indices of the FFT
Hn = H(dmt.iTonesTwoSided);

%% ICPD PSD
% Compute the ICPD PSD in order to account for it in the computation of
% the bit loading.
%
% Note: for the frequency and time-domain ICPD precoders (which should
% fully cancel the ICPD), the post-cursor ICPD is assumed to be null.
% In this case, only pre-cursor ICPD is taken into account.

switch (eqType)
    case {EQ_TEQ, EQ_NONE}
        % Compute the ISI matrices
        [Hisi, ~, ~, HpreIsi, ~] = ...
            dmtIsiIciMatrices(p, n0, nu, tau, Nfft, windowing);

        % Total ICPD PSD
        S_icpd = icpdPsdMtx(Hisi, HpreIsi, Rxx, Nfft);
    case {EQ_FREQ_PREC, EQ_TIME_PREC}
        % In this case, assume post-cursor ICPD is fully cancelled and
        % consider only the pre-cursor ICPD.

        % Compute the ISI matrices
        [~, ~, ~, HpreIsi, ~] = ...
            dmtIsiIciMatrices(p, n0, nu, tau, Nfft, windowing);

        % PSD given by pre-cursor ICPD
        S_icpd = icpdPsdMtx(zeros(Nfft), HpreIsi, Rxx, Nfft);
end

%% Gain to noise
switch (eqType)
    case EQ_TEQ
        % Notes:
        %   # 1) The water-filling solution assumes no ISI/ICI. Even though
        %   the TEQ constrains the pulse response energy to a portion that
        %   can be covered by the guard band (commonly referred to the
        %   "window" of the SIR), the out-of-window response of the
        %   shortened response may still be significant and introduce
        %   non-negligible ISI/ICI.
        %   # 2) Note the feed-foward TEQ at the receiver shapes the
        %   spectrum of the noise, so the noise PSD becomes |H_w|^2 * N0/2,
        %   where H_w is given below:
        H_w = fft(dmt.w, Nfft);
        %   # 3) Meanwhile, the transmit signal is subject to the
        %   compounded response of the channel + equalizer, namely the
        %   effective response p_eff, whose DFT is "H".
        %   # 4) Then, assuming that the ICPD is uncorrelated to the noise,
        %   the gain to noise ratio at the receiver becomes:
        gn = (abs(H).^2)./((N0_over_2 * abs(H_w).^2) + S_icpd.');
        %   Store only the used tones
        gn = gn(dmt.iTonesTwoSided);
        %   Note that if ICPD was not accounted, since H_eff = Hn .* H_w,
        %   the above gain-to-noise ratio would tend to be equivalent to
        %   (for all non-zero H_W):
        %       gn = (abs(Hn).^2) / N0_over_2;
        %   Ultimately, the gain-to-noise ratio is not being affected by
        %   the TEQ in the expression. This can be the fallacious in the
        %   model.
    otherwise
        gn = (abs(Hn).^2) ./...
            (N0_over_2 + S_icpd(dmt.iTonesTwoSided).');
end

end

