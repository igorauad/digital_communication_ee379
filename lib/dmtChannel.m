function [y] = dmtChannel(txSeq, p, dmt)
% DMT AWGN Channel
% [y] = dmtChannel(txSeq, p, dmt)
%
% Inputs:
%  txSeq   -> Transmit sequence
%  p       -> Channel pulse response
%  dmt     -> Struct with DMT parameters
%
% Output:
%  y       -> Channel output sequence

% Parameters
N0_over_2 = dmt.N0_over_2;

%% Main
y = conv(txSeq, p);
% Note:
%   In contrast to the derivation of Chapter 3, here the scaling of Ts
%   is not used. The reason is that "u" comes from an orthonormal
%   expansion (the normalized IFFT) and p also satisfies the inner
%   product invariance, due to the more involved explanation in the
%   sequel.
%
%   The text mentions that "the notational use of P for the channel
%   matrix suggests that any anti-alias analog filters at transmitter
%   and receiver have been convolved with the channel impulse response
%   h(t) and included in the discrete-time response of the matrix
%   channel". The model is given in (4.185).
%
%   So let us further interpret that: once we have a measurement of the
%   channel impulse response, we have samples and, therefore,
%   coefficients of an orthogonal expansion. Thus, convolution or inner
%   product (for energy computation) can not be applied directly,
%   only with a Ts factor in front. Nonetheless, in practice the
%   samples would be obtained by first passing the signal through an
%   anti-alias LPF, ideally with unitary-energy. Such filter
%   effectively scales each "pre-ADC" sample by sqrt(Ts). To understand
%   that, note a unitary energy anti-alias filter has continuous-time
%   response:
%
%       1/sqrt(Ts) * sinc(t/Ts)
%
%   Then, by sampling at t = kTs, the only non-zero value of the
%   sampled sequence is "1/sqrt(Ts)". Finally, by convolving the
%   "pre-ADC" samples (from an orthogonal expansion) with the receive
%   filter samples, one obtains:
%
%       h = Ts * conv(rx_filter, h)
%         = Ts * (1/sqrt(Ts)) * h
%         = sqrt(Ts) * h_pre_adc
%
%   Essentially, when we sample the channel, we get h, instead of
%   h_pre_adc. The channel impulse response energy can be computed by
%   either "Ts * norm(h_pre_adc)^2" or "norm(h)^2", because both are
%   equivalent. Similarly, the sampled response "h" can be used
%   directly in the convolution, without a Ts factor in front.

% Add noise
noise = (sqrt(N0_over_2) * randn(length(y),1));
y = y + noise;
% Important considerations:
%
% First, recall the noise continuous-time PSD coincides with the noise
% energy per dimension. Second, remember that the sinc functions in the
% orthogonal expansion of the sampling theorem have energy 1/2W, so the
% variance of each real and imaginary coefficient in the noise
% expansion must be scaled up by 2W from the noise energy N0/2 per
% degree of freedom. Since AWGN is flat, multiplication of N0/2
% (two-sided PSD) by 2W yields the total transmit power. Hence, if "y"
% consisted of samples, the target variance for the "randn" sequence
% would be the noise power N0_over_2 * 2W. However, the catch here is
% that y does not represent the samples
end