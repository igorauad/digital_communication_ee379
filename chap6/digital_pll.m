clearvars, close all, clc
% Digital PLL Implementation
%
% Implements a DSP-based Phase-Locked Loop featured with a PI controller.
%
% The loop has three main components:
%   1) Direct Digital Synthesizer (DDS): iteratively generates a a complex
%   sinusoid by computing exp(j*phi_loop), where "phi_loop" is its
%   instantaneous phase that continuously grows based on a configurable
%   increment value. The angle increment added to the DDS (accumulated in
%   "phi_loop", i.e. its "phase accumulator") at each clock cycle (sample
%   period) is the nominal increment 2*pi*f0/fs, where f0 is the nominal
%   frequency of the DDS, plus a varying correction term. The latter term
%   should vary considerably until the instant after which the DDS output
%   is synchronized to the input, namely until the loop "locks".
%       Note that the DDS is a numerically controlled oscillator (NCO).
%   Just like a VCO, the derivative of its output signal's phase is
%   proportional to the signal passed to it as input, or, equivalently, its
%   phase output is the integral of the input. The NCO outputs a complex IQ
%   sinusoid:
%
%           s_loop(n) = exp(j*(2*pi*(f0/fs)n + theta_hat(n))),     (1)
%
%   where "theta_hat" is a phase that should track the phase "theta" from
%   an input of the form:
%
%           s_in(n) = exp(j*(2*pi*(f0/fs)*n + theta(n))).          (2)
%
%   This phase "theta_hat" comes from a "phase accumulator" that
%   continuously integrates the NCO input signal (the filtered error
%   explained below). Hence,
%
%           theta_hat(n) = sum_{n=0}^{n} filtered_error(n).        (3)
%
%   Clearly, if the input signal s_in has a fixed frequency offset f_offset
%   with respect to the nominal frequency, then its "theta" term is:
%
%           theta(n) = 2*pi*(f_offset/fs)*n                        (4)
%
%   In this case, the filtered error has to converge to
%   "2*pi*(f_offset/fs)", such that the NCO output phase difference in (3)
%   can track the input phase difference term of (4).
%
%   2) Phase Detector: extracts the phase error from the input complex
%   exponential and the loop complex exponential
%
%   3) Loop Filter: filters the phase error such that, after sufficient
%   iterations, it is driven to zero (for null frequency offset and Ki=0 or
%   for non-null frequency offset, but also non-null Ki=0) or a constant
%   value (for constant frequency offset and Ki=0). The output of the loop
%   filter (called filtered phase error) is accumulated by the NCO (DDS),
%   as described in (3), and composes the correction term "theta_hat" of
%   Equation (1).
%
%   In the case when Ki > 0, namely the integral controller is enabled,
%   since the phase error from the phase detector is integrated, once the
%   loop converges to a "locked" state, the integrator has already
%   accumulated all the phase error required to output a constant value to
%   the NCO accumulator. In the previous example of a fixed frequency
%   offset, the integral filter output would rise and converge to
%   2*pi*(f_offset/fs) while the phase error (input to the integral filter)
%   would converge to zero.
%
% Suggested experiments:
%   1) Set Ki = 0 and configure the frequency offset to zero (y_ppm = 0).
%   See that the PLL can correct the constant phase error.
%   2) Set Ki = 0 and configure a large enough frequency offset. See that
%   the phase error converges to a constant value 2*pi*f_offset/fs.
%   3) With the same configuration, increase Kp and see that the
%   steady-state phase error to which the PLL converges becomes smaller,
%   at the expense of a noisier phase output due to the fact that the PLL
%   bandwidth increases with Kp.
%   4) Still with the same configuration (Ki = 0), reduce Kp such that the
%   phase error can exceed pi, case in wich the "atan" phase detector
%   suffers from ambiguity. In another words, use a Kp that allows a phase
%   error exceeding the pull-range of the PLL. Maybe increase the number of
%   iterations (nIterations) for better visualization.

nIterations = 2e4;
fs          = 1e6;   % Nominal sampling frequency
y_ppm       = 5e4;    % Fractional frequency offset in ppm
f0          = 1e3;   % Nominal clock frequency
pn_var      = 1e-9;  % Phase noise variance
Kp          = 0.0005;  % Proportional Constant
Ki          = 0;  % Integral Constant

%% Derived parameters
y        = y_ppm * 1e-6;  % Fractional frequency offset in ppm
f_offset = y * f0;        % Absolute freq error (in Hz)
F_offset = f_offset/fs;   % Normalized absolute frequency offset

% Nominal Phase Increment (used in the loop phase accumulator)
delta_phi_c = 2*pi*f0/fs;

% PLL pull range
% Normalize frequency offset must be smaller than Kp*pi, namely:
%
%   F_offset < Kp*pi
if (F_offset > Kp*pi)
    warning('Frequency offset is not within the pull range\n');
end

%% Generate input signal
% The input signal obviously does not depend on the loop processing, so we
% can generate it before simulating the loop.
%
% The input here is considered to be a pure complex sinusoid (exponential)
% slightly corrupted by phase noise and eventually with frequency offset.
% We can simulate such an exponential by considering that in each clock
% cycle (sample period), the exponential grows by the following increment:
delta_phi_in = 2*pi*(f0 + f_offset)/fs;
% which when combined to the following phase noise random sequence:
phase_noise  = sqrt(pn_var)*randn(nIterations, 1);
% results in the actual instantanous phase:
phi_in       = (0:nIterations-1).' * delta_phi_in + phase_noise;
% Finally, generate the exponential:
s_in         = exp(1j * phi_in);

%% Loop

% Preallocate
s_loop             = zeros(nIterations, 1); % DDS Complex Value
phi_loop           = zeros(nIterations, 1); % DDS Phase Accumulator
dds_mult           = zeros(nIterations, 1); % Conjugate Product
phi_error          = zeros(nIterations, 1); % Phase Error
phi_error_filtered = zeros(nIterations, 1); % Filtered Phase Error

% Initialize the loop DDS to a random phase
phi_loop(1)  = sqrt(pi)*randn;

% Initialize integral filter output
integral_out = 0;

for i = 1:nIterations

%% Loop DDS:
s_loop(i) = exp(1j*phi_loop(i));

%% Phase Detector
% Phase error is obtained by computing the angle of the conjugate product
%   Input                 : e^(j*phi_in)
%   Loop                  : e^(j*phi_loop)
%   Conjugate Product     : e^[j*(phi_in - phi_loop)]
%   ---------------------------------------------------
%   Conj. Product Angle   : phi_in - phi_loop

% Multiply the input complex exponential by the conjugate of the loop DDS:
dds_mult(i) = s_in(i) * conj(s_loop(i));

% Phase error:
phi_error(i) = angle(dds_mult(i));

%% Loop Filter
% The loop filter consists of a Proportional+Integral (PI) controller

% Proportional term
proportional_out = phi_error(i)*Kp;

% Integral term
integral_out     = phi_error(i)*Ki + integral_out;

% PI Filter output:
phi_error_filtered(i) = proportional_out + integral_out;

%% Update the phase accumulator of the loop DDS
% The angle of the DDS complex exponential in the next clock cycle results
% from the sum between the nominal phase increment and the filtered phase
% error term:
phi_loop(i+1) = phi_loop(i) + delta_phi_c + phi_error_filtered(i);

end

%% Expected filter steady-state value
% The values to which the phase error and the filtered phase error converge
% depend on the loop order. For a second-order loop (with a PI controller),
% assuming an input with fixed frequency offset, the phase error converges
% to 0. In contrast, for a first-order loop (with a Proportional controller
% only), the phase error converges to 2*pi*(f_offset/fs)/Kp. In both cases,
% the filter output converges to 2*pi*(f_offset/fs).
if (Ki == 0)
    phi_error_ss_expected = 2*pi * F_offset / Kp;
else
    phi_error_ss_expected = 0;
end

phi_error_filtered_ss_expected = 2*pi * F_offset;

%% Performance
% Input vs. Output Instantaneous Phase
figure
plot(phi_in)
hold on
plot(phi_loop, 'r')
title('Instantaneous Phase')
xlabel('Sample')
ylabel('Phase (rad)')
legend('Input', 'Output')

% Filtered Phase Error
figure
plot(phi_error)
hold on
plot(phi_error_filtered, 'g')
hold on
plot(phi_error_ss_expected * ones(nIterations,1), '-.k')
hold on
plot(phi_error_filtered_ss_expected *  ones(nIterations,1), '--r')
title('Loop Filter Output')
legend('Phase Error', ...
    'Filtered Phase Error', ...
    'Expected Phase Error Steady-State', ...
    'Expected Filtered Steady-State')
xlabel('Sample')
ylabel('Error (rad)')

% Input vs. Output Sinusoids
figure
subplot(121)
plot(real(s_in))
hold on
plot(real(s_loop), 'r')
ylabel('Amplitude')
xlabel('Sample')
title('Real part')
legend('Input', 'Output')

subplot(122)
plot(imag(s_in))
hold on
plot(imag(s_loop), 'r')
ylabel('Amplitude')
xlabel('Sample')
title('Imaginary part')
legend('Input', 'Output')
