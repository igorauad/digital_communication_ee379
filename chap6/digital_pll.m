clearvars clc
% Digital PLL Implementation. Main scheme is described in prof. Fred Harris
% material.
%
% Experiment with the loop filter constants Kp and Ki to verify classical
% observatinos about PLL, such that frequency offsets can not be corrected
% if Ki=0;

nIterations = 2e3;
fs          = 1e6;   % Nominal sampling frequency
y_ppm       = 1e5;   % Fractional frequency offset in ppm
fc          = 1e3;   % Nominal clock frequency
pn_var      = 1e-8; % Phase noise variance
Kp          = 0.01;  % Proportional Constant
Ki          = 0.01;  % Integral Constant

%% Derived parameters
y        = y_ppm * 1e-6;  % Fractional frequency offset in ppm
f_offset = y * fc;        % Absolute freq error (in Hz)

% Nominal Phase Increment (used in the loop phase accumulator)
delta_phi_c   = 2*pi*fc/fs;
% Input Phase Increment (eventually with frequency offset). This particular
% increment is used to simulate the input signal.
delta_phi_in   = 2*pi*(fc + f_offset)/fs;

%% Loop

% Preallocate
dds_in   = zeros(nIterations, 1);
dds_loop = zeros(nIterations, 1);
dds_mult = zeros(nIterations, 1);
phi_error          = zeros(nIterations, 1);
phi_error_filtered = zeros(nIterations, 1);
phi_in   = zeros(nIterations, 1);
phi_loop = zeros(nIterations, 1);

% Initialization values
phi_in(1) = 0;
phi_loop(1) = sqrt(pi)*randn;
integral_out_last = 0;

for i = 1:nIterations

% DDSs:
dds_in(i) = exp(1j*phi_in(i));
dds_loop(i) = exp(1j*phi_loop(i));

% Multiply the input DDS by the conjugate of the loop DDS
dds_mult(i) = dds_in(i) * conj(dds_loop(i));

% Compute the phase error by inspecting the angle of the product (atan of
% the complex number)
phi_error(i) = angle(dds_mult(i));

% Loop Filter
% Proportional term
proportional_out = phi_error(i)*Kp;
% Integral term
integral_out = phi_error(i)*Ki + integral_out_last;
integral_out_last = integral_out; 
% PI Filter output:
phi_error_filtered(i) = proportional_out + integral_out;

% Future values in the phase accumulators
phi_in(i+1) = phi_in(i) + delta_phi_in + sqrt(pn_var)*randn;
phi_loop(i+1) = phi_loop(i) + delta_phi_c + phi_error_filtered(i); 
end

%% Performance
figure
plot(phi_in)
hold on
plot(phi_loop, 'r')
title('Instantaneous Phase')
xlabel('Sample')
ylabel('Phase (rad)')
legend('Input', 'Output')

figure
plot(phi_error_filtered)
title('Loop Filter Output')
xlabel('Sample')
ylabel('Error (rad)')

figure
plot(real(dds_in))
hold on
plot(real(dds_loop), 'r')
ylabel('Amplitude')
xlabel('Sample')
title('Real part of DDS Sequences')
legend('Input', 'Output')

