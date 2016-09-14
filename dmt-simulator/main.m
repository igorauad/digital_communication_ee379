clearvars, clc
addpath(genpath('../lib'))

% Create an object that will contain all the DMT parameters
Dmt = [];

%% Debug levels
Dmt.debug.enable        = 1;  % Enable debug information
Dmt.debug.Pe            = 1;  % Debug error probabilities
Dmt.debug.loading       = 1;  % Debug bit loading
Dmt.debug.tx_energy     = 1;  % Debug transmit energy
Dmt.debug.teq           = 1;  % Debug TEQ design

%% Parameters
Dmt.alpha      = 1;       % Increase FFT size by this factor preserving Fs
% Note: this is useful to evaluate the DMT performance as N -> infty
Dmt.L          = 2;       % Oversampling (support only for integer values)
Dmt.Px         = 1e-3;    % Transmit Power (W)
Dmt.N0_over_2  = 1e-10;   % Noise PSD (W/Hz/dim) and variance per dimension
Dmt.N          = 128;     % FFT size and the number of used real dimensions
Dmt.nu         = 2;       % Cyclic Prefix Length
Dmt.gap_db     = 8.8;     % SNR gap to capacity (dB)
Dmt.delta_f    = 51.75e3; % Subchannel bandwidth
Dmt.nSymbols   = 1e3;     % Number of DMT symbols transmitted per iteration
Dmt.max_load   = inf;     % Maximum allowed bit load for each subchannel
Dmt.equalizer  = 1;       % 0 - None; 1) TEQ;
Dmt.dcNyquist  = 0;       % Flag to enable loading DC and Nyquist subchan
% MMSE-TEQ Parameters
Dmt.teqType    = 0;       % 0 - MMSE; 1 - SSNR; 2 - GeoSNR
% Monte-Carlo Parameters
Dmt.montecarlo = 1;       % 1) Enabled; 0) Disabled
Dmt.maxNumErrs = inf;
Dmt.maxIterations = 30;   % Number of iterations (nSymbols per iteration)

%% Run

[Pe_bar, Rb] = dmt(Dmt);
