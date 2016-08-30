clearvars, clc
addpath(genpath('../lib'))

% Create an object that will contain all the DMT parameters
Dmt = [];

%% Debug levels
Dmt.debug.enable        = 1;  % Enable debug information
Dmt.debug.constellation = 0;  % Debug a certain subchannel constellation
Dmt.debug.tone          = 16; % Tone whose constellation is debugged
Dmt.debug.Pe            = 1;  % Debug error probabilities
Dmt.debug.loading       = 0;  % Debug bit loading
Dmt.debug.tx_energy     = 0;  % Debug transmit energy
Dmt.debug.teq           = 1;  % Debug TEQ design

%% Parameters
Dmt.alpha      = 1;       % Increase FFT size by this factor preserving Fs
% Note: this is useful to evaluate the DMT performance as N -> infty
Dmt.L          = 2;       % Oversampling (support only for integer values)
Dmt.Px         = 1e-3;    % Transmit Power (W)
Dmt.N0_over_2  = 1e-15;   % Noise PSD (W/Hz/dim) and variance per dimension
Dmt.N          = 2048;    % FFT size and the number of used real dimensions
Dmt.nu         = 64;      % Cyclic Prefix Length
Dmt.tau        = 64;      % Cyclic Suffix
Dmt.windowing  = 1;       % Activate Lcs windowing + Overlap
Dmt.gap_db     = 8.8;     % SNR gap to capacity (dB)
Dmt.delta_f    = 51.75e3; % Subchannel bandwidth
Dmt.nSymbols   = 1e3;     % Number of DMT symbols transmitted per iteration
Dmt.max_load   = inf;     % Maximum allowed bit load for each subchannel
Dmt.equalizer  = 1;       % 0 - None; 1) TEQ; 2) Cheong; 3) Time Domain
Dmt.dcNyquist  = 0;       % Flag to enable loading DC and Nyquist subchan
% MMSE-TEQ Parameters
Dmt.teqType    = 0;       % 0 - MMSE; 1 - SSNR; 2 - GeoSNR
% Monte-Carlo Parameters
Dmt.montecarlo = 1;       % 1) Enabled; 0) Disabled
Dmt.maxNumErrs = inf;
Dmt.maxNumDmtSym = 1e12;
% Channel
channelChoice = 1;

%% Run
[Pe_bar, Rb] = dmt(Dmt, channelChoice);