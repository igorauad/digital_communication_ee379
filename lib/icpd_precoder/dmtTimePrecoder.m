function [ Precoder, Info ] = dmtTimePrecoder(p, n0, nu, tau, N, windowing)
%   Compute Time Domain Precoder and Equalizer Matrices
% ----------------------------------------------------------
%   [ Params ] = timeDomainPrecoderEqualizer(Params, Channel)
% Calculates the time domain ICI precoding matrix and time domain ISI
% equalizer matrix.
%
%   Input Parameters
% p                 Channel Pulse Response
% n0                CIR delay (index of its peak)
% nu                Cyclic Prefix Length
% tau               Cyclic Suffix length
% N                 FFT size
% windowing         When spectrum shaping windowing should be accounted
%
%   Output parameters
% This function outputs a structure containing the following information:
%
% Wt                Time domain precoding matrix
% wk                Precoder row norms
% Wisi              Time domain ISI equalizer
% sig_col           Significant columns in Wt
% sig_col_isi       Significant columns in Wisi
% sig_rows_isi      Significant rows in Wt

%% Initialization
TD_THRESHOLD = 1e-5;

% Preallocate metadata struct
Info.nColIci             = [];
Info.nRowIci             = [];
Info.nColIsi             = [];
Info.nRowIsi             = [];
Info.sizeWt              = [];
Info.sizeWisi            = [];
Info.complexity          = [];
Info.complexityNormFFT   = [];

% Normalized FFT Matrix
Q = (1/sqrt(N))*fft(eye(N));

%% Time-domain ISI, ICI and Circulant Matrices

[ Hisi, Hici, Hcirc ] = dmtIsiIciMatrices(p,...
            n0, nu, tau, N, windowing);

%% Precoder/Equalizer Design

% Time domain ICI precoder:
Wt = (Hcirc + Hici)\Hcirc;

% Time domain ISI generator:
if (sum(Hisi(:)) ~= 0) % if Hisi is not empty (when Post-cursor exists)
    Wisi = Hisi*circshift(Wt, [-tau 0]);
    % A circular rotation in the rows of Wt by -tau corresponds to a
    % rotation of the precoded transmit vector \tilde{x}^{i-1} by -tau.
    % Further details can be found in the paper.
else
    Wisi = zeros(N, N);
end

clear W_common

% Further reductions on Wt and Wisi dimensions
Wt(abs(Wt) < TD_THRESHOLD) = 0; % Round
% Use "full" Precoder matrix to compute the energy increase factors for
% each subchannel:
w_norm_n = sum(abs(Q * Wt * Q').^2,2).';
% Subtract the diagonal
Wt = Wt - eye(N);

%% Truncation for complexity reduction

% Extract only significant columns and rows in Wt and Wisi:
sig_col_ici = find(max(abs(Wt)) > TD_THRESHOLD);
Wt = Wt(:, sig_col_ici);
fprintf('Percent number of significant columns (ICI): %f\n', ...
    100*length(sig_col_ici)/N)

sig_row_ici = find(max(abs(Wt), [], 2) > TD_THRESHOLD);
Wt = Wt(sig_row_ici, :);
fprintf('Percent number of significant rows (ICI): %f\n', ...
    100*length(sig_row_ici)/N)

sig_col_isi = find(max(abs(Wisi)) > TD_THRESHOLD);
Wisi = Wisi(:, sig_col_isi);
fprintf('Percent number of significant columns (ISI): %f\n', ...
    100*length(sig_col_isi)/N)

sig_row_isi = find(max(abs(Wisi), [], 2) > TD_THRESHOLD);
Wisi = Wisi(sig_row_isi, :);
fprintf('Percent number of significant lines (ISI): %f\n', ...
    100*length(sig_row_isi)/N)

%% Summary of Results

MACs = (length(sig_row_ici)*length(sig_col_ici)) + ...
    (length(sig_col_isi)*length(sig_row_isi)) +  N*log2(N);

fprintf('Wt Dimensions: %d x %d\n', ...
    length(sig_col_ici), length(sig_row_ici));

fprintf('Wisi Dimensions: %d x %d\n', ...
    length(sig_col_isi), length(sig_row_isi));

fprintf('Complexity:     \t %d MACs\n', MACs);

fprintf('MACs per sample:\t %d \n', ...
    MACs/(N+nu));

% Copy to struct for external usage:
Precoder.ici.Wt              = Wt;
Precoder.ici.significantCols = sig_col_ici;
Precoder.ici.significantRows = sig_row_ici;
Precoder.ici.wk              = w_norm_n;

Precoder.isi.Wisi            = Wisi;
Precoder.isi.significantCols = sig_col_isi;
Precoder.isi.significantRows = sig_row_isi;

% Metadata:
Info.nColIci           = 100*length(sig_col_ici)/N;
Info.nRowIci           = 100*length(sig_row_ici)/N;
Info.nColIsi           = 100*length(sig_col_isi)/N;
Info.nRowIsi           = 100*length(sig_row_isi)/N;
Info.sizeWt            = [length(sig_col_ici) length(sig_row_ici)];
Info.sizeWisi          = [length(sig_col_isi) length(sig_row_isi)];
Info.complexity        = MACs;
Info.complexityNormFFT = MACs/(N*log2(N));

Precoder.Info = Info;

end

