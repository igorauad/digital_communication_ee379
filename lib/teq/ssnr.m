function [ ssnr ] = ssnr( w, p, delta, nu )
% [ ssnr ] = ssnr( w, p, delta, nu )
% -------------------------------------------------------------------------
%   Compute the "Shortening SNR" (SSNR), which measures the ratio between
%   the energy within the window of samples where MMSE should concentrate
%   the energy and outside the window.

% Effective pulse response
p_eff = conv(p,w);
p_eff_indexes = 1:length(p_eff);
% nu + 1 consecutive samples where energy should be concentrated

% Window indexes
win_indexes = delta+1:delta+nu+1;
p_win = p_eff(win_indexes);

% Remaining samples
wall_indexes = setdiff(p_eff_indexes, win_indexes);
p_wall = p_eff(wall_indexes);

% SSNR:
ssnr = (p_win * p_win') / (p_wall * p_wall');

end

