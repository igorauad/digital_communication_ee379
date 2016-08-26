function [dims_subch, dims_dft] = dmtSubchDimLookup(Nfft, iSubch)
% Generate look-up vectors with the number of dimensions in the DFT
%
% Inputs:
%   Nfft    -> DFT Length
%   iSubch  -> Indexes in the DFT corresponding to the subchannels
%
% Outputs:
%   dims_subch  -> Number of dimensions per subchannel
%   dims_dft    -> Number of dimensions per DFT tone

% Number of real dimensions in each tone of the DFT:
dims_dft = [1 2*ones(1, Nfft/2-1) 1 2*ones(1, Nfft/2-1)];

% Then store the number of dimensions corresponding to each used subchannel
dims_subch = dims_dft(iSubch);

end