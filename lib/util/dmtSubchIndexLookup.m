function [iTones, iTonesFull] = dmtSubchIndexLookup(N, Nfft, L, DcNyquist)
% Generate a look-up table with the DFT indexes of the subchannels
% [iTones, iTonesFull] = dmtSubchIndexLookup(N, Nfft, L, DcNyquist)
%
% Inputs:
%  N          -> Number of Subchannels
%  Nfft       -> DFT Size
%  L          -> Oversampling factor
%  DcNyquist  -> Whether DC and Nyquist 1-dim tones shall be used
%
% Outputs:
%
% iTonesFull -> Contains the indexes that can be used as subchannels from
%               the full Hermitian DFT. Bit loading will determine whether
%               used or not.
% iTones     -> Contains the indexes that can be used as subchannels from
%               the positive half of the DFT. Again, bitloading determines
%               whether used or not.

%% Main
% Vector of used subchannels among the DFT indices
% Assume DC is at tone 1 and Nyquist at Nfft/2 + 1. To compute the
% Hermitian symmetric of the k-th tone, we compute "Nfft + 2 - k".
if (L == 1)
    if (~DcNyquist)
        iTonesFull = [2:(N/2), (Nfft +2 - N/2):Nfft].';
        iTones      = 2:(N/2);
    else
        iTonesFull = [1:(N/2 +1), (Nfft +2 - N/2):Nfft].';
        iTones      = 1:(N/2 +1);
    end
else
    % When oversampling is used, it is easier to force the disabling of DC
    % and Nyquist tones. In this case, the bin that for L=1 (no
    % oversampling) represents Nyquist is loaded as a complex subchannel.
    warning('DC and Nyquist are tones are not loaded due to oversampling');
    iTonesFull = [2:(N/2 + 1), (Nfft +2 - N/2 - 1):Nfft].';
    % N/2 complex subcarriers at the positive half and the corresponding
    % N/2 complex conjugate subcarriers at the negative half.
    iTones      = 2:(N/2 + 1);
end

end