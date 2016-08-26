function [tx_data] = dmtRndData(dmt)
% Random DMT Data Generation
% [tx_data] = dmtRndData(dmt)
%
% Inputs:
%  dmt     -> Struct with DMT parameters
%
% Outputs:
%  tx_data -> Bits allocated per subchannel

% Preallocate
tx_data = zeros(dmt.N_subch, dmt.nSymbols);

% Iterate over the distinct modulators
for iModem = 1:length(dmt.modulator)
    % Modulation order
    M = dmt.modulator{iModem}.M;

    % Loaded subchannels
    iSubChs = (dmt.modem_n == iModem);

    % Generate random data
    tx_data(iSubChs, :) = randi(M, sum(iSubChs), dmt.nSymbols) - 1;
end

end