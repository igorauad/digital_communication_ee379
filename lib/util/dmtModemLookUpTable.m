function [ modem_n ] = dmtModemLookUpTable( M, dim )
% Creates a look-up table indicating the modem for each subchannel
% ----------------------------------------------------------------
%   It assumes the cell array of modems stores first the two-dimensional
%   modems, then, in its last elements, the one-dimensional modem objects.
%
%   Inputs:
%  M   -> Vector with the constellation order for each subchannel
%  dim -> Vector with the number of dimensions associated to each
%  subchannel

%% Infer constellation orders
oneDimModOrders     = M(dim == 1);
twoDimModOrders     = M(dim == 2);
oneDim_const_orders = unique(oneDimModOrders(oneDimModOrders~=1));
twoDim_const_orders = unique(twoDimModOrders(twoDimModOrders~=1));

%% Look-up generation:

% Preallocate
modem_n = zeros(length(M), 1);

% Iterate over subchannels
for k = 1:length(M)
    % For two-dimensional subchannels
    if (dim(k) == 2)
        iModem = find(twoDim_const_orders == M(k));
        if (iModem)
            modem_n(k) = iModem;
        end
    else % For one-dimensional subchannels
        iModem = find(oneDim_const_orders == M(k));
        if (iModem)
            modem_n(k) = length(twoDim_const_orders) + iModem;
        end
    end
end

end

