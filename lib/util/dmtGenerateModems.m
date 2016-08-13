function [modulator, demodulator] = dmtGenerateModems( M, dim )
% Instantiates the mod and demod objects for each bit loading
% ----------------------------------------------------------------
%   It outputs a cell array of modems, first the two-dimensional modems,
%   then, the one-dimensional modem objects.
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

%% Generate Modems

%Preallocate modems
modulator   = cell(length(twoDim_const_orders), 1);
demodulator = cell(length(twoDim_const_orders), 1);

% Configure 2-dimensional modems for each distinct bit loading:
for i = 1:length(twoDim_const_orders)
    M = twoDim_const_orders(i);

    if (mod(log2(M),2) ~= 0)
        modulator{i} = modem.genqammod('Constellation', ...
            qamHybridConstellation(M));
        demodulator{i} = modem.genqamdemod('Constellation', ...
            qamHybridConstellation(M));
    else
        modulator{i} = modem.qammod('M', M, 'SymbolOrder', 'Gray');
        demodulator{i} = modem.qamdemod('M', M, 'SymbolOrder', 'Gray');
    end
end

% Configure 1-dimensional modems for each distinct bit loading:
for l = 1:length(oneDim_const_orders)
    i = i + 1;
    M = oneDim_const_orders(l);
    modulator{i} = modem.pammod('M', M, 'SymbolOrder', 'Gray');
    demodulator{i} = modem.pamdemod('M', M, 'SymbolOrder', 'Gray');
end

end

