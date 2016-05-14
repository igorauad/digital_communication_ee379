function [h] = truncateCir(ht)
% Truncate the Channel Impulse Response
% ---------------------------------------
% function [h] = truncateCir(ht)
%   This function returns a truncated version of ht, which contains
%   99.9% of the energy.
%
% Input Parameters:
%   ht    => the original impulse response of channel
%
% Output Parameters:
%    h    => Truncated impulse response of channel with 99.9% of original
% energy

target = 0.999;

% Total Impulse Response energy:
energyHt = sum(ht.^2);

accEnergy = 0;

% Index in which the energy is reached:
iReached = 1;

% Cumulative energy at each index of the CIR:
for l = 1 : length(ht)
    accEnergy = accEnergy + ht(l).^2;
    if(accEnergy >= (target * energyHt))
        iReached = l;
        break;
    end
end

% Truncated Response:
h = ht(1 : iReached);

% THRESHOLD = 4e-4;
% STOP_NUMBER = 150;
% L = length(ht) - 1;
%
% diff_h = abs(diff(ht));
%
% stall_counter = 0;
%
% [~, iPeak] = max(ht);
%
% for i = 1:L
%     if(diff_h(i) > THRESHOLD)
%         stall_counter = 0;
%     else
%         stall_counter = stall_counter + 1;
%     end
%
%     % Check how long the derivative is stall and if the moving cursor is
%     % already beyond the CIR peak.
%     if(stall_counter > STOP_NUMBER && i > iPeak)
%        iReached = i;
%        break;
%     end
% end

% % Truncated Response:
% h = ht(1 : iReached);

return