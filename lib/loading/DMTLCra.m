function [En,bn] = DMTLCra(gn, Ex_bar, N, gap_db, max_load, noDcNyquist)
%
% EE379C 2008 Spring
%
% Levin Campello?s Method with DMT
%
% Inputs
% gn is the unitary-energy SNRs
% Ex_bar is the normalized energy
% N is the total number of real dimensions
% gap_db is the gap in dB
% max_load is the maximum allowed bit load for a subchannel
%
% Outputs
% En is the vector energy distribution (PAM or QAM) per subchannel
% bn is the vector bit distribution (PAM or QAM) per subchannel
%
% The first and last bins are PAM; the rest are QAM.
% dB into normal scale

if (nargin < 5)
    max_load = inf;
end
if (nargin < 6)
    noDcNyquist = 0;
end

gap = 10^(gap_db/10);

% initialization
En = zeros(1,N/2+1);
bn = zeros(1,N/2+1);
decision_table = zeros(1,N/2+1);

%debugging purpose
%plot(gn)
%%%%%%%%%%%%%%%%%%%%%%%
% Levin Campello Loading %
%%%%%%%%%%%%%%%%%%%%%%%
%initialization
%used energy so far
E_so_far=0;
%decision table - QAM and PAM
decision_table(2:N/2)=2*gap./gn(2:N/2);
% Gap formula incremental energies.
if (noDcNyquist)
    decision_table(1)=inf;
else
    decision_table(1)=3*gap/gn(1);
end
if (noDcNyquist)
    decision_table(N/2+1)=inf;
else
    decision_table(N/2+1)=3*gap/gn(N/2+1);
end
%decision_table: debugging purpose
while(1)
    [y,index] = min(decision_table);
    if (isinf(y))
        warning('Consider pouring the remaining energy over subchannels');
        error('All suchannels are maximally loaded');
    end
    E_so_far = E_so_far+y;
    if E_so_far > Ex_bar*N
        break;
    else
        En(index)=En(index)+y;
        bn(index)=bn(index)+1;
        % Prevent constellations above a certain maximum load
        if( (bn(index) == max_load) )
            decision_table(index) = inf;
        end
        % Update the decision table
        if (index ==1 || index == N/2+1)
            decision_table(index) = 4*decision_table(index);
        else
            % From 4.121, the incremental energy is 3dB for each bit.
            % However, note this only works when the gap for the two loads
            % (before and after the bit increment) are the same. Thus, care
            % must be taken for QAM, for which different gaps might be
            % applicable if QAM-Cross is used for odd number of bits and
            % QAM-SQ is used for even number of bits.
            decision_table(index) = 2*decision_table(index);
        end
    end
end

end