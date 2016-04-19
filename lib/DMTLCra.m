function [En,bn] = DMTLCra(gn,Ex_bar,N,gap)
%
% EE379C 2008 Spring
%
% Levin Campello?s Method with DMT
%
% Inputs
% gn is the unitary-energy SNRs
% Ex_bar is the normalized energy
% N is the total number of real dimensions
% gap is the gap in dB
%
% Outputs
% En is the vector energy distribution (PAM or QAM) per subchannel
% bn is the vector bit distribution (PAM or QAM) per subchannel
%
% The first and last bins are PAM; the rest are QAM.
% dB into normal scale

gap = 10^(gap/10);

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
if gn(1) ~= 0
    decision_table(1)=3*gap/gn(1);
else
    decision_table(1)=inf;
end
if gn(N/2+1) ~=0
    decision_table(N/2+1)=3*gap/gn(N/2+1);
else
    decision_table(N/2+1)=inf;
end
%decision_table: debugging purpose
while(1)
    [y,index]=min(decision_table);
    E_so_far=E_so_far+y;
    if E_so_far > Ex_bar*N
        break;
    else
        En(index)=En(index)+y;
        bn(index)=bn(index)+1;
        if (index ==1 || index == N/2+1)
            decision_table(index) = 4*decision_table(index);
        else
            decision_table(index) = 2*decision_table(index);
        end
    end
end

end