function [dfseSNR,w_t,opt_delay]=dfsecolorsnr(l,p,nff,nbb,delay,Ex,noise);
% --------------------------------------------------------------
% [dfseSNR,w_t,opt_delay] = dfecolor(l,p,nff,nbb,delay,Ex,noise);
% l     = oversampling factor
% p     = pulse response, oversampled at l (size)
% nff   = number of feedforward taps
% nbb   = number of feedback taps 
% delay = delay of system <= nff+length of p - 2
%         if delay = -1, then choose best delay
% Ex    = average energy of signals
% noise = noise autocorrelation vector (size l*nff)
% NOTE: noise is assumed to be stationary
%
% outputs:
% dfseSNR = equalizer SNR, unbiased in dB
% w_t = equalizer coefficients [w -b]
% opt_delay = optimal delay found if delay =-1 option used.
%             otherwise, returns delay value passed to function
% created 4/96;
% ---------------------------------------------------------------

size = length(p);
nu = ceil(size/l)-1;
p = [p zeros(1,(nu+1)*l-size)];

% error checks
if nff<=0
  error('number of feedforward taps must be > 0');
end
if nbb<0
  error('number of feedback taps must be >= 0');
end
if delay > (nff+nu-1)
  error('delay must be <=  (nff+(length of p)-2)');
elseif delay < -1 
  error('delay must be >= -1');
elseif delay == -1
    delay = 0:1:nff+nu-1;
end

%form ptmp = [p_0 p_1 ... p_nu] where p_i=[p(i*l) p(i*l-1)... p((i-1)*l+1)
ptmp(1:l,1) = [p(1); zeros(l-1,1)];
for i=1:nu
  ptmp(1:l,i+1) = conj((p(i*l+1:-1:(i-1)*l+2))');
end

dfseSNR = -100;
%loop over all possible delays
for d = delay,
    %% COMPUTE EXACT number of feedback taps WE NEED to 
    %% consider given delay d
    nbb_used = min(nff+nu-1-d,nbb);
    %form matrix P, vector channel matrix
    P = zeros(nff*l+nbb_used,nff+nu);
    for i=1:nff,
        P(((i-1)*l+1):(i*l),i:(i+nu)) = ptmp;
    end
    P(nff*l+1:nff*l+nbb_used,d+2:d+1+nbb_used) = eye(nbb_used);
    %compute Rn matrix
    Rn  = zeros(nff*l+nbb_used);
    Rn(1:nff*l,1:nff*l) = toeplitz(noise);
    temp = zeros(1,nff+nu);
    temp(d+1)=1;
    %construct matrices
    Ry  = Ex*P*P' + Rn;
    Rxy = Ex*temp*P';
    new_w_t = Rxy*inv(Ry);
    sigma_dfse = Ex - real(new_w_t*Rxy');
    new_dfseSNR = 10*log10(Ex/sigma_dfse - 1);
    %save setting of this delay if best performance so far
    if new_dfseSNR >= dfseSNR
        w_t = new_w_t;
        dfseSNR = new_dfseSNR;
        opt_delay = d; 
        nbb_final = nbb_used;
    end
end

if nbb_final < nbb
    warning(strcat('For OPTIMAL FB filter settings N_bb=',num2str(nbb_final),' was used instead of the N_bb=',num2str(nbb), ' inputed'));
end