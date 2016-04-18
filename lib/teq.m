function [P,w,b,SNRmfb,alpha] = teq(p,L,nu,delta,sigma,Ex_bar,filtertype)
% Programmed by Ryoulhee Kwak
%
% filtertype 1 = FIR , 2 = IIR(just for one pole filter)
% p = channel impulse response
% in case of FIR,  p = 1+0.5D+2.3D^2 -> p = [1 0.5 2.3]
% in case of IIR,  p = 1/(2.4-0.5D)  -> p = [2.4 -0.5]
% sigma => noise power in linear scale
% delta => delay
% L  => number of taps
% c  => w*p
% v  =>eigen  vectors
% d  =>eigent values
% nu => the dispersion of the desirable pulse response (b)

size_p = length(p);

if filtertype==1 % FIR
    
    norm_p = norm(p);
    
    % Channel matrix:
    P = toeplitz([p(1), zeros(1,L-1)]',[p,zeros(1,L-1)]);
    
    % Channel output autocorrelation for white noise:
    ryy = Ex_bar*P*P' + sigma*eye(L,L);
    % See (4.342)
    
    % Cross-correlation matrix:    
    rxy = [zeros(nu+1,delta) eye(nu+1) zeros(nu+1,L+size_p-2-nu-delta)]...
        * P' * Ex_bar;
    % See comments in (4.343). The above is rxy(delta), namely a function
    % of delta, the delay.
    
    % MMSE-LE error auto-correlation:
    rle = eye(nu+1)*Ex_bar - (rxy * inv(ryy) * rxy');
    % See (4.327)
    
    % The MMSE TEQ problem can be stated as:
    %   - Find a vector of coefficients "b" (target pulse response) and the
    %   optimal delay "delta" such that the minimum mean square error is
    %   obtained.
    
    % Solution is shown to be the eigenvector of R_LE that corresponds to
    % the minimum eigenvalue. Thus, find first eigenvectors and eigen
    % values:
    [v d] = eig(rle);
    % Find the minimum eigenvalue:
    d_temp = diag(d)';
    [lambda_min,n] = min(d_temp);
    % Find b as the scaled version (by ||p||) of the chosen eigenvector:
    b = norm_p * v(:,n)';
    
    % Then, obtain the feed-forward equalizer from (4.325):
    w = b * rxy * inv(ryy);
   
    % Find the MMSE:
    mmse=b*rle*b';
    
    % "The bias factor itself is obtained by first forming the
    % convolution":
    c = conv(w,p);
    alpha = c(delta+1)./b(1); % Bias. See (4.334)
    
    % Biased error energy, from (4.338):
    biased_error_energy = norm_p^2*(lambda_min - Ex_bar*(1-alpha)^2);
    
    % Bias removed by multiplication of the equalizer output by 1/alpha:
    unbiased_error_energy = biased_error_energy / alpha^2;
    
    % Finally, the SNRmfb can be computed, according to the definition of
    % (4.300):
    SNRmfb = norm_p^2*Ex_bar / unbiased_error_energy;
    
else %IIR
    %in case of IIR filter P matrix is infinite dimesion but there is a trick to get exact rxy, ryy
    norm_p = sqrt(1/(p(1)^2*(1-(p(2)/p(1))^2 )));
    
    ptemp=[(-p(2)/p(1)).^(0:1:L-1)]/p(1); %-1 factor!
    P=toeplitz([ptemp(1) zeros(1,L-1)]',ptemp);
    Ptemp=toeplitz(ptemp',ptemp);
    
    ryy=Ex_bar*Ptemp*norm_p^2+sigma*eye(L,L);
    rxy=[zeros(nu+1,delta) eye(nu+1) zeros(nu+1,L-1-nu-delta)]*P'*Ex_bar;
    rle=eye(nu+1)*Ex_bar-rxy*inv(ryy)*rxy';
    [v d]=eig(rle);
    d_temp=diag(d)';
    [lambda_min,n ]= min(d_temp);
    b=norm_p*v(:,n)';
    w=b*rxy*inv(ryy);
    c=conv(w,p);
    sum(conv(w,p).^2)-norm(b)^2*(w(1)/b(1))^2;
    mmse=b*rle*b';
    %by using easy way to get error energy
    alpha=c(delta+1)./b(1)
    biased_error_energy=norm_p^2*(lambda_min-Ex_bar*(1-alpha)^2)
    unbiased_error_energy=norm_p^2*(lambda_min-Ex_bar*(1-alpha)^2)/alpha^2
    
    SNRmfb=norm_p^2*Ex_bar/unbiased_error_energy;
    SNRmfb_in_dB=10*log10(SNRmfb);
end
