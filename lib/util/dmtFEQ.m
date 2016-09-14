function [ FEQn ] = dmtFEQ(p, dmt )
% Design the DMT Frequency Equalizer

% Parameters
Nfft      = dmt.Nfft;
n0        = dmt.n0;

% Frequency domain response (use the effective pulse response in the FEQ)
H = fft(p, Nfft);

% Store only the response at the used indices of the FFT
Hn = H(dmt.iTonesTwoSided);

% Corresponding phase shift due to cursor
phaseShift = exp(1j*2*pi*(n0/Nfft)*(dmt.iTonesTwoSided.' - 1));

FEQn    = (1 ./ (Hn .* phaseShift));

end

