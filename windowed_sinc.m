function [n,coeff,f,H] = windowed_sinc(L,window_type,fc,fs)
% L has to be odd
n = 0:L-1;
sinc_signal = sinc(2*fc/fs*(n-L/2));    % sinc(x) = sin(pi*x)/pi*x
% creating window
switch(window_type)
    case 'Rectangular'
        window = rectwin(L);
    case 'Blackman'
        window = blackman(L);
    case 'Chebyshev'
        window = chebwin(L);
    case 'Kaiser'
        window = kaiser(L,1);
end
coeff = sinc_signal.*window';           % multiplying sinc by window

[H, f] = freqz(coeff, 1);               % getting response in frequency domain

end