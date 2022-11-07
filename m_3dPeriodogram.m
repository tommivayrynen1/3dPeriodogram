function [p,f] = m_3dPeriodogram(X,taper,nfft,srate)

% Implement 3dPeriodogram function from AFNI

% X = row vector or 2D data array channel x time
% taper = fraction of data to taper 
% nfft = FFT length in samples
% sampling rate

    X=detrend(X',1)';           % remove linear trend
    npts = size(X,2);           % number of time points analyzed (<= nfft) (i.e., the length of the input dataset)
    ntaper = taper * npts / 2;  % = number of points to taper on each end
    ktop = npts - ntaper;
    phi = pi / ntaper;
    % k = 0 ... nfft-1
    %kp=[0:1:nfft-1]; -> needs to be : 0 .. npts-1
    kp = [0:1:npts-1];

    syms y(k) % Define fractional hamming window function.
    y(k) = piecewise(0 <= k < ntaper,0.54 - 0.46 * cos(k*phi),...
           ktop <= k < npts,0.54 + 0.46 * cos((k-ktop+1)*phi),...
           1);
       
    w = double(y(kp));
    
    [pxx,f] = periodogram(X',w,nfft,srate); % units in V^2 / Hz
    p = pxx*max(f);
    
%   p is in same units as in 3dPeriodogram
%   "The result is the squared magnitude of the FFT of w(k)*data(k),
%   divided by P.  This division makes the result be the 'power',
%   which is to say the data's sum-of-squares ('energy') per unit
%   time (in units of 1/TR, not 1/sec) ascribed to each FFT bin."

end
