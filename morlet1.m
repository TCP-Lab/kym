function mor = morlet1(n,Y,sigma2,lft,nvoice)

%
%---------------------------------------------------------------------------
% Morlet Wavelet ver.1
%---------------------------------------------------------------------------
%
%
% Function Definition
%
% mor = morlet1(n,Y,sigma2,lft,nvoice)
%
% INPUT       TYPE        MEANING
% -----       ----        -------
% n        -> scalar   -> Signal Length
% Y        -> array    -> Signal FFT
% sigma2   -> scalar   -> Time/Frequency Resolution Tradeoff
% lft      -> scalar   -> Low-Frequency Threshold
% nvoice   -> scalar   -> Lines of Pixel per Octave
%
% OUTPUT      TYPE        MEANING
% ------      ----        -------
% mor      -> matrix   -> Morlet Wavelet Transform
%

% Frequency Shift - Useful for Calibration
fs = 6.4;

omega = [(0:(n/2)),(((-n/2)+1):-1)]*(2*pi/n);
omega = omega(:);

noctave = floor(log2(n))-lft;
nscale = nvoice*noctave;
mor = zeros(n,nscale);

kscale = 1;
scale = 2;

for jo = 1:noctave
	
	for jv = 1:nvoice
		
		a = scale*(2^(jv/nvoice));
		freq = (fs - a*omega);
		% Time convolution as product in the transformed domain
		Psi = exp(-(sigma2/2)*freq.**2);
		Psi = Psi.*(omega > 0);
		mor(1:n,kscale) = ifft(Y.*Psi);
		kscale = kscale+1;
	
	endfor
	
	scale = scale*2;

endfor

% Normalization
mor = ((sigma2/pi)**(1/4))*mor;

% The matrix mor is ordered from high to low frequencies
mor = mor';
mor = flipud(mor);


%---------------------------------------------------------------------%
%                                                                     %
% A.A. 2009 / 2010                                                    %
% Original code by Federico Alessandro Ruffinatti                     %
% Universit� degli Studi di Torino - Italy - DBAU - Scienze MFN       %
% Scuola di Dottorato in Neuroscienze - XXV ciclo                     %
%                                                                     %
% Wavelet computation is regarded as a time convolution and it is     %
% implemented as a product in the Fourier transformed domain.         %
% A standard code for this algorithm can be found, for instance,      %
% in WaveLab850 - http://www-stat.stanford.edu/~wavelab/              %
%                                                                     %
% Peaks detection uses a technique that is based on images dilation.  %
% See, for instance, localMaximum.m m-file by Yonathan Nativ          %
% http://www.mathworks.com/matlabcentral/fileexchange/authors/26510/  %
%                                                                     %
%---------------------------------------------------------------------%