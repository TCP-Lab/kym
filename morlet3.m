function [mor,fac] = morlet3(Y,sigma2,lft,nvoice)

%
%--------------------------------------------------------------------------------
% Morlet Wavelet ver.3
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% [mor,fac] = morlet3(Y,sigma2,lft,nvoice)
%
% INPUT       TYPE        MEANING
% -----       ----        -------
% Y        -> array    -> Signal FFT
% sigma2   -> scalar   -> Time/Frequency Resolution Tradeoff
% lft      -> scalar   -> Low-Frequency Threshold
% nvoice   -> scalar   -> Amount of Inter-Octave Frequencies
%
% OUTPUT      TYPE        MEANING
% ------      ----        -------
% mor      -> matrix   -> Morlet Wavelet Transform
% fac      -> scalar   -> COI e-folding Factor
%

% Frequency shift - Useful for calibration
fs = 7.4;

n = length(Y);

omega = [(0:(n/2)),(((-n/2)+1):-1)]*(2*pi/n);
omega = omega(:);

noctave = floor(log2(n))-lft;
nscale = nvoice*noctave;
mor = zeros(n,nscale);

kscale = 1;
scale = 2^(lft-1);

for jo = 1:noctave
	
	for jv = 1:nvoice
		
		a = scale*(2^(jv/nvoice));
		freq = n*(omega/a);
		
		% Time convolution as product in the transformed domain
		Psi = exp(-(sigma2/2)*(freq - fs).^2) - exp(-(freq.^2 + fs.^2)/2);
		
		% Renormalization
		Psi = Psi / sqrt(a);
		
		mor(1:n,kscale) = ifft(Y.*Psi);
		kscale = kscale+1;
	
	end
	
	scale = scale*2;

end

% Normalization
mor = ((1+exp(-fs^2)-2*exp(-(3/4)*fs^2))^(-(1/2)))*((sigma2/pi)^(1/4))*mor;

% The matrix mor is ordered from low to high frequencies 
mor = mor';

% Cone of influence e-folding factor
fac = (sqrt(sigma2)*fs)/(sqrt(2)*pi);


%%------------------------------------------------------------------------------------------------------%%
%%------------------------------------------------------------------------------------------------------%%
%%                                                                                                      %%
%% KYM Project                                                                                          %%
%% -----------                                                                                          %%
%% First Released in 2010                                                                               %%
%% Original code by Federico Alessandro Ruffinatti                                                      %%
%%                                                                                                      %%
%% UNIVERSITY OF TORINO                                                                                 %%
%% DOCTORAL SCHOOL IN LIFE AND HEALTH SCIENCES                                                          %%
%% Neurosciences Ph.D. - Experimental Neurosciences - XXV Cycle                                         %%
%% Department of Life Sciences and Systems Biology                                                      %%
%% Laboratory of Cellular Neurophysiology                                                               %%
%% Via Accademia Albertina 13 10123 Torino                                                              %%
%%                                                                                                      %%
%% Acknowledgements:                                                                                    %%
%% -----------------                                                                                    %%
%% Wavelet Transform computation is here implemented as a product in the Fourier transformed domain.    %%
%% A standard code for this algorithm can be found, for instance, in WaveLab850.                        %%
%% http://www-stat.stanford.edu/~wavelab/                                                               %%
%%                                                                                                      %%
%% Peaks detection uses a technique that is based on images dilation.                                   %%
%% See, for instance, localMaximum.m m-file by Yonathan Nativ.                                          %%
%% http://www.mathworks.com/matlabcentral/fileexchange/authors/26510/                                   %%
%%                                                                                                      %%
%%------------------------------------------------------------------------------------------------------%%
%%------------------------------------------------------------------------------------------------------%%