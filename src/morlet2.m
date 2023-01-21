function [mor,fac] = morlet2(Y,sigma2,lft,nvoice)

%
%--------------------------------------------------------------------------------
% Morlet Wavelet ver.2
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% [mor,fac] = morlet2(Y,sigma2,lft,nvoice)
%
% INPUT       TYPE         MEANING
% -----       ----         -------
% Y        -> array     -> FFT Signal
% sigma2   -> scalar    -> Time/Frequency Resolution Tradeoff
% lft      -> scalar    -> Low-Frequency Threshold
% nvoice   -> scalar    -> Amount of Inter-Octave Frequencies
%
% OUTPUT      TYPE         MEANING
% ------      ----         -------
% mor      -> matrix    -> CWT with Morlet Mother Wavelet
% fac      -> structure -> Wavelet Characteristic Factors
% fac.coi  -> scalar    -> COI e-folding Factor
% fac.adm  -> scalar    -> Admissibility Condition Factor
%

% Variables Assignment

% Frequency shift - Useful for calibration
% This 's' value makes 1/a==ordinary_frequency
s = 2*pi;

n = length(Y);

noctave = floor(log2(n))-lft; % Number of octaves (Dyadic algorithm)
nscale = nvoice*noctave; % nvoice = Number of voices per octave
mor = zeros(n,nscale); % CWT transform empty matrix

stp = 1/n; % Frequency step

% Frequency domain [-1/2,+1/2-1/n], written as [0,+1/2-1/n]U[-1/2,-1/n] in accordance with FFT convention
% Simplify(1): freq == [-1/2,+1/2-1/n]*(1/dt) with dt==sampling interval -> 1/(2*dt)==f.Nyquist
freq = [0:stp:1/2 -1/2+stp:stp:-stp]';

% No need for negative scales since signal is real and wavelet is Hermitian
kscale = 1;
scale = 2;

for jo = 1:noctave % Octave index
	
	for jv = 0:nvoice-1 % Voice index
		
		% Dyadic algorithm: 'scale' ranges from 2 (finer scale corresponds to the highest frequency==f.Nyquist) to ~2^(noctave+1)
		% Simplify(2): a = scale*(2^(jv/nvoice))*dt == [1/Nyquist frequency,Signal_Lenght]
		a = scale*(2^(jv/nvoice)); 
		
		% Fourier transform of the Morlet wavelet (Exact form)
		Psi = exp(-(1/2)*(2*pi*freq*a - s).^2) - exp(-(1/2)*((2*pi*freq*a).^2 + s^2));
		
		% Energetic normalization (Unitary Transform)
		%Psi = Psi*sqrt(a*USdt);
		
		% Time convolution as product in the transformed domain
		mor(1:n,kscale) = ifft(Psi.*Y); % 1/dt factor is not needed because of the previous Y=fft(y);
		kscale = kscale+1;
	
	end
	
	scale = scale*2; % Dyadic algorithm: 'scale' ranges from 2 to 2^(noctave) 

end

% Normalization
mor = ((1+exp(-s^2)-2*exp(-(3/4)*s^2))^(-1/2))*((4*pi*sigma2)^(1/4))*mor;

% The matrix mor is ordered from high to low frequencies
mor = mor';
mor = flipud(mor);

% Cone of influence e-folding factor
fac.coi = (sqrt(sigma2)*s)/(sqrt(2)*pi);
% Admissibility condition factor (1/Cpsi)
fac.adm	= 1/1.013179900657;


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
%% Neurosciences PhD - Experimental Neurosciences - XXV Cycle                                           %%
%% Department of Life Sciences and Systems Biology                                                      %%
%% Laboratory of Cellular Neurophysiology                                                               %%
%% Via Accademia Albertina 13 10123 Torino                                                              %%
%%                                                                                                      %%
%% Acknowledgements:                                                                                    %%
%% -----------------                                                                                    %%
%% CWT convolution is implemented as a product in the Fourier transformed domain.                       %%
%% In particular, the code for CWT computation is a refinement of WaveLab850 dyadic algorithm.          %%
%% http://www-stat.stanford.edu/~wavelab/                                                               %%
%%                                                                                                      %%
%% A technique based on image dilation has been used for the detection of peaks and maxima.             %%
%% This idea comes from Yonathan Nativ's localMaximum.m m-file.                                         %%
%% http://www.mathworks.com/matlabcentral/fileexchange/authors/26510/                                   %%
%%                                                                                                      %%
%%------------------------------------------------------------------------------------------------------%%
%%------------------------------------------------------------------------------------------------------%%