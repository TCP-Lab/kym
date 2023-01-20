function [d1,fac] = dog1(Y,lft,nvoice)

%
%--------------------------------------------------------------------------------
% DOG-1 Wavelet (First Derivative of Gaussian)
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% [mor,fac] = dog1(Y,lft,nvoice)
%
% INPUT       TYPE         MEANING
% -----       ----         -------
% Y        -> array     -> FFT Signal
% lft      -> scalar    -> Low-Frequency Threshold
% nvoice   -> scalar    -> Amount of Inter-Octave Frequencies
%
% OUTPUT      TYPE         MEANING
% ------      ----         -------
% d1       -> matrix    -> CWT with DOG-1 Mother Wavelet
% fac      -> structure -> Wavelet Characteristic Factors
% fac.coi  -> scalar    -> COI e-folding Factor
% fac.adm  -> scalar    -> Admissibility Condition Factor
%

% Variables Assignment

% Frequency shift - Useful for calibration
% This 's' value makes 1/a==ordinary_frequency
s = 1/(2*pi);

n = length(Y);

noctave = floor(log2(n))-lft; % Number of octaves (Dyadic algorithm)
nscale = nvoice*noctave; % nvoice = Number of voices per octave
d1 = zeros(n,nscale); % CWT transform empty matrix

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
		asxi = a*s*2*pi*freq;
		
		% Fourier transform of (real) DOG-1 wavelet
		Psi = asxi.*exp(-0.5*asxi.^2);
		% Complex analytic wavelet
		Psi = Psi.*(freq >= 0);
		
		% Energetic normalization (Unitary Transform)
		%Psi = Psi*sqrt(a*USdt);
		
		% Time convolution as product in the transformed domain
		d1(1:n,kscale) = ifft(Psi.*Y); % 1/dt factor is not needed because of the previous Y=fft(y);
		kscale = kscale+1;
	
	end
	
	scale = scale*2; % Dyadic algorithm: 'scale' ranges from 2 to 2^(noctave) 

end

% Normalization
d1 = (2*sqrt(2*s)*pi^(1/4))*d1;

% The matrix d1 is ordered from high to low frequencies
d1 = d1';
d1 = flipud(d1);

% Cone of influence e-folding factor
fac.coi = 0.478147; %%%%%%%%%%%%%% NOOOOOOOOOOOOOOOOOOO DA CALCOLARE!!!!
% Admissibility condition factor (1/Cpsi)
fac.adm	= 1/(2/sqrt(pi));


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