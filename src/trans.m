function [wt,fac] = trans(y,lft,nvoice,wavelet)

%
%--------------------------------------------------------------------------------
% Wavelet Transform Computation
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% [wt,fac] = trans(y,lft,nvoice,wavelet)
%
% INPUT       TYPE         MEANING
% -----       ----         -------
% y        -> scalar    -> Signal to be Transformed
% lft      -> scalar    -> Low-Frequency Threshold
% nvoice   -> scalar    -> Amount of Inter-Octave Frequencies
% wavelet  -> string    -> Mother Wavelet Type
%
% OUTPUT      TYPE         MEANING
% ------      ----         -------
% wt       -> matrix    -> Wavelet Transform of y
% fac      -> structure -> Wavelet Characteristic Factors
%

% Variables Assignment
n = length(y);

% Fast Fourier Transform of the signal
Y = fft(y); % dt factor is not needed because of the following ifft(Psi.*Y);

switch (wavelet)
	
	case 'M1' % Morlet Wavelet ver.1
		
		[wt,fac] = morlet1(Y,1,lft,nvoice);
	
	case 'M2' % Morlet Wavelet ver.2
		
		[wt,fac] = morlet2(Y,1,lft,nvoice);
	
	case {'R','D2'} % Ricker Wavelet
		
		[wt,fac] = ricker(Y,lft,nvoice);
	
	case 'D1' % DOG-1 Wavelet
		
		[wt,fac] = dog1(Y,lft,nvoice);
	
end


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