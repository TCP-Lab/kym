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
% INPUT       TYPE        MEANING
% -----       ----        -------
% y        -> scalar   -> Signal to Be Transformed
% lft      -> scalar   -> Low-Frequency Threshold
% nvoice   -> scalar   -> Amount of Inter-Octave Frequencies
% wavelet  -> string   -> Mother Wavelet Type
%
% OUTPUT      TYPE         MEANING
% ------      ----         -------
% wt       -> matrix   -> Wavelet Transform of y
% fac      -> scalar   -> COI e-folding Factor
%

% Fourier Transform of the signal
Y = fft(y);

switch (wavelet)
	
	case 'M1' % Morlet Wavelet ver.1
		
		[wt,fac] = morlet1(Y,2,lft,nvoice);
	
	case 'M2' % Morlet Wavelet ver.2
		
		[wt,fac] = morlet2(Y,1.5,lft,nvoice);
	
	case 'M3' % Morlet Wavelet ver.3
		
		[wt,fac] = morlet3(Y,1.5,lft,nvoice);
		
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