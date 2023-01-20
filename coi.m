function coimask = coi(nvoice,nscale,n,dt,lowest,fac,cone)

%
%--------------------------------------------------------------------------------
% Cone of Influence Handling - COI
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% coimask = coi(nvoice,nscale,n,dt,lowest,fac,cone)
%
% INPUT       TYPE        MEANING
% -----       ----        -------
% nvoice   -> scalar   -> Amount of Inter-Octave Frequencies
% nscale   -> scalar   -> size(wt,1)
% n        -> scalar   -> size(wt,2)
% dt       -> scalar   -> Sampling Time Mode
% lowest   -> scalar   -> Frequency of the Lowest Octave
% fac      -> scalar   -> COI e-folding Factor
% cone     -> string   -> Cone of Influence Handling Method
%
% OUTPUT      TYPE        MEANING
% ------      ----        -------
% coimask  -> matrix   -> COI Mask
%

switch (cone)
	
	case 'COI' % Cone of influence mask
		
		coimask = ones(nscale,n);
		
		COI = fac./(lowest.*(2.^([1:1:nscale]./nvoice)));
		COI = (COI.*1000)./dt;
		COI(find(COI > n)) = n;
		
		% You can use a more severe COI to further reduce edge effects on WAI indexes (especially J)
		%COI = COI * 1.5;
		%COI(find(COI > n)) = n;
		
		for u = 1:nscale
			coimask(u,1:floor(COI(u))) = 0;
			coimask(u,n-floor(COI(u))+1:n) = 0;
		end
		
	case {'PAD','NONE'} % Do nothing
		
		coimask = ones(nscale,n);
		
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