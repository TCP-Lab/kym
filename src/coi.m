function coimask = coi(wtabs,nvoice,coifactor,cone,test2)

%
%--------------------------------------------------------------------------------
% Cone of Influence - COI
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% coimask = coi(wtabs,nvoice,coifactor,cone,test2)
%
% INPUT        TYPE        MEANING
% -----        ----        -------
% wtabs     -> matrix   -> CWT Modulus
% nvoice    -> scalar   -> Amount of Inter-Octave Frequencies
% coifactor -> scalar   -> COI e-folding Factor
% cone      -> string   -> Cone of Influence Handling Method
% test2     -> boolean  -> 1 if filename=='test2', 0 otherwise
%
% OUTPUT       TYPE        MEANING
% ------       ----        -------
% coimask   -> matrix   -> COI Mask
%

% Variables Assignment
[nscale,n] = size(wtabs);

switch (cone)
	
	case 'COI' % Cone of influence mask
		
		coimask = ones(nscale,n);
		
		freqVecComp = (1/2).*2.^(-([0:1:nscale-1])./nvoice); % Short-Circuit(1): f.Nyquist=1/(2*dt)
		freqVecComp = flipud(freqVecComp');
		
		COI = coifactor./freqVecComp; % Short-Circuit(2): COI=COI/dt;
		COI(find(COI > n)) = n;
		
		% You can use a more severe COI to further reduce edge effects on WAI indexes
		%COI = COI * 1.5;
		%COI(find(COI > n)) = n;
		
		for u = 1:nscale
			coimask(u,1:round(COI(u))) = 0;
			coimask(u,n-round(COI(u))+1:n) = 0;
		end
		
		% e-folding Test
		if test2
			r1 = wtabs(1:nscale,1)';
			for h = 1:nscale
				r2(h) = wtabs(h,min(find(coimask(h,:))));
			end
			fold = flipud((r1./r2)')
			mean(fold)
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