function fwt = FILT2(wt,varargin)

% 
%--------------------------------------------------------------------------------
% Filter 2 - Band Pass Filter
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% fwt = FILT2(wt,varargin)
%
% INPUT       TYPE          MEANING
% -----       ----          -------
% wt       -> matrix     -> 1st WT Output - Continuous Wavelet Transform
% varargin -> cell array -> Band Pass Arrays - [lower_voice : upper_voice]
%
% OUTPUT      TYPE          MEANING
% ------      ----          -------
% fwt      -> matrix     -> Filtered Wavelet Transform
%

mask = zeros(size(wt,1),size(wt,2));

for k = 1:length(varargin)
	
	mask(varargin{k},:) = 1;
	
end

fwt = wt .* mask;


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