function fwt = FILT1(wt,thr)

% 
%--------------------------------------------------------------------------------
% Filter 1 - Amplitude Filter
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% fwt = FILT1(wt,thr)
%
% INPUT       TYPE         MEANING
% -----       ----         -------
% wt       -> matrix    -> 1st WT Output - Continuous Wavelet Transform
% thr      -> scalar    -> Threshold Percentage
%
% OUTPUT      TYPE         MEANING
% ------      ----         -------
% fwt      -> matrix    -> Filtered Wavelet Transform
%

% Variables assignment
x = abs(wt);

% Threshold
x = (x-min(min(x)))/max(max(x-min(min(x))));
threshold = thr/100;

fwt = wt .* (x >= threshold);


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