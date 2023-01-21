function dataout = FILTsig(datain,denoisetype,L,fc,graph)

% 
%--------------------------------------------------------------------------------
% Digital FIR/IIR Filters for 1D signals
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% dataout = FILTsig(datain,denoisetype,L,fc,graph)
%
% INPUT           TYPE        MEANING
% -----           ----        -------
% datain       -> matrix   -> Matrix of the Input Signals to Be Filtered
% denoisetype  -> string   -> Name of the Chosen Filter
% L            -> scalar   -> Filter Order
% fc           -> scalar   -> Cutoff frequency (-3dB) in Normalized Units
% graph        -> boolean  -> Enable Filter Visualization Tool
%
% OUTPUT          TYPE         MEANING
% ------          ----         -------
% dataout      -> matrix   -> Matrix of the Filtered Output Signals
% -none-       -> plot     -> Filter Visualization Tool
%

switch (denoisetype)
	
	case 'none' % NO FILTER
	
		dataout = datain;
	
	case 'mean' % Mean Convolutive Filter
	
		% Window Width = Filter Order + 1 = L+1
		% Symmetric FIR -> Linear Phase
		% Lowpass with Stopband Ripples (but no overshoot or ringing)
		
		b = ones(1,L+1)/(L+1);
		a = 1;
		%dataout = filter(b,a,datain);
		dataout = conv(datain,b,'same'); % Alternative - to avoid the Delay introduced by a FIR Filter
	
	case 'gauss' % Gaussian Filter
	
		% Window Width = Filter Order + 1 = L+1
		% Symmetric FIR -> Linear Phase 
		% Lowpass with fewer Stopband Ripples (and no overshoot or ringing)
		
		b = gausswin(L+1)/sum(gausswin(L+1));
		a = 1;
		%dataout = filter(b,a,datain);
		dataout = conv(datain,b,'same'); % Alternative - to avoid the Delay introduced by a FIR Filter
	
	case 'wmean' % Minimal (Weighted-Mean) Convolutive Filter
		
		% Window Width = 2 -> Filter Order = 1
		% Asymmetric FIR -> Non-Linear Phase
		% Lowpass with no Ripples (and no overshoot or ringing)
		
		w = 0.1; % Weight 0 <= w <= 1
		
		b = [w 1]/(1+w);
		a = 1;
		%dataout = filter(b,a,datain);
		dataout = conv(datain,b,'same'); % Alternative - to avoid the Delay introduced by a FIR Filter
	
	case 'butterworth' % Lowpass Digital IIR Butterworth Filter (maximally flat magnitude filter)
		
		% Window Width = 2 -> Filter Order = 1
		% IIR Filter -> Non-Linear Phase
		% Lowpass with no Ripples (maximally flat magnitude filter)
		
		[b,a] = butter(L,fc);
		%dataout = filter(b,a,datain,[],1);
		dataout = filtfilt(b,a,datain); % It compensates for the delay introduced by the IIR filter and corrects for phase distortion
		% Using 'filtfilt' the actual filter transfer function equals the squared magnitude of the original filter transfer function
		% and the filter order is double the order of the original filter specified by b and a
	
	case 'FIRflat' % Lowpass Maximally Flat Digital FIR (Generalized Butterworth Filter)
		
		% Window Width = Filter Order + 1 = L+1
		% Symmetric FIR -> Linear Phase 
		% Lowpass with no Ripples (and little overshoot)
		
		b = maxflat(L,'sym',fc);
		a = 1;
		
		%dataout = filter(b,a,datain);
		% FIR filters with linear phase have no phase distortion, but just a pure delay that should equals half the filter order
		% You can check it through:
		%[Gd,W] = grpdelay(b,a);
		%delay = mean(Gd);
		% Compensates for the Delay Introduced by the FIR Filter
		%delay = L/2;
		%ntrace = size(datain,2);
		%dataout(end+1:end+delay,:) = zeros(delay,ntrace);
		%dataout(1:delay,:) = [];
		
		dataout = conv(datain,b,'same'); % Alternative - to avoid the Delay introduced by a FIR Filter
	
end

if (~strcmp(denoisetype,'none') && graph)
	fvtool(b,a)
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