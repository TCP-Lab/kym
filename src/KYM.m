function ratio = KYM(wt,par,sig,thr,output)

% 
%--------------------------------------------------------------------------------
% CWT Analysis and Peaks Detection
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% ratio = KYM(wt,par,sig,thr,output)
%
% INPUT       TYPE         MEANING
% -----       ----         -------
% wt       -> matrix    -> 1st WT Output - Continuous Wavelet Transform
% par      -> structure -> 2nd WT Output - Parameters
% sig      -> array     -> 3rd WT Output - Calcium Signal
% thr      -> scalar    -> Scalogram Threshold Percentage
% output   -> boolean   -> Print .eps Output Graphs
%
% OUTPUT      TYPE         MEANING
% ------      ----         -------
% ratio    -> array     -> post/pre Spectrum Ratio Vector
% -none-   -> plot      -> 8 Plots Resulting from Analysis
%

%---------------------------------------------------------------------------
% Graphic Parameters
%---------------------------------------------------------------------------

s1 = 16; % X-Y TickLabel size
s2 = 19; % X-Y Label and text size
s3 = 24; % Title size

%---------------------------------------------------------------------------
% Default Output
%---------------------------------------------------------------------------

ratio = [];

%---------------------------------------------------------------------------
% Variable Assignment
%---------------------------------------------------------------------------

coimask = par.k;
timeVec = par.x;
freqVec = par.y;
periodVec = par.z;
unit = par.u;
mark = par.m;
roi = par.r;
admis = par.c;
US = par.US;

x = abs(wt);
y = abs(real(wt));

%---------------------------------------------------------------------------
% Input Control
%---------------------------------------------------------------------------

% Print plots - Default value = false = Print nothing
if (nargin < 5)
	output = false;
end

% Thresholds - Default value = 0
if (nargin < 4)
	thr = 0;
end

%--------------------------------------------------------------------------------
% CWT Analysis: selectively enable/disable each one of the following features
%--------------------------------------------------------------------------------

ENERGY		= false;
FOURIER		= true;
PEAK		= false;
PATHS		= false;
SKELETON	= false;
HURST		= false;
HOLDER		= false;
WAI			= false;
VECTOR		= false;

%--------------------------------------------------------------------------------
% Wavelet Power Spectrum
%--------------------------------------------------------------------------------

if ENERGY
	figure
	
	energy(x,coimask,unit,US,timeVec,freqVec,mark,admis)
	
	% Print output graph
	if (output)
		print -depsc energyout.eps
	end
end

%--------------------------------------------------------------------------------
% Fourier Spectrum vs. Wavelet Energy Density
%--------------------------------------------------------------------------------

if FOURIER
	
	filter = 0;
	ylog = 1;
	
	kymfourier(x,coimask,sig,unit,US,freqVec,roi,admis,filter,ylog)
	
	% Print output graph
	if (output)
		print -depsc fourierout.eps
	end
end

%--------------------------------------------------------------------------------
% CWT Peaks Detection
%--------------------------------------------------------------------------------

if PEAK
	figure
	
	peak(x,coimask,thr,unit,US,timeVec,freqVec,mark,5,false)
	
	title([num2str(thr),'% Thresholded - NU CWT Peak Detection - ROI ',num2str(par.r)],'FontSize',s3)
	
	% Print output graph
	if (output)
		print -depsc maximaout.eps
	end
end

%--------------------------------------------------------------------------------
% CWT Frequency Paths
%--------------------------------------------------------------------------------

if PATHS
	figure
	
	maxima = paths(x,y,coimask,thr,unit,US,timeVec,freqVec,mark,5,false,true);
	
	title([num2str(thr),'% Thresholded - NU CWT Frequency Paths - ROI ',num2str(par.r)],'FontSize',s3)
	
	% Print output graph
	if (output)
		print -depsc pathout.eps
	end
end

%--------------------------------------------------------------------------------
% CWT Skeleton
%--------------------------------------------------------------------------------

if SKELETON
	figure
	
	[skex,skey] = skeleton(x,y,coimask,thr,unit,US,timeVec,freqVec,mark,5,false,true);
	
	title([num2str(thr),'% Thresholded - NU CWT Skeleton - ROI ',num2str(par.r)],'FontSize',s3)
	
	% Print output graph
	if (output)
		print -depsc pathout.eps
	end
end

%--------------------------------------------------------------------------------
% Hurst Exponent (H)
%--------------------------------------------------------------------------------

if HURST
	
	% Number of dyadic iteration
	iter = 0;
	
	hurst(x,sig,unit,US,timeVec,freqVec,mark,iter);
	
	% Print output graph
	if (output)
		print -depsc holderout.eps
	end
end

%--------------------------------------------------------------------------------
% Holder Exponent (alpha)
%--------------------------------------------------------------------------------

if HOLDER
	
	holder(x,sig,unit,US,timeVec,freqVec,periodVec,skex);
	
	% Print output graph
	if (output)
		print -depsc holderout.eps
	end
end

%--------------------------------------------------------------------------------
% CWT Wave-Activity Index - WAI
%--------------------------------------------------------------------------------

if WAI
	figure
	
	indextype = 'I';
	
	wai(x,coimask,thr,sig,timeVec,freqVec,mark,maxima,indextype,true);
	
	title([num2str(thr),'% Thresholded WAI'],'FontSize',s3)

	% Print output graph
	if (output)
		print -depsc indexout.eps
	end
end

%--------------------------------------------------------------------------------
% CWT Vector Approach
%--------------------------------------------------------------------------------

if VECTOR
	if (length(mark) > 0)
		
		ratio = vector(x,coimask,freqVec,periodVec,mark,par.r,output);
		
	end
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