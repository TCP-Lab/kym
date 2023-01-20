function ratio = PD(wt,par,sig,thr1,thr2,output)

% 
%--------------------------------------------------------------------------------
% CWT Analysis and Peaks Detection
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% ratio = PD(wt,par,sig,thr1,thr2,output)
%
% INPUT       TYPE         MEANING
% -----       ----         -------
% wt       -> matrix    -> 1st WT Output - Continuous Wavelet Transform
% par      -> structure -> 2nd WT Output - Parameters
% sig      -> array     -> 3rd WT Output - Calcium Signal
% thr1     -> scalar    -> Peaks Detection Threshold Percentage
% thr2     -> scalar    -> Frequency Paths Threshold Percentage
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
mark = par.m;

x = abs(wt);
y = abs(real(wt));

nscale = size(x,1);
nvoice = (size(x,1))/(length(freqVec)-1);
n = size(x,2);

%---------------------------------------------------------------------------
% Input Control
%---------------------------------------------------------------------------

% Print plots - Default value = false = Print nothing
if (nargin < 6)
	output = false;
end

% Thresholds - Default value = 0
if (nargin < 5)
	thr2 = 0;
end
if (nargin < 4)
	thr1 = 0;
end

%--------------------------------------------------------------------------------
% CWT Analysis: selectively enable/disable each one of the following features
%--------------------------------------------------------------------------------

%--------------------------------------------------------------------------------
% Wavelet Power Spectrum
%--------------------------------------------------------------------------------

if 1
	figure
	
	energy(x,coimask,timeVec,freqVec,nvoice,mark)
	
	% Print output graph
	if (output)
		print -depsc energyout.eps
	end
end

%--------------------------------------------------------------------------------
% Fourier Spectrum vs. Wavelet Energy Density
%--------------------------------------------------------------------------------

if 1
	figure
	
	fourier(x,coimask,sig,freqVec,1,1)
	
	title(['Fourier Spectrum vs. Wavelet Energy Density - ROI ',num2str(par.r)],'FontSize',s3)
	
	% Print output graph
	if (output)
		print -depsc fourierout.eps
	end
end

%--------------------------------------------------------------------------------
% CWT Peaks Detection
%--------------------------------------------------------------------------------

if 1
	figure
	
	peak(x,coimask,thr1,timeVec,freqVec,mark,5,false)
	
	title([num2str(thr1),'% Thresholded - Continuous Wavelet Transform Peak Detection - ROI ',num2str(par.r)],'FontSize',s3)
	
	% Print output graph
	if (output)
		print -depsc maximaout.eps
	end
end

%--------------------------------------------------------------------------------
% CWT Frequency Paths
%--------------------------------------------------------------------------------

if 1
	figure
	
	[maxima,maximaNT] = paths(x,y,coimask,thr2,timeVec,freqVec,mark,5,false,true);
	
	title([num2str(thr2),'% Thresholded - Continuous Wavelet Transform Frequency Paths - ROI ',num2str(par.r)],'FontSize',s3)
	
	% Print output graph
	if (output)
		print -depsc pathout.eps
	end
end

%--------------------------------------------------------------------------------
% CWT Wave-Activity Index - WAI
%--------------------------------------------------------------------------------

if 1
	figure
	
	% Thresholded index
	wai(x,timeVec,freqVec,nvoice,mark,maxima,0,true);
	
	title([num2str(thr2),'% Thresholded WAI'],'FontSize',s3)
	
	% NON Thresholded index
	wai(x,timeVec,freqVec,nvoice,mark,maximaNT,1,true);
	
	title('NO Thresholded WAI','FontSize',s3)
	
	% Print output graph
	if (output)
		print -depsc indexout.eps
	end
end

%--------------------------------------------------------------------------------
% CWT Complete Wave-Activity Index - Complete-WAI
%--------------------------------------------------------------------------------

if 1
	figure
	
	% Thresholded index
	wtabst = (x-min(min(x)))/max(max(x-min(min(x))));
	threshold = thr2/100;
	wtabst = (wtabst >= threshold);
	wtabs = x .* wtabst .* coimask;
	wtabsNT = x .* coimask;
	
	wai(wtabs,timeVec,freqVec,nvoice,mark,[],0,true);
	
	title([num2str(thr2),'% Thresholded Complete-WAI'],'FontSize',s3)
	
	% NON Thresholded index
	wai(wtabsNT,timeVec,freqVec,nvoice,mark,[],1,true);
	
	title('NO Thresholded Complete-WAI','FontSize',s3)
	
	% Print output graph
	if (output)
		print -depsc indexcompout.eps
	end
end

%--------------------------------------------------------------------------------
% CWT Vector Approach
%--------------------------------------------------------------------------------

if 1
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