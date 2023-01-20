function ratio = PD(wt,par,sig,thr1,thr2,output)

%
%---------------------------------------------------------------------------
% Morlet Wavelet Transform Analysis and Peaks Detection
%---------------------------------------------------------------------------
%
%
% Function Definition
%
% PD(wt,par,thr1,thr2,output)
%
% INPUT       TYPE         MEANING
% -----       ----         -------
% wt       -> matrix    -> 1st WT Output - Morlet Wavelet Transform
% par      -> structure -> 2nd WT Output - Parameters
% thr1     -> scalar    -> Peaks Detection Threshold Percentage
% thr2     -> scalar    -> Frequency Paths Threshold Percentage
% output   -> boolean   -> Print .eps Output Graphs
%
% OUTPUT      TYPE         MEANING
% ------      ----         -------
% ratio    -> array     -> post/pre Ratio Vector
% -none-   -> plot      -> 6 Plots Resulting from Analysis
%

% Variables Assignment
timeVec = par.x;
freqVec = par.y;
periodVec = par.z;
lowest = freqVec(end)/(2**(length(freqVec)-1));
mark = par.m;

x = abs(wt);
y = abs(real(wt));

nscale = size(x)(1);
nvoice = (size(x)(1))/(length(freqVec)-1);
n = size(x)(2);

% Output Print control - Default value = false
if (nargin < 6)
	output = false;
endif
% Threshold control - Default value = 0
if (nargin < 5)
	thr2 = 0;
endif
if (nargin < 4)
	thr1 = 0;
endif


%--------------------------------------------------------------------------------
% Wavelet Power Spectrum
%--------------------------------------------------------------------------------

figure

energy(x,timeVec,freqVec,nvoice,mark)

% Print output graph
if (output)
	print -depsc energyout.eps
endif

%--------------------------------------------------------------------------------
% Fourier Spectrum vs. Wavelet Energy Density
%--------------------------------------------------------------------------------

figure

fourier(x,sig,freqVec,1)

title(['Fourier Spectrum vs. Wavelet Energy Density - ROI ',num2str(par.r)],'FontSize',18)

% Print output graph
if (output)
	print -depsc fourierout.eps
endif

%--------------------------------------------------------------------------------
% Morlet Wavelet Transform Peaks Detection
%--------------------------------------------------------------------------------

figure

peak(x,thr1,timeVec,freqVec,mark,5,false)

title([num2str(thr1),'% Thresholded - Morlet Wavelet Transform Peaks Detection - ROI ',num2str(par.r)],'FontSize',18)

% Print output graph
if (output)
	print -depsc maximaout.eps
endif

%--------------------------------------------------------------------------------
% Morlet Wavelet Transform Frequency Paths
%--------------------------------------------------------------------------------

figure

[maxmap1,maxmap2,maxima,maximaNT] = paths(x,y,thr2,mark,5,false);

subplot(2,1,1), hold on
imagesc(timeVec,[],maxmap1)
set(gca,'XTick',[get(gca,'XTick'),floor(timeVec(end))]);
set(gca,'YDir','normal');
set(gca,'YTick',[0:nvoice:nscale]); % Maximum number of displayable ticks on ordinates
set(gca,'YTick',[get(gca,'YTick'),1]); % In order to display freqVec(1) value
set(gca,'YTickLabel',freqVec);
ylabel('Frequency (mHz)','FontSize',18)
title([num2str(thr2),'% Thresholded - Morlet Wavelet Transform Frequency Paths - ROI ',num2str(par.r)],'FontSize',18)

subplot(2,1,2), hold on
imagesc(timeVec,[],maxmap2)
set(gca,'XTick',[get(gca,'XTick'),floor(timeVec(end))]);
set(gca,'YDir','normal');
set(gca,'YTick',[0:nvoice:nscale]); % Maximum number of displayable ticks on ordinates
set(gca,'YTick',[get(gca,'YTick'),1]); % In order to display freqVec(1) value
set(gca,'YTickLabel',freqVec);
ylabel('Frequency (mHz)','FontSize',18)
xlabel('Time (s)','FontSize',18)

% Print output graph
if (output)
	print -depsc pathout.eps
endif

%--------------------------------------------------------------------------------
% Morlet Wavelet Transform Wave-Activity Index (WAI)
%--------------------------------------------------------------------------------

figure

% l = Regularization Window Width = Mean Convolutive Filter
% l value must be odd! - l=1 means NO FILTER
l = 7;

% Threshold Index
wai(x,timeVec,lowest,nvoice,mark,maxima,l,0)

title([num2str(thr2),'% Thresholded WAI'],'FontSize',18)

% NO Threshold Index
wai(x,timeVec,lowest,nvoice,mark,maximaNT,l,1)

title('NO Thresholded WAI','FontSize',18)

% Print output graph
if (output)
	print -depsc indexout.eps
endif

%--------------------------------------------------------------------------------
% Morlet Wavelet Transform Complete-Wave-Activity Index (WAI)
%--------------------------------------------------------------------------------

figure

% Threshold Index

minimox = min(min(x));
xt = x - minimox;
massimox = max(max(xt));
xt = xt / massimox;

threshold = thr2/100;

xt = (xt >= threshold);
xt = x .* xt;

waiComp(xt,timeVec,freqVec,nvoice,mark,0)

title([num2str(thr2),'% Thresholded Complete-WAI'],'FontSize',18)

% NO Threshold Index

waiComp(x,timeVec,freqVec,nvoice,mark,1)

title('NO Thresholded Complete-WAI','FontSize',18)

% Print output graph
if (output)
	print -depsc indexcompout.eps
endif

%--------------------------------------------------------------------------------
% Morlet Wavelet Transform Vector Approach
%--------------------------------------------------------------------------------

if (length(mark) < 3)
	
	figure
	
	ratio = vector(x,freqVec,periodVec,mark);
	
	title(['Morlet Wavelet Transform Vector Analysis - ROI ',num2str(par.r)],'FontSize',18)
	
	% Print output graph
	if (output)
		print -depsc vectorout.eps
	endif
	
endif


%---------------------------------------------------------------------%
%                                                                     %
% A.A. 2009 / 2010                                                    %
% Original code by Federico Alessandro Ruffinatti                     %
% Università degli Studi di Torino - Italy - DBAU - Scienze MFN       %
% Scuola di Dottorato in Neuroscienze - XXV ciclo                     %
%                                                                     %
% Wavelet computation is regarded as a time convolution and it is     %
% implemented as a product in the Fourier transformed domain.         %
% A standard code for this algorithm can be found, for instance,      %
% in WaveLab850 - http://www-stat.stanford.edu/~wavelab/              %
%                                                                     %
% Peaks detection uses a technique that is based on images dilation.  %
% See, for instance, localMaximum.m m-file by Yonathan Nativ          %
% http://www.mathworks.com/matlabcentral/fileexchange/authors/26510/  %
%                                                                     %
%---------------------------------------------------------------------%