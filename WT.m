function [wt,par,sig] = WT(filename,unit,mark,roi,lft,nvoice,wavelet,cone,output)

%
%--------------------------------------------------------------------------------
% Continuous Wavelet Transform - CWT
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% [wt,par,sig] = WT(filename,unit,mark,roi,lft,nvoice,wavelet,cone,output)
%
% INPUT       TYPE       EXAMPLE           MEANING
% -----       ----       -------           -------
% filename -> string  -> 'filename.csv' -> File Name
% unit     -> string  -> 's'            -> Time Unit: 's' or 'ms'
% mark     -> array   -> [300,450,500]  -> Time Marker Set (Step Number)
% roi      -> array   -> [2,3,9]        -> ROI Set ([x:y] = from x to y)
% lft      -> scalar  -> 4              -> Low-Frequency Threshold
% nvoice   -> scalar  -> 48             -> Amount of Inter-Octave Frequencies
% wavelet  -> string  -> 'M1'           -> Mother Wavelet Type
% cone     -> string  -> 'PAD'          -> Cone of Influence Handling Method
% output   -> boolean -> true           -> Print .eps Output Graph
%
% OUTPUT      TYPE                         MEANING
% ------      ----                         -------
% wt       -> matrix                    -> Last ROI Wavelet Transform
% par      -> structure                 -> Last ROI Parameters
% sig      -> array                     -> Last ROI Calcium Signal
% -none-   -> plot                      -> N Plots Resulting from Analysis
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

wt = [];
par = [];
sig = [];

%---------------------------------------------------------------------------
% Input Control
%---------------------------------------------------------------------------

% Print plots - Default value = false = Print nothing
if (nargin < 9)
	output = false;
end

% Cone Of Influence (COI) handling method - Default value = 'NONE'
if (nargin < 8)
  cone = 'NONE';
end
if ~(strcmp(cone,'NONE') | strcmp(cone,'PAD') | strcmp(cone,'COI'))
	fprintf('\n\nWARNING: Invalid COI Handling Method\n');
	fprintf('\n\n');
	return
end

% Wavelet version - Default value = 'M1'
if (nargin < 7)
  wavelet = 'M1';
end
if ~(strcmp(wavelet,'M1') | strcmp(wavelet,'M2') | strcmp(wavelet,'M3'))
	fprintf('\n\nWARNING: Invalid Mother Wavelet Type\n');
	fprintf('\n\n');
	return
end

% Voices number - Default value = 24 voices per octave
if (nargin < 6)
  nvoice = 24;
end

% Low-Frequency Threshold - Default value = 3
% This constant sets the amount of discarded octaves starting from the lowest frequency
% Lowest displayable frequency (lft=0) turns out to be f.0 = 1/(2*N*dt),
% where N is the largest power of 2 not greater than n: N = 2^(floor(log2(n))),
% even if the lowest significant frequency has always to be considered equal to f.low = 1/(n*dt)
% High-Frequency threshold turns always out to be = f.Nyquist = 1/(2*dt) automatically
if (nargin < 5)
	lft = 3;
end

% ROI control - Default value = All traces
if (nargin < 4)
	roi = [];
end

% Marker array - Default value = [] = No markers
if (nargin < 3)
	mark = [];
end

% Time unit - Default value = 's'
if (nargin < 2)
	unit = 's';
end
if ~(strcmp(unit,'s') | strcmp(unit,'ms'))
	fprintf('\n\nWARNING: Invalid Time Unit\n');
	fprintf('\n\n');
	return
end

%---------------------------------------------------------------------------
% Data Loading
%---------------------------------------------------------------------------

switch (filename)
	
	case 'test1' % Plane wave test (for Calibration): data = [time_dt=1s wave2_500mHz wave3_250mHz wave4_125mHz wave5_62.5mHz ...]
		
		OctFreq = 1./(2.^[1:10]); % Octave frequency (in Hz); OctFreq(1)=1/2*dt according to Nyquist
		data = (2*pi*OctFreq)'*[1:1/OctFreq(end)]; % (2*pi*nu*t); timeVec(end)=1/OctFreq(end)
		data = data';
		data = exp(i*data);
		%data = real(data); % Cosine
		%data = imag(data); % Sine
		data = [[1:1/OctFreq(end)]',data];
		
	case 'test2' % Edge effect masking test (with a linear trend)
		
		data = [1:1024]';
		data = [data (2/data(end))*data-1];
		
		lft = 4;
		cone = 'COI';
		
	otherwise % Open empirical data file
		
		data = dlmread(filename);
		
end

% M-Downsampling
%M = 2;
%data = data([1:M:size(data,1)],:);

% Samples must be even
if (mod(size(data,1),2) ~= 0)
	data = data(1:size(data,1)-1,:);
end

% Extract time vector
timeVec = data(:,1);
timeVec = timeVec - timeVec(1); % Start from t=0s
if (strcmp(unit,'ms')) % Measured in s
	timeVec = timeVec/1000;
end
dt = mode(diff(timeVec)); % Sampling time mode

% Data matrix
data(:,1) = [];
n = size(data,1);

% Select data of interest, keeping into account ROI 1 is the background
if (length(roi) == 0)
	roi = [2:size(data,2)+1];
end
data = data(:,roi(1:end)-1);

%---------------------------------------------------------------------------
% Characteristic parameters
%---------------------------------------------------------------------------

noctave = floor(log2(n))-lft;
nscale = nvoice*noctave;
scaleVec = 2.^[1:noctave+1];
scaleVec = fliplr(scaleVec);
% Express frequency in mHz, with 1 decimal digit; freqVec(end) = f.Nyquist = 1/(2*dt)
freqVec = floor((1000./(scaleVec*dt))*10)/10; % In order to get only 1 decimal digit
lowest = 1000/(2*dt*2^(length(freqVec)-1)); % To preserve an accurate value for freqVec(1)
% Express period in s, with 1 decimal digit - periodVec(end) = 1/f.Nyquist = 2*dt
periodVec = floor((scaleVec*dt)*10)/10; % In order to get only 1 decimal digit

%---------------------------------------------------------------------------
% Transform and Plot
%---------------------------------------------------------------------------

for j = 1:size(data,2)
	
	% Original trace
	yraw = data(:,j);
	y = yraw;
	
	%---------------------------------------------------------------------------
	% Denoise
	%---------------------------------------------------------------------------
	
	% Trace denoising
	% l = Regularization window width = Mean convolutive filter
	% l value must be odd! - l=1 means NO FILTER
	l = 3;
	
	for h = 1:n-l+1
		y(h+(l-1)/2) = (1./l)*sum(y(h:h+l-1));
	end
	sig = y;
	
	%---------------------------------------------------------------------------
	% Detrend and Pad
	%---------------------------------------------------------------------------
	
	if (strcmp(cone,'PAD'))
		% Trace detrending
		y1=y(1);
		x1=timeVec(1);
		y2=y(end);
		x2=timeVec(end);
		trend = (y2-y1)/(x2-x1).*(timeVec-x1)+y1;
		y = y-trend;
		
		sig = y;
		
		% Padding with n zeros
		y = [y;zeros(n,1)];
	end
	
	%---------------------------------------------------------------------------
	% Wavelet Transform
	%---------------------------------------------------------------------------
	
	[wt,fac] = trans(y,lft,nvoice,wavelet);
	
	%---------------------------------------------------------------------------
	% Cone of Influence (COI)
	%---------------------------------------------------------------------------
	
	coimask = coi(nvoice,nscale,n,dt,lowest,fac,cone);
	
	%---------------------------------------------------------------------------
	% Rescale to [0,1] and Convert to RGB
	%---------------------------------------------------------------------------
	
	if (strcmp(cone,'PAD')) % Cut the "zero zone": every point greater than n and the lowest octave
		wt = wt(nvoice+1:end,1:n);
	end
	
	wtabs = abs(wt);
	wtreal = abs(real(wt));
	
	if strcmp(filename,'test2') % e-folding Test
		r1 = wtabs(1:nscale,1)';
		for h = 1:nscale
			r2(h) = wtabs(h,min(find(coimask(h,:))));
		end
		fold = (r1./r2)'
		mean(fold)
	end
	
	massimo = max([max(max(wtabs)),max(max(wtreal))]);
	minimo = min([min(min(wtabs)),min(min(wtreal))]);
	
	wtabs = (wtabs-minimo);
	wtreal = (wtreal-minimo);
	numax = max([max(max(wtabs)),max(max(wtreal))]);
	wtabs = wtabs/numax;
	wtreal = wtreal/numax;
	
	[wtabsind,b] = gray2ind(wtabs,128); % gray2ind: floating point images may only contain values between 0 and 1
	[wtrealind,b] = gray2ind(wtreal,128);
	
	wtabs = ind2rgb(wtabsind,jet(128));
	wtreal = ind2rgb(wtrealind,jet(128));
	
	%---------------------------------------------------------------------------
	% Apply COI Mask
	%---------------------------------------------------------------------------
	
	% Gray scale scaleograms
	gswtabs = ind2rgb(wtabsind,gray(128));
	gswtreal = ind2rgb(wtrealind,gray(128));
	
	% Shade masked regions gray
	index = find(coimask == 0);
	
	% Gray scale - If length(index)==0 this does nothing
	wtabs(index) = gswtabs(index);
	wtabs(index + n*nscale) = gswtabs(index + n*nscale);
	wtabs(index + 2*n*nscale) = gswtabs(index + 2*n*nscale);
	
	% Gray scale - If length(index)==0 this does nothing
	wtreal(index) = gswtreal(index);
	wtreal(index + n*nscale) = gswtreal(index + n*nscale);
	wtreal(index + 2*n*nscale) = gswtreal(index + 2*n*nscale);
	
	%---------------------------------------------------------------------------
	% Plot Time Course
	%---------------------------------------------------------------------------
	
	figure
	
	subplot(3,2,1), hold on
		
		plot(timeVec,yraw,'-b','LineWidth',1)
		
		if (strcmp(cone,'PAD'))
			plot(timeVec,trend,'-g','LineWidth',1)
			plot(timeVec,sig+mean(yraw),'-r','LineWidth',1)
		else
			plot(timeVec,sig,'-r','LineWidth',1)
		end
		
		xlim([timeVec(1),timeVec(end)])
		set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
		ylabel('Ratio','FontSize',s2)
		title('Continuous Wavelet Transform - Modulus and Real Part','FontSize',s3)
		text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['ROI ',num2str(roi(j))],'FontSize',s2,'Color','k')
		
		for k = 1:length(mark) % If length(mark)==0 this does nothing
			plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
		end
	
	%---------------------------------------------------------------------------
	% Plot Modulus and Real Part
	%---------------------------------------------------------------------------
	
	subplot(3,2,3), hold on
		
		image(timeVec(1):timeVec(end),1:size(wtabs,1),wtabs)
		xlim([timeVec(1),timeVec(end)])
		ylim([0.5,size(wtabs,1)+0.5])
		set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
		set(gca,'YDir','normal','YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
		set(gca,'YTickLabel',num2str(freqVec'));
		ylabel('Frequency (mHz)','FontSize',s2)
		
	subplot(3,2,5), hold on
		
		image(timeVec(1):timeVec(end),1:size(wtreal,1),wtreal)
		xlim([timeVec(1),timeVec(end)])
		ylim([0.5,size(wtreal,1)+0.5])
		set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
		set(gca,'YDir','normal','YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
		set(gca,'YTickLabel',num2str(freqVec'));
		ylabel('Frequency (mHz)','FontSize',s2)
		xlabel('Time (s)','FontSize',s2)
		
	%---------------------------------------------------------------------------
	% Add Color Bar
	%---------------------------------------------------------------------------
	
	subplot(3,2,[2 4 6]), hold on
	
		t = [0:1:100];
		imagesc(1:2,1:101,t') % X-support [1:2] prevents a "division by zero" in OCTAVE
		xlim([0.5,2.5])
		ylim([0.5,101.5])
		set(gca,'FontSize',s1,'XTick',[0]);
		set(gca,'XTickLabel','');
		set(gca,'YDir','normal','YTick',[1:10:101]);
		magnVec = floor([minimo:(massimo-minimo)/10:massimo]*1000*10)/10; % In order to get only 1 decimal digit
		set(gca,'YTickLabel',num2str(magnVec'));
		xlabel('10^3 Amp','FontSize',s2)
		
	%---------------------------------------------------------------------------
	% Resize - Sintax template: set(gca,'Position',[left bottom width height])
	%---------------------------------------------------------------------------
	
	set(gca,'Position',[0.96,0.11,0.02,0.815])
	
	subplot(3,2,1)
	set(gca,'Position',[0.12,0.70926,0.77,0.21574])
	
	subplot(3,2,3)
	set(gca,'Position',[0.12,0.40963,0.77,0.21574])
	
	subplot(3,2,5)
	set(gca,'Position',[0.12,0.11,0.77,0.21574])
	
	%---------------------------------------------------------------------------
	% Print output graph
	%---------------------------------------------------------------------------
	
	if (output)
		print -depsc transformout.eps
	end
	
end

%---------------------------------------------------------------------------
% Return last ROI parameters as output - Using data structures
%---------------------------------------------------------------------------

par.k = coimask;
par.x = timeVec;
par.y = freqVec;
par.z = periodVec;
par.r = roi(end);
par.m = mark;


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