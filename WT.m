function [wt,par,sig] = WT(filename,unit,mark,roi,lft,nvoice,wavelet,output)

%
%---------------------------------------------------------------------------
% Morlet Wavelet Transform
%---------------------------------------------------------------------------
%
%
% Function Definition
%
% [wt,par] = WT(filename,unit,mark,roi,lft,nvoice,wavelet,output)
%
% INPUT       TYPE       EXAMPLE           MEANING
% -----       ----       -------           -------
% filename -> string  -> 'filename.csv' -> File Name
% unit     -> string  -> 's'            -> Time Unit: 's' or 'ms'
% mark     -> array   -> [300,450,500]  -> Time Marker Set (Step Number)
% roi      -> array   -> [2,3,9]        -> ROI Set
% lft      -> scalar  -> 4              -> Low-Frequency Threshold
% nvoice   -> scalar  -> 48             -> Lines of Pixel per Octave
% wavelet  -> string  -> 'M1'           -> Morlet Wavelet Version
% output   -> boolean -> true           -> Print .eps Output Graph
%
% OUTPUT      TYPE                         MEANING
% ------      ----                         -------
% wt       -> matrix                    -> Last ROI Morlet Wavelet Transform
% par      -> structure                 -> Last ROI Parameters
% -none-   -> plot                      -> 1 Plot Resulting from Analysis
%

% Open file
data = dlmread(filename);

% Samples must be even
if (mod((size(data))(1),2) != 0)
	data = data(1:(size(data))(1)-1,:);
endif

% Time Unit of Measure control - Default value = 's'
if (nargin < 2)
	unit = 's';
endif

% Extract time vector - Measured in s - Starting from t = 0s
timeVec = data(:,1);
timeVec = timeVec - timeVec(1);
if (strcmp(unit,'ms'))
	timeVec = timeVec/1000;
endif
dt = mean(diff(timeVec)); % Sampling Rate: df = 1/dt;

% Data Matrix
data(:,1) = [];
[n,noTraces] = size(data);

% Output Print control - Default value = false
if (nargin < 8)
	output = false;
endif
% Wavelet Version control - Default value = 'M1'
if (nargin < 7)
  wavelet = 'M1';
endif
% Voices Number control - Default value = 24 voices
if (nargin < 6)
  nvoice = 24;
endif
% Low-Frequency Threshold control - Default value = 4
% This constant sets the amount of discarded octaves starting from the lowest frequency = 1/(n*dt)
% High-Frequency Threshold turns always out to be = f.Nyquist = 1/(2*dt) automatically
if (nargin < 5)
	lft = 4;
endif
% ROIs control - Default value = All traces
% Keeping into account ROI 1 is the background
if (nargin < 4)
	roi = [2:(size(data))(2)+1];
endif
% Marker Array control - Default value = No markers
if (nargin < 3)
	mark = [];
endif

% Select data of interest
% Keeping into account ROI 1 is the background
data = data(:,roi(1:end)-1);

% Characteristic Parameters
noctave = floor(log2(n))-lft;
nscale = nvoice*noctave;
scaleVec = 2.**[1:noctave+1];
scaleVec = fliplr(scaleVec);
% Express frequency in mHz, with 1 decimal digit - freqVec(end) = f.Nyquist = 1/(2*dt)
freqVec = floor((1 ./ (scaleVec*dt))*1000*10)/10; % In order to have only 1 decimal digit
% Express period in s, with 1 decimal digit - periodVec(end) = 1/f.Nyquist = 2*dt
periodVec = floor((scaleVec*dt)*10)/10; % In order to have only 1 decimal digit

% Compute and Plot
for j = 1:(size(data))(2)
	
	% Original trace
	yraw = data(:,j);
	
	% Trace Denoising
	% l = Regularization Window Width = Mean Convolutive Filter
	% l value must be odd! - l=1 means NO FILTER
	y = yraw;
	l = 3;
	for h = 1:n-l+1
		y(h+(l-1)/2) = (1./l)*sum(y(h:h+l-1));
	endfor
	
	% Trace Detrending
	%y = detrend(y,2); % Less effective alternative
	y1=mean(y(1:10));
	x1=mean(timeVec(1:10));
	y2=mean(y(end+1-10:end));
	x2=mean(timeVec(end+1-10:end));
	trend = (y2-y1)/(x2-x1).*(timeVec-x1)+y1;
	y = y-trend;
	sig = y;
	
	% Fourier Transform of signal
	Y = fft(y);
	
	figure
	
	%---------------------------------------------------------------------------
	% Original Trace
	%---------------------------------------------------------------------------
	
	subplot(3,2,1), hold on
	
	plot(timeVec,yraw,'-b','LineWidth',1)
	plot(timeVec,y+mean(yraw),'-r','LineWidth',1)
	plot(timeVec,trend,'-g','LineWidth',1)
	for k = 1:length(mark)
		plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
	endfor
	
	axis([timeVec(1),timeVec(end)])
	set(gca,'XTick',[get(gca,'XTick'),floor(timeVec(end))]);
	text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),[" ROI ", num2str(roi(j))],'FontSize',18,'Color','k')
	ylabel('Ratio','FontSize',18)
	title('Morlet Wavelet Transform - Modulus and Real Part','FontSize',18)
	
	if (strcmp(wavelet,'M1'))
		
		%---------------------------------------------------------------------------
		% Morlet Wavelet ver.1
		%---------------------------------------------------------------------------
		
		wt = morlet1(n,Y,2,lft,nvoice);
		
	elseif (strcmp(wavelet,'M2'))
		
		%---------------------------------------------------------------------------
		% Morlet Wavelet ver.2
		%---------------------------------------------------------------------------
		
		wt = morlet2(n,Y,1.5,lft,nvoice);
		
	elseif (strcmp(wavelet,'M3'))
		
		%---------------------------------------------------------------------------
		% Morlet Wavelet ver.3
		%---------------------------------------------------------------------------
		
		wt = morlet3(n,Y,1.5,lft,nvoice);
	
	endif
	
	%---------------------------------------------------------------------------
	% Plotting Modulus and Real Part
	%---------------------------------------------------------------------------
	
	subplot(3,2,3), hold on
		imagesc(timeVec,[],abs(wt))
		set(gca,'XTick',[get(gca,'XTick'),floor(timeVec(end))]);
		set(gca,'YDir','normal');
		set(gca,'YTick',[0:nvoice:nscale]); % Maximum number of displayable ticks on ordinates
		set(gca,'YTick',[get(gca,'YTick'),1]); % In order to display freqVec(1) value
		set(gca,'YTickLabel',freqVec);
		ylabel('Frequency (mHz)','FontSize',18)
		
	subplot(3,2,5), hold on
		imagesc(timeVec,[],abs(real(wt)))
		set(gca,'XTick',[get(gca,'XTick'),floor(timeVec(end))]);
		set(gca,'YDir','normal');
		set(gca,'YTick',[0:nvoice:nscale]); % Maximum number of displayable ticks on ordinates
		set(gca,'YTick',[get(gca,'YTick'),1]); % In order to display freqVec(1) value
		set(gca,'YTickLabel',freqVec);
		ylabel('Frequency (mHz)','FontSize',18)
		xlabel('Time (s)','FontSize',18)
	
	%---------------------------------------------------------------------------
	% Resize and Add Color Bar
	%---------------------------------------------------------------------------
	
	% Non so perchè, ma l'ordine con cui si ridimensionano i subplot deve essere questo...
	% ...altrimenti qualcuno non viene visualizzato.
	
	% Sintax Legend: set(gca,'position',[left bottom width height])
	
	subplot(3,2,[2 4 6]), hold on
	t = [1:1:100];
	t = [t;t;t;t;t];
	imagesc(t')
	set(gca,'XTick',[]);
	set(gca,'YDir','normal');
	set(gca,'YTick',[0:10:100]);
	set(gca,'YTick',[get(gca,'YTick'),1]); % In order to display magnVec(1) value
	massimo = max(max(abs(wt)));
	minimo = min(min(abs(wt)));
	magnVec = floor([minimo:(massimo-minimo)/10:massimo]*1000*10)/10; % In order to have only 1 decimal digit
	set(gca,'YTickLabel',magnVec);
	xlabel('Amp x 1000','FontSize',8)
	set(gca,'position',[0.96,0.11,0.02,0.815])
	
	subplot(3,2,1)
	set(gca,'position',[0.12,0.70926,0.77,0.21574])
	
	subplot(3,2,3)
	set(gca,'position',[0.12,0.40963,0.77,0.21574])
	
	subplot(3,2,5)
	set(gca,'position',[0.12,0.11,0.77,0.21574])
	
	% Print output graph
	if (output)
		print -depsc transformout.eps
	endif
	
endfor

% Return last ROI as output - Using Data Structures
par.x = timeVec;
par.y = freqVec;
par.z = periodVec;
par.r = roi(end);
par.m = mark;


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