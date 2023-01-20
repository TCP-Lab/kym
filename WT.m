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
% -none-   -> plot                      -> 1 Plot Resulting from Analysis
%

% Open file
data = dlmread(filename);

% Plane Wave Test 
%data = [1:1:1024]';
%data = [data,e.**(i*2*pi*7.8125/1000*data),sin(2*pi*7.8125/1000*data),cos(2*pi*7.8125/1000*data)];

% M-Downsampling
%M = 2;
%data = data([1:M:size(data)(1)],:);

% Samples must be even
if (mod((size(data))(1),2) != 0)
	data = data(1:(size(data))(1)-1,:);
endif

% Time Unit control - Default value = 's'
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
if (nargin < 9)
	output = false;
endif
% Cone of Influence Handling Method - Default value = 'NONE'
if (nargin < 8)
  cone = 'NONE';
endif
% Wavelet Version control - Default value = 'M1'
if (nargin < 7)
  wavelet = 'M1';
endif
% Voices Number control - Default value = 24 voices
if (nargin < 6)
  nvoice = 24;
endif
% Low-Frequency Threshold control - Default value = 3
% This constant sets the amount of discarded octaves starting from the lowest frequency
% Lowest Displayable Frequency (lft=0) turns out to be f.0 = 1/(2*N*dt),
% where N is the largest power of two not greater than n: N = 2**(floor(log2(n))),
% even if the lowest significant frequency has always to be considered equal to f.low = 1/(n*dt)
% High-Frequency Threshold turns always out to be = f.Nyquist = 1/(2*dt) automatically
if (nargin < 5)
	lft = 3;
endif
% ROIs control - Default value = All traces
% Keeping into account ROI 1 is the background
if (nargin < 4)
	roi = [2:(size(data))(2)+1];
elseif (length(roi) == 0)
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
lowest = 1000/((2*dt)*2**(length(freqVec)-1)); % To preserve an accurate value for freqVec(1)
% Express period in s, with 1 decimal digit - periodVec(end) = 1/f.Nyquist = 2*dt
periodVec = floor((scaleVec*dt)*10)/10; % In order to have only 1 decimal digit

% Trace Denoising
% l = Regularization Window Width = Mean Convolutive Filter
% l value must be odd! - l=1 means NO FILTER
l = 3;

%---------------------------------------------------------------------------
% Transform and Plot
%---------------------------------------------------------------------------

for j = 1:(size(data))(2)
	
	% Original trace
	yraw = data(:,j);
	y = yraw;
	
	%---------------------------------------------------------------------------
	% Denoising, Padding and Fourier Transform
	%---------------------------------------------------------------------------
	
	% Trace Denoising
	for h = 1:n-l+1
		y(h+(l-1)/2) = (1./l)*sum(y(h:h+l-1));
	endfor
	sig = y;
	
	if (strcmp(cone,'PAD'))
		
		% Trace Detrending
		%y = detrend(y,2); % Less effective alternative
		y1=y(1);
		x1=timeVec(1);
		y2=y(end);
		x2=timeVec(end);
		trend = (y2-y1)/(x2-x1).*(timeVec-x1)+y1;
		y = y-trend;
		
		sig = y;
		
		% Padding with n zeros
		y = [y;zeros(n,1)];
		
	endif
	
	% Fourier Transform of the Signal
	Y = fft(y);
	
	%---------------------------------------------------------------------------
	% Plotting Original Trace
	%---------------------------------------------------------------------------
	
	figure
	
	subplot(3,2,1), hold on
	
	plot(timeVec,yraw,'-b','LineWidth',1)
	
	if (strcmp(cone,'PAD'))
		plot(timeVec,trend,'-g','LineWidth',1)
		plot(timeVec,y(1:n)+mean(yraw),'-r','LineWidth',1)
	else
		plot(timeVec,y,'-r','LineWidth',1)
	endif
	
	for k = 1:length(mark)
		plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
	endfor
	
	axis([timeVec(1),timeVec(end)])
	set(gca,'XTick',[get(gca,'XTick'),floor(timeVec(end))]);
	text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),[" ROI ", num2str(roi(j))],'FontSize',18,'Color','k')
	ylabel('Ratio','FontSize',18)
	title('Continuous Wavelet Transform - Modulus and Real Part','FontSize',18)
	
	if (strcmp(wavelet,'M1'))
		
		%---------------------------------------------------------------------------
		% Morlet Wavelet ver.1
		%---------------------------------------------------------------------------
		
		[wt,fac] = morlet1(Y,2,lft,nvoice);
		
	elseif (strcmp(wavelet,'M2'))
		
		%---------------------------------------------------------------------------
		% Morlet Wavelet ver.2
		%---------------------------------------------------------------------------
		
		[wt,fac] = morlet2(Y,1.5,lft,nvoice);
		
	elseif (strcmp(wavelet,'M3'))
		
		%---------------------------------------------------------------------------
		% Morlet Wavelet ver.3
		%---------------------------------------------------------------------------
		
		[wt,fac] = morlet3(Y,1.5,lft,nvoice);
	
	endif
	
	%---------------------------------------------------------------------------
	% Cone of Influence Handling
	%---------------------------------------------------------------------------
	
	if (strcmp(cone,'COI'))
		
		coimask = ones(size(wt)(1),size(wt)(2));
		
		COI = fac./(lowest.*(2.**([1:1:size(wt)(1)]./nvoice)));
		COI = (COI.*1000)./dt;
		COI(find(COI > n)) = n;
		
		% e-folding Test (with Lin1.csv and Lin2.csv)
		%r1=(abs(wt)(1:size(wt)(1),1))';
		%for c = 2*nvoice:size(wt)(1)
		%	r2(c)=abs(wt)(c,floor(COI(c)));
		%endfor
		%r1./r2
		%mean(r1./r2)
		
		% You can use a more severe COI to further reduce edge effects on WAI indexes (especially J)
		%COI = COI * 1.5;
		%COI(find(COI > n)) = n;
		
		for u = 1:size(wt)(1)
			coimask(u,1:floor(COI(u))) = 0;
			coimask(u,n-floor(COI(u))+1:n) = 0;
		endfor
		
		wtgraph = wt .* coimask;
		
	elseif (strcmp(cone,'PAD'))
		
		% Cut the "zeros zone": every point greater than n and the lower octave
		wt = wt(nvoice+1:end,1:n);
		coimask = ones(size(wt)(1),size(wt)(2));
		wtgraph = wt;
		
	elseif (strcmp(cone,'NONE'))
		
		% Do Nothing
		coimask = ones(size(wt)(1),size(wt)(2));
		wtgraph = wt;
		
	endif
	
	%---------------------------------------------------------------------------
	% Plotting Modulus and Real Part
	%---------------------------------------------------------------------------
	
	subplot(3,2,3), hold on
		imagesc(timeVec,[],abs(wtgraph))
		%%%imagesc(timeVec,[],real(wtgraph))
		set(gca,'XTick',[get(gca,'XTick'),floor(timeVec(end))]);
		set(gca,'YDir','normal');
		set(gca,'YTick',[0:nvoice:nscale]); % Maximum number of displayable ticks on ordinates
		set(gca,'YTick',[get(gca,'YTick'),1]); % In order to display freqVec(1) value
		set(gca,'YTickLabel',freqVec);
		ylabel('Frequency (mHz)','FontSize',18)
		
	subplot(3,2,5), hold on
		imagesc(timeVec,[],abs(real(wtgraph)))
		%%%imagesc(timeVec,[],imag(wtgraph))
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
par.k = coimask;
par.x = timeVec;
par.y = freqVec;
par.z = periodVec;
par.r = roi(end);
par.m = mark;


%---------------------------------------------------------------------%
%                                                                     %
% A.A. 2009/2010 - 2010/2011                                          %
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