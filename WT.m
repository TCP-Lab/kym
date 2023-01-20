function [wt,par,sig] = WT(filename,unit,mark,roi,lft,nvoice,wavelet,cone,addoct,output)

%
%--------------------------------------------------------------------------------
% Non Unitary Continuous Wavelet Transform - NU CWT
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% [wt,par,sig] = WT(filename,unit,mark,roi,lft,nvoice,wavelet,cone,addoct,output)
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
% addoct   -> scalar  -> 2              -> Number of Octaves Added by Upsampling
% output   -> boolean -> true           -> Print .eps Output Graph
%
% OUTPUT      TYPE                         MEANING
% ------      ----                         -------
% wt       -> matrix                    -> Last ROI Wavelet Transform
% par      -> structure                 -> Last ROI Parameters
% sig      -> array                     -> Last ROI Signal
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

% Print .eps plots - Default value = false = Print nothing
if (nargin < 10)
	output = false;
end

% Number of Octaves Added by Upsampling - Default value = 1 -> 2^1 Upsamplig Factor
if (nargin < 9)
	addoct = 1;
end

% Cone Of Influence (COI) handling method - Default value = 'NONE'
if (nargin < 8)
  cone = 'NONE';
end
if ~(strcmp(cone,'NONE') || strcmp(cone,'PAD') || strcmp(cone,'COI'))
	fprintf('\n\nWARNING: Invalid COI Handling Method\n');
	fprintf('\n\n');
	return
end

% Wavelet version - Default value = 'M1'
if (nargin < 7)
  wavelet = 'M1';
end
if ~(strcmp(wavelet,'M1') || strcmp(wavelet,'M2') || strcmp(wavelet,'R') || strcmp(wavelet,'D2') || strcmp(wavelet,'D1'))
	fprintf('\n\nWARNING: Invalid Mother Wavelet Type\n');
	fprintf('\n\n');
	return
end

% Voices number - Default value = 24 voices per octave
if (nargin < 6)
  nvoice = 24;
end

% Low-Frequency Threshold (lft) - Default value = 3
% This constant sets the amount of discarded octaves starting from the lowest frequency.
% Lowest displayable frequency (lft=0) is f.0==1/(2*N*dt),
% where N is the largest power of 2 not greater than n: N = 2^(floor(log2(n))),
% even if the lowest significant frequency (apart from the DC component) is f.low==1/(n*dt).
% High-Frequency limit is always == f.Nyquist==1/(2*dt) by construction.
if (nargin < 5)
	lft = 3;
end

% ROI control - Default value = All traces
if (nargin < 4)
	roi = [];
end
if ~isempty(find(roi < 2))
	fprintf('\n\nWARNING: Invalid ROI Number\n');
	fprintf('It must be positive and greater than 2 (ROI 1 is the background noise)\n');
	fprintf('\n\n');
	return
end

% Marker array - Default value = [] = No markers
if (nargin < 3)
	mark = [];
end

% Time unit - Default value = 's'
if (nargin < 2)
	unit = 's';
end
if ~(strcmp(unit,'s') || strcmp(unit,'ms'))
	fprintf('\n\nWARNING: Invalid Time Unit\n');
	fprintf('\n\n');
	return
end

%---------------------------------------------------------------------------
% Data Loading and Trace Selection
%---------------------------------------------------------------------------

DS = 1; %Downsampling Rate

[data,timeVec,dt,roi] = loader(filename,unit,roi,DS);
if (isempty(timeVec))
	return
end

% Data size
[n,ntrace] = size(data);

%---------------------------------------------------------------------------
% Bridge Detrending
%---------------------------------------------------------------------------

y1 = data(1,:);
y2 = data(end,:);
x2 = timeVec(end);
trend = timeVec*((y2-y1)/x2) + ones(n,1)*y1;
if strcmp(filename,'test1')
	trend = zeros(n,ntrace);
end
yraw = data;
data = data - trend;

%---------------------------------------------------------------------------
% Denoising
%---------------------------------------------------------------------------

% Filter Type
denoisetype = 'FIRflat';

% Filter Order
% For symmetric FIR filters, only even orders are allowed (if denoisetype='wmean' -> L=1 anyway)
L = 8;

% Only for 'butterworth' and 'FIRflat' - Cutoff frequency for the 3dB point below the passband value
% It is specified in normalized frequency units, must be 0 < Wn < 1, with 1 corresponding to half the sample rate (Nyquist frequency)
fc = 0.5;

data = FILTsig(data,denoisetype,L,fc,true);

%---------------------------------------------------------------------------
% Upsampling
%---------------------------------------------------------------------------

% US-Upsampling to avoid edge-effects around Nyquist frequency (and ensuing artifacts in energy evaluation)
% US factor must be a power of 2 to be compliant with CWT dyadic algorithm - US=2^0=1 means no Upsampling
US = 2^addoct;

USn = US*n;
USdt = dt/US;
USmark = US*(mark-1)+1;
UStimeVec = ([0:1:US*n-1]')*USdt;
UStrend = UStimeVec*((y2-y1)/x2) + ones(USn,1)*y1;
if strcmp(filename,'test1')
	UStrend = zeros(USn,ntrace);
end

% If US=1 this does nothing
% Previous detrending avoids "resample" edge artifacts
% In case of UPsampling "resample" and "interp" are almost the same in MATLAB but not in Octave!
for j = 1:ntrace
	USdata(:,j) = resample(data(:,j),US,1);
end

% Bridge Retrending
if ~(strcmp(cone,'PAD'))
	USdata = USdata + UStrend;
end

%---------------------------------------------------------------------------
% Frequency Characteristic Parameters
%---------------------------------------------------------------------------

noctave = floor(log2(USn))-lft;
nscale = nvoice*noctave;

scaleVec = 2.^[1:noctave+1];
scaleVec = fliplr(scaleVec);

if (strcmp(unit,'ms'))
	freqVec = 1./(scaleVec*USdt); % Frequency in Hz - freqVec(end) = f.Nyquist = 1/(2*dt)
	periodVec = 1000*scaleVec*USdt; % Period in ms - periodVec(end) = 1/f.Nyquist = 2*dt
else
	freqVec = 1000./(scaleVec*USdt); % Frequency in mHz - freqVec(end) = f.Nyquist = 1/(2*dt)
	periodVec = scaleVec*USdt; % Period in s - periodVec(end) = 1/f.Nyquist = 2*dt
end

%---------------------------------------------------------------------------
% Print General Information
%---------------------------------------------------------------------------

fprintf('\nUpsampling Factor = %d\n',US);
fprintf('\nData Size:');
fprintf('\nNumber of Samples = %d -> brought to %d by Upsampling',n,USn);
fprintf('\nNumber of Traces  = %d\n',ntrace);
if (strcmp(unit,'ms'))
	fprintf('\nSampling Time (Mode) = %.6f ms -> brought to %.6f ms by Upsampling\n',1000*dt,1000*USdt);
else
	fprintf('\nSampling Time (Mode) = %.4f s -> brought to %.4f s by Upsampling\n',dt,USdt);
end
if (strcmp(unit,'ms'))
	fprintf('\nNyquist Frequency = %.2f Hz -> brought to %.2f Hz by Upsampling\n',1/(2*dt),1/(2*USdt));
else
	fprintf('\nNyquist Frequency = %.2f mHz -> brought to %.2f mHz by Upsampling\n',1000/(2*dt),1000/(2*USdt));
end
fprintf('\nLowest Frequency = %.2f mHz\n',1000/timeVec(end));
fprintf('\nOctaves Number = %d -> brought to %d by Upsampling\n',noctave-log2(US),noctave);
fprintf('\n\n');

%---------------------------------------------------------------------------
% Transform and Plot
%---------------------------------------------------------------------------

for j = 1:ntrace
	
	y = USdata(:,j);
	sig = y;
	
	%---------------------------------------------------------------------------
	% Padding with Zeros
	%---------------------------------------------------------------------------
	
	if (strcmp(cone,'PAD'))
		y = [y;zeros(USn,1)];
	end
	
	%---------------------------------------------------------------------------
	% Wavelet Transform
	%---------------------------------------------------------------------------
	
	[wt,fac] = trans(y,lft,nvoice,wavelet);
	
	% Cut the "zero zone": every point greater than USn and the lowest octave
	if (strcmp(cone,'PAD'))
		wt = wt(nvoice+1:end,1:USn);
	end
	
	wtabs = abs(wt);
	wtreal = abs(real(wt));
	
	%---------------------------------------------------------------------------
	% Cone of Influence (COI)
	%---------------------------------------------------------------------------
	
	coimask = coi(wtabs,nvoice,fac.coi,cone,strcmp(filename,'test3'));
	
	%---------------------------------------------------------------------------
	% Rescale to [0,1] and Convert to RGB
	%---------------------------------------------------------------------------
	
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
	% Apply Shadows (ultra-Nyquist zone and COI mask)
	%---------------------------------------------------------------------------
	
	% Gray scale scalograms
	gswtabs = ind2rgb(wtabsind,gray(128));
	gswtreal = ind2rgb(wtrealind,gray(128));
	
	% Shade the "ultra-Nyquist zone": the log2(US) highest octaves
	wtabs(end-log2(US)*nvoice+1:end,:,:) = gswtabs(end-log2(US)*nvoice+1:end,:,:);
	wtreal(end-log2(US)*nvoice+1:end,:,:) = gswtreal(end-log2(US)*nvoice+1:end,:,:);
	
	% Shade masked regions gray
	index = find(coimask == 0);
	
	% Gray scale - If length(index)==0 this does nothing
	wtabs(index) = gswtabs(index);
	wtabs(index + USn*nscale) = gswtabs(index);
	wtabs(index + 2*USn*nscale) = gswtabs(index);
	
	% Gray scale - If length(index)==0 this does nothing
	wtreal(index) = gswtreal(index);
	wtreal(index + USn*nscale) = gswtreal(index);
	wtreal(index + 2*USn*nscale) = gswtreal(index);
	
	%---------------------------------------------------------------------------
	% Plot Time Course
	%---------------------------------------------------------------------------
	
	figure
	
	subplot(3,2,1), hold on
		
		if (strcmp(cone,'PAD'))
			plot(UStimeVec,UStrend(:,j)-mean(yraw(:,j)),'-g','LineWidth',1)
			yraw(:,j) = yraw(:,j) - mean(yraw(:,j));
		end
		
		plot(timeVec,yraw(:,j),'-b','LineWidth',1)
		plot(UStimeVec,sig,'-r','LineWidth',1)		
		
		xlim([UStimeVec(1),UStimeVec(end)])
		set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(UStimeVec(end))]));
		ylabel('(nominal Volt)','FontSize',s2)
		title('NU CWT (W/sHz)^{1/2} - Modulus and Real Part','FontSize',s3)
		text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['ROI ',num2str(roi(j))],'FontSize',s2,'Color','k')
		
		for k = 1:length(mark) % If length(mark)==0 this does nothing
			plot([UStimeVec(USmark(k)),UStimeVec(USmark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
		end
	
	%---------------------------------------------------------------------------
	% Plot Modulus and Real Part
	%---------------------------------------------------------------------------
	
	subplot(3,2,3), hold on
		
		image(UStimeVec,1:nscale,wtabs)
		xlim([UStimeVec(1),UStimeVec(end)])
		ylim([0.5,nscale+0.5])
		set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(UStimeVec(end))]));
		set(gca,'YDir','normal','YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
		set(gca,'YTickLabel',num2str(freqVec',4));
		if (strcmp(unit,'ms'))
			ylabel('Frequency (Hz)','FontSize',s2)
		else
			ylabel('Frequency (mHz)','FontSize',s2)
		end
		
	subplot(3,2,5), hold on
		
		image(UStimeVec,1:nscale,wtreal)
		xlim([UStimeVec(1),UStimeVec(end)])
		ylim([0.5,nscale+0.5])
		set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(UStimeVec(end))]));
		set(gca,'YDir','normal','YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
		set(gca,'YTickLabel',num2str(freqVec',4));
		if (strcmp(unit,'ms'))
			ylabel('Frequency (Hz)','FontSize',s2)
		else
			ylabel('Frequency (mHz)','FontSize',s2)
		end
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
		magnVec = [minimo:(massimo-minimo)/10:massimo];
		set(gca,'YTickLabel',num2str(magnVec',2));
		xlabel('CWT','FontSize',s2)
		
	%---------------------------------------------------------------------------
	% Resize - Syntax template: set(gca,'Position',[left bottom width height])
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
par.x = UStimeVec;
par.y = freqVec;
par.z = periodVec;
par.u = unit;
par.m = USmark;
par.r = roi(end);
par.c = fac.adm;
par.US = US;


%%------------------------------------------------------------------------------------------------------%%
%%------------------------------------------------------------------------------------------------------%%
%%                                                                                                      %%
%% KYM Project                                                                                          %%
%% -----------                                                                                          %%
%% First Release in 2010                                                                                %%
%% Original code by Federico Alessandro Ruffinatti                                                      %%
%%                                                                                                      %%
%% UNIVERSITY OF TORINO                                                                                 %%
%% DOCTORAL SCHOOL IN LIFE AND HEALTH SCIENCES                                                          %%
%% Neuroscience PhD - Experimental Neuroscience                                                         %%
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