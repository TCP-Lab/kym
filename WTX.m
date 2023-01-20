function WTX(ratiofile,unit,mark,roi,lft,nvoice,wavelet,cone,index,singlefile,output)

%
%--------------------------------------------------------------------------------
% CWT Wave-Activity Index as a Function of Time and Space - ATX
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% WTX(ratiofile,unit,mark,roi,lft,nvoice,wavelet,cone,index,singlefile,output)
%
% INPUT          TYPE       EXAMPLE              MEANING
% -----          ----       -------              -------
% ratiofile   -> string  -> 'file_ratio.csv'  -> Calcium File (340/380 ratio)
% unit        -> string  -> 's'               -> Time Unit: 's' or 'ms'
% mark        -> array   -> [300,450,500]     -> Time Marker Set (Step Number)
% roi         -> array   -> [2,3,9]           -> ROI Set ([x:y] = from x to y)
% lft         -> scalar  -> 4                 -> Low-Frequency Threshold
% nvoice      -> scalar  -> 48                -> Lines of Pixel per Octave
% wavelet     -> string  -> 'M1'              -> Mother Wavelet Type
% cone        -> string  -> 'PAD'             -> Cone of Influence Handling Method
% index       -> string  -> 'COMP'            -> Activity Index Computation Method
% singlefile  -> string  -> 'file_380.csv'    -> Single Wavelength File (380 nm)
% output      -> boolean -> true              -> Print .eps Output Graph
%
% OUTPUT         TYPE                            MEANING
% ------         ----                            -------
% -none-      -> plot                         -> 5 Plots Resulting from Analysis
%

%---------------------------------------------------------------------------
% Graphic Parameters
%---------------------------------------------------------------------------

s1 = 16; % X-Y TickLabel size
s2 = 19; % X-Y Label and text size
s3 = 24; % Title size

%---------------------------------------------------------------------------
% Input Control
%---------------------------------------------------------------------------

% Print plots - Default value = false = Print nothing
if (nargin < 11)
	output = false;
end

% Single wavelength file (Cell thickness) - Default value = 'NONE'
if (nargin < 10)
	singlefile = 'NONE';
end

% Activity Index Type (Complete WAI or Maxima WAI) - Default value = 'COMP'
if (nargin < 9)
  index = 'COMP';
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

% Wavelet Version - Default value = 'M1'
if (nargin < 7)
  wavelet = 'M1';
end
if ~(strcmp(wavelet,'M1') | strcmp(wavelet,'M2') | strcmp(wavelet,'M3'))
	fprintf('\n\nWARNING: Invalid Mother Wavelet Type\n');
	fprintf('\n\n');
	return
end

% Voices Number - Default value = 24 voices per octave
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

% Marker Array - Default value = [] = No markers
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

data = dlmread(ratiofile);

% M-Downsampling
%M = 2;
%data = data([1:M:size(data,1)],:);

% Samples must be even
if (mod(size(data,1),2) ~= 0)
	data = data(1:size(data,1)-1,:);
end

% Extract Time Vector
timeVec = data(:,1);
timeVec = timeVec - timeVec(1); % Start from t=0s
if (strcmp(unit,'ms')) % Measured in s
	timeVec = timeVec/1000;
end
dt = mode(diff(timeVec)); % Sampling time mode

% Data Matrix
data(:,1) = [];
n = size(data,1);

% Open Single Wavelength File (Cell Thickness)
if ~(strcmp(singlefile,'NONE'))
	thick = dlmread(singlefile);
	thick(:,1) = [];
	if (size(thick,2) ~= size(data,2))
		fprintf('\n\nWARNING: ''Activity function'' and ''Surface to Volume profile'' have different domains\n');
		return
	end
end

% Select data of interest, keeping into account ROI 1 is the background
if (length(roi) == 0)
	roi = [2:size(data,2)+1];
end
data = data(:,roi(1:end)-1);
if ~(strcmp(singlefile,'NONE'))
	thick = thick(:,roi(1:end)-1);
end

%---------------------------------------------------------------------------
% Characteristic parameters
%---------------------------------------------------------------------------

noctave = floor(log2(n))-lft;
nscale = nvoice*noctave;
scaleVec = 2.^[1:noctave+1];
scaleVec = fliplr(scaleVec);
% Express frequency in mHz, with 1 decimal digit; freqVec(end) = f.Nyquist = 1/(2*dt)
freqVec = round((1000./(scaleVec*dt))*10)/10; % In order to get only 1 decimal digit
lowest = 1000/(2*dt*2^(length(freqVec)-1)); % To preserve an accurate value for freqVec(1)
% Express period in s, with 1 decimal digit - periodVec(end) = 1/f.Nyquist = 2*dt
periodVec = round((scaleVec*dt)*10)/10; % In order to get only 1 decimal digit

%---------------------------------------------------------------------------
% Transform and Plot
%---------------------------------------------------------------------------

for j = 1:size(data,2)
	
	% Original trace
	y = data(:,j);
	
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
		
		% Padding with n zeros
		y = [y;zeros(n,1)];
	end
	
	%---------------------------------------------------------------------------
	% Wavelet Transform
	%---------------------------------------------------------------------------
	
	[wt,fac] = trans(y,lft,nvoice,wavelet);
	
	if (strcmp(cone,'PAD')) % Cut the "zero zone": every point greater than n and the lowest octave
		wt = wt(nvoice+1:end,1:n);
	end
	
	wtabs = abs(wt);
	wtreal = abs(real(wt));
	
	%---------------------------------------------------------------------------
	% Cone of Influence (COI)
	%---------------------------------------------------------------------------
	
	coimask = coi(nvoice,nscale,n,dt,lowest,fac,cone);
	
	%--------------------------------------------------------------------------------
	% CWT WAI and ATX
	%--------------------------------------------------------------------------------
	
	% CWT Threshold Percentage for Maxima WAI and Complete WAI computation
	thr = 0;
	
	switch (index)
		
		case 'MAX' % Maxima WAI
			
			[maxima,maximaNT] = paths(wtabs,wtreal,coimask,thr,timeVec,freqVec,mark,5,false,false);
			
			[WAI,MEAN] = wai(wtabs,timeVec,freqVec,nvoice,mark,maxima,0,false);
			
			wtabs = wtabs .* coimask;
		
		case 'COMP' % Complete WAI
			
			if (thr ~= 0)
				wtabst = (wtabs-min(min(wtabs)))/max(max(wtabs-min(min(wtabs))));
				threshold = thr/100;
				wtabst = (wtabst >= threshold);
				wtabs = wtabs .* wtabst;
			end
			
			wtabs = wtabs .* coimask;
			
			[WAI,MEAN] = wai(wtabs,timeVec,freqVec,nvoice,mark,[],0,false);
		
	end
	
	% Choose the Index (1=I 2=Y 3=J) and assemble functions ATX and WAIX
	ATX(j,:) = WAI(3,:); % Activity Index as a Function of Time and Space
	WAIX(j,:) = MEAN(3,:); % Mean Activity Indexes as a Function of Space
	
	%--------------------------------------------------------------------------------
	% VfX - Frequency Spectrum as a Function of Space
	%--------------------------------------------------------------------------------
	
	if (length(mark) > 0)
		v(1,:) = 1000*(sum(wtabs(:,1:mark(1)),2)/(mark(1)))';
		if (length(mark) > 1)
			for h = 2:length(mark)
				v(h,:) = 1000*(sum(wtabs(:,mark(h-1):mark(h)),2)/(mark(h)-mark(h-1)+1))';
			end
		end
		v(length(mark)+1,:) = 1000*(sum(wtabs(:,mark(end):end),2)/(n-mark(end)+1))';
		
		vfirst(j,:) = v(1,:);
		vsecond(j,:) = v(2,:);
		vratio(j,:) = v(2,:)./v(1,:);
		vratio(j,find(isnan(vratio(j,:)) | isinf(vratio(j,:)))) = 0; % Prevent to NaN and Inf from being plotted
	end
	
	VfX(j,:) = 1000*(sum(wtabs,2)/n)';
	
end

%---------------------------------------------------------------------------
% Plot ATX (CWT WAI as a Function of Time and Space)
%---------------------------------------------------------------------------

figure

subplot(2,2,1), hold on
	
	for k = 1:size(data,2)
		hold on
		p = plot(timeVec,ATX(k,:),'-b','LineWidth',1);
		% Blue->Red color gradient
		set(p,'Color',[(k-1)/(size(data,2)-1),0,(size(data,2)-k)/(size(data,2)-1)])
	end
	
	xlim([timeVec(1),timeVec(end)])
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	ylabel('Activity Index','FontSize',s2)
	title('WAI as a Function of Time and Space - ATX','FontSize',s3)
	text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['ROI from ',num2str(roi(1)),' to ',num2str(roi(end))],'FontSize',s2,'Color','k')
	
	for k = 1:length(mark) % If length(mark)==0 this does nothing
		plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
	end
	
subplot(2,2,3), hold on
	
	imagesc(timeVec(1):timeVec(end),1:size(ATX,1),ATX)
	xlim([timeVec(1),timeVec(end)])
	ylim([0.5,size(ATX,1)+0.5])
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	
	if (size(data,2) >= 5)
		set(gca,'YDir','normal','YTick',unique([1:floor(size(data,2)/5):size(data,2),size(data,2)]));
		set(gca,'YTickLabel',num2str([roi(1:floor(size(data,2)/5):end),roi(end)]'));
	else
		set(gca,'YDir','normal','YTick',[1:1:size(data,2)]);
		set(gca,'YTickLabel',num2str(roi'));
	end
	
	ylabel('ROI number','FontSize',s2)
	xlabel('Time (s)','FontSize',s2)
	
% Add Color Bar
subplot(2,2,[2 4]), hold on
	
	t = [0:1:100];
	imagesc([1:2],[1:101],t') % X-support [1:2] prevents a "division by zero" in OCTAVE
	xlim([0.5,2.5])
	ylim([0.5,101.5])
	set(gca,'FontSize',s1,'XTick',[0]);
	set(gca,'XTickLabel','');
	set(gca,'YDir','normal','YTick',[1:10:101]);
	massimo = max(max(ATX));
	minimo = min(min(ATX));
	magnVec = round([minimo:(massimo-minimo)/10:massimo]*100)/100; % In order to get only 2 decimal digits
	set(gca,'YTickLabel',num2str(magnVec'));
	xlabel('Ratio','FontSize',s2)
	
	% Resize -> Syntax Template: set(gca,'Position',[left bottom width height])
	set(gca,'Position',[0.96,0.11,0.02,0.815])
	
	subplot(2,2,1)
	set(gca,'Position',[0.12,0.58384,0.77,0.34116])
	
	subplot(2,2,3)
	set(gca,'Position',[0.12,0.11,0.77,0.34116])
	
% Print output graph
if (output)
	print -depsc atxout.eps
end

%---------------------------------------------------------------------------
% Plot Mean WAIs and r as a Function of Space
%---------------------------------------------------------------------------

figure

if (length(mark) > 0)
	
	subplot(2,1,1), hold on
		
		for k = 1:length(mark)+1
			hold on
			p = plot([1:1:size(data,2)],WAIX(:,k),'-b','LineWidth',1);
			% Blue->Red color gradient
			set(p,'Color',[(k-1)/length(mark),0,(length(mark)+1-k)/length(mark)])
			lgn{k} = ['W_',num2str(k)]; % Cell Arrays
		end
		
		xlim([1,size(data,2)])
		
		if (size(data,2) >= 10)
			set(gca,'FontSize',s1,'XTick',unique([1:floor(size(data,2)/10):size(data,2),size(data,2)]));
			set(gca,'XTickLabel',num2str([roi(1:floor(size(data,2)/10):end),roi(end)]'));
		else
			set(gca,'FontSize',s1,'XTick',[1:1:size(data,2)]);
			set(gca,'XTickLabel',num2str(roi'));
		end
		
		ylabel('Mean WAI','FontSize',s2)
		mylegend = legend(lgn,'Location','North');
		set(mylegend,'FontSize',s1);
		title('Mean WAI and r as a Function of Space','FontSize',s3)
		
	subplot(2,1,2), hold on
		
		plot([1:1:size(data,2)],WAIX(:,2)./WAIX(:,1),'-b','LineWidth',1), hold on
		plot([1:1:size(data,2)],ones(1,size(data,2)),'-r','LineWidth',1)
		
		xlim([1,size(data,2)])
		
		if (size(data,2) >= 10)
			set(gca,'FontSize',s1,'XTick',unique([1:floor(size(data,2)/10):size(data,2),size(data,2)]));
			set(gca,'XTickLabel',num2str([roi(1:floor(size(data,2)/10):end),roi(end)]'));
		else
			set(gca,'FontSize',s1,'XTick',[1:1:size(data,2)]);
			set(gca,'XTickLabel',num2str(roi'));
		end
		
		xlabel('ROI number','FontSize',s2)
		ylabel('r','FontSize',s2)
		
else
	
	plot([1:1:size(data,2)],WAIX(:,1),'-b','LineWidth',1);
	
	xlim([1,size(data,2)])
	
	if (size(data,2) >= 10)
		set(gca,'FontSize',s1,'XTick',unique([1:floor(size(data,2)/10):size(data,2),size(data,2)]));
		set(gca,'XTickLabel',num2str([roi(1:floor(size(data,2)/10):end),roi(end)]'));
	else
		set(gca,'FontSize',s1,'XTick',[1:1:size(data,2)]);
		set(gca,'XTickLabel',num2str(roi'));
	end
	
	xlabel('ROI number','FontSize',s2)
	ylabel('Mean WAI','FontSize',s2)
	title('Mean WAI as a Function of Space','FontSize',s3)
	
end

% Print output graph
if (output)
	print -depsc rxout.eps
end

%---------------------------------------------------------------------------
% Plot Mean WAI vs. Surface to Volume Ratio
%---------------------------------------------------------------------------

if ~(strcmp(singlefile,'NONE'))
	
	thick = mean(thick,1);
	
	figure
	
	subplot(2,1,1), hold on
	
		plot([1:1:size(data,2)],thick,'-g','LineWidth',1);
		
		xlim([1,size(data,2)])
		
		if (size(data,2) >= 10)
			set(gca,'FontSize',s1,'XTick',unique([1:floor(size(data,2)/10):size(data,2),size(data,2)]));
			set(gca,'XTickLabel',num2str([roi(1:floor(size(data,2)/10):end),roi(end)]'));
		else
			set(gca,'FontSize',s1,'XTick',[1:1:size(data,2)]);
			set(gca,'XTickLabel',num2str(roi'));
		end
		
		ylabel('Cellular Thickness (AU)','FontSize',s2)
		title('Mean WAI vs. Surface to Volume Ratio','FontSize',s3)
	
	subplot(2,1,2), hold on
		
		stvr = thick.^(-1); % Surface to Volume (StV) ratio
		stvr = stvr/max(stvr); % Rescaling StV ratio between 0 and 1
		
		for k = 1:length(mark)+1
			WAIX(:,k) = WAIX(:,k)/max(WAIX(:,k)); % Rescaling Mean WAIs between 0 and 1
		end
		
		gro = corr(stvr',WAIX); % Global Correlation (Pearson product-moment correlation coefficient)
		
		% Local (Point By Point) Correlation
		%%l3=5;
		%%% Standard
		%%for h = 1:size(data,2)-l3+1
		%%	lro(h+(l3-1)/2,:) = corr((stvr(h:h+l3-1))',(WAIX(h:h+l3-1,:))');
		%%end
		%%% Edges
		%%for h = 1:(l3-1)/2
		%%	lro(h,:) = corr((stvr(1:((l3-1)/2)+h))',(WAIX(1:((l3-1)/2)+h,:))');
		%%	lro(size(data,2)-h+1,:) = corr((stvr(size(data,2)-(l3-1)/2-h+1:size(data,2)))',(WAIX(size(data,2)-(l3-1)/2-h+1:size(data,2),:))');
		%%end
		
		if (length(mark) > 0)
			
			plot([1:1:size(data,2)],stvr,'-g','LineWidth',1);
			for k = 1:length(mark)+1
				p = plot([1:1:size(data,2)],WAIX(:,k),'-b','LineWidth',1);
				col(k,:) = [(k-1)/length(mark),0,(length(mark)+1-k)/length(mark)]; % Color code
				set(p,'Color',col(k,:))
			end
			
			xlim([1,size(data,2)])
			
			if (size(data,2) >= 10)
				set(gca,'FontSize',s1,'XTick',unique([1:floor(size(data,2)/10):size(data,2),size(data,2)]));
				set(gca,'XTickLabel',num2str([roi(1:floor(size(data,2)/10):end),roi(end)]'));
			else
				set(gca,'FontSize',s1,'XTick',[1:1:size(data,2)]);
				set(gca,'XTickLabel',num2str(roi'));
			end
			
			xlabel('ROI number','FontSize',s2)
			ylabel('Mean WAI vs. StV Ratio (AU)','FontSize',s2)
			for k = 1:length(mark)+1
				text(max(xlim)*(50/100),max(ylim)-(max(ylim)-min(ylim))*((5+10*k)/100),['\rho_',num2str(k),' = ',num2str(gro(k))],'FontSize',18,'Color',col(k,:))
			end
			
		else
			
			plot([1:1:size(data,2)],stvr,'-g','LineWidth',1);
			plot([1:1:size(data,2)],WAIX(:,1),'-b','LineWidth',1);
			
			xlim([1,size(data,2)])
			
			if (size(data,2) >= 10)
				set(gca,'FontSize',s1,'XTick',unique([1:floor(size(data,2)/10):size(data,2),size(data,2)]));
				set(gca,'XTickLabel',num2str([roi(1:floor(size(data,2)/10):end),roi(end)]'));
			else
				set(gca,'FontSize',s1,'XTick',[1:1:size(data,2)]);
				set(gca,'XTickLabel',num2str(roi'));
			end
			
			xlabel('ROI number','FontSize',s2)
			ylabel('Mean WAI vs. StV Ratio (AU)','FontSize',s2)
			text(max(xlim)*(50/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['\rho = ',num2str(gro)],'FontSize',s2,'Color','k')
			
		end
	
	% Print output graph
	if (output)
		print -depsc stvout.eps
	end

end

%---------------------------------------------------------------------------
% Plot VfX
%---------------------------------------------------------------------------

figure

subplot(2,2,1), hold on
	
	for k = 1:size(data,2)
		hold on
		p = plot([1:nscale],VfX(k,:),'-b','LineWidth',1);
		% Blue->Red color gradient for superimposition
		set(p,'Color',[(k-1)/(size(data,2)-1),0,(size(data,2)-k)/(size(data,2)-1)])
	end
	
	xlim([1,nscale])
	set(gca,'FontSize',s1,'XTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
	set(gca,'XTickLabel',num2str(freqVec'));
	ylabel('10^3 V Amplitude','FontSize',s2)
	title('Frequency Spectrum as a Function of Space - VfX','FontSize',s3)
	text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['ROI from ',num2str(roi(1)),' to ',num2str(roi(end))],'FontSize',s2,'Color','k')
	
subplot(2,2,3), hold on
	
	imagesc([1:nscale],[1:size(data,2)],VfX)
	xlim([1,nscale])
	ylim([0.5,size(data,2)+0.5])
	set(gca,'FontSize',s1,'XTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
	set(gca,'XTickLabel',num2str(freqVec'));
	
	if (size(data,2) >= 5)
		set(gca,'YDir','normal','YTick',unique([1:floor(size(data,2)/5):size(data,2),size(data,2)]));
		set(gca,'YTickLabel',num2str([roi(1:floor(size(data,2)/5):end),roi(end)]'));
	else
		set(gca,'YDir','normal','YTick',[1:1:size(data,2)]);
		set(gca,'YTickLabel',num2str(roi'));
	end
	
	ylabel('ROI number','FontSize',s2)
	xlabel('Frequency (mHz)','FontSize',s2)
	
% Add Color Bar
subplot(2,2,[2 4]), hold on
	
	t = [0:1:100];
	imagesc([1:2],[1:101],t') % X-support [1:2] prevents a "division by zero" in OCTAVE
	xlim([0.5,2.5])
	ylim([0.5,101.5])
	set(gca,'FontSize',s1,'XTick',[0]);
	set(gca,'XTickLabel','');
	set(gca,'YDir','normal','YTick',[1:10:101]);
	massimo = max(max(VfX));
	minimo = min(min(VfX));
	magnVec = round([minimo:(massimo-minimo)/10:massimo]*100)/100; % In order to get only 2 decimal digits
	set(gca,'YTickLabel',num2str(magnVec'));
	xlabel('10^3 V Amp','FontSize',s2)
	
	% Resize -> Syntax Template: set(gca,'Position',[left bottom width height])
	set(gca,'Position',[0.96,0.11,0.02,0.815])
	
	subplot(2,2,1)
	set(gca,'Position',[0.12,0.58384,0.77,0.34116])
	
	subplot(2,2,3)
	set(gca,'Position',[0.12,0.11,0.77,0.34116])
	
% Print output graph
if (output)
	print -depsc vfxout.eps
end

if (length(mark) > 0)
	
	% Rescale to [0,1] and Convert to RGB
	
	massimo = max([max(max(vfirst)),max(max(vsecond))]);
	minimo = min([min(min(vfirst)),min(min(vsecond))]);
	
	vfirst = (vfirst-minimo);
	vsecond = (vsecond-minimo);
	numax = max([max(max(vfirst)),max(max(vsecond))]);
	vfirst = vfirst/numax;
	vsecond = vsecond/numax;
	
	[vfirst,b] = gray2ind(vfirst,128);
	[vsecond,b] = gray2ind(vsecond,128);
	
	vfirst = ind2rgb(vfirst,jet(128));
	vsecond = ind2rgb(vsecond,jet(128));
	
	figure
	
	subplot(3,2,1), hold on
		
		image([1:nscale],[1:size(data,2)],vfirst)
		xlim([1,nscale])
		ylim([0.5,size(data,2)+0.5])
		set(gca,'FontSize',s1,'XTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
		set(gca,'XTickLabel',num2str(freqVec'));
		
		if (size(data,2) >= 5)
			set(gca,'YDir','normal','YTick',unique([1:floor(size(data,2)/5):size(data,2),size(data,2)]));
			set(gca,'YTickLabel',num2str([roi(1:floor(size(data,2)/5):end),roi(end)]'));
		else
			set(gca,'YDir','normal','YTick',[1:1:size(data,2)]);
			set(gca,'YTickLabel',num2str(roi'));
		end
		
		ylabel('W_1 - ROIs','FontSize',s2)
		title('Second Window Spectra and Ratio as a Function of Space','FontSize',s3)
	
	subplot(3,2,3), hold on
		
		image([1:nscale],[1:size(data,2)],vsecond)
		xlim([1,nscale])
		ylim([0.5,size(data,2)+0.5])
		set(gca,'FontSize',s1,'XTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
		set(gca,'XTickLabel',num2str(freqVec'));
		
		if (size(data,2) >= 5)
			set(gca,'YDir','normal','YTick',unique([1:floor(size(data,2)/5):size(data,2),size(data,2)]));
			set(gca,'YTickLabel',num2str([roi(1:floor(size(data,2)/5):end),roi(end)]'));
		else
			set(gca,'YDir','normal','YTick',[1:1:size(data,2)]);
			set(gca,'YTickLabel',num2str(roi'));
		end
		
		ylabel('W_2 - ROIs','FontSize',s2)
		
	subplot(3,2,5), hold on
		
		imagesc([1:nscale],[1:size(data,2)],vratio)
		xlim([1,nscale])
		ylim([0.5,size(data,2)+0.5])
		set(gca,'FontSize',s1,'XTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
		set(gca,'XTickLabel',num2str(freqVec'));
		
		if (size(data,2) >= 5)
			set(gca,'YDir','normal','YTick',unique([1:floor(size(data,2)/5):size(data,2),size(data,2)]));
			set(gca,'YTickLabel',num2str([roi(1:floor(size(data,2)/5):end),roi(end)]'));
		else
			set(gca,'YDir','normal','YTick',[1:1:size(data,2)]);
			set(gca,'YTickLabel',num2str(roi'));
		end
		
		ylabel('W_2 / W_1 - ROIs','FontSize',s2)
		xlabel('Frequency (mHz)','FontSize',s2)
	
	% Add first color bar
	subplot(3,2,[2 4]), hold on
		
		t = [0:1:100];
		imagesc([1:2],[1:101],t') % X-support [1:2] prevents a "division by zero" in OCTAVE
		xlim([0.5,2.5])
		ylim([0.5,101.5])
		set(gca,'FontSize',s1,'XTick',[0]);
		set(gca,'XTickLabel','');
		set(gca,'YDir','normal','YTick',[1:10:101]);
		magnVec = round([minimo:(massimo-minimo)/10:massimo]*100)/100; % In order to get only 2 decimal digits
		set(gca,'YTickLabel',num2str(magnVec'));
		xlabel('10^3 V Amp','FontSize',s2)
		
		% Resize -> Syntax Template: set(gca,'Position',[left bottom width height])
		set(gca,'Position',[0.96,0.40963,0.02,0.517])
		
	% Add second color bar
	subplot(3,2,6), hold on
		
		t = [0:1:100];
		imagesc([1:2],[1:101],t') % X-support [1:2] prevents a "division by zero" in OCTAVE
		xlim([0.5,2.5])
		ylim([0.5,101.5])
		set(gca,'FontSize',s1,'XTick',[0]);
		set(gca,'XTickLabel','');
		set(gca,'YDir','normal','YTick',[1:20:101]);
		massimo = max(max(vratio));
		minimo = min(min(vratio));
		magnVec = round([minimo:(massimo-minimo)/5:massimo]*100)/100; % In order to get only 2 decimal digits
		set(gca,'YTickLabel',num2str(magnVec'));
		xlabel('fold','FontSize',s2)
		
		% Resize -> Syntax Template: set(gca,'Position',[left bottom width height])
		set(gca,'Position',[0.96,0.11,0.02,0.216])
		
		subplot(3,2,1)
		set(gca,'Position',[0.12,0.70926,0.77,0.21574])
		
		subplot(3,2,3)
		set(gca,'Position',[0.12,0.40963,0.77,0.21574])
		
		subplot(3,2,5)
		set(gca,'Position',[0.12,0.11,0.77,0.21574])
		
	% Print output graph
	if (output)
		print -depsc vratioout.eps
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