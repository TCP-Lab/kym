function energy(x,timeVec,freqVec,nvoice,mark)

%
%---------------------------------------------------------------------------
% Wavelet Power Spectrum
%---------------------------------------------------------------------------
%
%
% Function Definition
%
% energy(x,timeVec,freqVec,nvoice,mark)
%
% INPUT       TYPE        MEANING
% -----       ----        -------
% x        -> matrix   -> Morlet Wavelet Modulus
% timeVec  -> array    -> Time Vector
% freqVec  -> array    -> Frequency Vector
% nvoice   -> scalar   -> Lines of Pixel per Octave
% mark     -> array    -> Time Marker Set (Step Number)
%
% OUTPUT      TYPE        MEANING
% ------      ----        -------
% -none-   -> plot     -> Plot Resulting from Analysis
%

nscale = size(x)(1);
n = size(x)(2);

PWR = (10**5)*sum(x.**2,1)/(nvoice);
NRG = (10**3)*sum(x.**2,2);

% Power Mean Values
MEAN(1) = (1./mark(1))*sum(PWR(1:mark(1)));
if (length(mark) > 1)
	for h = 2:length(mark)
		MEAN(h) = (1./(mark(h)-mark(h-1)+1))*sum(PWR(mark(h-1):mark(h)));
	endfor
endif
MEAN = [MEAN,(1./(length(timeVec)-mark(end)+1))*sum(PWR(mark(end):end))];
d = MEAN(2)/MEAN(1);
d = floor(d*100)/100; % In order to have only 2 decimal digits
MEAN = floor(MEAN*10)/10; % In order to have only 1 decimal digits

% Plotting
subplot(2,2,2), hold on
	plot(timeVec,PWR)
	for k = 1:length(mark)
			plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
	endfor
	axis([timeVec(1),timeVec(end)])
	set(gca,'XTick',[get(gca,'XTick'),timeVec(end)]);
	set(gca,'XTickLabel',[]);
	ylabel('Power P','FontSize',18)
	title('Wavelet Power Spectrum','FontSize',18)
	text(timeVec(15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN(1))],'FontSize',18,'Color','k')
	if (length(mark) > 1)
		for h = 2:length(mark)
			text(timeVec(mark(h-1)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN(h))],'FontSize',18,'Color','k')
		endfor
	endif
	text(timeVec(mark(end)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN(end))],'FontSize',18,'Color','k')
	text(min(xlim)-(max(xlim)-min(xlim))*(35/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(d)],'FontSize',18,'Color','r')

subplot(2,2,3), hold on
	plot(NRG,[1:nscale])
	axis([0,max(NRG),1,nscale])
	set(gca,'XDir','reverse');
	set(gca,'XTick',[0:floor(max(NRG))/4:max(NRG)]);
	set(gca,'XTick',[get(gca,'XTick'),floor(max(NRG))]);
	set(gca,'YTick',[0:nvoice:nscale]); % Maximum number of displayable ticks on ordinates
	set(gca,'YTick',[get(gca,'YTick'),1]); % In order to display freqVec(1) value
	set(gca,'YTickLabel',freqVec);
	ylabel('Frequency (mHz)','FontSize',18)
	xlabel('Energy Density E','FontSize',18)

subplot(2,2,4), hold on
	imagesc(timeVec,[],x)
	set(gca,'XTick',[get(gca,'XTick'),floor(timeVec(end))]);
	set(gca,'YDir','normal');
	set(gca,'YTick',[0:nvoice:nscale]); % Maximum number of displayable ticks on ordinates
	set(gca,'YTickLabel',[]);
	xlabel('Time (s)','FontSize',18)
	
% Resize Plot
subplot(2,2,2)
set(gca,'position',[0.35,0.7,0.62,0.2])

subplot(2,2,3)
set(gca,'position',[0.12,0.11,0.2,0.55])

subplot(2,2,4)
set(gca,'position',[0.35,0.11,0.62,0.55])


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