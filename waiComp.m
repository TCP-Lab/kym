function waiComp(x,coimask,timeVec,freqVec,nvoice,mark,flag)

%
%--------------------------------------------------------------------------------
% CWT Complete Wave-Activity Index - CompWAI
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% waiComp(x,coimask,timeVec,freqVec,nvoice,mark,flag)
%
% INPUT       TYPE        MEANING
% -----       ----        -------
% x        -> matrix   -> Continuous Wavelet Modulus
% coimask  -> matrix   -> COI Mask
% timeVec  -> array    -> Time Vector
% freqVec  -> array    -> Frequency Vector
% nvoice   -> scalar   -> Lines of Pixel per Octave
% mark     -> array    -> Time Marker Set (Step Number)
% flag     -> scalar   -> 0 = Thresholded WAI / 1 = NO Threshold
%
% OUTPUT      TYPE        MEANING
% ------      ----        -------
% -none-   -> plot     -> Plot Resulting from Analysis
%

nscale = size(x)(1);
n = size(x)(2);

x = x .* coimask;

WAI1=[];
WAI2=[];
WAI3=[];

% Frequency Matrix for the Complete WAI
freqVecComp = freqVec(end).*2.**(-([0:1:nscale-1])./nvoice);
freqVecComp = freqVecComp';
freqVecComp = flipud(freqVecComp);
freqVecComp = freqVecComp*ones(1,size(x)(2));

% Indexes
WAI1 = (10**5)*sum((x.**2),1)/(nvoice);
WAI2 = (10**4)*sum((log2(freqVecComp).*(x.**2)),1)/(nvoice);
WAI3 = (10**3)*sum((freqVecComp.*(x.**2)),1)/(nvoice);

% Mean Values
MEAN1(1) = (1./mark(1))*sum(WAI1(1:mark(1)));
if (length(mark) > 1)
	for h = 2:length(mark)
		MEAN1(h) = (1./(mark(h)-mark(h-1)+1))*sum(WAI1(mark(h-1):mark(h)));
	endfor
endif
MEAN1 = [MEAN1,(1./(length(timeVec)-mark(end)+1))*sum(WAI1(mark(end):end))];
d1 = MEAN1(2)/MEAN1(1);
d1 = floor(d1*100)/100; % In order to have only 2 decimal digits
%MEAN1 = floor(MEAN1*10)/10; % In order to have only 1 decimal digits
MEAN1 = floor(MEAN1*100)/100; % In order to have only 2 decimal digits

MEAN2(1) = (1./mark(1))*sum(WAI2(1:mark(1)));
if (length(mark) > 1)
	for h = 2:length(mark)
		MEAN2(h) = (1./(mark(h)-mark(h-1)+1))*sum(WAI2(mark(h-1):mark(h)));
	endfor
endif
MEAN2 = [MEAN2,(1./(length(timeVec)-mark(end)+1))*sum(WAI2(mark(end):end))];
d2 = MEAN2(2)/MEAN2(1);
d2 = floor(d2*100)/100; % In order to have only 2 decimal digits
%MEAN2 = floor(MEAN2*10)/10; % In order to have only 1 decimal digits
MEAN2 = floor(MEAN2*100)/100; % In order to have only 2 decimal digits

MEAN3(1) = (1./mark(1))*sum(WAI3(1:mark(1)));
if (length(mark) > 1)
	for h = 2:length(mark)
		MEAN3(h) = (1./(mark(h)-mark(h-1)+1))*sum(WAI3(mark(h-1):mark(h)));
	endfor
endif
MEAN3 = [MEAN3,(1./(length(timeVec)-mark(end)+1))*sum(WAI3(mark(end):end))];
d3 = MEAN3(2)/MEAN3(1);
d3 = floor(d3*100)/100; % In order to have only 2 decimal digits
%MEAN3 = floor(MEAN3*10)/10; % In order to have only 1 decimal digits
MEAN3 = floor(MEAN3*100)/100; % In order to have only 2 decimal digits

% Plot Indexes
subplot(3,2,1+flag), hold on
plot(timeVec,WAI1)
for k = 1:length(mark)
		plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
endfor
axis([timeVec(1),timeVec(end)])
set(gca,'XTick',[floor(timeVec(1)),get(gca,'XTick'),floor(timeVec(end))]);
if (flag == 0)
	ylabel('Index I','FontSize',18)
endif
text(timeVec(15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN1(1))],'FontSize',18,'Color','k')
if (length(mark) > 1)
	for h = 2:length(mark)
		text(timeVec(mark(h-1)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN1(h))],'FontSize',18,'Color','k')
	endfor
endif
text(timeVec(mark(end)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN1(end))],'FontSize',18,'Color','k')
text(max(xlim)+15,max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(d1)],'FontSize',18,'Color','r')

subplot(3,2,3+flag), hold on
plot(timeVec,WAI2)
for k = 1:length(mark)
		plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
endfor
axis([timeVec(1),timeVec(end)])
set(gca,'XTick',[floor(timeVec(1)),get(gca,'XTick'),floor(timeVec(end))]);
if (flag == 0)
	ylabel('Index Y','FontSize',18)
endif
text(timeVec(15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN2(1))],'FontSize',18,'Color','k')
if (length(mark) > 1)
	for h = 2:length(mark)
		text(timeVec(mark(h-1)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN2(h))],'FontSize',18,'Color','k')
	endfor
endif
text(timeVec(mark(end)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN2(end))],'FontSize',18,'Color','k')
text(max(xlim)+15,max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(d2)],'FontSize',18,'Color','r')

subplot(3,2,5+flag), hold on
plot(timeVec,WAI3)
for k = 1:length(mark)
		plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
endfor
axis([timeVec(1),timeVec(end)])
set(gca,'XTick',[floor(timeVec(1)),get(gca,'XTick'),floor(timeVec(end))]);
if (flag == 0)
	ylabel('Index J','FontSize',18)
endif
xlabel('Time (s)','FontSize',18)
text(timeVec(15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN3(1))],'FontSize',18,'Color','k')
if (length(mark) > 1)
	for h = 2:length(mark)
		text(timeVec(mark(h-1)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN3(h))],'FontSize',18,'Color','k')
	endfor
endif
text(timeVec(mark(end)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN3(end))],'FontSize',18,'Color','k')
text(max(xlim)+15,max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(d3)],'FontSize',18,'Color','r')

subplot(3,2,1+flag), hold on


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