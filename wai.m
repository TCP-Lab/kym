function wai(x,timeVec,lowest,nvoice,mark,maxima,l,flag)

%
%---------------------------------------------------------------------------
% Morlet Wavelet Transform Wave-Activity Index (WAI)
%---------------------------------------------------------------------------
%
%
% Function Definition
%
% wai(x,timeVec,lowest,nvoice,mark,maxima,l)
%
% INPUT       TYPE        MEANING
% -----       ----        -------
% x        -> matrix   -> Morlet Wavelet Modulus
% timeVec  -> array    -> Time Vector
% lowest   -> scalar   -> Lowest Frequency Taken into Account
% nvoice   -> scalar   -> Lines of Pixel per Octave
% mark     -> array    -> Time Marker Set (Step Number)
% maxima   -> matrix   -> 3rd-4th paths Output - Maxima Coordinate
% l        -> scalar   -> Mean Filter Window Width
% flag     -> scalar   -> 0 = Thresholded WAI / 1 = NO Threshold
%
% OUTPUT      TYPE        MEANING
% ------      ----        -------
% -none-   -> plot     -> Plot Resulting from Analysis
%

% l = Regularization Window Width = Mean Convolutive Filter
% l value must be odd! - l=1 means NO FILTER

n = size(x)(2);

WAI1=[];
WAI2=[];
WAI3=[];
WAI4=[];

% Input Control
if (mark(1) <= (l-1)/2)
	printf("\nWARNING!\nFirst mark before Regularization-Window-Half-Width\n\n");
	return
endif

% Indexes
for h = 1:n
	m = find(maxima(:,h));
	WAI1(h) = sum(x(m,h).**2);
	WAI2(h) = sum(log2(lowest*2.**(m/nvoice)).*(x(m,h).**2));
	WAI3(h) = sum((lowest*2.**(m/nvoice)).*(x(m,h).**2));
	WAI4(h) = sum((lowest*2.**(m/nvoice)).*(x(m,h).**2)) / sum((lowest*2.**(m/nvoice)));
endfor

% Convolution and Scale Factor
for h = 1:n-l+1
	WAI1(h) = (1./l)*sum(WAI1(h:h+l-1))*(10**5);
	WAI2(h) = (1./l)*sum(WAI2(h:h+l-1))*(10**4);
	WAI3(h) = (1./l)*sum(WAI3(h:h+l-1))*(10**3);
	WAI4(h) = (1./l)*sum(WAI4(h:h+l-1))*(10**6);
endfor

WAI1 = WAI1(1:n-l+1);
WAI2 = WAI2(1:n-l+1);
WAI3 = WAI3(1:n-l+1);
WAI4 = WAI4(1:n-l+1);

% Mean Values
c{1} = WAI1(1:mark(1)-(l-1)/2); % Cell Arrays
MEAN1(1) = (1./(mark(1)-(l-1)/2))*sum(c{1});
if (length(mark) > 1)
	for h = 2:length(mark)
		c{h} = WAI1(mark(h-1)-(l-1)/2:mark(h)-(l-1)/2);
		MEAN1(h) = (1./(mark(h)-mark(h-1)+1))*sum(c{h});
	endfor
endif
c{length(c)+1} = WAI1(mark(end)-(l-1)/2:end);
MEAN1 = [MEAN1,(1./(length(timeVec)-(l-1)-(mark(end)-(l-1)/2)+1))*sum(c{end})];
if (length(mark) < 3)
	[pval,t,dof] = t_test_2(c{1},c{2});
	d1 = MEAN1(2)/MEAN1(1);
	d1 = floor(d1*100)/100; % In order to have only 2 decimal digits
	printf(["\nr-I = ",num2str(d1)]);
	%printf(["\nDelta-I Significance Level = ",num2str(pval*100)]);
	printf(["\n"]);
endif
MEAN1 = floor(MEAN1*10)/10; % In order to have only 1 decimal digits

c{1} = WAI2(1:mark(1)-(l-1)/2); % Cell Arrays
MEAN2(1) = (1./(mark(1)-(l-1)/2))*sum(c{1});
if (length(mark) > 1)
	for h = 2:length(mark)
		c{h} = WAI2(mark(h-1)-(l-1)/2:mark(h)-(l-1)/2);
		MEAN2(h) = (1./(mark(h)-mark(h-1)+1))*sum(c{h});
	endfor
endif
c{length(c)+1} = WAI2(mark(end)-(l-1)/2:end);
MEAN2 = [MEAN2,(1./(length(timeVec)-(l-1)-(mark(end)-(l-1)/2)+1))*sum(c{end})];
if (length(mark) < 3)
	[pval,t,dof] = t_test_2(c{1},c{2});
	d2 = MEAN2(2)/MEAN2(1);
	d2 = floor(d2*100)/100; % In order to have only 2 decimal digits
	printf(["\nr-Y = ",num2str(d2)]);
	%printf(["\nDelta-Y Significance Level = ",num2str(pval*100)]);
	printf(["\n"]);
endif
MEAN2 = floor(MEAN2*10)/10; % In order to have only 1 decimal digits

c{1} = WAI3(1:mark(1)-(l-1)/2); % Cell Arrays
MEAN3(1) = (1./(mark(1)-(l-1)/2))*sum(c{1});
if (length(mark) > 1)
	for h = 2:length(mark)
		c{h} = WAI3(mark(h-1)-(l-1)/2:mark(h)-(l-1)/2);
		MEAN3(h) = (1./(mark(h)-mark(h-1)+1))*sum(c{h});
	endfor
endif
c{length(c)+1} = WAI3(mark(end)-(l-1)/2:end);
MEAN3 = [MEAN3,(1./(length(timeVec)-(l-1)-(mark(end)-(l-1)/2)+1))*sum(c{end})];
if (length(mark) < 3)
	[pval,t,dof] = t_test_2(c{1},c{2});
	d3 = MEAN3(2)/MEAN3(1);
	d3 = floor(d3*100)/100; % In order to have only 2 decimal digits
	printf(["\nr-J = ",num2str(d3)]);
	%printf(["\nDelta-J Significance Level = ",num2str(pval*100)]);
	printf(["\n"]);
endif
MEAN3 = floor(MEAN3*10)/10; % In order to have only 1 decimal digits

c{1} = WAI4(1:mark(1)-(l-1)/2); % Cell Arrays
MEAN4(1) = (1./(mark(1)-(l-1)/2))*sum(c{1});
if (length(mark) > 1)
	for h = 2:length(mark)
		c{h} = WAI4(mark(h-1)-(l-1)/2:mark(h)-(l-1)/2);
		MEAN4(h) = (1./(mark(h)-mark(h-1)+1))*sum(c{h});
	endfor
endif
c{length(c)+1} = WAI4(mark(end)-(l-1)/2:end);
MEAN4 = [MEAN4,(1./(length(timeVec)-(l-1)-(mark(end)-(l-1)/2)+1))*sum(c{end})];
if (length(mark) < 3)
	[pval,t,dof] = t_test_2(c{1},c{2});
	d4 = MEAN4(2)/MEAN4(1);
	d4 = floor(d4*100)/100; % In order to have only 2 decimal digits
	printf(["\nr-K = ",num2str(d4)]);
	%printf(["\nDelta-K Significance Level = ",num2str(pval*100)]);
	printf(["\n\n"]);
endif
MEAN4 = floor(MEAN4*10)/10; % In order to have only 1 decimal digits

% Plot Indexes
timeVec2 = timeVec(((l-1)/2)+1:end-((l-1)/2));

subplot(4,2,1+flag), hold on
plot(timeVec2,WAI1)
for k = 1:length(mark)
		plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'--r','LineWidth',1)
endfor
axis([timeVec2(1),timeVec2(end)])
set(gca,'XTick',[floor(timeVec2(2)),get(gca,'XTick'),floor(timeVec2(end))]); % floor(timeVec2(1)) results sometimes to small to be displayed
if (flag == 0)
	ylabel('Index I','FontSize',18)
endif
text(timeVec((l-1)/2+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN1(1))],'FontSize',18,'Color','k')
if (length(mark) > 1)
	for h = 2:length(mark)
		text(timeVec(mark(h-1)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN1(h))],'FontSize',18,'Color','k')
	endfor
endif
text(timeVec(mark(end)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN1(end))],'FontSize',18,'Color','k')
text(max(xlim)+15,max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(d1)],'FontSize',18,'Color','r')

subplot(4,2,3+flag), hold on
plot(timeVec2,WAI2)
for k = 1:length(mark)
		plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'--r','LineWidth',1)
endfor
axis([timeVec2(1),timeVec2(end)])
set(gca,'XTick',[floor(timeVec2(2)),get(gca,'XTick'),floor(timeVec2(end))]);
if (flag == 0)
	ylabel('Index Y','FontSize',18)
endif
text(timeVec((l-1)/2+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN2(1))],'FontSize',18,'Color','k')
if (length(mark) > 1)
	for h = 2:length(mark)
		text(timeVec(mark(h-1)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN2(h))],'FontSize',18,'Color','k')
	endfor
endif
text(timeVec(mark(end)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN2(end))],'FontSize',18,'Color','k')
text(max(xlim)+15,max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(d2)],'FontSize',18,'Color','r')

subplot(4,2,5+flag), hold on
plot(timeVec2,WAI3)
for k = 1:length(mark)
		plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'--r','LineWidth',1)
endfor
axis([timeVec2(1),timeVec2(end)])
set(gca,'XTick',[floor(timeVec2(2)),get(gca,'XTick'),floor(timeVec2(end))]);
if (flag == 0)
	ylabel('Index J','FontSize',18)
endif
text(timeVec((l-1)/2+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN3(1))],'FontSize',18,'Color','k')
if (length(mark) > 1)
	for h = 2:length(mark)
		text(timeVec(mark(h-1)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN3(h))],'FontSize',18,'Color','k')
	endfor
endif
text(timeVec(mark(end)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN3(end))],'FontSize',18,'Color','k')
text(max(xlim)+15,max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(d3)],'FontSize',18,'Color','r')

subplot(4,2,7+flag), hold on
plot(timeVec2,WAI4)
for k = 1:length(mark)
		plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'--r','LineWidth',1)
endfor
axis([timeVec2(1),timeVec2(end)])
set(gca,'XTick',[floor(timeVec2(2)),get(gca,'XTick'),floor(timeVec2(end))]);
if (flag == 0)
	ylabel('Index K','FontSize',18)
endif
xlabel('Time (s)','FontSize',18)	
text(timeVec((l-1)/2+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN4(1))],'FontSize',18,'Color','k')
if (length(mark) > 1)
	for h = 2:length(mark)
		text(timeVec(mark(h-1)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN4(h))],'FontSize',18,'Color','k')
	endfor
endif
text(timeVec(mark(end)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN4(end))],'FontSize',18,'Color','k')
text(max(xlim)+15,max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(d4)],'FontSize',18,'Color','r')

subplot(4,2,1+flag), hold on


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