function ratio = vector(x,freqVec,periodVec,mark)

%
%---------------------------------------------------------------------------
% Morlet Wavelet Transform Vector Approach
%---------------------------------------------------------------------------
%
%
% Function Definition
%
% ratio = vector(x,freqVec,periodVec,mark)
%
% INPUT        TYPE        MEANING
% -----        ----        -------
% x         -> matrix   -> Morlet Wavelet Modulus
% freqVec   -> array    -> Frequency Vector
% periodVec -> array    -> Period Vector
% mark      -> array    -> Time Marker Set (Step Number)
%
% OUTPUT       TYPE        MEANING
% ------       ----        -------
% ratio     -> array    -> post/pre Ratio Vector
% -none-    -> plot     -> Plot Resulting from Analysis
%

n = size(x)(2);
nscale = size(x)(1);
nvoice = (size(x)(1))/(length(freqVec)-1);

pre = sum(x(:,1:mark(1)),2);
pre = 1000*pre/(mark(1));

if (length(mark) > 1)
	
	mid = sum(x(:,mark(1):mark(2)),2);
	mid = 1000*mid/(mark(2)-mark(1)+1);
	
	dist1 = sqrt(sum(abs(mid-pre).**2)/(nvoice));
	normpre = sqrt(sum(abs(pre).**2)/(nvoice));
	normmid = sqrt(sum(abs(mid).**2)/(nvoice));
	delta1 = (normmid-normpre);
	theta1 = acos((sum(mid.*pre)/(nvoice))/(normmid*normpre));
	theta1 = (theta1/(pi))*180; % Radians -> Degrees conversion
	
	dist1 = floor(dist1*100)/100; % In order to have only 2 decimal digits
	delta1 = floor(delta1*100)/100; % In order to have only 2 decimal digits
	theta1 = floor(theta1*100)/100; % In order to have only 2 decimal digits

endif

post = sum(x(:,mark(end):end),2);
post = 1000*post/(n-mark(end)+1);

if (length(mark) > 1)
	ratio = mid./pre;
	ratio = ratio(:)';
elseif (length(mark) == 1)
	ratio = post./pre;
	ratio = ratio(:)';
endif

dist2 = sqrt(sum(abs(post-pre).**2)/(nvoice));
normpre = sqrt(sum(abs(pre).**2)/(nvoice));
normpost = sqrt(sum(abs(post).**2)/(nvoice));
delta2 = (normpost-normpre);
theta2 = acos((sum(post.*pre)/(nvoice))/(normpost*normpre));
theta2 = (theta2/(pi))*180; % Radians -> Degrees conversion

dist2 = floor(dist2*100)/100; % In order to have only 2 decimal digits
delta2 = floor(delta2*100)/100; % In order to have only 2 decimal digits
theta2 = floor(theta2*100)/100; % In order to have only 2 decimal digits

if (length(mark) > 1)
	
	dist3 = sqrt(sum(abs(post-mid).**2)/(nvoice));
	delta3 = (normpost-normmid);
	theta3 = acos((sum(post.*mid)/(nvoice))/(normpost*normmid));
	theta3 = (theta3/(pi))*180; % Radians -> Degrees conversion
	
	dist3 = floor(dist3*100)/100; % In order to have only 2 decimal digits
	delta3 = floor(delta3*100)/100; % In order to have only 2 decimal digits
	theta3 = floor(theta3*100)/100; % In order to have only 2 decimal digits

endif

plot(pre,'-b'), hold on
if (length(mark) > 1)
	plot(mid,'-g')
endif
plot(post,'-r')
axis([0,nscale])
set(gca,'XTick',[0:nvoice:nscale]);
set(gca,'XTickLabel',freqVec);
xlabel('Frequency (mHz)','FontSize',18)
ylabel('W Amplitude x 1000','FontSize',18)
set(gca,'position',[0.130,0.110,0.775,0.68])

if (length(mark) > 1)
	legend('pre','mid','post','Location','East')
elseif (length(mark) == 1)
	legend('pre','post','Location','East')
endif

text(max(xlim)*(80/100),max(ylim)-(max(ylim)-min(ylim))*(5/100),['pre vs. post'],'FontSize',16,'Color','r')
text(max(xlim)*(80/100),max(ylim)-(max(ylim)-min(ylim))*(10/100),['D = ', num2str(dist2)],'FontSize',16,'Color','k')
text(max(xlim)*(80/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['Delta = ', num2str(delta2)],'FontSize',16,'Color','k')
text(max(xlim)*(80/100),max(ylim)-(max(ylim)-min(ylim))*(20/100),['Theta = ', num2str(theta2), '°'],'FontSize',16,'Color','k')

if (length(mark) > 1)
	
	text(max(xlim)*(60/100),max(ylim)-(max(ylim)-min(ylim))*(5/100),['mid vs. post'],'FontSize',16,'Color','r')
	text(max(xlim)*(60/100),max(ylim)-(max(ylim)-min(ylim))*(10/100),['D = ', num2str(dist3)],'FontSize',16,'Color','k')
	text(max(xlim)*(60/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['Delta = ', num2str(delta3)],'FontSize',16,'Color','k')
	text(max(xlim)*(60/100),max(ylim)-(max(ylim)-min(ylim))*(20/100),['Theta = ', num2str(theta3), '°'],'FontSize',16,'Color','k')
	
	text(max(xlim)*(40/100),max(ylim)-(max(ylim)-min(ylim))*(5/100),['pre vs. mid'],'FontSize',16,'Color','r')
	text(max(xlim)*(40/100),max(ylim)-(max(ylim)-min(ylim))*(10/100),['D = ', num2str(dist1)],'FontSize',16,'Color','k')
	text(max(xlim)*(40/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['Delta = ', num2str(delta1)],'FontSize',16,'Color','k')
	text(max(xlim)*(40/100),max(ylim)-(max(ylim)-min(ylim))*(20/100),['Theta = ', num2str(theta1), '°'],'FontSize',16,'Color','k')
	
endif

axes('XAxisLocation','top');
set(gca,'YTick',[]);
axis([0,nscale])
set(gca,'XTick',[0:nvoice:nscale]);
set(gca,'XTickLabel',periodVec);
xlabel('Period (s)','FontSize',18)
set(gca,'position',[0.130,0.110,0.775,0.68])


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