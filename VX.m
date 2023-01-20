function VX(filename,unit,mark,roi)

%
%---------------------------------------------------------------------------
% Traces Visualization
%---------------------------------------------------------------------------
%
%
% Function Definition
%
% VX(filename,unit,mark,roi)
%
% INPUT       TYPE      EXAMPLE           MEANING
% -----       ----      -------           -------
% filename -> string -> 'filename.csv' -> File Name
% unit     -> string -> 's'            -> Time Unit: 's' or 'ms'
% mark     -> array  -> [300,450,500]  -> Time Marker Set (Step Number)
% roi      -> array  -> [2,9]          -> ROI Interval
%
% OUTPUT      TYPE                        MEANING
% ------      ----                        -------
% -none-   -> plot                     -> 1 Plot Sampling Quality Control
% -none-   -> plot                     -> N Plot Showing the Traces
%

% Set the accepted time error in seconds for each sample
% This variable checks that the time steps are evenly spaced
timeError = 0.25;

% Marker Array control - Default value = No markers
if (nargin < 3)
	mark = [];
endif

% Time Unit of Measure control - Default value = 's'
if (nargin < 2)
	unit = 's';
endif

% Open file
A = dlmread(filename);

% Extract time vector - measured in s - starting from t = 0s
timeVec = A(:,1);
timeVec = timeVec - timeVec(1);
if (strcmp(unit,'ms'))
	timeVec = timeVec/1000;
endif
dt = mean(diff(timeVec));

% Check that time steps are equal (evenly spaced)
% timeError is the accepted time difference
if !(isempty(find(abs(diff(timeVec) - dt) > timeError)))
	printf("\nWARNING!\nNot equal time steps\n");
endif
printf("\n(Mean) Sampling Time = %f s\n\n", dt);

% Sampling quality control
figure
plot(diff(timeVec),'--b','LineWidth',1), hold on

for j = 1:length(mark)
		plot([mark(j),mark(j)],[min(ylim),max(ylim)],'--r','LineWidth',1)
endfor

axis([1,length(timeVec)])
set(gca,'XTick',[get(gca,'XTick'),length(timeVec)]);

xlabel('Time Step')
ylabel('dt')
text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),["Sampling Quality Control"],'FontSize',14,'Color','r')

% Data Matrix
A(:,1) = [];

% ROIs control - Default value = All traces
% Keeping into account ROI 1 is the background
if (nargin < 4)
	roi = [2,(size(A))(2)+1];
endif

% Select data of interest
% Keeping into account ROI 1 is the background
A = A(:,roi(1)-1:roi(2)-1);

% Check user input
if ((size(A))(2) < 1)
	printf("\nWARNING!\nThis file appears not to have any calcium recordings\n\n");
	return; %Function breaks
endif

% Plot original traces calcium recording
for i = 1:(size(A))(2)
	
	y = A(:,i);
	
	if (mod(i,3)==1)
		figure
		count = 1;
	endif
	
	subplot(3,1,count), hold on
	
	plot(timeVec,y,'--b','LineWidth',1)
	
	for j = 1:length(mark)
		plot([timeVec(mark(j)),timeVec(mark(j))],[min(ylim),max(ylim)],'--r','LineWidth',1)
	endfor
	
	axis([timeVec(1),timeVec(end)])
	set(gca,'XTick',[get(gca,'XTick'),timeVec(end)]);
	
	if (count==3)
		xlabel('Time (s)')
	endif
	
	ylabel('Ratio')
	text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),[" ROI ", num2str(roi(1)-1+i)],'FontSize',14,'Color','r')
	
	count = count+1;
	
endfor


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