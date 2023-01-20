function VX(sub,filename,unit,mark,roi,output)

%
%--------------------------------------------------------------------------------
% Traces Visualization
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% VX(sub,filename,unit,mark,roi,output)
%
% INPUT       TYPE       EXAMPLE           MEANING
% -----       ----       -------           -------
% sub      -> scalar  -> 3              -> Traces per Plot (Integer or 'a' or 'b')
% filename -> string  -> 'filename.csv' -> File Name
% unit     -> string  -> 's'            -> Time Unit: 's' or 'ms'
% mark     -> array   -> [300,450,500]  -> Time Marker Set (Step Number)
% roi      -> array   -> [2,3,9]        -> ROI Set ([x:y] = from x to y)
% output   -> boolean -> true           -> Print .eps Output Graphs
%
% OUTPUT      TYPE                        MEANING
% ------      ----                        -------
% -none-   -> plot                     -> 1 Plot Sampling Quality Control
% -none-   -> plot                     -> N Plots Showing the Traces
%

% Accepted time error in seconds for each sample
% in order to check that the time steps are evenly spaced
timeError = 0.5;

% Output Print control - Default value = false
if (nargin < 6)
	output = false;
endif

% Marker Array control - Default value = No markers
if (nargin < 4)
	mark = [];
endif

% Time Unit of Measure control - Default value = 's'
if (nargin < 3)
	unit = 's';
endif

% Open file
data = dlmread(filename);

% Extract Time Vector - measured in s - starting from t = 0s
timeVec = data(:,1);
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
plot(diff(timeVec),'-b','LineWidth',1), hold on

for k = 1:length(mark)
		plot([mark(k),mark(k)],[min(ylim),max(ylim)],'-r','LineWidth',1)
endfor

axis([1,length(timeVec)])
set(gca,'XTick',[1,get(gca,'XTick'),length(timeVec)]);
xlabel('Time Step','FontSize',18)
ylabel('dt (s)','FontSize',18)
title('Sampling Quality Control','FontSize',18)
text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),[" (Mean) Sampling Time = ", num2str(dt)],'FontSize',14,'Color','k')

% Print other generic information
printf("\nNyquist Frequency = %f mHz\n\n", 1000/(2*dt));
printf("\nTime Vector Length = %f\n\n", length(timeVec));
printf("\nLowest Frequency = %f mHz\n\n", 1000/timeVec(end));
printf("\nOctaves Number = %f\n\n", floor(log2(length(timeVec))));

% Data Matrix
data(:,1) = [];

% ROIs control - Default value = All traces
% Keeping into account ROI 1 is the background
if (nargin < 5)
	roi = [2:(size(data))(2)+1];
elseif (length(roi) == 0)
	roi = [2:(size(data))(2)+1];
endif

% Select data of interest
% keeping into account ROI 1 is the background
data = data(:,roi(1:end)-1);

% Check user input
if ((size(data))(2) < 1)
	printf("\nWARNING!\nThis file appears not to have any calcium recordings\n\n");
	return; %Function breaks
endif

% Plot original traces calcium recording
if (ischar(sub) == 0)

	if (sub == 1)
	
		for k = 1:(size(data))(2)
			
			y = data(:,k);
			
			figure
			plot(timeVec,y,'-b','LineWidth',1), hold on
			
			for l = 1:length(mark)
				plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
			endfor
			
			axis([timeVec(1),timeVec(end)])
			set(gca,'XTick',[get(gca,'XTick'),floor(timeVec(end))]);
			xlabel('Time (s)','FontSize',18)
			ylabel('Ratio','FontSize',18)
			title('Trace Visualization','FontSize',18)
			
			text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),[" ROI ", num2str(roi(k))],'FontSize',18,'Color','k')
			
			% Print output graph
			if (output)
				print(['traceout',num2str(k),'.eps'],'-depsc') %Progressive enumeration of output files
			endif
		
		endfor
		
	else
	
		for k = 1:(size(data))(2)
		
			y = data(:,k);
			
			if (mod(k,sub) == 1)
				figure
				count = 1;
			endif
			
			subplot(sub,1,count), hold on
			plot(timeVec,y,'-b','LineWidth',1)
			
			for l = 1:length(mark)
				plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
			endfor
			
			axis([timeVec(1),timeVec(end)])
			set(gca,'XTick',[get(gca,'XTick'),floor(timeVec(end))]);
			ylabel('Ratio','FontSize',18)
			text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),[" ROI ", num2str(roi(k))],'FontSize',18,'Color','k')
			
			if (count == sub | k == (size(data))(2))
				xlabel('Time (s)','FontSize',18)
				subplot(sub,1,1)
				title('Traces Visualization','FontSize',18)
				
				% Print output graph
				if (output)
					print(['traceout',num2str(floor(k/sub)+(mod(k,sub)!=0)),'.eps'],'-depsc') %Progressive enumeration of output files
				endif
			endif
			
			count = count+1;
			
		endfor
		
	endif

else

	switch (sub)
	
		case 'a'
		
			figure
			
			for k = 1:(size(data))(2)
			
				y = data(:,k);
				
				hold on
				p = plot(timeVec,y,'-b','LineWidth',1);
				
				% Blue->Red color gradient in case of superimposition (ischar(sub) == 1)
				set(p,'Color',[(k-1)/((size(data))(2)-1),0,((size(data))(2)-k)/((size(data))(2)-1)])
				
			endfor
			
			for l = 1:length(mark)
				plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
			endfor
			
			axis([timeVec(1),timeVec(end)])
			set(gca,'XTick',[get(gca,'XTick'),floor(timeVec(end))]);
			xlabel('Time (s)','FontSize',18)
			ylabel('Ratio','FontSize',18)
			title('Traces Visualization','FontSize',18)
			
			text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),[" ROI from ", num2str(roi(1)), " to " , num2str(roi(end))],...
			'FontSize',18,'Color','k')
			
			% Print output graph
			if (output)
				print -depsc traceout.eps
			endif
			
		case 'b'
		
			figure
			subplot(2,2,1), hold on
			
			for k = 1:(size(data))(2)
			
				y = data(:,k);
				
				hold on
				p = plot(timeVec,y,'-b','LineWidth',1);
				
				% Blue->Red color gradient in case of superimposition (ischar(sub) == 1)
				set(p,'Color',[(k-1)/((size(data))(2)-1),0,((size(data))(2)-k)/((size(data))(2)-1)])
				
			endfor
			
			for l = 1:length(mark)
				plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
			endfor
			
			axis([timeVec(1),timeVec(end)])
			set(gca,'XTick',[get(gca,'XTick'),floor(timeVec(end))]);
			ylabel('Ratio','FontSize',18)
			title('Traces Visualization','FontSize',18)
			
			text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),[" ROI from ", num2str(roi(1)), " to " , num2str(roi(end))],...
			'FontSize',18,'Color','k')
			
			subplot(2,2,3), hold on
			
			imagesc(timeVec,[],data')
			set(gca,'XTick',[get(gca,'XTick'),floor(timeVec(end))]);
			%set(gca,'YDir','normal');
			if ((size(data))(2) >= 5)
				set(gca,'YTick',[1:floor((size(data))(2)/5):(size(data))(2),(size(data))(2)]);
				set(gca,'YTickLabel',[roi(1:floor((size(data))(2)/5):end),roi(end)]);
			else
				set(gca,'YTick',[1:1:(size(data))(2)]);
				set(gca,'YTickLabel',roi);
			endif
			ylabel('ROI number','FontSize',18)
			xlabel('Time (s)','FontSize',18)
			
			% Resize and Add Color Bar
			% Sintax Legend: set(gca,'position',[left bottom width height])
			
			subplot(2,2,[2 4]), hold on
			t = [1:1:100];
			t = [t;t;t;t;t];
			imagesc(t')
			set(gca,'XTick',[]);
			set(gca,'YDir','normal');
			set(gca,'YTick',[0:10:100]);
			set(gca,'YTick',[get(gca,'YTick'),1]); % In order to display magnVec(1) value
			massimo = max(max(data));
			minimo = min(min(data));
			magnVec = floor([minimo:(massimo-minimo)/10:massimo]*100)/100; % In order to have only 2 decimal digit
			set(gca,'YTickLabel',magnVec);
			xlabel('Ratio','FontSize',8)
			set(gca,'position',[0.96,0.11,0.02,0.815])
			
			subplot(2,2,1)
			set(gca,'position',[0.12,0.58384,0.77,0.34116])
			
			subplot(2,2,3)
			set(gca,'position',[0.12,0.11,0.77,0.34116])
			
			% Print output graph
			if (output)
				print -depsc traceout.eps
			endif
		
		otherwise
		
			printf("\nWARNING!\nInvalid first parameter\n\n");
			
	endswitch
	
endif


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